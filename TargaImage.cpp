///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <map>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial

bool compare_by_value(const pair<float, int>& p1, const pair<float, int>& p2) {
    return p1.first < p2.first;
}

float distance_formula(float x1, float y1, float z1, float x2, float y2, float z2) {
    return sqrt( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2) );
}

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (data) {
        for (int i = 0; i < height * width; i++) {
            int gray_value = data[i * 4] * 0.299 + data[i * 4 + 1] * 0.587 + data[i * 4 + 2] * 0.114;
            data[i * 4 + 2] = data[i * 4 + 1] = data[i * 4] = gray_value; //;
        }
        return true;
    }
    else {
        ClearToBlack();
        return false;
    }
}// To_Grayscale

vector<pair<float, int>> Normalize_Grayscale(unsigned char* data, int pixels) {
    vector<pair<float, int>> v;
    for (int i = 0; i < pixels; i++) {
        v.push_back(make_pair(data[4 * i] / (float)256, i));
    }
    return v;
}


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if (data) {
        for (int i = 0; i < height * width; i++) {
            data[i * 4] = (data[i * 4] / 32) * 32;
            data[i * 4 + 1] = (data[i * 4 + 1] / 32) * 32;
            data[i * 4 + 2] = (data[i * 4 + 2] / 64) * 64;
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    if (data) {
        // Step down to 32 colors R G B
        for (int i = 0; i < height * width; i++) {
            for (int j = 0; j < 3; j++) {
                data[i * 4 + j] = (data[i * 4 + j] / 8) * 8;
            }
        }
        // Create a map of colors to counts
        map<int, int> hist;

        for (int i = 0; i < height * width; i++) {
            int key = data[i * 4] * 1000000 + data[i * 4 + 1] * 1000 + data[i * 4 + 2];
            if ( hist.find(key) == hist.end() ) {
                hist[key] = 1;
             }
            else {
                hist[key]++;
            }
        }

        // Find most popular colors from the map
        // create an empty vector of pairs
        vector<pair<int, int>> hist_vec;

        // copy key-value pairs from the map to the vector
        copy(hist.begin(), hist.end(),
             back_inserter<std::vector<pair<int,int>>>(hist_vec)
        );

        // sort the vector by increasing the order of its pair's second value
        // if the second value is equal, order by the pair's first value
        sort(hist_vec.begin(), hist_vec.end(),
            [](const pair<int,int>& l, const pair<int,int>& r)
            {
                if (l.second != r.second) {
                    return l.second > r.second;
                }

                return l.first > r.first;
            });

        int colors[256];

        for (int i = 0; i < 256; i++) {
            colors[i] = hist_vec[i].first;
        }

        for (int i = 0; i < height * width; i++) {
            int r = data[i * 4], g = data[i * 4 + 1], b = data[i * 4 + 2], min_r, min_g, min_b;

            float current_min_distance = 100000.;

            for (int j = 0; j < 256; j++) {

                int t_r = colors[j] / 1000000;
                int t_g = (colors[j] % 1000000) / 1000;
                int t_b = (colors[j] % 1000000) % 1000;
                float distance = distance_formula(r, g, b, t_r, t_g, t_b);

                if (distance < current_min_distance) {
                    min_r = t_r;
                    min_g = t_g;
                    min_b = t_b;
                    current_min_distance = distance;
                }
            }

            data[i * 4] = min_r;
            data[i * 4 + 1] = min_g;
            data[i * 4 + 2] = min_b;
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (data) {
        To_Grayscale();
        vector<pair<float, int>> norm_data = Normalize_Grayscale(data, width * height);
        for (int i = 0; i < norm_data.size(); i++) {
            if (norm_data[i].first < 0.5) {
                data[4*i] = data[4*i + 1] = data[4*i + 2] = 0;
            }
            else {
                data[4 * i] = data[4 * i + 1] = data[4 * i + 2] = 255;
                } 
        }
    }
    else {
        ClearToBlack();
        return false;
    }

}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//      
//      Gaussian noise is good for commercial applications. That will help get a better result.
//      
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    if (data) {
        To_Grayscale();
        vector<pair<float, int>> norm_data = Normalize_Grayscale(data, width * height);
        for (int i = 0; i < norm_data.size(); i++) {
            float rand_val = 0.4 * (((float)rand() / RAND_MAX) - 0.5);
            if (( norm_data[i].first + rand_val ) < 0.5) {
                data[4 * i] = data[4 * i + 1] = data[4 * i + 2] = 0;
            }
            else {
                data[4 * i] = data[4 * i + 1] = data[4 * i + 2] = 255;
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    ClearToBlack();
    return false;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation. 
//  
//  Normalize vector 
//
//  Need to sort the vector and assign the bottom 40% as 0 and top 60% as 1
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    if (data) {
        To_Grayscale();
        vector<pair<float, int>> norm_data = Normalize_Grayscale(data, width * height);

        float running_sum = 0, average_intensity;
        for (auto x : norm_data) {
            running_sum += x.first;
        }
        average_intensity = running_sum / norm_data.size();

        // Sort pixels and populate the top average_intensity % with 1
        sort(norm_data.begin(), norm_data.end(), &compare_by_value);

        for (int i = 0; i < norm_data.size(); i++) {
            if (i <= (1 - average_intensity) * norm_data.size()) {
                data[4 * norm_data[i].second] = data[4 * norm_data[i].second + 1] = data[4 * norm_data[i].second + 2] = 0;
            }
            else {
                data[4 * norm_data[i].second] = data[4 * norm_data[i].second + 1] = data[4 * norm_data[i].second + 2] = 255;
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{

    float mask[4][4] = {
        {0.7500, 0.3750, 0.6250, 0.2500 },
        {0.0625, 1.0000, 0.8750, 0.4375 },
        {0.5000, 0.8125, 0.9375, 0.1250 },
        {0.1875, 0.5625, 0.3125, 0.6875 },
    };
    if (data) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if ( data[(i * width + j) * 4] < mask[i % 4][j % 4] * 255.) {
                    data[(i * width + j) * 4] = data[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 2] = 0;
                }
                else {
                    data[(i * width + j) * 4] = data[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 2] = 255;
                }
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    if (data) {
        for (int i = 0; i < height * width; i++) {
            float alpha_f = data[i * 4 + 3] / 255.;
            for (int j = 0; j < 4; j++) {
                data[i * 4 + j] = (data[i * 4 + j] + (1.0 - alpha_f) * pImage->data[i * 4 + j]);
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }


}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    if (data) {
        for (int i = 0; i < height * width; i++) {
            float alpha_g = pImage->data[i * 4 + 3] / 255.;
            for (int j = 0; j < 4; j++) {
                data[i * 4 + j] = alpha_g*data[i * 4 + j];
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    if (data) {
        for (int i = 0; i < height * width; i++) {
            float alpha_g = pImage->data[i * 4 + 3] / 255.;
            for (int j = 0; j < 4; j++) {
                data[i * 4 + j] = (1 - alpha_g) * data[i * 4 + j];
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }

}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    if (data) {
        for (int i = 0; i < height * width; i++) {
            float alpha_f = data[i * 4 + 3] / 255.;
            float alpha_g = pImage->data[i * 4 + 3] / 255.;
            for (int j = 0; j < 4; j++) {
                data[i * 4 + j] = alpha_g * data[i * 4 + j] + (1.0 - alpha_f) * pImage->data[i * 4 + j];
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    if (data) {
        for (int i = 0; i < height * width; i++) {
            float alpha_f = data[i * 4 + 3] / 255.;
            float alpha_g = pImage->data[i * 4 + 3] / 255.;
            for (int j = 0; j < 4; j++) {
                data[i * 4 + j] = (1.0 - alpha_g) * data[i * 4 + j] + (1.0 - alpha_f) * pImage->data[i * 4 + j];
            }
        }
    }
    else {
        ClearToBlack();
        return false;
    }


}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference

bool TargaImage::Filter_5x5(float filter[5][5]){

    unsigned char* tmp = new unsigned char[width * height * 4];

    if (data) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                float sum_r = 0., sum_g = 0., sum_b = 0.;
                // k represents moving through the matrix
                for (int k = -2; k < 3; k++) {
                    for (int l = -2; l < 3; l++) {
                        int height_idx;
                        int width_idx;
                        // Handling edge cases with if statements
                        if (i + k < 0) {
                            // In this case k must be <= -1
                            height_idx = -1 * k;
                        }
                        else if (i + k >= height) {
                            // In this case k must be >= 1
                            height_idx = height - k;
                        }
                        else {
                            height_idx = i + k;
                        }

                        if (j + l < 0) {
                            // In this case l must be <= -1
                            width_idx = -1 * l;
                        }
                        else if (j + l >= width) {
                            // In this case l but >= 1
                            width_idx = width - l;
                        }
                        else {
                            width_idx = j + l;
                        }

                        sum_r += data[(height_idx * width + width_idx) * 4] * filter[k + 2][l + 2];
                        sum_g += data[(height_idx * width + width_idx) * 4 + 1] * filter[k + 2][l + 2];
                        sum_b += data[(height_idx * width + width_idx) * 4 + 2] * filter[k + 2][l + 2];
                    }
                }
                tmp[(i * width + j) * 4] = (unsigned char) round(sum_r);
                tmp[(i * width + j) * 4 + 1] = (unsigned char) round(sum_g);
                tmp[(i * width + j) * 4 + 2] = (unsigned char) round(sum_b);
            }
        }
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                    data[(i * width + j) * 4] = tmp[(i * width + j) * 4];
                    data[(i * width + j) * 4 + 1] = tmp[(i * width + j) * 4 + 1];
                    data[(i * width + j) * 4 + 2] = tmp[(i * width + j) * 4 + 2];
            }
        }
        delete[] tmp;

        // Clamping values
        for (int i = 0; i < width * height; i++) {
            for (int j = 0; j <= 2; j++) {
                if (data[4 * i + j] > 255) {
                    data[4 * i + j] = 255;
                }
                else if (data[4 * i + j] < 0) {
                    data[4 * i + j] = 0;
                }
            }
        }
        return true;
    }
    else {
        ClearToBlack();
        return false;
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    float filter[5][5] = {
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04}
    };

    if (data) {
        return TargaImage::Filter_5x5(filter);
    }
    else {
        ClearToBlack();
        return false;
    }
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    float filter[5][5] = {
        {1. / 81., 2. / 81., 3. / 81., 2. / 81., 1. / 81.},
        {2. / 81., 4. / 81., 6. / 81., 4. / 81., 2. / 81.},
        {3. / 81., 6. / 81., 9. / 81., 6. / 81., 3. / 81.},
        {2. / 81., 4. / 81., 6. / 81., 4. / 81., 2. / 81.},
        {1. / 81., 2. / 81., 3. / 81., 2. / 81., 1. / 81.}
    };


    if (data) {
        return TargaImage::Filter_5x5(filter);
    }
    else {
        ClearToBlack();
        return false;
    }
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    float filter[5][5] = {
        {1. / 256., 4. / 256., 6. / 256., 4. / 256., 1. / 256.},
        {4. / 256., 16. / 256., 24. / 256., 16. / 256., 4. / 256.},
        {6. / 256., 24. / 256., 36. / 256., 24. / 256., 6. / 256.},
        {4. / 256., 16. / 256., 24. / 256., 16. / 256., 4. / 256.},
        {1. / 256., 4. / 256., 6. / 256., 4. / 256., 1. / 256.},
    };

    if (data) {
        return TargaImage::Filter_5x5(filter);
    }
    else {
        ClearToBlack();
        return false;
    }
   
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
        
        if (data) {

            float filter[5][5] = {
            {-1. / 256., -4. / 256., -6. / 256., -4. / 256., -1. / 256.},
            {-4. / 256., -16. / 256., -24. / 256., -16. / 256., -4. / 256.},
            {-6. / 256., -24. / 256., 220. / 256., -24. / 256., -6. / 256.},
            {-4. / 256., -16. / 256., -24. / 256., -16. / 256., -4. / 256.},
            {-1. / 256., -4. / 256., -6. / 256., -4. / 256., -1. / 256.},
            };

            unsigned char* original_data = new unsigned char[width * height * 4];
            memcpy(original_data, data, sizeof(unsigned char) * width * height * 4);

            TargaImage::Filter_5x5(filter);

            // Clamping values
            for (int i = 0; i < width * height; i++) {
                for (int j = 0; j <= 2; j++) {
                    if (data[4 * i + j] > 255) {
                        data[4 * i + j] = 255;
                    }
                    else if (data[4 * i + j] < 0) {
                        data[4 * i + j] = 0;
                    }
                }
            }
            return true;
        }
        else {
            ClearToBlack();
            return false;
        }

}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

