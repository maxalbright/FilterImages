#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;


void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);

void convolve(unsigned char out[][SIZE][3], unsigned char in[][SIZE][3], int N,
              double kernal[][11]);

void sobel(unsigned char out[][SIZE][3], unsigned char in[][SIZE][3]);

void gaussian(double k[][11], int N, double sigma);

void gaussian_filter(unsigned char output[][SIZE][3], 
                     unsigned char input[][SIZE][3], int N, double sigma);

void unsharp(unsigned char output[][SIZE][3], unsigned char input[][SIZE][3],
             int N, double sigma, double alpha);





#ifndef AUTOTEST

int main(int argc, char* argv[])
{
   //First check argc
  if(argc < 3)
    {
      //we need at least ./filter <input file> <filter name> to continue
      cout << "usage: ./filter <input file> <filter name> <filter parameters>";
      cout << " <output file name>" << endl;
      return -1;
    }
   //then check to see if we can open the input file
   unsigned char input[SIZE][SIZE][RGB];
   unsigned char output[SIZE][SIZE][RGB];
   char* outfile;
   int N;
   double sigma, alpha;
   //double kernel[11][11];

   // read file contents into input array
   int status = readRGBBMP(argv[1], input); 
   if(status != 0)
   {
      cout << "unable to open " << argv[1] << " for input." << endl;
      return -1;
   }
   //Input file is good, now look at next argument
   if( strcmp("sobel", argv[2]) == 0)
   {
     sobel(output, input);
     outfile = argv[3];
   }
  
   else if( strcmp("blur", argv[2]) == 0)
   {
     if(argc < 6)
       {
	 cout << "not enough arguments for blur." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     outfile = argv[5];
     gaussian_filter(output, input, N, sigma);
   }
   else if( strcmp("unsharp", argv[2]) == 0)
   {
     if(argc < 7)
       {
	 cout << "not enough arguments for unsharp." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     alpha = atof(argv[5]);
     outfile = argv[6];
     unsharp(output, input, N, sigma, alpha);

   }
   
   else if( strcmp("dummy", argv[2]) == 0)
   {
     //do dummy
     dummy(output, input);
     outfile = argv[3];
   }
   else
   {
      cout << "unknown filter type." << endl;
      return -1;
   }

   if(writeRGBBMP(outfile, output) != 0)
   {
      cout << "error writing file " << argv[3] << endl;
   }
}   

#endif 



// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()

void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   for (int i = 0; i < 3; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         k[i][j] = 0;
      }
   }
   k[1][1] = 1;
   convolve(out, in, 3, k);
}


// Convolves an input image with an NxN kernel to produce the output kernel
// You will need to complete this function by following the 
//  instructions in the comments
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
	      int N, double kernel[][11])
{
 
   int padded[SIZE+10][SIZE+10][RGB];  // Use for input image with appropriate 
                                       // padding
   int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel 
                                       // values then copy from temp to out, 
                                       // applying clamping 

   //first set all of padded to 0 (black)
  for (int row=0; row <= (SIZE + 10 - 1); row++) {
    for (int column=0; column <= (SIZE + 10 - 1); column++) {
        padded[row][column][0] = 0;
        padded[row][column][1] = 0;
        padded[row][column][2] = 0;
    }
  }

   //now copy input into padding to appropriate locations
  for (int row=0; row < SIZE; row++) {
    for (int column=0; column < SIZE; column++) {
      padded[row + N/2][column + N/2][0] = in[row][column][0];
      padded[row + N/2][column + N/2][1] = in[row][column][1];
      padded[row + N/2][column + N/2][2] = in[row][column][2];
    }
  }
   //initialize temp pixels to 0 (black)
  for (int row=0; row < SIZE; row++) {
    for (int column=0; column < SIZE; column++) {
      temp[row][column][0] = 0;
      temp[row][column][1] = 0;
      temp[row][column][2] = 0;
    }
  }
  //now perform convolve (using convolution equation on each pixel of the 
  // actual image) placing the results in temp (i.e. unclamped result)
  for (int row = 0; row <= 255; row++) {
    for (int column = 0; column <= 255; column++) {
      for (int color=0; color < 3; color++) {
        for (int i = -1 * N/2; i <= N/2; i++) {
          for (int j = - 1 * N/2; j <= N/2; j++) {
            temp[row][column][color] += padded[i + row + N/2][j + column + N/2]
              [color] * kernel[i+N/2][j+N/2];
          }
        }
      }
    }
  }


   //now clamp and copy to output 
  for (int row = 0; row <= 255; row++) {
    for (int column = 0; column <= 255; column++) {
      for (int color=0; color < 3; color++) {
        if (temp[row][column][color] < 0) {
          temp[row][column][color] = 0;
        }
        else if (temp[row][column][color] > 255) {
          temp[row][column][color] = 255;
        }
        out[row][column][color] =  temp[row][column][color];
      }
    }
  }
}


// This is a certain type of filter that is used to highlight edges

void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   double s_h1[3][3] = { {-1, 0, 1}, 
                         {-2, 0, 2}, 
                         {-1, 0, 1} };
   double s_h2[3][3] = { {1, 0, -1}, 
                         {2, 0, -2}, 
                         {1, 0, -1} };
   
   unsigned char h1_sobel[SIZE][SIZE][RGB]; //hold intemediate images
   unsigned char h2_sobel[SIZE][SIZE][RGB]; 

   for (int i = 0; i < 11; i++)
   {
      for(int j=0; j < 11; j++)
      {
         k[i][j] = 0;
      }
   }


   // Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      k[i][j] = s_h1[i][j];
    }
  }
   // Call convolve to apply horizontal sobel kernel with result in h1_sobel
  convolve(out, in, 3, k);
  for (int row=0; row<SIZE; row++) {
    for (int column=0; column<SIZE; column++) {
      for (int color=0; color < 3; color++) {
        h1_sobel[row][column][color] = out[row][column][color];
      }
    }
  }
   // Copy in 2nd 3x3 horizontal sobel kernel (i.e. copy s_h2 into k)
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      k[i][j] = s_h2[i][j];
    }
  }

   // Call convolve to apply horizontal sobel kernel with result in h2_sobel
  convolve(out, in, 3, k);
  for (int row=0; row<SIZE; row++) {
    for (int column=0; column<SIZE; column++) {
      for (int color=0; color < 3; color++) {
        h2_sobel[row][column][color] = out[row][column][color];
      }
    }
  }
  
   // Add the two results (applying clamping) to produce the final output 
  for (int row=0; row<SIZE; row++) {
    for (int column=0; column<SIZE; column++) {
      for (int color=0; color < 3; color++) {
        out[row][column][color] = h1_sobel[row][column][color]
          + h2_sobel[row][column][color];
        if (out[row][column][color] > 255) {
          out[row][column][color] = 255;
        }
        else if (out[row][column][color] < 0) {
          out[row][column][color] = 0;
        }
      }
    }
  }
}
// Creates Kernal for Gaussian Blur
void gaussian(double k[][11], int N, double sigma) {
  
  double sum=0;
  // fill in raw values
  for (int y = 0; y<N; y++) {
    for (int x = 0; x<N; x++) {
      double value = (x - N/2) * (x - N/2) / (2 * sigma * sigma) + (y - N/2) * 
        (y - N/2) / (2 * sigma * sigma); 
      k[y][x] = exp(-1 * value);
      sum += k[y][x];
    }
  }
  // normalize
  for (int y = 0; y<N; y++) {
    for (int x = 0; x<N; x++) {
      k[y][x] /= sum;
    }
  }
}


//Applies Gaussian Blur filter by convolving with a new kernal
void gaussian_filter(unsigned char output[][SIZE][3], 
                     unsigned char input[][SIZE][3], int N, double sigma) {
  double k [11][11];
  // calls gaussian to get kernel values 
  gaussian(k, N, sigma);
  // convolves and clamps
  convolve(output, input, N, k);
}

// Applies a sharpened filter by subtracting the color values of a blurred image
// from the original image
void unsharp(unsigned char output[][SIZE][3], unsigned char input[][SIZE][3],
             int N, double sigma, double alpha) {
  // blurs the whole picture and stores it in output
  gaussian_filter(output, input, N, sigma);
  double detail;
  // for each color value of each pixel, the detail is calculated 
  for (int row=0; row < SIZE; row++) {
    for (int column=0; column < SIZE; column++) {
      for (int color=0; color < 3; color++) {
        detail = alpha * (input[row][column][color] - 
                                 output[row][column][color]);
        // clamps and changes output to original value + detail
        if (input[row][column][color] + detail > 255) {
          output[row][column][color] = 255;
        }
        else if (input[row][column][color] + detail < 0) {
          output[row][column][color] = 0;
        }        
        else {
          output[row][column][color] = (int)(input[row][column][color]
            + detail);
        }
      }
    }
  } 
}