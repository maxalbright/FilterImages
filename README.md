This program applies a filter to a given 255x255 image, like those in the repository, given certain parameters.
It uses kernals and kernal math to apply the filters, padding the edges temporarily and clamping for RGB values that exceed 255. 

To compile, just compile filter.cpp and bmp.cpp

To apply the dummy filter (which does nothing):
./executableName inputFilename dummy outputFilename

To apply the sobel filter (which highlights the edges of a picture:
./executableName inputFilename sobel outputFilename

To apply a Gaussian blur:
./executableName inputFilename blur N sigma outputFilename
Where N is the width of the Kernal that the filter uses and sigma is taken from the following equation:
G(x,y)={{1} / {2 * pi * sigma ^{2}}} * e^{-{{x^{2}+y^{2}} / {2 * sigma ^{2}}}}
where x is the horizontal distace from the center of the kernal, y is the vertical, and G is the kernal value 
before being normalized. 

To apply a sharpening filter:
./executableName inputFilename unsharp N sigma alpha outputFilename
Where N is the width of the Kernal, sigma is the value from the guassian equation, and alpha is what percent of
the Gaussian Blur is subtracted from the original image.
This filter works by creating a blurred version of the original image with the gaussian function, and
subtracting some fraction of those RGB values from the original image.
