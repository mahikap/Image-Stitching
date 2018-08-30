# Image-Stitching
This algorithm takes two images and uses feature extraction + homography estimation to stitch both images, forming a panoramic effect.

Descriptor vectors are extracted using the VLFeat library: http://www.vlfeat.org/

This is an example of 4 images being stitched when this algorithm was called multiple times. The original images can be found in ./TestImages
![Panoroma with 4 images](/TestImages/stitchedimages1234.jpg)
