lib3dreconstruct
================

This library contains methods used in 3D reconstructions from photographies.Actually, it contains an interest point detector based on Harris corner detector, a basic method to find stereo correspondences in two photographies and method to find 3D point from stereo correspondences.

The stero matches algorithm can be summarized as follow : first we look for interest points in each photography using Harris corner detection. This give us two list of point that could possibly match. Then, for each point of the first list, we look for the point in the second list that maximize the Normalized Crossed Correlation (NCC). If we found a pair of point with NCC > threshold, we consider that the pair is a stereo correspondence. This method will be improved in a further work.

Concerning the 3D reconstruction, the algorithm is summarized as follows : Given a pair of photographies, the associated projection matrices and a pair of stereo correspondences, we use a method based on DLT to find the point in space that originated the pair of stereo correspondences.
