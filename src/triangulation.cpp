/*
 *  triangulate.cpp
 *
 *  Copyright (C) 2014 Foxel www.foxel.ch
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Feb 20, 2014
 *      Authors: s.flotron@foxel.ch, luc.deschenaux@foxel.ch
 */
#include "triangulation.hpp"

#define M 4
#define N 4

// convert unsigned char to double (used in automatic matching)
double charToDouble(unsigned char* c){
     return (double)(*c) ;
}

// extract image basename
std::string path(const std::string &filename) {
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot);
}

/*--------------------------------------------------------------------------
 *
 *  callSVD : compute a svd decomposition of matrix A using svd routine from lapack
 *
 *  inputs  : p  : number of row of matrix A
 *            q  : number of cols of matrix A
 *            a  : matrix A^T, stored in an array. Must be transposed for fortran linkage !
 *            jobu = 'A'. Return matrix U.
 *            jobu = 'S'. U contain only left singular vectors
 *            jobu = 'N'  U is not computed
 *            jobv = 'A'. Return matrix Vt.
 *            jobv = 'S'. Vt contain only right singular vectors
 *            jobv = 'N'  Vt is not computed
 *
 *  outputs :  s : array of singular values
 *             u : left singular vectors
 *             vt : right singular vectors
 *
 *  warning : the (row,col) storage in fortran is the transposed of c++. So
 *            you must give the transpose of A as input ! Moreover the outputs
 *            arrays u and vt are the transpose of the matrix U, Vt !
 *            Be aware of that in your computations !!!
 *
 * --------------------------------------------------------------------------
 */

void callSVD(int* p, int* q, double* a, double* s, double* u, double* vt, char jobu, char jobv){
	// loading and creating paramter for svd call (see doc of lapack/dgesvd for more details)
     int m = *p, n = *q, lda = *p, ldu = *p, ldvt = *q;
	 int    lwork(-1);  // size of the array work, which is a debbuging array.
	 double wkopt;      // debugging array
     double* work;
	 int info;          // return value of the routine. If =0 success, >0 did not converge, <0 big problems.

     /* Executable statements */
 //    printf( " DGESVD Example Program Results\n" );
     /* Query and allocate the optimal workspace */

     dgesvd_(&jobu,&jobv, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork,&info);
     lwork = (int)wkopt;
     work = (double*)malloc( lwork*sizeof(double) );
     /* Compute SVD */
     dgesvd_(&jobu,&jobv, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info );

     /* Check for convergence */
     if( info > 0 ) {
             printf( "The algorithm computing SVD failed to converge.\n" );
             exit( 1 );
     }
#if 0
     /* Print singular values */
     print_matrix( "Singular values", 1, min(n,m), s, 1 );
     cout << " smallest eigenvalue = " << s[min(n,m)-1] << endl;
     /* Print right singular vectors */
     print_matrix( "Right singular vectors (stored rowwise)", n, n, Vt, ldvt );
     /* Free workspace */
#endif
     free( (void*)work );

}

/* ---------------------------------------------------------------------------------------
   find_match_stat (double (&p0)[2], double (&p1)[2])

   Given a pixel (u,v) in first image, find the pixel (um,vm) in the second image
   using statistical comparison.

   inputs : p0=(u,v)      : pixels in first image
            imageZeroPath : path of first image
            imageOnePath  : path of second image

   output : p1=(um,vm)   : pixels in second image that match the pixel (u,v) of first image.
------------------------------------------------------------------------------------------
*/

void findMatchStat(double (&p0)[3], double (&p1)[3], const char *imageZeroPath, const char *imageOnePath, double& corr){
	// boolean variable used to check if a correspondence is found
	bool found_match(0);

	// load images
    CImg<unsigned char> img0;
    CImg<unsigned char> img1;

    std::string image0path;
    image0path=path(imageZeroPath)+".tiff";

    std::string image1path;
    image1path=path(imageOnePath)+".tiff";

    img0=CImg<unsigned char>(image0path.c_str());
    img1=CImg<unsigned char>(image1path.c_str());

    int width=img1.width(); int height=img1.height();
    int neighSize = 5; // size of surrounding patch (in pixels)
    corr=-1.0;  // correlation of the points. If close to one : good.

    // compute mean and RMS of first patch
    double  p0_mean[3]={0,0,0};
    double  p0_rms(0.0);

    // compute average
    for(int i(-neighSize); i<neighSize+1; ++i){
    	for(int j(-neighSize); j<neighSize+1; ++j){
    		if((p0[0]+i>=0) && (p0[0]+i <= width) && (p0[1]+j>=0) && (p0[1]+j< height)){
				double r=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,0));
				double g=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,1));
				double b=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,2));

				p0_mean[0] += r/(2.0*neighSize+1)/(2.0*neighSize+1) ;
				p0_mean[1] += g/(2.0*neighSize+1)/(2.0*neighSize+1);
				p0_mean[2] += b/(2.0*neighSize+1)/(2.0*neighSize+1);
    		}
    	}
    }

    // compute rms
    for(int i(-neighSize); i<neighSize+1; ++i){
    	for(int j(-neighSize); j<neighSize+1; ++j){
    		if((p0[0]+i>=0) && (p0[0]+i <= width) && (p0[1]+j>=0) && (p0[1]+j< height)){
				double r=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,0));
				double g=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,1));
				double b=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,2));

				p0_rms += (r-p0_mean[0])*(r-p0_mean[0])/((2.0*neighSize+1)*(2.0*neighSize+1));
				p0_rms += (g-p0_mean[1])*(g-p0_mean[1])/((2.0*neighSize+1)*(2.0*neighSize+1));
				p0_rms += (b-p0_mean[2])*(b-p0_mean[2])/((2.0*neighSize+1)*(2.0*neighSize+1));
    		}
    	}
    }
    double rms0 = sqrt(p0_rms);
    int um=0; int vm=0;

    // check if root mean square is not too small
    if( rms0 < 1e-5){
    	cout << " Variance of ring too small, return default pixel (0,0) " << endl;
    	p1[0] = 0; p1[1] = 0;
    }
    // rms is big enough to make a search in second image
    else{
        // look for a matching pixel in second image using Normalized Cross Correlation (NCC)
    	int u(p1[0]), v(p1[1]);
//		for(int u(0); u< width; ++u){
//			for(int v(0); v < height; ++v){
				double p1_mean[3] = {0,0,0};
				double p1_rms = 0.0;

				// compute average
				for(int i(-neighSize); i<neighSize+1; ++i){
					for(int j(-neighSize); j<neighSize+1; ++j){
						if( (u+i>=0) && (u+i < width) && (v+j>=0) && (v+j< height) ){
							double r=charToDouble(img1.data(u+i,v+j,0,0));
							double g=charToDouble(img1.data(u+i,v+j,0,1));
							double b=charToDouble(img1.data(u+i,v+j,0,2));

							p1_mean[0] += r/(2.0*neighSize+1)/(2.0*neighSize+1) ;
							p1_mean[1] += g/(2.0*neighSize+1)/(2.0*neighSize+1);
							p1_mean[2] += b/(2.0*neighSize+1)/(2.0*neighSize+1);
						};
					};
				}

				// compute rms
				for(int i(-neighSize); i<neighSize+1; ++i){
					for(int j(-neighSize); j<neighSize+1; ++j){
						if((u+i>=0) && (u+i < width) && (v+j>=0) && (v+j< height)){
							double r=charToDouble(img1.data(u+i,v+j,0,0));
							double g=charToDouble(img1.data(u+i,v+j,0,1));
							double b=charToDouble(img1.data(u+i,v+j,0,2));

							p1_rms += (r-p1_mean[0])*(r-p1_mean[0])/((2.0*neighSize+1)*(2.0*neighSize+1));
							p1_rms += (g-p1_mean[1])*(g-p1_mean[1])/((2.0*neighSize+1)*(2.0*neighSize+1));
							p1_rms += (b-p1_mean[2])*(b-p1_mean[2])/((2.0*neighSize+1)*(2.0*neighSize+1));
						}
					}
				}
				double rms1 = sqrt(p1_rms);

				double crossCorr=0.0;

				// if rms1 is big enough, compute cross correlation associated to pixel (u,v)
				if(rms1 > 1e-5){
					found_match=1;

					// compute cross correlation
					for(int i(-neighSize); i<neighSize+1; ++i){
						for(int j(-neighSize); j<neighSize+1; ++j){
							if((u+i>=0) && (u+i <= width) && (v+j>=0) && (v+j< height)){
								if((p0[0]+i>=0) && (p0[0]+i <= width) && (p0[1]+j>=0) && (p0[1]+j< height)){
									double r0=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,0));
									double g0=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,1));
									double b0=charToDouble(img0.data(p0[0]+i,p0[1]+j,0,2));

									double r1=charToDouble(img1.data(u+i,v+j,0,0));
									double g1=charToDouble(img1.data(u+i,v+j,0,1));
									double b1=charToDouble(img1.data(u+i,v+j,0,2));

									crossCorr += (r0-p0_mean[0])*(r1-p1_mean[0])/((2.0*neighSize+1)*(2.0*neighSize+1))
											   + (g0-p0_mean[1])*(g1-p1_mean[1])/((2.0*neighSize+1)*(2.0*neighSize+1))
											   + (b0-p0_mean[2])*(b1-p1_mean[2])/((2.0*neighSize+1)*(2.0*neighSize+1));
								}
							}
						}
					} // end cross corr

					crossCorr /= rms1*rms0;

					// update output values
					if(crossCorr > corr){
						corr = crossCorr;
						um = u; vm=v;
					}
				}
//			}
		}
//    }

    // print a warning if we don't found a match
    if(found_match==0){
 //   	cout << " found no match in second picture, return default value (0,0)" << endl;
    	p1[0] = 0; p1[1]=0;
    } else {
        p1[0] = um; p1[1]=vm;
    }
}

/* ---------------------------------------------------------------------------------------
   find_match(int u, int v, int um, int vm, double (&P0)[3][4], double (&P1)[3][4],
		double (&K0)[9], double (&K1)[9], double (&F)[9], double (&R0)[3][3], double (&R1)[3][3])

   Given a pixel (u,v) in first image, find the pixel (um,vm) in the second image
   such that the reprojection error is minimized. Equally check that the pixel (um,vm) that we found
   is the right one using the fundamental matrix constraint.

   inputs : (u,v)     : pixels in first image
             pixWidth : pixels width in rectified second image
             P0, P1   : projection matrices of camera 0 and 1 (expressed in cam0 referential)
             K0, K1   : camera matrices of camera 0 and 1
             F        : fundamental matrix
             R0, R1   : rotation matrices used for rectification of image 0 and 1 resp.

   output : (um,vm)   : pixels in second image that match the pixel (u,v) of first image.

------------------------------------------------------------------------------------------
*/

void findMatch(double (&pm)[2], double u, double v, int pixWidth, double (&P0)[3][4], double (&P1)[3][4],
		double (&K0)[9], double (&K1)[9], double (&F)[9], double (&R0)[3][3], double (&R1)[3][3],
		const char *imageZeroPath, const char *imageOnePath){

	// create point p0 used in DLT triangulation
	double p0[3] = {u,v,1.0};

	// load some inputs used in rectification
    double ful  = K0[0]; double fvl  = K0[4];
    double uol  = K0[2]; double vol  = K0[5];

    double fur  = K1[0]; double fvr  = K1[4];
    double uor  = K1[2]; double vor  = K1[5];

    double fnu = ful; double fnv = fvl;
    double up0 = uol; double vp0 = vol;

    // load rotation matrices
    Matx<double,3,3> RT0(
         R0[0][0], R0[1][0], R0[2][0],
         R0[0][1], R0[1][1], R0[2][1],
         R0[0][2], R0[1][2], R0[2][2]
     );

    // transpose of second rotation matrice
    Matx<double,3,3> Rot1(
        R1[0][0], R1[0][1], R1[0][2],
        R1[1][0], R1[1][1], R1[1][2],
        R1[2][0], R1[2][1], R1[2][2]
     );

	// in rectified image, we only search matching pixel on a horizontal line
	// so we compute the rectified pixel corresponding to (u,v)
    Matx<double,3,1> P ((u-uol)/ful, (v-vol)/fvl, 1);
    Matx<double,3,1> PR (RT0*P);

    double uo = up0 + fnu*PR.val[0] / PR.val[2];
    double vo = vp0 + fnv*PR.val[1] / PR.val[2];

    // look for point on epipolar line that minimize reprojection error
    double error(1.0e10);
    double crossCorr(-1.0);

    for(int s(0); s < 1000*pixWidth; ++s){
    	for(int k(0); k<1; ++k){
    	// go into pixel of second image
        Matx<double,3,1> Q ((s/1000.0-up0)/fnu, (vo+k/10.0-vp0)/fnv, 1);
        Matx<double,3,1> QR (Rot1*Q);

        double ur = uor + fur*QR.val[0] / QR.val[2];
        double vr = vor + fvr*QR.val[1] / QR.val[2];

        if( (ur >= 0) && (ur < 2592) && (vr >=0) && (vr < 1934)){
        // compute point in space with pixels ur,vr
        double X[4]  = {0,0,0,0};
        double p1[3] = {ur,vr,1.0};
        linearDLTTriangulation(X,p0,p1,P0,P1);

        // compute reprojection
        Matx<double,3,4> Pr0(
          P0[0][0], P0[0][1], P0[0][2], P0[0][3],
          P0[1][0], P0[1][1], P0[1][2], P0[1][3],
          P0[2][0], P0[2][1], P0[2][2], P0[2][3]);

        Matx<double,3,4> Pr1(
          P1[0][0], P1[0][1], P1[0][2], P1[0][3],
          P1[1][0], P1[1][1], P1[1][2], P1[1][3],
          P1[2][0], P1[2][1], P1[2][2], P1[2][3]);

        Matx<double,4,1> X0(X);

        Matx<double,3,1> xop(Pr0*X0);
        Matx<double,3,1> x1p(Pr1*X0);

        // compute reprojection error
        double  errProj0 = sqrt((xop.val[0]/xop.val[2]-u) *(xop.val[0]/xop.val[2]-u)
        		              + (xop.val[1]/xop.val[2]-v) *(xop.val[1]/xop.val[2]-v));
        double  errProj1 = sqrt((x1p.val[0]/x1p.val[2]-ur)*(x1p.val[0]/x1p.val[2]-ur)
        		              + (x1p.val[1]/x1p.val[2]-vr)*(x1p.val[1]/x1p.val[2]-vr));

        // compute correlation
        double corr=-1.0;
       // findMatchStat(p0,p1,imageZeroPath, imageOnePath, corr);

        if(((errProj0+errProj1) < error) || (crossCorr < corr)){
              pm[0]=ur; pm[1]=vr; error = (errProj0+errProj1);
              crossCorr=corr;
           }
        }
    }
    }

}

/* ---------------------------------------------------------------------------------------
void linearDLTTriangulation(double (&X)[4], double (&p0)[3], double (&p1)[3],
              double (&P0)[3][4], double (&P1)[3][4]);

   Given a pixel (u,v) in first image and the matching pixel (um,vm) in the second image, we find
   the original point X in space such that (u,v,1) = P0.X and (um,vm,1) = P1.X. This is performed
   using the DLT algorithm (see. Hartley, Zisserman, multiple view geometry in computer vision, ch. 12)

   inputs :  p0 = (u,v,1)    : pixel in first image
	         p1 = (um,vm,1)  : pixel in 2nd image
	         P0, P1          : projection matrices of camera 0 and 1

   output :  X point in R^3 in homogenous coordinate such that (u,v,1)=  P0.X and
             (um,vm,1) = P1.X

------------------------------------------------------------------------------------------
*/

void linearDLTTriangulation(double (&X)[4], double (&p0)[3],
		double (&p1)[3], double (&P0)[3][4], double (&P1)[3][4]){

	// loading and creating paramter for svd call (see doc of lapack/dgesvd for more details)
     int m = 4, n = 4;
	 double S[min(n,m)];    // singular values of A  of size = min(n,m)
	 double U[m*m];   // the matrix U
	 double Vt[n*n]; // the matrix V^T

#if 1
	 double rmsp0 = 933.585; double rmsp1=933.585;
	 double mean_p0[2]={1296,967};
	 double mean_p1[2]={1296,967};
#else
	 // do not normalize input... Give not good stability results !!!
	 double rmsp0 = sqrt(2); double rmsp1=sqrt(2);
	 double mean_p0[2]={0,0};
	 double mean_p1[2]={0,0};
#endif

	 double x=(p0[0]-mean_p0[0])/rmsp0*sqrt(2); 	 double y=(p0[1]-mean_p0[1])/rmsp0*sqrt(2);
	 double xp=(p1[0]-mean_p1[0])/rmsp1*sqrt(2);   double yp=(p1[1]-mean_p1[1])/rmsp1*sqrt(2);

	 // normalization matrices
     Matx<double,3,3> T0(
        sqrt(2)/rmsp0,  0           ,-sqrt(2)*mean_p0[0]/rmsp0,
        0,             sqrt(2)/rmsp0,-sqrt(2)*mean_p0[1]/rmsp0,
        0,              0           ,1.0);

     Matx<double,3,3> T1(
    	        sqrt(2)/rmsp1,  0           ,-sqrt(2)*mean_p1[0]/rmsp1,
    	        0,             sqrt(2)/rmsp1,-sqrt(2)*mean_p1[1]/rmsp1,
    	        0,              0           ,1.0);

     // normalized projection matrices
     Matx<double,3,4> Pun0(
    	 P0[0][0], P0[0][1], P0[0][2], P0[0][3],
    	 P0[1][0], P0[1][1], P0[1][2], P0[1][3],
    	 P0[2][0], P0[2][1], P0[2][2], P0[2][3]);

     Matx<double,3,4> Pun1(
    	 P1[0][0], P1[0][1], P1[0][2], P1[0][3],
    	 P1[1][0], P1[1][1], P1[1][2], P1[1][3],
    	 P1[2][0], P1[2][1], P1[2][2], P1[2][3]);

     Matx<double,3,4> PN0(T0*Pun0);
     Matx<double,3,4> PN1(T1*Pun1);

     double A[M*N] = {
         x*PN0(2,0) - PN0(0,0), y*PN0(2,0) - PN0(1,0),  xp*PN1(2,0) - PN1(0,0), yp*PN1(2,0) - PN1(1,0),
         x*PN0(2,1) - PN0(0,1), y*PN0(2,1) - PN0(1,1),  xp*PN1(2,1) - PN1(0,1), yp*PN1(2,1) - PN1(1,1),
         x*PN0(2,2) - PN0(0,2), y*PN0(2,2) - PN0(1,2),  xp*PN1(2,2) - PN1(0,2), yp*PN1(2,2) - PN1(1,2),
         x*PN0(2,3) - PN0(0,3), y*PN0(2,3) - PN0(1,3),  xp*PN1(2,3) - PN1(0,3), yp*PN1(2,3) - PN1(1,3),
     };

     callSVD(&m, &n, A, S, U, Vt, 'N', 'S');

     // compute output
     for(int i(0); i<n; ++i){
        X[i] = Vt[(n-1)+i*n] / Vt[n*n-1];
     }
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}
