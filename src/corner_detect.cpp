/*
 * corner_detect.cpp
 *
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
 *  Created on: Apr 22, 2014
 *      Authors: s.flotron@foxel.ch
 */

#define DEBUG 1

#include "corner_detect.hpp"

// merged two list sorted in decreasing order...
vector< vector <double> > merge(vector < vector <double > >& p, int P, vector < vector <double > >& q, int Q){

	int i(0); int j(0);

	vector <vector <double> > m;

	while ( i < P && j < Q){
		if( p[i][0] < q[j][0] ){
			m.push_back(q[j]); ++j;
		}
		else{
			m.push_back(p[i]);
			++i;
			if(P == 1){
				while( j < Q){
					m.push_back(q[j]); ++j;
				}
			}
		}
	}

	if( i < P){
		while(i < P){
		   m.push_back(p[i]); ++i;
		}
	}

	if( j < Q){
		while(j < Q){
			m.push_back(q[j]); ++j;
		}
	}

	return  m;
}

/* ---------------------------------------------------------------------------------------
   mergeSort (vector < vector <double> >& v,  const int& n)

   Sort an unsorted list v by decreasing first component

   inputs :  v     : list to sort
             n     : size of list

   output :  v ordered by decreasing first component
------------------------------------------------------------------------------------------
*/

void mergeSort( vector < vector <double> >& v,  const int& n){
	int Q = n/2;
	int P = n-Q;
	vector < vector <double> > q;
	vector < vector <double> > p;

	for(int i(Q); i < n ; ++i)
		p.push_back(v[i]);

	for(int i(0); i < Q; ++i)
		q.push_back(v[i]);

	if(P > 1)
		mergeSort(p,P);
    if(Q >1 )
		mergeSort(q,Q);

	v = merge(p,P,q,Q);
}

/* ---------------------------------------------------------------------------------------
   harrisANMS(const char *grayPath, const char *corners, vector <vector <int> >& points,
		double k, double tol, int nip)

   Given an image, find interest point using Harris corner detector and Adaptive Non Max
   suppression (ANMS) algorithm.

   inputs :  grayPath     : path to grayscale image used to find interest point.
             double k     : Harris response parameters (0.04 is the standard value)
             double tol   : tolerance used for harris response function (3e4 is a good one)
             int nip      : number of interest point to keep ( >1000 in general)

   output :   corners     : image path of image containing interest point (in DEBUG mode)
              points      : array containing interest points
------------------------------------------------------------------------------------------
*/

void harrisANMS(const char *grayPath, const char *corners, vector <vector <int> >& points,
		double k, double tol, int nip){

	   //load image and create output
	    CImg<unsigned char> gray;

	    std::string image0path;
	    image0path=path(grayPath)+".tiff";

	    std::string cornerpath;
	    cornerpath=path(corners)+".tiff";

	    gray=CImg<unsigned char>(image0path.c_str());

	#if DEBUG
	    // create output bitmap
	    CImg<unsigned char> gradX;
	    CImg<unsigned char> gradY;
	    gradX.assign(gray.width(),gray.height(),1,1,0);
	    gradY.assign(gray.width(),gray.height(),1,1,0);
	#endif

	    // gaussian filter for sigma=1 (could automatize this step!!)
	    double gaussFilt[5][5] = {
	    		1/273.0,  4/273.0,  7/273.0,  4/273.0, 1/273.0,
	    		4/273.0, 16/273.0, 26/273.0, 16/273.0, 4/273.0,
	    		7/273.0, 26/273.0, 41/273.0, 26/273.0, 7/273.0,
	    		4/273.0, 16/273.0, 26/273.0, 16/273.0, 4/273.0,
	    		1/273.0,  4/273.0,  7/273.0,  4/273.0, 1/273.0
	    };

	    // smooth intensity before computing derivatives
	    for(int u(0); u< gray.width(); ++u)
	    	for(int v(0); v <gray.height(); ++v){
	    		double smooth = 0.0;

	    		if( (u>1) && (u < gray.width()-2) && (v>1) && (v<gray.height()-2)){
					for(int i(-2); i<3; ++i)
						for(int j(-2); j<3; ++j)
						   smooth  += gaussFilt[2+i][2+j]*( (double) gray(u,v,0,0));
	    		}

				gray(u,v,0,0) = round(smooth);
	    	}

	    // create table product of derivatives
	    vector< vector<double> > gxx;
	    vector< vector<double> > gxy;
	    vector< vector<double> > gyy;
	    vector< vector<double> > gradx;
	    vector< vector<double> > grady;

	    double maxgx(-1.0e10); double mingx(1.0e10);
	    double maxgy(-1.0e10); double mingy(1.0e10);

	    // compute derivatives and product of it (for Harris corner detection)
	    for(int u(0); u < gray.width(); ++u){
	    	vector<double > tmpx;
			vector<double > tmpy;

	    	for(int v(0); v < gray.height(); ++v){
	    		int gx = 0;
	    		int gy = 0;
	    		// compute derivative using finite difference
	    		if( (u>0) && (u < gray.width()-1) && (v>0) && (v<gray.height()-1)){
					gx = ((double) gray(u+1,v,0,0)) -((double) gray(u-1,v,0,0));
					gy = ((double) gray(u,v+1,0,0)) -((double) gray(u,v-1,0,0));
	    		}

	    		if( gx > maxgx) maxgx = gx;
	    		if( gx < mingx) mingx = gx;
	    		if( gy > maxgy) maxgy = gy;
	    		if( gy < mingy) mingy = gy;

	    		// scale
	#if DEBUG
	    		gradX(u,v,0,0) = round(gx);
	    		gradY(u,v,0,0) = round(gy);
	#endif
	    		tmpx.push_back(gx);
	    		tmpy.push_back(gy);

	    	}
	    	gradx.push_back(tmpx);
	    	grady.push_back(tmpy);
	    }

	    // normalize derivative
	    for(int u(0); u < gray.width(); ++u){

	    	vector<double > tmpxx;
			vector<double > tmpxy;
			vector<double > tmpyy;

	    	for(int v(0); v < gray.height(); ++v){
	    		double gx, gy;

	    		gx = gradx[u][v]*255.0/(maxgx-mingx);
	    		gy = grady[u][v]*255.0/(maxgy-mingy);

	    		// scale
	    		//if( gx > 255) gx = 255;
	    		//if( gy > 255) gy = 255;
	#if DEBUG
	    		gradX(u,v,0,0) = round(fabs(gx));
	    		gradY(u,v,0,0) = round(fabs(gy));
	#endif
	    		tmpxx.push_back(gx*gx);
	    		tmpxy.push_back(gx*gy);
	    		tmpyy.push_back(gy*gy);
	    	}
	    	gxx.push_back(tmpxx);
	    	gxy.push_back(tmpxy);
	    	gyy.push_back(tmpyy);
	    }

	#if DEBUG
	    // save image gradient
	    gradX.save_tiff("./gradX_norm.tiff");
	    gradY.save_tiff("./gradY_norm.tiff");
	#endif

	    vector <vector <double > > response;
	    vector <vector <double > > intPoint;

	    // compute Harris response function
	    for(int u(0); u<gray.width(); ++u){
	    	vector <double > tmpresp;

	    	for(int v(0); v<gray.height(); ++v){
	        	vector <double > tmpip;
	    		double gxxfilter=0.0;
	    		double gyyfilter=0.0;
	    		double gxyfilter=0.0;

	    		if( (u>1) && (u < gray.width()-2) && (v>1) && (v<gray.height()-2)){
					// compute convolution
					for(int i(-2); i<3; ++i)
						for(int j(-2); j<3; ++j){
						   gxxfilter  += gaussFilt[2+i][2+j]*gxx[u+i][v+j];
						   gyyfilter  += gaussFilt[2+i][2+j]*gxy[u+i][v+j];
						   gxyfilter  += gaussFilt[2+i][2+j]*gyy[u+i][v+j];
						}
	    		}

	    		double harris_resp = gxxfilter*gyyfilter-gxyfilter*gxyfilter -k*(gxxfilter+gyyfilter)*(gxxfilter+gyyfilter);

	    		// response = |M| -k*Tr(M)^2 where M = ( gxfilter, gxyfilter ; gxyfilter, gyfilter)
	    		tmpresp.push_back(harris_resp);
	    		if( harris_resp < -tol ){
	    			tmpip.push_back( -harris_resp );
	        		tmpip.push_back( (double) u);
	        		tmpip.push_back( (double) v);
	            	intPoint.push_back(tmpip);
	    		}
	    	}
	    	response.push_back(tmpresp);
	    }

	    // Keep non local maximas using ANMS
	    // step one : sort interest point by decreasing harris function
	    mergeSort(intPoint, intPoint.size());

	    // keep only local maximums using adaptive non-maximal suppression algorithm
	    vector < vector <double> > radius;
	    vector < double > first;
	    first.push_back(1e10); first.push_back(intPoint[0][1]); first.push_back(intPoint[0][2]);
	    radius.push_back(first);

	    for(int i(1); i< (int) intPoint.size(); ++i){
	    	double ri(1e10);
	    	int j(0);
	    	double u0 = intPoint[i][1]; double v0 = intPoint[i][2];

	    	while(0.90*intPoint[j][0] > intPoint[i][0] && (j < (int) intPoint.size() ) ){
	        	double u1 = intPoint[j][1]; double v1 = intPoint[j][2];
	        	double dist = sqrt((u1-u0)*(u1-u0)+(v1-v0)*(v1-v0));
	        	if(dist < ri)
	        		ri = dist;
	        	j++;
	    	}
	    	vector <double> tmp;

	    	tmp.push_back(ri); tmp.push_back(u0); tmp.push_back(v0);
	    	radius.push_back(tmp);
	    }

	    // order by decreasing radius in order to keep only first biggest radius
	    mergeSort(radius, radius.size());

	    // export only first nip point radius = interest point
	    points.clear();

	    for(int i(0); i < min(nip, (int) intPoint.size()); ++i){
	    	vector <int> pixels;
	    	if( radius[i][1] > 40 && (radius[i][1] < 2550) && (radius[i][2] > 40) && (radius[i][2] < 1900) ){
	    		pixels.push_back( (int) radius[i][1]); pixels.push_back( (int) radius[i][2]);
	    		points.push_back(pixels);
	    	}
	    }

#if DEBUG
	    // export response function
	    CImg<unsigned char> harris;
	    harris.assign(gray.width(),gray.height(),1,3,0);

	    for(int u(0); u<gray.width(); ++u){
	        for(int v(0); v<gray.height(); ++v){
	        		harris(u,v,0,0) = gray(u,v,0,0);
	        		harris(u,v,0,1) = gray(u,v,0,0);
	        		harris(u,v,0,2) = gray(u,v,0,0);
	        }
	    }

	    for(int i(0);  i< (int) points.size(); ++i){
	    	int u = points[i][0];
	    	int v = points[i][1];
	    	for(int i(-2); i<3; ++i)
	    		for(int j(-2); j<3; ++j){
					harris(u+i,v+j,0,0) = 255;
					harris(u+i,v+j,0,1) = 0;
					harris(u+i,v+j,0,2) = 0;
	    		}
	    }

	    // save image
	    harris.save_tiff(cornerpath.c_str());
#endif
}

/* ---------------------------------------------------------------------------------------
  harrisMatch(vector < vector <double> >&  pixels, const char *grayPath0, const char *grayPath1,
		vector <vector <int> >& ip0, vector <vector <int> >& ip1, double (&F)[9], double NCCtol,
		int neighSize, int dist, double Ftol)

   Given two images and two lists of interest point, find stereo matches between theses lists
   using Normalized Cross Correlation (NCC). We recall that NCC is a number in [-1,1] and
   that if NCC is close to one, we have a good match.

   inputs :  grayPath0    : path to first grayscale image
             grayPath1    : path to second grayscale image
             ip0          : list of interest point in first image
             ip1          : list of interest point in 2nd image
             F            : fundamental matrix (used to remove outliers)
             Ncctol       : tolerance for Normalized Cross Correlation. A good value is
                            NCCtol = 0.9.
             neighSize    : size of patch around pixels. A good value of this parameter is 20.
                            we consider patch the patch [u-neighSize, u+neighSize]x[v-neighSize, v+neighSize]
                            to compute the NCC.
             dist         : distance in pixels between the match. We made the assumptions that
                            the matching pixels are in the same region in the two images.
                            Good value is 200.
             Ftol         : Let x = (u0,v0,1) and x'=(u1,v1,1) be a possible stereo match,
                            where (u0,v0) is a pixel in first image, and (u1,v1) a pixel in second image.
                            if |x'Fx / F[3][3] | < Ftol, we consider that we have matching pixels.
                            Good value of Ftol is 0.15

   output :   pixels      : list of matches. We have pixel[i] = ( u0, v0, u1, v1) where (u0,v0) is
                            a pixel in first image, and (u1,v1) is the corresponding one in
                            the second image.

------------------------------------------------------------------------------------------
*/

void harrisMatch(vector < vector <double> >&  pixels, const char *grayPath0, const char *grayPath1,
		vector <vector <int> >& ip0, vector <vector <int> >& ip1, double (&F)[9], double NCCtol,
		int neighSize, int dist, double Ftol) {

	// clean inputs
	pixels.clear();

	// load images
    CImg<unsigned char> img0;
    CImg<unsigned char> img1;
    CImg<unsigned char> match0;
    CImg<unsigned char> match1;

    std::string image0path;
    image0path=path(grayPath0)+".tiff";

    std::string image1path;
    image1path=path(grayPath1)+".tiff";

    img0=CImg<unsigned char>(image0path.c_str());
    img1=CImg<unsigned char>(image1path.c_str());

    int width=img1.width(); int height=img1.height();
    match0.assign(width,height,1,3,0);
    match1.assign(width,height,1,3,0);

	vector < vector <double> > imgInt0;
	vector < vector <double> > imgInt1;

	// store intensity for cpu time optimization
    for(int u(0); u<width; ++u){
		vector <double> tmpI0;
		vector <double> tmpI1;
    	for(int v(0); v < height; ++v){
    		tmpI0.push_back((double) img0(u,v,0,0));
    		tmpI1.push_back((double) img1(u,v,0,0));
    	}
    	imgInt0.push_back(tmpI0);
    	imgInt1.push_back(tmpI1);
    }

    vector < vector <double> > matches;

    // create a list of matche status
    vector < bool >  ip1match;

    for(int q(0); q < (int) ip1.size(); ++q)
    	ip1match.push_back(false);

	// create a list of pixel matches using NCC
	for(int p(0); p< (int) ip0.size(); ++p){
	    // compute mean and RMS of first patch
	    int u0 = ip0[p][0]; int v0 = ip0[p][1];
	    double  p0_mean=0.0;
	    double  p0_rms(0.0);

	    // compute average
	    for(int i(-neighSize); i<neighSize+1; ++i){
	    	for(int j(-neighSize); j<neighSize+1; ++j){
	    		if((u0+i>=0) && (u0+i < width) && (v0+j>=0) && (v0+j< height)){
					double r=imgInt0[u0+i][v0+j];

					p0_mean += r/(2.0*neighSize+1)/(2.0*neighSize+1) ;
	    		}
	    	}
	    }

	    // compute rms
	    for(int i(-neighSize); i<neighSize+1; ++i){
	    	for(int j(-neighSize); j<neighSize+1; ++j){
	    		if((u0+i>=0) && (u0+i < width) && (v0+j>=0) && (v0+j< height)){
					double r=imgInt0[u0+i][v0+j];

					p0_rms += (r-p0_mean)*(r-p0_mean)/((2.0*neighSize+1)*(2.0*neighSize+1));
	    		}
	    	}
	    }
	    double rms0 = sqrt(p0_rms);
    	double corr(-1.0);
    	int um(0), vm(0), qMatch=0;;

	    // check if root mean square is not too small
	    if( rms0 < 1e-5){
	    	cout << " Variance of neigbourgh hood too small, return default pixel (0,0) " << endl;
	    }
	    // rms is big enough to make a search in second image
	    else{
		for(int q(0); q< (int) ip1.size(); ++q){
			int u(ip1[q][0]), v(ip1[q][1]);
				// look for a matching pixel in second image using Normalized Cross Correlation (NCC)
				double p1_mean = 0.0;
				double p1_rms = 0.0;

				// compute average
				for(int i(-neighSize); i<neighSize+1; ++i){
					for(int j(-neighSize); j<neighSize+1; ++j){
						if( (u+i>=0) && (u+i < width) && (v+j>=0) && (v+j< height) ){
							double r=imgInt1[u+i][v+j];
							p1_mean += r/(2.0*neighSize+1)/(2.0*neighSize+1);
						};
					};
				}

				// compute rms
				for(int i(-neighSize); i<neighSize+1; ++i){
					for(int j(-neighSize); j<neighSize+1; ++j){
						if((u+i>=0) && (u+i < width) && (v+j>=0) && (v+j< height)){
							double r=imgInt1[u+i][v+j];
							p1_rms += (r-p1_mean)*(r-p1_mean)/((2.0*neighSize+1)*(2.0*neighSize+1));
						}
					}
				}

				double rms1 = sqrt(p1_rms);
				double crossCorr=0.0;

				// if rms1 is big enough, compute cross correlation associated to pixel (u,v)
				if(rms1 > 1e-5){
					// compute cross correlation
					for(int i(-neighSize); i<neighSize+1; ++i){
						for(int j(-neighSize); j<neighSize+1; ++j){
							if((u+i>=0) && (u+i < width) && (v+j>=0) && (v+j< height)){
								if((u0+i>=0) && (u0+i < width) && (v0+j>=0) && (v0+j< height)){

									double r0=imgInt0[u0+i][v0+j];
									double r1=imgInt1[u+i][v+j];
									crossCorr += (r0-p0_mean)*(r1-p1_mean)/((2.0*neighSize+1)*(2.0*neighSize+1));

								}
							}
						}
					} // end cross corr

					crossCorr /= rms1*rms0;

					// update output values
					if(crossCorr > corr && crossCorr > NCCtol){
						um = u; vm=v;
						corr = crossCorr;
						qMatch=q;
					}
				}
			}
		}

		 // export value
		if(corr > NCCtol){ // if we found a matching point...
			double x1Fx= 0.0;
			int p1[3] = { um, vm, 1};
			int p0[3] = { u0, v0, 1};

			for(int i(0); i<3; ++i){
				for(int j(0); j<3;++j)
					 x1Fx += F[3*i+j]*p0[j]*p1[i]/F[8];
			}

			// be sure that (um,vm) is not too far from (u0,v0) and it satisfies the fundamental matrix condition.
			if(sqrt((u0-um)*(u0-um) + (v0-vm)*(v0-vm) ) < dist && fabs(x1Fx) < Ftol && ip1match[qMatch]==false){
				  vector < double >  stereo_match;
				   stereo_match.push_back(u0);
				   stereo_match.push_back(v0);
				   stereo_match.push_back(um);
				   stereo_match.push_back(vm);
				   matches.push_back(stereo_match);

				   ip1match[qMatch]=true;
			}
	    
	      }
	}

	// eliminate multiple matches (i.e. matches that are the same at most two pixels diff per compounds)
	mergeSort(matches, (int) matches.size());

   	for(int n(0); n < (int) matches.size(); ++n){

    		if(n ==0 ){
    			pixels.push_back(matches[n]);
    		}
    		else{
			double crit = fabs(matches[n][0]-matches[n-1][0])
						  + fabs(matches[n][1]-matches[n-1][1])
						  + fabs(matches[n][2]-matches[n-1][2])
						  + fabs(matches[n][3]-matches[n-1][3]);
			if( crit > 8 )
				pixels.push_back(matches[n]);
    	}
    }

#if DEBUG
    cout << " orginal number of matches " << matches.size() << endl;
    cout << " Eliminated " << (int) matches.size() - (int) pixels.size() << " correspondences " << endl;
	cout << " found " << pixels.size() << "  correspondances " << endl;

	// export matches on original image
	for(int u(0); u<width; ++u){
		for(int v(0); v<height; ++v){
			match0(u,v,0,0) = img0(u,v,0,0);
			match0(u,v,0,1) = img0(u,v,0,0);
			match0(u,v,0,2) = img0(u,v,0,0);

			match1(u,v,0,0) = img1(u,v,0,0);
			match1(u,v,0,1) = img1(u,v,0,0);
			match1(u,v,0,2) = img1(u,v,0,0);
		}
	}

	// visualize matches
    int size = (int) pixels.size();
	for(int n(0); n < size; ++n){
		int u0=pixels[n][0]; int v0 =pixels[n][1];
		int u1=pixels[n][2]; int v1 =pixels[n][3];

		int r,g,b;

		r=0; b=0; g=0;

		int half = size/ 2.0;

		if( n < half)
		   r = floor(255*(1.0 - 2.0*n/( (double) size ) ) );

		if( n < half )
		    g = floor(255*(2.0*n/( (double) size ) ) );
		else
		    g = floor(255*(1.0-2.0*(n-half)/( (double) size ) ) );

		if( n >= half )
		    b = floor(255*(2.0*(n-half)/( (double) size ) ) );

		for(int i(-5); i<6; ++i){
			for(int j(-5); j<6; ++j){
				if((u0+i >= 0) && (u0+i < width) && (v0+j >=0 ) && (v0+j < height) ){
					match0(u0+i,v0+j,0,0)=r;
					match0(u0+i,v0+j,0,1)=g;
					match0(u0+i,v0+j,0,2)=b;
				}

				if((u1+i >= 0) && (u1+i < width) && (v1+j >=0 ) && (v1+j < height) ){
					match1(u1+i,v1+j,0,0)=r;
					match1(u1+i,v1+j,0,1)=g;
					match1(u1+i,v1+j,0,2)=b;
				}
			}
		}
	}

	match0.save_tiff("./match0.tiff");
	match1.save_tiff("./match1.tiff");
#endif
}

// convert a RGB image into a grayscale one
void rgb2gray(const char *RGBPath, const char *grayPath){
    //load image and create output
    CImg<unsigned char> rgb;
    CImg<unsigned char> gray;

    std::string image0path;
    image0path=path(RGBPath)+".tiff";

    std::string image1path;
    image1path=path(grayPath)+".tiff";

    rgb=CImg<unsigned char>(image0path.c_str());

    // create output bitmap
    gray.assign(rgb.width(),rgb.height(),1,1,0);

    // convert image (blue is omitted, because it blurs image)
    for(int u(0); u < rgb.width(); ++u)
    	for(int v(0); v < rgb.height(); ++v){
    		gray(u,v,0,0) = round( 0.229*((double) rgb(u,v,0,0)) +0.771*( (double) rgb(u,v,0,1)) );
    	}

    // save image
    gray.save_tiff(image1path.c_str());
};

