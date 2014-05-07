/*
 *  triangulation.hpp
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
 *  Created on: Apr 1, 2014
 *      Author: s.flotron@foxel.ch
 */

#ifndef TRIANGULATION_HPP_
#define TRIANGULATION_HPP_

#include "common.hpp"

using namespace std;
using namespace cv;
using namespace cimg_library;

// extern library functions (in lapack)
extern "C" {
     void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a,
             int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
             double* work, int* lwork, int* info );
}

// functions used in triangulation

void print_matrix( char* desc, int m, int n, double* a, int lda );
void callSVD(int* m, int* n, double* a, double* s, double* u, double* vt, char jobu, char jobv);
// triangulation using DLT method
void linearDLTTriangulation(double (&X)[4], double (&p0)[3], double (&p1)[3], double (&P0)[3][4], double (&P1)[3][4]);
// given a pixel in an image, find the matching pixel in the second image minimizing reprojection error.
void findMatch(double (&pm)[2], double u, double v, int pixWidth, double (&P0)[3][4], double (&P1)[3][4],
		double (&K0)[9], double (&K1)[9], double (&F)[9], double (&R0)[3][3], double (&R1)[3][3],
		const char *imageZeroPath, const char *imageOnePath);
// given a pixel in first image, find a matching pixel in second image using statistical comparison
void findMatchStat(double (&p0)[3], double (&p1)[3], const char *imageZeroPath, const char *imageOnePath, double& corr);
// convert unsigned char to double
double charToDouble(unsigned char* c);
// extract basename of image
std::string path(const std::string &filename);

#endif /* TRIANGULATION_HPP_ */
