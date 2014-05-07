/*
 * corner_detect.hpp
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
 *      Author: s.flotron@foxel.ch
 */

#ifndef CORNER_DETECT_HPP_
#define CORNER_DETECT_HPP_

#include "common.hpp"
#include "triangulation.hpp"

void rgb2gray(const char *RBGPath, const char *grayPath); // convert rgb image into grayscale
// merge sort
void mergeSort(vector < vector <double > >& p, const int& n);
vector< vector <double> > merge(vector < vector <double > >& p, int P, vector < vector <double > >& q, int Q);

// find interest point using harris corner detector
void harrisANMS(const char *grayPath, const char *corners, vector <vector <int> >& points,
		double k, double tol, int nip);

// find matches using harris interest points and fundamental matrix
void harrisMatch(vector < vector <double> >&  pixels, const char *grayPath0, const char *grayPath1,
		vector <vector <int> >& ip0, vector <vector <int> >& ip1, double (&F)[9], double NCCtol,
		int neighSize, int dist, double Ftol);

#endif /* CORNER_DETECT_HPP_ */
