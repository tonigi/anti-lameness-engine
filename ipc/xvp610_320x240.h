// Copyright 2003 David Hilvert <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __xvp610_320x240_h__
#define __xvp610_320x240_h__

#include <assert.h>
#include "../point.h"

/*
 * This is a first hack at a response function for an IBM XVP610 webcam.  
 *
 * The camera driver used was the ibmcam.o module for Linux 2.5.69, and images
 * were captured by Xawtv with the following settings: hue=20508/65535,
 * color=31653/65535, contrast=44582/65535.  sharpness was set to 6/6 at driver
 * initialization.  
 *
 * The script ale-calibrate was used for calibration.  The calibration measured
 * output after 3 iterations.  Only one device was used in calibration.  A
 * single set of images was used in calibration (36 images, each about 14
 * pixels by 22 pixels), and the output of processing this set was compared
 * against a cropped and transformed close-up of the same scene (captured by
 * the same device) for the purposes of calibration.
 *
 * The results of calibration were multiplied by 10 to obtain values in the
 * domain of the C signed character type.
 */

#if 1

/*
 * This filter may offer better quality than the filter below, but may also
 * be more computationally expensive.  Calibrated at 3 iterations.
 */

#define FILTER_DIM 5
#define RADIUS 1.5

static const signed char calibrated_response[FILTER_DIM][FILTER_DIM][3] = {
	{
		{ 15,  0,  2 },
		{  4,  0,  1 },
		{ -5,  0,  5 },
		{  8,  0,  4 },
		{ -4,  2,  1 },
	}, {
		{  4, 13,  0 },
		{  4,  1, 11 },
		{ 12, -6, 30 },
		{ -4, -3,  0 },
		{ 21, 26,  0 },
	}, {
		{  0,  0,  0 },
		{ 13, 30,  0 },
		{ 47, 58, 28 },
		{ 19, 31, 15 },
		{  0,  0,  0 },
	}, {
		{  0,  0,  0 },
		{  1,  0,  0 },
		{ 36, 56,  1 },
		{ 13,  4,  0 },
		{  0,  0,  0 },
	}, {
		{  1,  1,  1 },
		{  1,  1,  1 },
		{  1,  1,  1 },
		{  1,  1,  0 },
		{  1,  0,  0 },
	}
};

#else

/*
 * This filter may be less computationally expensive and offer lower quality than
 * the filter above.  Calibrated at 2 iterations.
 */

#define FILTER_DIM 3
#define RADIUS 1.1

static const signed char calibrated_response[FILTER_DIM][FILTER_DIM][3] = {
	{
		{  6, -2,  6 },
		{  1,  2,124 },
		{ -2, -6, 10 },
	}, {
		{ -2, 10,  0 },
		{ 40, 24, 11 },
		{ 18, 18, 12 },
	}, {
		{  0,  1,  0 },
		{ 37, 30,  1 },
		{ -5,  1,  0 },
	}
};

#endif

class xvp610_320x240 {
public:
	/*
	 * Result type.
	 */
	class ipc_result {
		friend class xvp610_320x240;
		double response[3];
	public:
		/*
		 * Response intensity on channel K2 resulting from
		 * stimulus on channel K1.  Channel R=0, G=1, B=2.
		 */
		double operator()(int k1, int k2) {
			if (k1 == k2)
				return response[k1];
			else
				return 0;
		}
	};

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	float min_i() { return -RADIUS; }
	float max_i() { return  RADIUS; }
	float min_j() { return -RADIUS; }
	float max_j() { return  RADIUS; }

	/*
	 * Get the response to the quadrilateral (p[0], p[1], p[2], p[3]).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  The position (i, j) of the responding pixel is
	 * provided, in case response is not uniform for all pixels (e.g. some
	 * sensor arrays stagger red, green, and blue sensors).
	 */
	struct ipc_result operator()(const point p[4], int i, int j) {
		struct ipc_result result;

		result.response[0] = 0;
		result.response[1] = 0;
		result.response[2] = 0;

		/* 
		 * XXX: We make some unsafe assumptions here.
		 */

		float top = p[0][0];
		float bot = p[2][0];
		float lef = p[0][1];
		float rig = p[1][1];

		if (top < min_i())
			top = min_i();
		if (bot > max_i())
			bot = max_i();
		if (lef < min_j())
			lef = min_j();
		if (rig > max_j())
			rig = max_j();

		int il = (int) floor((top - min_i()) / (max_i() - min_i()) * FILTER_DIM);
		int ih = (int) floor((bot - min_i()) / (max_i() - min_i()) * (FILTER_DIM - 0.001));
		int jl = (int) floor((lef - min_j()) / (max_j() - min_j()) * FILTER_DIM);
		int jh = (int) floor((rig - min_j()) / (max_j() - min_j()) * (FILTER_DIM - 0.001));

		for (int ii = il; ii <= ih; ii++)
		for (int jj = jl; jj <= jh; jj++) {

			float ltop = ((float) ii) / FILTER_DIM * (max_i() - min_i()) + min_i();
			float lbot = ((float) ii + 1) / FILTER_DIM * (max_i() - min_i()) + min_i();
			float llef = ((float) jj) / FILTER_DIM * (max_j() - min_j()) + min_j();
			float lrig = ((float) jj + 1) / FILTER_DIM * (max_j() - min_j()) + min_j();

			if (ltop < top)
				ltop = top;
			if (lbot > bot)
				lbot = bot;
			if (llef < lef)
				llef = lef;
			if (lrig > rig)
				lrig = rig;

			assert (ii >= 0);
			assert (ii < FILTER_DIM);
			assert (jj >= 0);
			assert (jj < FILTER_DIM);

			for (int k = 0; k < 3; k++)
				result.response[k] += ((lbot - ltop) * (lrig - llef)
					      * calibrated_response[ii][jj][k]);
		}

		return result;
	}
};

#endif
