// Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __hf_filter_h__
#define __hf_filter_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "image.h"

#define FILTER_SCALE 2.5

/* 
 * High frequency filter function
 */
static inline int hf_filter(double factor, int i, int j, int k, image *im) {

	double filter_radius; 
	double adj_factor;
	double filter_dimension[4];
	int fr_ceil[4];
	double fr_frac[4];
	double result;
	int ii, jj;
	int d;

	filter_radius = ((FILTER_SCALE * factor) - 1) / 2;

	filter_dimension[0] = (i - filter_radius < 0) ? i : filter_radius;
	filter_dimension[1] = (j - filter_radius < 0) ? j : filter_radius;
	filter_dimension[2] = (i + filter_radius > height(im) - 1)
			    ? (height(im) - i - 1) : filter_radius;
	filter_dimension[3] = (j + filter_radius > width(im) - 1)
			    ? (width(im) - j - 1) : filter_radius;

	for (d = 0; d < 4; d++) {
		fr_ceil[d] = ceil(filter_dimension[d]);
		fr_frac[d] = 1 - (fr_ceil[d] - filter_dimension[d]);
	}

	adj_factor = (filter_dimension[0] + filter_dimension[2] + 1)
		   * (filter_dimension[1] + filter_dimension[3] + 1);

	result = adj_factor * get_pixel_component(im, i, j, k);

	result -= fr_frac[0] * fr_frac[1] 
		* get_pixel_component(im, i - fr_ceil[0], j - fr_ceil[1], k);
	result -= fr_frac[2] * fr_frac[1] 
		* get_pixel_component(im, i + fr_ceil[2], j - fr_ceil[1], k);
	result -= fr_frac[2] * fr_frac[3] 
		* get_pixel_component(im, i + fr_ceil[2], j + fr_ceil[3], k);
	result -= fr_frac[0] * fr_frac[3] 
		* get_pixel_component(im, i - fr_ceil[0], j + fr_ceil[3], k);

	for (jj = j - fr_ceil[1] + 1; jj <= j + fr_ceil[3] - 1; jj++) {
		result -= fr_frac[0] * get_pixel_component(im, i - fr_ceil[0], jj, k);
		result -= fr_frac[2] * get_pixel_component(im, i + fr_ceil[2], jj, k);
	}

	for (ii = i - fr_ceil[0] + 1; ii <= i + fr_ceil[2] - 1; ii++) {
		result -= fr_frac[1] * get_pixel_component(im, ii, j - fr_ceil[1], k);
		result -= fr_frac[3] * get_pixel_component(im, ii, j + fr_ceil[3], k);
	}

	for (ii = i - fr_ceil[0] + 1; ii <= i + fr_ceil[2] - 1; ii++)
		for (jj = j - fr_ceil[1] + 1; jj <= j + fr_ceil[3] - 1; jj++)
			result -= get_pixel_component(im, ii, jj, k);

	return result / adj_factor;
}

#endif
