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

/* 
 * High frequency filter function
 */
static inline int hf_filter(double factor, int i, int j, int k, image *im) {

	double filter_radius = (factor - 1) / 2;
	double fr_ceil = ceil(filter_radius);
	double fr_frac = 1 - (fr_ceil - filter_radius);
	double result;
	int ii, jj;

	result = get_pixel_component(im, i, j, k) * factor * factor;

	if (i - fr_ceil < 0
         || j - fr_ceil < 0
	 || i + fr_ceil >= height(im)
	 || j + fr_ceil >= width(im))
		return 0;

	result -= fr_frac * fr_frac * get_pixel_component(im, i - fr_ceil, j - fr_ceil, k);
	result -= fr_frac * fr_frac * get_pixel_component(im, i + fr_ceil, j - fr_ceil, k);
	result -= fr_frac * fr_frac * get_pixel_component(im, i + fr_ceil, j + fr_ceil, k);
	result -= fr_frac * fr_frac * get_pixel_component(im, i - fr_ceil, j + fr_ceil, k);

	for (jj = j - fr_ceil + 1; jj <= j + fr_ceil - 1; jj++) {
		result -= fr_frac * get_pixel_component(im, i - fr_ceil, jj, k);
		result -= fr_frac * get_pixel_component(im, i + fr_ceil, jj, k);
	}

	for (ii = i - fr_ceil + 1; ii <= i + fr_ceil - 1; ii++) {
		result -= fr_frac * get_pixel_component(im, ii, j - fr_ceil, k);
		result -= fr_frac * get_pixel_component(im, ii, j + fr_ceil, k);
	}

	for (ii = i - fr_ceil + 1; ii <= i + fr_ceil - 1; ii++)
		for (jj = j - fr_ceil + 1; jj <= j + fr_ceil - 1; jj++)
			result -= get_pixel_component(im, ii, jj, k);

// 	if (result > 255)
// 		return 255;
// 	if (result < -255)
// 		return -255;

	return result;
}

#endif
