// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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

/*
 * zero_image.h: Image that is zero everywhere.
 */

#ifndef __image_zero_h__
#define __image_zero_h__

#include "point.h"
#include "pixel.h"
#include "exposure/exposure.h"

class image_zero : public image_weighted_avg {
public:

	pixel get_pixel(unsigned int y, unsigned int x) const {
		return pixel::zero();
	}

	ale_real maxval() const {
		return 0;
	}

	ale_real minval() const {
		return 0;
	}

	/*
	 * Get a color value at a given position using bilinear interpolation between the
	 * four nearest pixels.  Result values:
	 *
	 * result[0] == pixel value
	 * result[1] == pixel confidence
	 */
	void get_bl(point x, pixel result[2]) const {
		result[0] = pixel::zero();
		result[1] = pixel::zero();
	}

	pixel get_bl(point x) const {
		pixel result[2];

		get_bl(x, result);

		return result[0];
	}

	pixel get_scaled_bl(point x, ale_pos f) const {
		return pixel::zero();
	}

		
	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, char *name) const {

		image *is = new image_zero(height, width, depth, name);

		assert(is);

		return is;
	}

	/*
	 * Return an image scaled by some factor >= 1.0
	 */
	image *scale(ale_pos f, char *name) const {

		image *is = new image_zero(
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(is);

		return is;
	}

	/*
	 * Scale by half.  We use the following filter:
	 *
	 * 	1/16	1/8	1/16
	 * 	1/8	1/4	1/8
	 * 	1/16	1/8	1/16
	 *
	 * At the edges, these values are normalized so that the sum of
	 * contributing pixels is 1.
	 */
	image *scale_by_half(char *name) const {
		ale_pos f = 0.5;

		image *result = new image_zero(
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(result);

		return result;
	}

	/*
	 * Scale an image definition array by 1/2.
	 *
	 * ALE considers an image definition array as a special kind of image
	 * weight array (typedefs of which should appear below the definition
	 * of this class).  ALE uses nonzero pixel values to mean 'defined' and
	 * zero values to mean 'undefined'.  Through this interpretation, the
	 * image weight array implementation that ALE uses allows image weight
	 * arrays to also serve as image definition arrays.
	 *
	 * Whereas scaling of image weight arrays is not generally obvious in
	 * either purpose or method, ALE requires that image definition arrays
	 * be scalable, and the method we implement here is a fairly obvious
	 * one.  In particular, if any source pixel contributing to the value of 
	 * a scaled target pixel has an undefined value, then the scaled target
	 * pixel is undefined (zero).  Otherwise, it is defined (non-zero).
	 * 
	 * Since there are many possible ways of implementing this function, we
	 * choose an easy way and simply multiply the numerical values of the
	 * source pixels to obtain the value of the scaled target pixel.
	 *
	 * XXX: we consider all pixels within a 9-pixel square to contribute.
	 * There are other approaches.  For example, edge pixels could be
	 * considered to have six contributing pixels and corner pixels four
	 * contributing pixels.  To use this convention, change the ': 0' text
	 * in the code below to ': 1'.
	 */

	image *defined_scale_by_half(char *name) const {
		ale_pos f = 0.5;

		image *result = new image_zero(
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(result);

		return result;
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.
	 */
	virtual void extend(int top, int bottom, int left, int right) {
		_dimy += top + bottom;
		_dimx += left + right;
		_offset[0] -= top;
		_offset[1] -= left;
	}

	/*
	 * Clone 
	 */
	image *clone(char *name) const {
		return new image_zero(_dimy, _dimx, _depth, name);
	}

	/*
	 * Calculate the average (mean) clamped magnitude of a channel across
	 * all pixels in an image.  The magnitude is clamped to the range of 
	 * real inputs.
	 */
	ale_real avg_channel_clamped_magnitude(unsigned int k) const {
		return 0;
	}

	pixel avg_channel_clamped_magnitude() const {
		return pixel::zero();
	}

	/*
	 * Calculate the average (mean) magnitude of a channel across all
	 * pixels in an image.
	 */
	ale_real avg_channel_magnitude(unsigned int k) const {
		return 0;
	}

	pixel avg_channel_magnitude() const {
		return pixel::zero();
	}

	/*
	 * Calculate the average (mean) magnitude of a pixel (where magnitude
	 * is defined as the mean of the channel values).
	 */
	ale_real avg_pixel_magnitude() const {
		return 0;
	}

	image_zero(unsigned int dimy, unsigned int dimx, unsigned int depth,
			char *name = "anonymous") : image_weighted_avg(dimy, dimx, depth, name) {
	}

	int accumulate_norender(int i, int j) {
		return 1;
	}

	void accumulate(int i, int j, int f, pixel new_value, pixel new_weight) {
		assert(0);
	}

	image *get_colors() {
		return this;
	}

	image *get_weights() {
		return this;
	}
};

#endif
