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

/*
 * image.h: The internal representation of images used by ALE.
 */

#ifndef __image_h__
#define __image_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "point.h"

template <class colorT>
class image_base {
private:
	unsigned int _dimx, _dimy, _depth;
	point _offset;
	colorT *_p;
	int _apm_memo;
	double _apm;
public:
	image_base (unsigned int dimy, unsigned int dimx, unsigned int depth) {
		_dimx = dimx;
		_dimy = dimy;
		_offset = point(0, 0);
		_depth = depth;
		_apm_memo = 0;

		_p = (colorT *) calloc(dimx * dimy * depth, sizeof(colorT));

		assert (_p);

		if (!_p) {
			fprintf(stderr, "Could not allocate memory for image data.\n");
			exit(1);
		}
	}

	~image_base() {
		free(_p);
	}

	point offset() const {
		return _offset;
	}

	unsigned int width() const {
		return _dimx;
	}

	unsigned int height() const {
		return _dimy;
	}

	unsigned int depth() const {
		return _depth;
	}

	colorT get_pixel_component(unsigned int y, unsigned int x, 
			unsigned int color) const {

		assert (x < _dimx);
		assert (y < _dimy);
		assert (color < _depth);

		return _p[y * _dimx * _depth + x * _depth + color];
	}

	void set_pixel_component(unsigned int y, unsigned int x, 
			unsigned int color, colorT value) {

		assert (x < _dimx);
		assert (y < _dimy);
		assert (color < _depth);

		/* XXX: We could do much better than this. */
		_apm_memo = 0;

		_p[y * _dimx * _depth + x * _depth + color] = value;
	}

#if 0

	/* Tested for ALE 0.4.5 -- offers an increase in speed but breaks encapsulation. */

	unsigned char *get_pixel_array() {
		_apm_memo = 0;
		return _p;
	}
#endif

	/*
	 * Get a color value at a given position using bilinear interpolation between the
	 * four nearest pixels.
	 */
	colorT get_bl_component(double y, double x, unsigned int color) const {

		assert (y >= 0);
		assert (x >= 0);
		assert (y <= _dimy - 1);
		assert (x <= _dimx - 1);
		assert (color < _depth);

		int lx = (int) floor(x);
		int hx = (int) floor(x) + 1;
		int ly = (int) floor(y);
		int hy = (int) floor(y) + 1;

		return 
		(colorT) ((hx - x) 
			* (hy - y)
			* get_pixel_component(ly, lx, color)

			+ (hx - x) 
			* (y - ly)
			* get_pixel_component(hy % _dimy, lx, color)

			+ (x - lx) 
			* (y - ly)
			* get_pixel_component(hy % _dimy, hx % _dimx, color)

			+ (x - lx) 
			* (hy - y)
			* get_pixel_component(ly, hx % _dimx, color));

	}

	/*
	 * Get a color value at a given position using bilinear interpolation between the
	 * four nearest pixels.
	 */
	colorT get_scaled_bl_component(double y, double x, unsigned int color, double f) const {
		return get_bl_component(
			(y/f <= height() - 1) 
				? (y/f) 
				: (height() - 1), 
			(x/f <= width() - 1)
				? (x/f) 
				: (width() - 1), color);
	}

	/*
	 * Scale by some factor >= 1.0.
	 */
	void scale(double f) {

		if (f == 1.0)
			return;

		/*
		 * Assertion:
		 *
		 * This is not a good way to scale down images.  
		 */
		assert (f > 1.0);

		image_base<colorT> *is = new image_base<colorT> (
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(is);

		unsigned int i, j, k;

		for (i = 0; i < is->height(); i++)
		for (j = 0; j < is->width(); j++) 
		for (k = 0; k < is->depth(); k++)
			is->set_pixel_component(i, j, k,
				get_scaled_bl_component(i, j, k, f));

		free(_p);

		_p = is->_p;
		_dimx = is->_dimx;
		_dimy = is->_dimy;
		_depth = is->_depth;
		_apm_memo = 0;
		_offset = point(_offset[0] * f, _offset[1] * f);

		is->_p = NULL;

		delete is;
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
	image_base<colorT> *scale_by_half() const {
		double f = 0.5;

		image_base<colorT> *is = new image_base<colorT> (
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(is);

		unsigned int i, j, k;

		for (i = 0; i < is->height(); i++)
		for (j = 0; j < is->width(); j++) 
		for (k = 0; k < is->depth(); k++)

			is->set_pixel_component(i, j, k, (int) 

			      ( ( (i > 0 && j > 0) 
				    ? get_pixel_component(2 * i - 1, 2 * j - 1, k) * 0.0625
				    : 0
				+ (i > 0)
				    ? get_pixel_component(2 * i - 1, 2 * j, k) * 0.125
				    : 0
				+ (i > 0 && j < is->width() - 1)
				    ? get_pixel_component(2 * i - 1, 2 * j + 1, k) * 0.0625
				    : 0
				+ (j > 0)
				    ? get_pixel_component(2 * i, 2 * j - 1, k) * 0.125
				    : 0
				+ get_pixel_component(2 * i, 2 * j, k) * 0.25
				+ (j < is->width() - 1)
				    ? get_pixel_component(2 * i, 2 * j + 1, k) * 0.125
				    : 0
				+ (i < is->height() - 1 && j > 0)
				    ? get_pixel_component(2 * i + 1, 2 * j - 1, k) * 0.0625
				    : 0
				+ (i < is->height() - 1)
				    ? get_pixel_component(2 * i + 1, 2 * j, k) * 0.125
				    : 0
				+ (i < is->height() && j < is->width() - 1)
				    ? get_pixel_component(2 * i + 1, 2 * j + 1, k) * 0.0625 
				    : 0 ) 

			     /

				( (i > 0 && j > 0) 
				    ? 0.0625
				    : 0
				+ (i > 0)
				    ? 0.125
				    : 0
				+ (i > 0 && j < is->width() - 1)
				    ? 0.0625
				    : 0
				+ (j > 0)
				    ? 0.125
				    : 0
				+ 0.25
				+ (j < is->width() - 1)
				    ? 0.125
				    : 0
				+ (i < is->height() - 1 && j > 0)
				    ? 0.0625
				    : 0
				+ (i < is->height() - 1)
				    ? 0.125
				    : 0
				+ (i < is->height() && j < is->width() - 1)
				    ? 0.0625
				    : 0 ) ) );

		is->_offset = point(_offset[0] * f, _offset[1] * f);

		return is;
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

	image_base<colorT> *defined_scale_by_half() const {
		double f = 0.5;

		image_base<colorT> *is = new image_base<colorT> (
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(is);

		unsigned int i, j, k;

		for (i = 0; i < is->height(); i++)
		for (j = 0; j < is->width(); j++) 
		for (k = 0; k < is->depth(); k++)

			is->set_pixel_component(i, j, k,

				( (i > 0 && j > 0) 
				    ? get_pixel_component(2 * i - 1, 2 * j - 1, k)
				    : 0
				* (i > 0)
				    ? get_pixel_component(2 * i - 1, 2 * j, k)
				    : 0
				* (i > 0 && j < is->width() - 1)
				    ? get_pixel_component(2 * i - 1, 2 * j + 1, k)
				    : 0
				* (j > 0)
				    ? get_pixel_component(2 * i, 2 * j - 1, k)
				    : 0
				* get_pixel_component(2 * i, 2 * j, k)
				* (j < is->width() - 1)
				    ? get_pixel_component(2 * i, 2 * j + 1, k)
				    : 0
				* (i < is->height() - 1 && j > 0)
				    ? get_pixel_component(2 * i + 1, 2 * j - 1, k)
				    : 0
				* (i < is->height() - 1)
				    ? get_pixel_component(2 * i + 1, 2 * j, k)
				    : 0
				* (i < is->height() && j < is->width() - 1)
				    ? get_pixel_component(2 * i + 1, 2 * j + 1, k)
				    : 0 ) );

		is->_offset = point(_offset[0] * f, _offset[1] * f);

		return is;
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.
	 */
	void extend(unsigned int top, unsigned int bottom, 
			unsigned int left, unsigned int right) {

		image_base<colorT> *is = new image_base<colorT> (
			height() + top  + bottom,
			 width() + left + right , depth());

		assert(is);

		unsigned int i, j, k;

		for (i = 0; i < height(); i++)
		for (j = 0; j <  width(); j++) 
		for (k = 0; k <  depth(); k++)
			is->set_pixel_component(i + top, j + left, k,
				get_pixel_component(i, j, k));

		free(_p);

		_p = is->_p;
		_dimx = is->_dimx;
		_dimy = is->_dimy;
		_depth = is->_depth;
		_apm_memo = 0;

		_offset[0] -= top;
		_offset[1] -= left;

		is->_p = NULL;

		delete is;
	}


	/*
	 * Clone 
	 */
	image_base<colorT> *clone() const {
		image_base<colorT> *ic = new image_base<colorT>(
			height(), width(), depth());

		assert(ic);

#if 0
		/*
		 * Tested for 0.5.0.  No obvious increase in speed.
		 */

		memcpy(ic->_p, _p, height() * width() * depth() * sizeof(colorT));
#else
		unsigned int i, j, k;

		for (i = 0; i < height(); i++)
		for (j = 0; j < width();  j++)
		for (k = 0; k < depth();  k++)
			ic->set_pixel_component(i, j, k,
				get_pixel_component(i, j, k));
#endif

		ic->_apm_memo = _apm_memo;
		ic->_apm = _apm;
		ic->_offset = _offset;

		return ic;
	}


	/*
	 * Calculate the average (mean) magnitude of a pixel (where magnitude
	 * is defined as the mean of the channel values).
	 */
	double avg_pixel_magnitude() const {
		unsigned int i, j, k;

		image_base<colorT> *rw_ptr;
		
		rw_ptr = (image_base<colorT> *) this;

		if (_apm_memo)
			return _apm;

		rw_ptr->_apm_memo = 1;
		rw_ptr->_apm = 0;

		for (i = 0; i < _dimy;  i++)
		for (j = 0; j < _dimx;  j++)
		for (k = 0; k < _depth; k++)
			rw_ptr->_apm += get_pixel_component(i, j, k);

		rw_ptr->_apm /= (_dimy * _dimx * _depth);

		return _apm;
	}

};

typedef image_base<unsigned char> image;

static inline image *new_image(unsigned int dimy, unsigned int dimx, 
		unsigned int depth) {
	image *result = new image(dimy, dimx, depth);
	
	assert(result);

	if (result == NULL) {
		fprintf(stderr, "Unable to allocate memory for image.\n");
		exit(1);
	}

	return result;
}

typedef image_base<double> image_weights;

static inline image_weights *new_image_weights(unsigned int dimy, unsigned int dimx, 
		unsigned int depth) {
	image_weights *result = new image_weights(dimy, dimx, depth);
	
	assert(result);

	if (result == NULL) {
		fprintf(stderr, "Unable to allocate memory for image.\n");
		exit(1);
	}

	return result;
}

#endif
