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

#ifndef __image_h__
#define __image_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

template <class colorT>
class image_base {
private:
	unsigned int _dimx, _dimy, _depth;
	colorT *_p;
	int _apm_memo;
	double _apm;
public:
	image_base (unsigned int dimy, unsigned int dimx, unsigned int depth) {
		_dimx = dimx;
		_dimy = dimy;
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

	unsigned int width() {
		return _dimx;
	}

	unsigned int height() {
		return _dimy;
	}

	unsigned int depth() {
		return _depth;
	}

	colorT get_pixel_component(unsigned int y, unsigned int x, 
			unsigned int color) {

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

	/*
	 * Get a color value at a given position using bilinear interpolation between the
	 * four nearest pixels.
	 */
	colorT get_bl_component(double y, double x, unsigned int color) {

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
	 * Scale by some factor
	 */
	void scale(double f) {

		image_base<colorT> *is = new image_base<colorT> (
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth());

		assert(is);

		unsigned int i, j, k;

		for (i = 0; i < is->height(); i++)
		for (j = 0; j < is->width(); j++) 
		for (k = 0; k < is->depth(); k++)
			is->set_pixel_component(i, j, k,
				get_bl_component(
					(i/f <= height() - 1) 
						? (i/f) 
					        : (height() - 1), 
					(j/f <= width() - 1)
						? (j/f) 
						: (width() - 1), k));

		free(_p);

		_p = is->_p;
		_dimx = is->_dimx;
		_dimy = is->_dimy;
		_depth = is->_depth;
		_apm_memo = 0;

		is->_p = NULL;

		delete is;
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

		is->_p = NULL;

		delete is;
	}


	/*
	 * Clone 
	 */
	image_base<colorT> *clone() {
		image_base<colorT> *ic = new image_base<colorT>(
			height(), width(), depth());

		assert(ic);

		unsigned int i, j, k;

		for (i = 0; i < height(); i++)
		for (j = 0; j < width();  j++)
		for (k = 0; k < depth();  k++)
			ic->set_pixel_component(i, j, k,
				get_pixel_component(i, j, k));

		ic->_apm_memo = _apm_memo;

		return ic;
	}


	/*
	 * Calculate the average (mean) magnitude of a pixel (where magnitude
	 * is defined as the mean of the channel values).
	 */
	double avg_pixel_magnitude() {
		unsigned int i, j, k;

		if (_apm_memo)
			return _apm;

		_apm_memo = 1;
		_apm = 0;

		for (i = 0; i < _dimy;  i++)
		for (j = 0; j < _dimx;  j++)
		for (k = 0; k < _depth; k++)
			_apm += get_pixel_component(i, j, k);

		_apm /= (_dimy * _dimx * _depth);

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
