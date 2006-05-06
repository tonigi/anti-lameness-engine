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
 * image.h: Abstract base class for the internal representations of images used
 * by ALE.
 */

#ifndef __image_h__
#define __image_h__

#include "point.h"
#include "pixel.h"
#include "exposure/exposure.h"

#define IMAGE_BAYER_NONE 0

/*
 * This constant indicates that some other default value should be filled in.
 */

#define IMAGE_BAYER_DEFAULT 0x8

/*
 * Do not change these values without inspecting
 * image_bayer_ale_real::r_*_offset().
 */
#define IMAGE_BAYER_RGBG 0x4 /* binary 100 */
#define IMAGE_BAYER_GBGR 0x5 /* binary 101 */
#define IMAGE_BAYER_GRGB 0x6 /* binary 110 */
#define IMAGE_BAYER_BGRG 0x7 /* binary 111 */

class image : protected exposure::listener {
protected:
	unsigned int _dimx, _dimy, _depth;
	point _offset;
	char *name;
	mutable exposure *_exp;
	unsigned int bayer;
private:
	/*
	 * Memoized function variables.  We may want to change these even when
	 * *this is constant.
	 */
	mutable int _apm_memo;
	mutable ale_real _apm;
	mutable int _accm_memo;
	mutable pixel _accm;
	mutable int _acm_memo;
	mutable pixel _acm;

	void avg_channel_clamped_magnitude_memo() const {
		unsigned int i, j, k;
		pixel_accum accumulator;
		pixel_accum divisor;

		if (_accm_memo)
			return;

		_accm_memo = 1;

		accumulator = pixel_accum(0, 0, 0);

		for (i = 0; i < _dimy;  i++)
		for (j = 0; j < _dimx;  j++) {
			pixel value = get_pixel(i, j);

			for (k = 0; k < _depth; k++)
				if (finite(value[k])) {
					if (value[k] > 1)
						value[k] = 1;
					if (value[k] < 0)
						value[k] = 0;
					accumulator[k] += value[k];
					divisor[k] += 1;
				}
		}

		accumulator /= divisor;

		_accm = accumulator;
	}

	void avg_channel_magnitude_memo() const {
		unsigned int i, j, k;
		pixel_accum accumulator;
		pixel_accum divisor;

		if (_acm_memo)
			return;

		_acm_memo = 1;

		accumulator = pixel_accum(0, 0, 0);

		for (i = 0; i < _dimy;  i++)
		for (j = 0; j < _dimx;  j++) {
			pixel value = get_pixel(i, j);

			for (k = 0; k < _depth; k++) 
				if (finite(value[k])) {
					accumulator[k] += value[k];
					divisor[k] += 1;
				}
		}

		accumulator /= divisor;

		_acm = accumulator;
	}

protected:
	void image_updated() {
		_apm_memo = 0;
		_acm_memo = 0;
		_accm_memo = 0;
	}

public:
	image (unsigned int dimy, unsigned int dimx, unsigned int depth, 
			char *name = "anonymous", exposure *_exp = NULL,
			unsigned int bayer = IMAGE_BAYER_NONE) {

		assert (depth == 3);
		_depth = 3;

		_dimx = dimx;
		_dimy = dimy;
		_offset = point(0, 0);
		_apm_memo = 0;
		_acm_memo = 0;
		_accm_memo = 0;
		this->name = name;
		this->_exp = _exp;
		this->bayer = bayer;

		if (_exp != NULL)
			_exp->add_listener(this, name);
	}

	unsigned int get_bayer() const {
		return bayer;
	}

	exposure &exp() const {
		return *_exp;
	}

	point offset() const {
		return _offset;
	}

	void set_offset(int i, int j) {
		_offset[0] = i;
		_offset[1] = j;
	}

	void set_offset(point p) {
		_offset = p;
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

	virtual spixel &pix(unsigned int y, unsigned int x) = 0;

	virtual const spixel &pix(unsigned int y, unsigned int x) const = 0;

	virtual ale_real &chan(unsigned int y, unsigned int x, unsigned int k) = 0;

	virtual const ale_real &chan(unsigned int y, unsigned int x, unsigned int k) const = 0;

	virtual void set_pixel(unsigned int y, unsigned int x, spixel p) {
		pix(y, x) = p;
	}
	
	virtual pixel get_pixel(unsigned int y, unsigned int x) const {
		return ((const image *)this)->pix(y, x);
	}

	ale_real maxval() const {
		ale_real result = get_pixel(0, 0)[0];

		for (unsigned int i = 0; i < _dimy; i++)
		for (unsigned int j = 0; j < _dimx; j++) {
			pixel p = get_pixel(i, j);

			for (unsigned int k = 0; k < _depth; k++)
				if (p[k] > result || !finite(result))
					result = p[k];
		}

		return result;
	}

	ale_real minval() const {
		ale_real result = get_pixel(0, 0)[0];

		for (unsigned int i = 0; i < _dimy; i++)
		for (unsigned int j = 0; j < _dimx; j++) {
			pixel p = get_pixel(i, j);

			for (unsigned int k = 0; k < _depth; k++)
				if (p[k] < result || !finite(result))
					result = p[k];
		}

		return result;
	}

	/*
	 * Get the maximum difference among adjacent pixels.
	 */

	pixel get_max_diff(unsigned int i, unsigned int j) const {
		assert(i >= 0);
		assert(j >= 0);
		assert(i <= _dimy - 1);
		assert(j <= _dimx - 1);

		pixel max = get_pixel(i, j), min = get_pixel(i, j);

		for (int ii = -1; ii <= 1; ii++)
		for (int jj = -1; jj <= 1; jj++) {
			int iii = i + ii;
			int jjj = j + jj;

			if (iii < 0)
				continue;
			if (jjj < 0)
				continue;
			if ((unsigned int) iii > _dimy - 1)
				continue;
			if ((unsigned int) jjj > _dimx - 1)
				continue;

			pixel p = get_pixel(iii, jjj);

			for (int d = 0; d < 3; d++) {
				if (p[d] > max[d])
					max[d] = p[d];
				if (p[d] < min[d])
					min[d] = p[d];
			}
		}

		return max - min;
	}

	pixel get_max_diff(point x) const {
		assert (x[0] >= 0);
		assert (x[1] >= 0);
		assert (x[0] <= _dimy - 1);
		assert (x[1] <= _dimx - 1);

		unsigned int i = (unsigned int) round(x[0]);
		unsigned int j = (unsigned int) round(x[1]);

		return get_max_diff(i, j);
	}

	/*
	 * Get a color value at a given position using bilinear interpolation between the
	 * four nearest pixels.  Result values:
	 *
	 * result[0] == pixel value
	 * result[1] == pixel confidence
	 */
	void get_bl(point x, pixel result[2]) const {

		assert (x[0] >= 0);
		assert (x[1] >= 0);
		assert (x[0] <= _dimy - 1);
		assert (x[1] <= _dimx - 1);

		int lx = (int) floor(x[1]);
		int hx = (int) floor(x[1]) + 1;
		int ly = (int) floor(x[0]);
		int hy = (int) floor(x[0]) + 1;

		pixel neighbor[4];
		ale_pos factor[4];

		neighbor[0] = get_pixel(ly, lx);
		neighbor[1] = get_pixel(hy % _dimy, lx);
		neighbor[2] = get_pixel(hy % _dimy, hx % _dimx);
		neighbor[3] = get_pixel(ly, hx % _dimx);

		factor[0] = (hx - x[1]) * (hy - x[0]);
		factor[1] = (hx - x[1]) * (x[0] - ly);
		factor[2] = (x[1] - lx) * (x[0] - ly);
		factor[3] = (x[1] - lx) * (hy - x[0]);

		/*
		 * Use bilinear interpolation for pixel value
		 */

		result[0] = pixel(0, 0, 0);

		for (int n = 0; n < 4; n++)
			result[0] += factor[n] * neighbor[n];

#if 0
		/*
		 * Take the minimum confidence
		 */

		if (_exp) {
			result[1] = _exp->confidence(neighbor[0]);

			for (int n = 1; n < 4; n++)
				if (factor[n] > 0) {
					pixel confidence 
						= _exp->confidence(neighbor[n]);

					for (int k = 0; k < 3; k++)
						if (confidence[k] < result[1][k])
							result[1][k] = confidence[k];
				}
		} else {
			result[1] = pixel(1, 1, 1);
		}
#else
		/*
		 * Use bilinear interpolation for confidence
		 */

		if  (_exp) {
			result[1] = pixel(0, 0, 0);
			for (int n = 0; n < 4; n++)
				result[1] += factor[n] * _exp->confidence(neighbor[n]);
		} else {
			result[1] = pixel(1, 1, 1);
		}
#endif
	}

	int in_bounds(point x) const {
		if (x[0] < 0
		 || x[1] < 0
		 || x[0] > height() - 1
		 || x[1] > width() - 1)
			return 0;
		if (!x.defined())
			return 0;
		return 1;
	}

	pixel get_bl(point x) const {
		pixel result[2];

		get_bl(x, result);

		return result[0];
	}

	void get_scaled_bl(point x, ale_pos f, pixel result[2]) const {
		point scaled(
				x[0]/f <= height() - 1
			               ? (x[0]/f)
				       : (height() - 1),
				x[1]/f <= width() - 1
				       ? (x[1]/f)
				       : (width() - 1));

		get_bl(scaled, result);
	}

	pixel get_scaled_bl(point x, ale_pos f) const {
		pixel result[2];

		get_scaled_bl(x, f, result);

		return result[0];
	}

		
	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, char *name) const = 0;

	/*
	 * Generate an image of medians within a given radius
	 */

	image *medians(int radius) const {

		assert (radius >= 0);
		
		image *is = scale_generator(height(), width(), depth(), "median");
		assert(is);

		for (unsigned int i = 0; i < height(); i++)
		for (unsigned int j = 0; j < width(); j++) {

			std::vector<ale_real> p[3];

			for (int ii = -radius; ii <= radius; ii++)
			for (int jj = -radius; jj <= radius; jj++) {
				int iii = i + ii;
				int jjj = j + jj;

				if (in_bounds(point(iii, jjj)))
					for (int k = 0; k < 3; k++)
						if (finite(get_pixel(iii, jjj)[k]))
							p[k].push_back(get_pixel(iii, jjj)[k]);
			}

			is->pix(i, j) = d2::pixel::undefined();

			for (int k = 0; k < 3; k++) {
				std::sort(p[k].begin(), p[k].end());

				unsigned int pkc = p[k].size();

				if (pkc == 0)
					continue;

				if (pkc % 2 == 0) 
					is->chan(i, j, k) = (p[k][pkc / 2] + p[k][pkc / 2 - 1]) / 2;
				else
					is->chan(i, j, k) = p[k][pkc / 2];
			}
		}

		return is;
	}

	/*
	 * Generate an image of differences of the first channel.  The first
	 * coordinate differences are stored in the first channel, second in the
	 * second channel.
	 */

	image *fcdiffs() const {
		image *is = scale_generator(height(), width(), depth(), "diff");

		assert(is);

		for (unsigned int i = 0; i < height(); i++)
		for (unsigned int j = 0; j < width(); j++) {

			if (i + 1 < height()
			 && i > 0
			 && !finite(chan(i, j, 0))) {

				is->chan(i, j, 0) = (chan(i + 1, j, 0) - chan(i - 1, j, 0)) / 2;

			} else if (i + 1 < height()
			        && i > 0
			        && finite(chan(i + 1, j, 0)
			        && finite(chan(i - 1, j, 0)) {

				is->chan(i, j, 0) = ((chan(i, j, 0) - chan(i - 1, j, 0))
						   + (chan(i + 1, j, 0) - chan(i, j, 0))) / 2;

			} else if (i + 1 < height()
				&& finite(chan(i + 1, j, 0))) {

				is->chan(i, j, 0) = chan(i + 1, j, 0) - chan(i, j, 0);

			} else if (i > 0
				&& finite(chan(i - 1, j, 0))) {
				 
				is->chan(i, j, 0) = chan(i, j, 0) - chan(i - 1, j, 0);

			} else {
				is->chan(i, j, 0) = 0;
			}

			if (j + 1 < width()
			 && j > 0
			 && !finite(chan(i, j, 0))) {

				is->chan(i, j, 1) = (chan(i, j + 1, 0) - chan(i, j - 1, 0)) / 2;

			} else if (j + 1 < width()
			        && j > 0
			        && finite(chan(i, j + 1, 0)
			        && finite(chan(i, j - 1, 0)) {

				is->chan(i, j, 1) = ((chan(i, j, 0) - chan(i, j - 1, 0))
						   + (chan(i, j + 1, 0) - chan(i, j, 0))) / 2;

			} else if (j + 1 < width()
				&& finite(chan(i, j + 1, 0))) {

				is->chan(i, j, 1) = chan(i, j + 1, 0) - chan(i, j, 0);

			} else if (j > 0
				&& finite(chan(i, j - 1, 0))) {
				 
				is->chan(i, j, 1) = chan(i, j, 0) - chan(i, j - 1, 0);

			} else {
				is->chan(i, j, 1) = 0;
			}
		}

		return is;
	}

	/*
	 * Generate an image of median (within a given radius) difference of the
	 * first channel.
	 */

	image *fcdiff_median(int radius) const {
		image *diff = fcdiffs();

		assert(diff);

		image *median = diff->medians(radius);

		assert(median);

		delete diff;

		return median;
	}

	/*
	 * Scale by half.  We use the following filter:
	 *
	 * 	1/16	1/8	1/16
	 * 	1/8	1/4	1/8
	 * 	1/16	1/8	1/16
	 *
	 * At the edges, these values are normalized so that the sum of the
	 * weights of contributing pixels is 1.
	 */
	image *scale_by_half(char *name) const {
		ale_pos f = 0.5;

		image *is = scale_generator(
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth(), name);

		assert(is);

		for (unsigned int i = 0; i < is->height(); i++)
		for (unsigned int j = 0; j < is->width(); j++) 

			is->set_pixel(i, j, 

			      ( ( ((i > 0 && j > 0) 
				    ? get_pixel(2 * i - 1, 2 * j - 1) * (ale_real) 0.0625
				    : pixel(0, 0, 0))
				+ ((i > 0)
				    ? get_pixel(2 * i - 1, 2 * j) * 0.125
				    : pixel(0, 0, 0))
				+ ((i > 0 && j < is->width() - 1)
				    ? get_pixel(2 * i - 1, 2 * j + 1) * 0.0625
				    : pixel(0, 0, 0))
				+ ((j > 0)
				    ? get_pixel(2 * i, 2 * j - 1) * 0.125
				    : pixel(0, 0, 0))
				+ get_pixel(2 * i, 2 * j) * 0.25
				+ ((j < is->width() - 1)
				    ? get_pixel(2 * i, 2 * j + 1) * 0.125
				    : pixel(0, 0, 0))
				+ ((i < is->height() - 1 && j > 0)
				    ? get_pixel(2 * i + 1, 2 * j - 1) * 0.0625
				    : pixel(0, 0, 0))
				+ ((i < is->height() - 1)
				    ? get_pixel(2 * i + 1, 2 * j) * 0.125
				    : pixel(0, 0, 0))
				+ ((i < is->height() && j < is->width() - 1)
				    ? get_pixel(2 * i + 1, 2 * j + 1) * 0.0625 
				    : pixel(0, 0, 0)))

			     /

				( ((i > 0 && j > 0) 
				    ? 0.0625
				    : 0)
				+ ((i > 0)
				    ? 0.125
				    : 0)
				+ ((i > 0 && j < is->width() - 1)
				    ? 0.0625
				    : 0)
				+ ((j > 0)
				    ? 0.125
				    : 0)
				+ 0.25
				+ ((j < is->width() - 1)
				    ? 0.125
				    : 0)
				+ ((i < is->height() - 1 && j > 0)
				    ? 0.0625
				    : 0)
				+ ((i < is->height() - 1)
				    ? 0.125
				    : 0)
				+ ((i < is->height() && j < is->width() - 1)
				    ? 0.0625
				    : 0) ) ) );

		is->_offset = point(_offset[0] * f, _offset[1] * f);

		return is;
	}

	/*
	 * Scale by half.  This function uses externally-provided weights,
	 * multiplied by the following filter:
	 *
	 * 	1/16	1/8	1/16
	 * 	1/8	1/4	1/8
	 * 	1/16	1/8	1/16
	 *
	 * Values are normalized so that the sum of the weights of contributing
	 * pixels is 1.
	 */
	image *scale_by_half(const image *weights, char *name) const {

		if (weights == NULL)
			return scale_by_half(name);

		ale_pos f = 0.5;

		image *is = scale_generator(
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth(), name);

		assert(is);

		for (unsigned int i = 0; i < is->height(); i++)
		for (unsigned int j = 0; j < is->width(); j++) {

			pixel value = pixel

			      ( ( ((i > 0 && j > 0) 
				    ? get_pixel(2 * i - 1, 2 * j - 1) 
				    * weights->get_pixel(2 * i - 1, 2 * j - 1) 
				    * (ale_real) 0.0625
				    : pixel(0, 0, 0))
				+ ((i > 0)
				    ? get_pixel(2 * i - 1, 2 * j) 
				    * weights->get_pixel(2 * i - 1, 2 * j)
				    * 0.125
				    : pixel(0, 0, 0))
				+ ((i > 0 && j < is->width() - 1)
				    ? get_pixel(2 * i - 1, 2 * j + 1) 
				    * weights->get_pixel(2 * i - 1, 2 * j + 1)
				    * 0.0625
				    : pixel(0, 0, 0))
				+ ((j > 0)
				    ? get_pixel(2 * i, 2 * j - 1) 
				    * weights->get_pixel(2 * i, 2 * j - 1)
				    * 0.125
				    : pixel(0, 0, 0))
				+ get_pixel(2 * i, 2 * j) 
				    * weights->get_pixel(2 * i, 2 * j)
				    * 0.25
				+ ((j < is->width() - 1)
				    ? get_pixel(2 * i, 2 * j + 1) 
				    * weights->get_pixel(2 * i, 2 * j + 1)
				    * 0.125
				    : pixel(0, 0, 0))
				+ ((i < is->height() - 1 && j > 0)
				    ? get_pixel(2 * i + 1, 2 * j - 1) 
				    * weights->get_pixel(2 * i + 1, 2 * j - 1)
				    * 0.0625
				    : pixel(0, 0, 0))
				+ ((i < is->height() - 1)
				    ? get_pixel(2 * i + 1, 2 * j) 
				    * weights->get_pixel(2 * i + 1, 2 * j)
				    * 0.125
				    : pixel(0, 0, 0))
				+ ((i < is->height() && j < is->width() - 1)
				    ? get_pixel(2 * i + 1, 2 * j + 1) 
				    * weights->get_pixel(2 * i + 1, 2 * j + 1)
				    * 0.0625 
				    : pixel(0, 0, 0)))

			     /

				( ((i > 0 && j > 0) 
				    ? weights->get_pixel(2 * i - 1, 2 * j - 1)
				    * 0.0625
				    : pixel(0, 0, 0))
				+ ((i > 0)
				    ? weights->get_pixel(2 * i - 1, 2 * j)
				    * 0.125
				    : pixel(0, 0, 0))
				+ ((i > 0 && j < is->width() - 1)
				    ? weights->get_pixel(2 * i - 1, 2 * j + 1)
				    * 0.0625
				    : pixel(0, 0, 0))
				+ ((j > 0)
				    ? weights->get_pixel(2 * i, 2 * j - 1)
				    * 0.125
				    : pixel(0, 0, 0))
				+ weights->get_pixel(2 * i, 2 * j) 
				* 0.25
				+ ((j < is->width() - 1)
				    ? weights->get_pixel(2 * i, 2 * j + 1)
				    * 0.125
				    : pixel(0, 0, 0))
				+ ((i < is->height() - 1 && j > 0)
				    ? weights->get_pixel(2 * i + 1, 2 * j - 1)
				    * 0.0625
				    : pixel(0, 0, 0))
				+ ((i < is->height() - 1)
				    ? weights->get_pixel(2 * i + 1, 2 * j)
				    * 0.125
				    : pixel(0, 0, 0))
				+ ((i < is->height() && j < is->width() - 1)
				    ? weights->get_pixel(2 * i + 1, 2 * j + 1)
				    * 0.0625
				    : pixel(0, 0, 0)) ) );

			for (int k = 0; k < 3; k++)
				if (!finite(value[k]))
					value[k] = 0;

			is->set_pixel(i, j, value);
		}

		is->_offset = point(_offset[0] * f, _offset[1] * f);

		return is;
	}

	/*
	 * Return an image scaled by some factor != 1.0, using bilinear
	 * interpolation.
	 */
	image *scale(ale_pos f, char *name) const {

		/*
		 * We probably don't want to scale images by a factor of 1.0,
		 * or by non-positive values.
		 */
		assert (f != 1.0 && f > 0);

		if (f > 1.0) {
			image *is = scale_generator(
				(int) floor(height() * f), 
				(int) floor(width()  * f), depth(), name);

			assert(is);

			unsigned int i, j, k;

			for (i = 0; i < is->height(); i++)
			for (j = 0; j < is->width(); j++) 
			for (k = 0; k < is->depth(); k++)
				is->set_pixel(i, j,
					get_scaled_bl(point(i, j), f));

			return is;
		} else if (f == 0.5) {
			return scale_by_half(name);
		} else {
			image *is = scale(2*f, name);
			image *result = is->scale(0.5, name);
			delete is;
			return result;
		}

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

		image *is = scale_generator(
			(int) floor(height() * f), 
			(int) floor(width()  * f), depth(), name);

		assert(is);

		for (unsigned int i = 0; i < is->height(); i++)
		for (unsigned int j = 0; j < is->width(); j++) 

			is->set_pixel(i, j, 

				( ((i > 0 && j > 0) 
				    ? get_pixel(2 * i - 1, 2 * j - 1)
				    : pixel())
				* ((i > 0)
				    ? get_pixel(2 * i - 1, 2 * j)
				    : pixel())
				* ((i > 0 && j < is->width() - 1)
				    ? get_pixel(2 * i - 1, 2 * j + 1)
				    : pixel())
				* ((j > 0)
				    ? get_pixel(2 * i, 2 * j - 1)
				    : pixel())
				* get_pixel(2 * i, 2 * j)
				* ((j < is->width() - 1)
				    ? get_pixel(2 * i, 2 * j + 1)
				    : pixel())
				* ((i < is->height() - 1 && j > 0)
				    ? get_pixel(2 * i + 1, 2 * j - 1)
				    : pixel())
				* ((i < is->height() - 1)
				    ? get_pixel(2 * i + 1, 2 * j)
				    : pixel())
				* ((i < is->height() && j < is->width() - 1)
				    ? get_pixel(2 * i + 1, 2 * j + 1)
				    : pixel())));

		is->_offset = point(_offset[0] * f, _offset[1] * f);

		return is;
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.  Negative values
	 * shrink the image.
	 */
	virtual void extend(int top, int bottom, int left, int right) = 0;

	/*
	 * Clone 
	 */
	image *clone(char *name) const {
		image *ic = scale_generator(
			height(), width(), depth(), name);

		assert(ic);

		for (unsigned int i = 0; i < height(); i++)
		for (unsigned int j = 0; j < width();  j++)
			ic->set_pixel(i, j, 
				get_pixel(i, j));


		ic->_offset = _offset;

		ic->_apm_memo = _apm_memo;
		ic->_acm_memo = _acm_memo;
		ic->_accm_memo = _accm_memo;
		ic->_apm = _apm;
		ic->_acm = _acm;
		ic->_accm = _accm;

		return ic;
	}

	/*
	 * Calculate the average (mean) clamped magnitude of a channel across
	 * all pixels in an image.  The magnitude is clamped to the range of 
	 * real inputs.
	 */
	ale_real avg_channel_clamped_magnitude(unsigned int k) const {

		/*
		 * This is a memoized function
		 */

		assert (k < _depth);

		avg_channel_clamped_magnitude_memo();
		return _accm[k];
	}

	pixel avg_channel_clamped_magnitude() const {
		avg_channel_clamped_magnitude_memo();
		return _accm;
	}

	/*
	 * Calculate the average (mean) magnitude of a channel across all
	 * pixels in an image.
	 */
	ale_real avg_channel_magnitude(unsigned int k) const {

		/*
		 * This is a memoized function
		 */

		assert (k < _depth);

		avg_channel_magnitude_memo();
		return _acm[k];
	}

	pixel avg_channel_magnitude() const {
		avg_channel_magnitude_memo();
		return _acm;
	}

	/*
	 * Calculate the average (mean) magnitude of a pixel (where magnitude
	 * is defined as the mean of the channel values).
	 */
	ale_real avg_pixel_magnitude() const {
		unsigned int i, j, k;

		ale_accum accumulator;
		ale_accum divisor = 0;

		if (_apm_memo)
			return _apm;

		_apm_memo = 1;
		accumulator = 0;

		for (i = 0; i < _dimy;  i++)
		for (j = 0; j < _dimx;  j++) {
			pixel value = get_pixel(i, j);

			for (k = 0; k < _depth; k++) 
				if (finite(value[k])) {
					accumulator += value[k];
					divisor++;
				}
		}

		accumulator /= divisor;

		_apm = accumulator;

		return _apm;
	}

	virtual ~image() {
	}
};

#endif
