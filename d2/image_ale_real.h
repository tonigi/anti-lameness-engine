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
 * image_ale_real.h: Image represented by an array of ale_reals
 */

#ifndef __image_ale_real_h__
#define __image_ale_real_h__

#include "exposure/exposure.h"
#include "point.h"
#include "image.h"

class image_ale_real : public image {
private:
	spixel *_p;

private:
	void trigger(pixel multiplier) {
		for (unsigned int i = 0; i < _dimx * _dimy; i++)
			_p[i] *= multiplier;
	}

public:
	image_ale_real (unsigned int dimy, unsigned int dimx, unsigned int
			depth, char *name = "anonymous", exposure *exp = NULL) 
			: image(dimy, dimx, depth, name, exp) {

		_p = new spixel[dimx * dimy];

		assert (_p);

		if (!_p) {
			fprintf(stderr, "Could not allocate memory for image data.\n");
			exit(1);
		}
	}

	virtual ~image_ale_real() {
		delete[] _p;
	}

	spixel &pix(unsigned int y, unsigned int x) {
		assert (x < _dimx);
		assert (y < _dimy);

		image_updated();

		return _p[y * _dimx + x];
	}

	const spixel &pix(unsigned int y, unsigned int x) const {
		assert (x < _dimx);
		assert (y < _dimy);

		return _p[y * _dimx + x];
	}

	ale_real &chan(unsigned int y, unsigned int x, unsigned int k) {
		return pix(y, x)[k];
	}

	const ale_real &chan(unsigned int y, unsigned int x, unsigned int k) const {
		return pix(y, x)[k];
	}

	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, char *name) const {
		return new image_ale_real(height, width, depth, name, _exp);
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.
	 */
	void extend(unsigned int top, unsigned int bottom, 
			unsigned int left, unsigned int right) {

		image_ale_real *is = new image_ale_real (
			height() + top  + bottom,
			 width() + left + right , depth(), name, _exp);

		assert(is);

		for (unsigned int i = 0; i < height(); i++)
		for (unsigned int j = 0; j <  width(); j++) 
			is->pix(i + top, j + left)
				= get_pixel(i, j);

		delete[] _p;

		_p = is->_p;
		_dimx = is->_dimx;
		_dimy = is->_dimy;
		_depth = is->_depth;
		_offset[0] -= top;
		_offset[1] -= left;

		image_updated();

		is->_p = NULL;

		delete is;
	}

};

#endif
