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
 * image_bayer_ale_real.h: Bayer-patterned image represented by an array of ale_reals
 */

#ifndef __image_bayer_ale_real_h__
#define __image_bayer_ale_real_h__

#include "exposure/exposure.h"
#include "point.h"
#include "image.h"
#include "image_ale_real.h"

class image_bayer_ale_real : public image {
private:
	ale_real *_p;

private:
	/*
	 * X offset of 'R' element
	 */
	unsigned int r_x_offset() const {
		return bayer & 0x1;
	}

	/*
	 * Y offset of 'R' element
	 */
	unsigned int r_y_offset() const {
		return (bayer & 0x2) >> 1;
	}

public:
	/*
	 * Return the color of a given pixel.
	 */
	unsigned int bayer_color(unsigned int i, unsigned int j) const {
		return (i + r_y_offset()) % 2 + (j + r_x_offset()) % 2;
	}

private:
	void trigger(pixel multiplier) {
		for (unsigned int i = 0; i < _dimx; i++)
		for (unsigned int j = 0; j < _dimy; j++)
			_p[i * _dimy + j] *= multiplier[bayer_color(i, j)];
	}

public:
	image_bayer_ale_real (unsigned int dimy, unsigned int dimx, unsigned int depth,
			unsigned int bayer, char *name = "anonymous", exposure *exp = NULL) 
			: image(dimy, dimx, depth, name, exp, bayer) {

		assert (bayer == IMAGE_BAYER_BGRG
		     || bayer == IMAGE_BAYER_GBGR
		     || bayer == IMAGE_BAYER_GRGB
		     || bayer == IMAGE_BAYER_RGBG);

		_p = new ale_real[dimx * dimy];

		assert (_p);

		if (!_p) {
			fprintf(stderr, "Could not allocate memory for image data.\n");
			exit(1);
		}
	}

	virtual ~image_bayer_ale_real() {
		delete[] _p;
	}

	ale_real &chan(unsigned int y, unsigned int x, unsigned int k) {
		assert (k == bayer_color(y, x));
		return _p[y * _dimx + x];
	}

	const ale_real &chan(unsigned int y, unsigned int x, unsigned int k) const {
		assert (k == bayer_color(y, x));
		return _p[y * _dimx + x];
	}

	spixel &pix(unsigned int y, unsigned int x) {
		assert(0);

		static spixel foo;
		return foo;
	}

	const spixel &pix(unsigned int y, unsigned int x) const {
		static spixel foo = get_pixel(y, x);
		return foo;
	}

	/*
	 * This method throws away data not stored at this pixel
	 * position.
	 */
	void set_pixel(unsigned int y, unsigned int x, spixel p) {
		chan(y, x, bayer_color(y, x)) = p[bayer_color(y, x)];
	}

	/*
	 * This method uses bilinear interpolation.
	 */
	pixel get_pixel(unsigned int y, unsigned int x) const {
		pixel result;
		unsigned int k = bayer_color(y, x);
		ale_real sum;
		unsigned int num;

		result[k] = chan(y, x, k);

		if (k == 1) {
			unsigned int k1 = bayer_color(y + 1, x);
			unsigned int k2 = 2 - k1;

			sum = 0; num = 0;
			if (y > 0) {
				sum += chan(y - 1, x, k1);
				num++;
			}
			if (y < _dimy - 1) {
				sum += chan(y + 1, x, k1);
				num++;
			}
			assert (num > 0);
			result[k1] = sum / num;

			sum = 0; num = 0;
			if (x > 0) {
				sum += chan(y, x - 1, k2);
				num++;
			}
			if (x < _dimx - 1) {
				sum += chan(y, x + 1, k2);
				num++;
			}
			assert (num > 0);
			result[k2] = sum / num;

			return result;
		}

		sum = 0; num = 0;
		if (y > 0) {
			sum += chan(y - 1, x, 1);
			num++;
		}
		if (x > 0) {
			sum += chan(y, x - 1, 1);
			num++;
		}
		if (y < _dimy - 1) {
			sum += chan(y + 1, x, 1);
			num++;
		}
		if (x < _dimx - 1) {
			sum += chan(y, x + 1, 1);
			num++;
		}
		assert (num > 0);
		result[1] = sum / num;

		sum = 0; num = 0;
		if (y > 0 && x > 0) {
			sum += chan(y - 1, x - 1, 2 - k);
			num++;
		} 
		if (y > 0 && x < _dimx - 1) {
			sum += chan(y - 1, x + 1, 2 - k);
			num++;
		}
		if (y < _dimy - 1 && x > 0) {
			sum += chan(y + 1, x - 1, 2 - k);
			num++;
		}
		if (y < _dimy - 1 && x < _dimx - 1) {
			sum += chan(y + 1, x + 1, 2 - k);
			num++;
		}
		result[2 - k] = sum/num;

		return result;
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
	void extend(int top, int bottom, int left, int right) {
		/*
		 * Bayer-patterned images should always represent inputs,
		 * which should not ever be extended.
		 */
		assert(0);
	}

};

#endif
