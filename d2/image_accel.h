// Copyright 2002, 2003, 2004, 2008 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                                <dhilvert@ugcs.caltech.edu>,
//                                                <dhilvert@gmail.com>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
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
 * image_accel.h: Image represented by a (possibly accelerated) opaque type.
 */

#ifndef __image_accel_h__
#define __image_accel_h__

#include "exposure/exposure.h"
#include "point.h"
#include "image.h"

class image_accel : public image {
private:
	/*
	 * Libale opaque type.
	 */

	void *im;

public:

	image_accel (unsigned int dimy, unsigned int dimx, unsigned int
			depth, unsigned int bayer, const char *name = "anonymous", 
			exposure *exp = NULL) 
			: image(dimy, dimx, depth, name, exp, bayer) {

		assert (depth == 3);

		int libale_bayer;

		switch (bayer) {
		case IMAGE_BAYER_NONE:
			libale_bayer = ALE_FORMAT_RGB;
			break;
		case IMAGE_BAYER_RGBG:
			libale_bayer = ALE_FORMAT_BAYER_RGBG;
			break;
		case IMAGE_BAYER_GBGR:
			libale_bayer = ALE_FORMAT_BAYER_GBGR;
			break;
		case IMAGE_BAYER_GRGB:
			libale_bayer = ALE_FORMAT_BAYER_GRGB;
			break;
		case IMAGE_BAYER_BGRG:
			libale_bayer = ALE_FORMAT_BAYER_BGRG;
			break;
		default:
			assert(0);
		}

		im = ale_new_domain_2d(accel::context(), libale_bayer, ALE_TYPE_FLOAT_32, dimy, dimx);

		assert (im);

		if (!im) {
			fprintf(stderr, "Could not allocate image data.\n");
			exit(1);
		}
	}

	image_accel (const image *source) : image(source, 0) {

		assert (_depth == 3);

		int libale_bayer;

		switch (bayer) {
		case IMAGE_BAYER_NONE:
			libale_bayer = ALE_FORMAT_RGB;
			break;
		case IMAGE_BAYER_RGBG:
			libale_bayer = ALE_FORMAT_BAYER_RGBG;
			break;
		case IMAGE_BAYER_GBGR:
			libale_bayer = ALE_FORMAT_BAYER_GBGR;
			break;
		case IMAGE_BAYER_GRGB:
			libale_bayer = ALE_FORMAT_BAYER_GRGB;
			break;
		case IMAGE_BAYER_BGRG:
			libale_bayer = ALE_FORMAT_BAYER_BGRG;
			break;
		default:
			assert(0);
		}

		im = ale_new_domain_2d(accel::context(), libale_bayer, ALE_TYPE_FLOAT_32, _dimy, _dimx);

		if (!im) {
			fprintf(stderr, "Could not allocate Libale domain.\n");
			exit(1);
		}

		/*
		 * Copy image data
		 */

		float *data = NULL;

		if (bayer == IMAGE_BAYER_NONE) {
			data = (float *) malloc(_dimy * _dimx * _depth * sizeof(float));

			if (!data) {
				fprintf(stderr, "Could not allocate image data.\n");
				exit(1);
			}

			for (unsigned int i = 0; i < _dimy; i++)
			for (unsigned int j = 0; j < _dimx; j++)
			for (unsigned int k = 0; k < _depth; k++)
				data[i * _dimx * _depth + j * _depth + k] = source->get_chan(i, j, k);

		} else {
			data = (float *) malloc(_dimy * _dimx * sizeof(float));

			if (!data) {
				fprintf(stderr, "Could not allocate image data.\n");
				exit(1);
			}

			for (unsigned int i = 0; i < _dimy; i++)
			for (unsigned int j = 0; j < _dimx; j++)
				data[i * _dimx + j] = source->get_chan(i, j, source->bayer_color(i, j));
		}

		ale_load_into_domain(im, data);

		free(data);
	}

	virtual ~image_accel() {
		ale_delete_domain_2d(im);
	}

	virtual int accel_type() {
		return 1;
	}

	virtual image *unaccel_equiv() {
		/*
		 * XXX: An unaccelerated equivalent should be created here.
		 */

		assert(0);

		return NULL;
	}

	spixel get_pixel(unsigned int y, unsigned int x) const {
		assert(0);
	}

	void set_pixel(unsigned int y, unsigned int x, spixel p) {
		assert(0);
	}

	void mul_pixel(unsigned int y, unsigned int x, spixel p) {
		assert(0);
	}

	void add_pixel(unsigned int y, unsigned int x, pixel p) {
		assert(0);
	}

	ale_sreal get_chan(unsigned int y, unsigned int x, unsigned int k) const {
		assert(0);
	}

	void set_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert(0);
	}

	void div_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert(0);
	}

	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, const char *name) const {
		return new image_accel(height, width, depth, IMAGE_BAYER_NONE, name, _exp);
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.  Negative values
	 * shrink the image.
	 */
	image *_extend(int top, int bottom, int left, int right) {

		image *is = new image_accel (
			height() + top  + bottom,
			 width() + left + right , depth(), bayer, name, _exp);

		assert(is);

		unsigned int min_i = (-top > 0)
			      ? -top
			      : 0;

		unsigned int min_j = (-left > 0)
			      ? -left
			      : 0;

		unsigned int max_i = (height() < is->height() - top)
			      ? height()
			      : is->height() - top;

		unsigned int max_j = (width() < is->width() - left)
			      ? width()
			      : is->width() - left;

		for (unsigned int i = min_i; i < max_i; i++)
		for (unsigned int j = min_j; j < max_j; j++) 
			is->set_pixel(i + top, j + left, get_pixel(i, j));

		is->set_offset(_offset[0] - top, _offset[1] - left);

		return is;
	}

private:
	void trigger(pixel multiplier) {
		for (unsigned int i = 0; i < height(); i++)
		for (unsigned int j = 0; j < width(); j++) {
			mul_pixel(i, j, multiplier);
		}
	}
};

#endif
