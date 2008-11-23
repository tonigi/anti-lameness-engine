// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
 * image_ale_real.h: Image represented by an array of ale_reals
 */

#ifndef __image_ale_real_h__
#define __image_ale_real_h__

#include "exposure/exposure.h"
#include "point.h"
#include "image.h"

#define RESIDENT_DIVISIONS 200


template <int disk_support>
class image_ale_real : public image {
private:
	/*
	 * Data structures without file support.
	 */

	spixel *_p;

	/*
	 * Data structures for file support.
	 */

	FILE *support;
	mutable spixel *_p_segments[RESIDENT_DIVISIONS];
	mutable int dirty_segments[RESIDENT_DIVISIONS];
	mutable int resident_list[RESIDENT_DIVISIONS];
	mutable int resident_next;
	int resident_max;
	int rows_per_segment;
	mutable thread::rwlock_t rwlock;

public:

	/*
	 * Wrapper encapsulating details of the separation between the
	 * resident-checking implementation and non-checking.
	 */

	static image *new_image_ale_real(unsigned int dimy,
					 unsigned int dimx,
					 unsigned int depth,
					 const char *name = "anonymous",
					 exposure *exp = NULL) {

		double resident = image::get_resident();

		if (resident == 0 || resident * 1000000 >= dimy * dimx)
			return new image_ale_real<0>(dimy, dimx, depth, name, exp);

		return new image_ale_real<1>(dimy, dimx, depth, name, exp);

	}

	image_ale_real (unsigned int dimy, unsigned int dimx, unsigned int
			depth, const char *name = "anonymous", exposure *exp = NULL) 
			: image(dimy, dimx, depth, name, exp) {

		if (disk_support == 0) {
			_p = new spixel[dimx * dimy];

			assert (_p);

			if (!_p) {
				fprintf(stderr, "Could not allocate memory for image data.\n");
				exit(1);
			}
		} else {
			rows_per_segment = (int) ceil((double) dimy / (double) RESIDENT_DIVISIONS);

			assert (rows_per_segment > 0);

			for (int i = 0; i < RESIDENT_DIVISIONS; i++) {
				_p_segments[i] = NULL;
				dirty_segments[i] = 0;
				resident_list[i] = -1;
			}

			resident_max = (unsigned int) floor((image::get_resident() * 1000000)
				                          / (rows_per_segment * dimx));

			assert (resident_max <= RESIDENT_DIVISIONS);

			if (resident_max == 0) {
				ui::get()->error_hint(
					"No segments resident in image array.",
					"Try recompiling with more RESIDENT_DIVISIONS");
			}

			resident_next = 0;

			support = tmpfile();

			if (!support) {
				ui::get()->error_hint(
					"Unable to create temporary file to support image array.",
					"Set --resident 0, or Win32/64 users might run as admin.");
			}

			spixel *zero = new spixel[dimx];

			assert(zero);

			for (unsigned int i = 0; i < dimy; i++) {
				unsigned int c = fwrite(zero, sizeof(spixel), dimx, support);
				if (c < dimx)
					ui::get()->error_hint("Image array support file error.",
					                      "Submit a bug report.");
			}

			delete[] zero;
		}
	}

	virtual ~image_ale_real() {
		if (disk_support == 0) {
			delete[] _p;
		} else {
			for (int i = 0; i < RESIDENT_DIVISIONS; i++) {
				if (_p_segments[i])
					delete[] _p_segments[i];
			}

			fclose(support);
		}
	}

	void resident_begin(unsigned int segment) const {
		rwlock.rdlock();
		if (_p_segments[segment])
			return;
		rwlock.unlock();

		rwlock.wrlock();

		if (_p_segments[segment])
			return;

		if (resident_list[resident_next] >= 0) {
			/*
			 * Eject a segment
			 */

			if (dirty_segments[resident_list[resident_next]]) {
				fseek(support, rows_per_segment * _dimx * sizeof(spixel) 
				             * resident_list[resident_next], 
					       SEEK_SET);
				assert(_p_segments[resident_list[resident_next]]);
				size_t fwrite_result = fwrite(_p_segments[resident_list[resident_next]], 
					sizeof(spixel), rows_per_segment * _dimx, support);

				assert(fwrite_result == rows_per_segment * _dimx);

				dirty_segments[resident_list[resident_next]] = 0;
			}

			delete[] _p_segments[resident_list[resident_next]];
			_p_segments[resident_list[resident_next]] = NULL;
		}

		resident_list[resident_next] = segment;

		_p_segments[segment] = new spixel[_dimx * rows_per_segment];

		assert (_p_segments[segment]);

		fseek(support, rows_per_segment * _dimx * sizeof(spixel)
		             * segment,
			     SEEK_SET);

		size_t fread_result = fread(_p_segments[segment], sizeof(spixel), rows_per_segment * _dimx, support);

		assert(fread_result == rows_per_segment * _dimx);

		/*
		 * Update the next ejection candidate.
		 */
		resident_next++;
		if (resident_next >= resident_max)
			resident_next = 0;
	}

	void resident_end(unsigned int segment) const {
		rwlock.unlock();
	}

	spixel get_pixel(unsigned int y, unsigned int x) const {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			return _p[y * _dimx + x];
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);

			spixel result = _p_segments[segment][(y % rows_per_segment) * _dimx + x];
				
			resident_end(segment);

			return result;
		}
	}

	void set_pixel(unsigned int y, unsigned int x, spixel p) {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			_p[y * _dimx + x] = p;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x] = p;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	void mul_pixel(unsigned int y, unsigned int x, spixel p) {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			_p[y * _dimx + x] *= p;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x] *= p;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	void add_pixel(unsigned int y, unsigned int x, pixel p) {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			_p[y * _dimx + x] += p;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x] += p;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	ale_sreal get_chan(unsigned int y, unsigned int x, unsigned int k) const {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			return _p[y * _dimx + x][k];
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);

			ale_sreal result = _p_segments[segment]
			                              [(y % rows_per_segment) * _dimx + x][k];
				
			resident_end(segment);

			return result;
		}
	}

	void set_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			_p[y * _dimx + x][k] = c;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x][k] = c;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	void div_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert (x < _dimx);
		assert (y < _dimy);

		if (disk_support == 0) {
			_p[y * _dimx + x][k] /= c;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x][k] /= c;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, const char *name) const {
		return new_image_ale_real(height, width, depth, name, _exp);
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.  Negative values
	 * shrink the image.
	 */
	image *_extend(int top, int bottom, int left, int right) {

		image *is = new_image_ale_real (
			height() + top  + bottom,
			 width() + left + right , depth(), name, _exp);

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

/*
 * Wrapper encapsulating details of the separation between the
 * resident-checking implementation and non-checking.
 */
static inline image *new_image_ale_real(unsigned int dimy,
                                 unsigned int dimx,
				 unsigned int depth,
				 const char *name = "anonymous",
				 exposure *exp = NULL) {

	return image_ale_real<0>::new_image_ale_real(dimy, dimx, depth, name, exp);
}


#endif
