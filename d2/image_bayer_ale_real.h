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
 * image_bayer_ale_real.h: Bayer-patterned image represented by an array of ale_reals
 */

#ifndef __image_bayer_ale_real_h__
#define __image_bayer_ale_real_h__

#include "exposure/exposure.h"
#include "point.h"
#include "image.h"
#include "image_ale_real.h"

template <int disk_support>
class image_bayer_ale_real : public image {
private:
	/*
	 * Data structures without file support.
	 */

	ale_sreal *_p;

	/*
	 * Data structures for file support.
	 */

	FILE *support;
	mutable ale_sreal *_p_segments[RESIDENT_DIVISIONS];
	mutable int dirty_segments[RESIDENT_DIVISIONS];
	mutable int resident_list[RESIDENT_DIVISIONS];
	mutable int resident_next;
	int resident_max;
	int rows_per_segment;
	mutable thread::rwlock_t rwlock;

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

	char get_channels(int i, int j) const {
		return (1 << bayer_color(i, j));
	}

public:

	/*
	 * Wrapper encapsulating details of the separation between the
	 * resident-checking implementation and non-checking.
	 */
	static inline image *new_image_bayer_ale_real(unsigned int dimy,
					 unsigned int dimx,
					 unsigned int depth,
					 unsigned int bayer,
					 const char *name = "anonymous",
					 exposure *exp = NULL) {

		unsigned int resident = image::get_resident();

		if (resident == 0 || resident * 1000000 > dimy * dimx)
			return new image_bayer_ale_real<0>(dimy, dimx, depth, bayer, name, exp);

		return new image_bayer_ale_real<1>(dimy, dimx, depth, bayer, name, exp);

	}

	image_bayer_ale_real (unsigned int dimy, unsigned int dimx, unsigned int depth,
			unsigned int bayer, const char *name = "anonymous", exposure *exp = NULL) 
			: image(dimy, dimx, depth, name, exp, bayer) {

		assert (bayer == IMAGE_BAYER_BGRG
		     || bayer == IMAGE_BAYER_GBGR
		     || bayer == IMAGE_BAYER_GRGB
		     || bayer == IMAGE_BAYER_RGBG);

		if (disk_support == 0) {
			_p = new ale_sreal[dimx * dimy];

			assert (_p);

			if (!_p) {
				fprintf(stderr, "Could not allocate memory for image data.\n");
				exit(1);
			}

			for (unsigned int i = 0; i < dimx * dimy; i++) {
				_p[i] = 0;
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

			ale_sreal *zero = new ale_sreal[dimx];

			assert(zero);

			for (unsigned int i = 0; i < dimx; i++) 
				zero[i] = 0;

			for (unsigned int i = 0; i < dimy; i++) {
				unsigned int c = fwrite(zero, sizeof(ale_sreal), dimx, support);
				if (c < dimx)
					ui::get()->error_hint("Image array support file error.",
					                      "Submit a bug report.");
			}

			delete[] zero;
		}
	}

	virtual ~image_bayer_ale_real() {
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
				fseek(support, rows_per_segment * _dimx * sizeof(ale_sreal) 
				             * resident_list[resident_next], 
					       SEEK_SET);
				assert(_p_segments[resident_list[resident_next]]);
				size_t fwrite_result = fwrite(_p_segments[resident_list[resident_next]], 
					sizeof(ale_sreal), rows_per_segment * _dimx, support);

				assert(fwrite_result == rows_per_segment * _dimx);

				dirty_segments[resident_list[resident_next]] = 0;
			}

			delete[] _p_segments[resident_list[resident_next]];
			_p_segments[resident_list[resident_next]] = NULL;
		}

		resident_list[resident_next] = segment;

		_p_segments[segment] = new ale_sreal[_dimx * rows_per_segment];

		assert (_p_segments[segment]);

		fseek(support, rows_per_segment * _dimx * sizeof(ale_sreal)
		             * segment,
			     SEEK_SET);

		size_t fread_result = fread(_p_segments[segment], sizeof(ale_sreal), rows_per_segment * _dimx, support);

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

	void set_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert (k == bayer_color(y, x));
		if (disk_support == 0) {
			_p[y * _dimx + x] = c;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x] = c;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	void add_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert (k == bayer_color(y, x));
		if (disk_support == 0) {
			_p[y * _dimx + x] += c;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x] += c;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	void div_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert (k == bayer_color(y, x));
		if (disk_support == 0) {
			_p[y * _dimx + x] /= c;
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);
				
			_p_segments[segment][(y % rows_per_segment) * _dimx + x] /= c;
			dirty_segments[segment] = 1;
				
			resident_end(segment);
		}
	}

	ale_sreal get_chan(unsigned int y, unsigned int x, unsigned int k) const {
#if 0
		/*
		 * This may be expensive.
		 */
		assert (k == bayer_color(y, x));
#endif
		if (disk_support == 0) {
			return _p[y * _dimx + x];
		} else {
			int segment = y / rows_per_segment;
			assert (segment < RESIDENT_DIVISIONS);

			resident_begin(segment);

			ale_sreal result = _p_segments[segment]
			                              [(y % rows_per_segment) * _dimx + x];
				
			resident_end(segment);

			return result;
		}
	}

	/*
	 * This method throws away data not stored at this pixel
	 * position.
	 */
	void set_pixel(unsigned int y, unsigned int x, spixel p) {
		set_chan(y, x, bayer_color(y, x), p[bayer_color(y, x)]);
	}

	/*
	 * This method uses bilinear interpolation.
	 */
	spixel get_pixel(unsigned int y, unsigned int x) const {
		pixel result;
		unsigned int k = bayer_color(y, x);
		ale_real sum;
		unsigned int num;

		result[k] = get_chan(y, x, k);

		if (k == 1) {
			unsigned int k1 = bayer_color(y + 1, x);
			unsigned int k2 = 2 - k1;

			sum = 0; num = 0;
			if (y > 0) {
				sum += get_chan(y - 1, x, k1);
				num++;
			}
			if (y < _dimy - 1) {
				sum += get_chan(y + 1, x, k1);
				num++;
			}
			assert (num > 0);
			result[k1] = sum / num;

			sum = 0; num = 0;
			if (x > 0) {
				sum += get_chan(y, x - 1, k2);
				num++;
			}
			if (x < _dimx - 1) {
				sum += get_chan(y, x + 1, k2);
				num++;
			}
			assert (num > 0);
			result[k2] = sum / num;

			return result;
		}

		sum = 0; num = 0;
		if (y > 0) {
			sum += get_chan(y - 1, x, 1);
			num++;
		}
		if (x > 0) {
			sum += get_chan(y, x - 1, 1);
			num++;
		}
		if (y < _dimy - 1) {
			sum += get_chan(y + 1, x, 1);
			num++;
		}
		if (x < _dimx - 1) {
			sum += get_chan(y, x + 1, 1);
			num++;
		}
		assert (num > 0);
		result[1] = sum / num;

		sum = 0; num = 0;
		if (y > 0 && x > 0) {
			sum += get_chan(y - 1, x - 1, 2 - k);
			num++;
		} 
		if (y > 0 && x < _dimx - 1) {
			sum += get_chan(y - 1, x + 1, 2 - k);
			num++;
		}
		if (y < _dimy - 1 && x > 0) {
			sum += get_chan(y + 1, x - 1, 2 - k);
			num++;
		}
		if (y < _dimy - 1 && x < _dimx - 1) {
			sum += get_chan(y + 1, x + 1, 2 - k);
			num++;
		}
		result[2 - k] = sum/num;

		return result;
	}

	spixel get_raw_pixel(unsigned int y, unsigned int x) const {
		pixel result;
		int k = bayer_color(y, x);

		result[k] = get_chan(y, x, k);

		return result;
	}

	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, const char *name) const {
		return new_image_ale_real(height, width, depth, name, _exp);
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.
	 */
	image *_extend(int top, int bottom, int left, int right) {
		/*
		 * Bayer-patterned images should always represent inputs,
		 * which should not ever be extended.
		 */
		assert(0);

		return NULL;
	}

private:
	void trigger(pixel multiplier) {
		for (unsigned int i = 0; i < _dimy; i++)
		for (unsigned int j = 0; j < _dimx; j++) {
			unsigned int k = bayer_color(i, j);
			set_chan(i, j, k, get_chan(i, j, k) * multiplier[k]); 
		}
	}

};

/*
 * Wrapper encapsulating details of the separation between the
 * resident-checking implementation and non-checking.
 */
static inline image *new_image_bayer_ale_real(unsigned int dimy,
                                 unsigned int dimx,
				 unsigned int depth,
				 unsigned int bayer,
				 const char *name = "anonymous",
				 exposure *exp = NULL) {

	return image_bayer_ale_real<0>::new_image_bayer_ale_real(dimy, dimx, depth, bayer, name, exp);
}

#endif
