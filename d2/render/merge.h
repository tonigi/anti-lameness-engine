// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

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
 * merge.h: A render subclass implementing the default renderer used by ALE,
 * also used as the reference rendering for alignment.  A template value indicates
 * whether the merging weights should be altered so that a replacement algorithm is
 * used instead.
 *
 * Template value:
 * 	
 * 	0	Implement standard merging algorithm
 * 	1	Implement pixel replacement variant of the merging algorithm
 */

#ifndef __merge_h__
#define __merge_h__

#include "../gpt.h"
#include "../image.h"
#include "../point.h"

template <int replace>
class merge : public render {
private:
	image *accum_image;
	image *accum_weight;
	int step;
	int extend;
	ale_pos scale_factor;

	/*
	 * Increase extents of image and weight according to a new image to be
	 * merged.
	 */

	void increase_extents(transformation t) {

		assert (accum_weight != NULL);
		assert (accum_image  != NULL);

		ale_pos extend_offset_i = accum_image->offset()[0];
		ale_pos extend_offset_j = accum_image->offset()[1];

		assert (accum_weight->offset()[0] == extend_offset_i);
		assert (accum_weight->offset()[1] == extend_offset_j);

		int extend_top = 0;
		int extend_bottom = 0;
		int extend_left = 0;
		int extend_right = 0;

		ale_pos height = t.height(), width = t.width();

		point p[4];

		p[0] = t(point(   0  ,   0  ));
		p[1] = t(point(height,   0  ));
		p[2] = t(point(height, width));
		p[3] = t(point(   0  , width));

		for (int n = 0; n < 4; n++) {

			if (p[n][0] < extend_offset_i - extend_top) {
				extend_top = (int) ceil(-p[n][0] + extend_offset_i);
			}
			if (p[n][0] > accum_image->height() + extend_offset_i + extend_bottom) {
				extend_bottom = (int) ceil(p[n][0] - accum_image->height() -
						extend_offset_i);
			}
			if (p[n][1] < extend_offset_j - extend_left) {
				extend_left = (int) ceil(-p[n][1] + extend_offset_j);
			}
			if (p[n][1] > accum_image->width() + extend_offset_j + extend_right) {
				extend_right = (int) ceil(p[n][1] - accum_image->width() - 
						extend_offset_j);
			}
		}

		extend_offset_i -= extend_top;
		extend_offset_j -= extend_left;

		accum_image->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
		accum_weight->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
	}	

	/*
	 * Merge part of a delta frame with part of the accumulated image using
	 * the specified transformation.
	 */
	void
	_merge(const image *delta, transformation t) {

		assert (accum_image != NULL);
		assert (delta != NULL);
		assert (accum_weight != NULL);

		for (unsigned int i = 0; i < accum_image->height(); i++)
		for (unsigned int j = 0; j < accum_image->width(); j++) {

			/*
			 * Transform
			 */

			struct point q;
			struct point offset = accum_image->offset();

			q = t.inverse_transform(
				point(i + offset[0], j + offset[1]));	

			ale_pos ti = q[0];
			ale_pos tj = q[1];

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of the delta frame.
			 */

			if (ti >= 0
			 && ti <= delta->height() * scale_factor - 1
			 && tj >= 0
			 && tj <= delta->width() * scale_factor - 1){ 

				/*
				 * Pixel value to be merged, and the associated
				 * confidence
				 */

				pixel result[2];

				delta->get_scaled_bl(point(ti, tj), scale_factor, result);

				pixel value = result[0];
				pixel confidence = result[1];

				/*
				 * Determine and update merging weight at this pixel
				 */

				pixel old_weight, new_weight;

				if (replace) {
					old_weight = pixel(0, 0, 0);
				} else {
					old_weight = accum_weight->get_pixel(i, j);
				}

				new_weight = confidence;

				/*
				 * Update the weight
				 */

				accum_weight->set_pixel(i, j, old_weight + new_weight);

				/*
				 * Update the channel
				 */
				
				accum_image->set_pixel(i, j,
					(old_weight * accum_image->get_pixel(i, j)
				     +   new_weight * value)
				     / (old_weight + new_weight));
			}
		}
	}

public:

	/*
	 * Constructor
	 */
	merge(int extend, ale_pos scale_factor) {
		step = -1;
		this->extend = extend;
		this->scale_factor = scale_factor;
		accum_weight = NULL;
		accum_image = NULL;
	}

	/*
	 * Result of rendering.
	 */

	virtual const image *get_image() {
		assert (step >= 0);
		return accum_image;
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() {
		assert (step >= 0);
		return accum_weight;
	}

	/*
	 * Perform rendering steps requiring no frames beyond frame N.
	 */

	virtual void sync(int n) {
		assert (step >= -1);
		for (int i = step + 1; i <= n; i++) {
			if (i == 0) {
				const image *im = image_rw::open(0);

				accum_image = new image_ale_real((int) floor(im->height() * scale_factor),
						                 (int) floor(im->width() * scale_factor), 3);

				accum_weight = new image_ale_real(accum_image->height(),
						accum_image->width(), accum_image->depth());

				_merge(im, transformation::eu_identity(im, scale_factor));

				image_rw::close(0);
			} else if (align::match(i)) {
				transformation t = align::of(i);
				if (extend)
					increase_extents(t);
				const image *im = image_rw::open(i);
				_merge(im, t);
				image_rw::close(i);
			}
			step = i;
		}
	}
};

#endif
