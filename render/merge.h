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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../gpt.h"
#include "../image.h"
#include "../point.h"

template <int replace>
class merge : public render {
private:
	image *accum_image;
	image_weights *accum_weight;
	int step;
	int extend;
	double scale_factor;

	/*
	 * Increase extents of image and weight according to a new image to be
	 * merged.
	 */

	void increase_extents(transformation t) {

		assert (accum_weight != NULL);
		assert (accum_image  != NULL);

		double extend_offset_i = accum_image->offset()[0];
		double extend_offset_j = accum_image->offset()[1];

		assert (accum_weight->offset()[0] == extend_offset_i);
		assert (accum_weight->offset()[1] == extend_offset_j);

		int extend_top = 0;
		int extend_bottom = 0;
		int extend_left = 0;
		int extend_right = 0;

		my_real height = t.height(), width = t.width();

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
		int i, j, k;

		assert (accum_image != NULL);
		assert (delta != NULL);
		assert (accum_weight != NULL);

		for (i = 0; (unsigned int) i < accum_image->height(); i++)
			for (j = 0; (unsigned int) j < accum_image->width(); j++) {

				/*
				 * Transform
				 */

				struct point q;
				struct point offset = accum_image->offset();

				q = t.inverse_transform(
					point(i + offset[0], j + offset[1]));	

				my_real ti = q[0];
				my_real tj = q[1];

				/*
				 * Check that the transformed coordinates are within
				 * the boundaries of the delta frame.
				 */

				if (ti >= 0
				 && ti <= delta->height() * scale_factor - 1
				 && tj >= 0
				 && tj <= delta->width() * scale_factor - 1){ 

					/*
					 * Determine and update merging weight at this pixel
					 */

					double weight;

					if (replace) {
						weight = 0;
					} else {
						weight = accum_weight->get_pixel_component(i, j, 0);
					}

					accum_weight->set_pixel_component(i, j, 0, weight + 1);

					/*
					 * Update each channel
					 */
					
					for (k = 0; k < 3; k++)
						accum_image->set_pixel_component(i, j, k, (unsigned char)
							((weight * accum_image->get_pixel_component(i, j, k)
						     +   delta->get_scaled_bl_component(
							     ti, tj, k, scale_factor)) 
							 / (weight + 1)));
				}
			}
	}

public:

	/*
	 * Constructor
	 */
	merge(int extend, double scale_factor) {
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

	virtual const image_weights *get_defined() {
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

				accum_image = image_rw::copy(0);
				accum_image->scale(scale_factor);

				accum_weight = new_image_weights(accum_image->height(),
						accum_image->width(), 1);

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
