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
 * drizzle.h: A renderer that implements a drizzling algorithm based on research
 * by Richard Hook and Andrew Fruchter.  For more details, see:
 *
 * 	http://www.cv.nrao.edu/adass/adassVI/hookr.html
 *
 * A template value indicates whether the drizzle weights should be altered
 * so that a replacement algorithm is used instead.
 */

#ifndef __drizzle_h__
#define __drizzle_h__

#include "../gpt.h"
#include "../image.h"
#include "../image_ale_real.h"
#include "../point.h"

template <int replace>
class drizzle : public render {
private:
	image *drizzle_image;
	image *drizzle_weight;
	int step;
	int extend;
	ale_pos radius;
	ale_pos scale_factor;

	/*
	 * Increase extents of image and weight according to a new image to be
	 * drizzled.
	 */

	void increase_extents(transformation t) {

		assert (drizzle_weight != NULL);
		assert (drizzle_image  != NULL);

		ale_pos extend_offset_i = drizzle_image->offset()[0];
		ale_pos extend_offset_j = drizzle_image->offset()[1];

		assert (drizzle_weight->offset()[0] == extend_offset_i);
		assert (drizzle_weight->offset()[1] == extend_offset_j);

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
			if (p[n][0] > drizzle_image->height() + extend_offset_i + extend_bottom) {
				extend_bottom = (int) ceil(p[n][0] - drizzle_image->height() -
						extend_offset_i);
			}
			if (p[n][1] < extend_offset_j - extend_left) {
				extend_left = (int) ceil(-p[n][1] + extend_offset_j);
			}
			if (p[n][1] > drizzle_image->width() + extend_offset_j + extend_right) {
				extend_right = (int) ceil(p[n][1] - drizzle_image->width() - 
						extend_offset_j);
			}
		}

		extend_offset_i -= extend_top;
		extend_offset_j -= extend_left;

		drizzle_image->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
		drizzle_weight->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
	}	

	/*
	 * Drizzle part of a delta frame onto part of a target frame using the
	 * specified transformation.  Following the conventions of the
	 * technical description, the target frame is also called the
	 * 'accumulated image' for the drizzle renderer.
	 */
	void _drizzle(const image *delta, transformation t) {

		assert (drizzle_image != NULL);
		assert (delta != NULL);
		assert (drizzle_weight != NULL);

		/*
		 * Iterate over all pixels in the accumulated image.
		 */

		for (unsigned int i = 0; i < drizzle_image->height(); i++)
		for (unsigned int j = 0; j < drizzle_image->width();  j++) {

			/*
			 * Obtain the position Q and dimensions D of
			 * accumulated image pixel (i, j) in the coordinate
			 * system of the delta frame.
			 */

                        point p = point(i + drizzle_image->offset()[0], j + drizzle_image->offset()[1]);
			point q;
			ale_pos d[2];

			t.unscaled_map_area_inverse(p, &q, d);

			/*
			 * Iterate over all delta frame pixels contributing to
			 * the accumulated image pixel (i, j).
			 */

			for (int ii = (int) floor(q[0] - d[0] - radius); 
				ii <= ceil(q[0] + d[0] + radius); ii++)
			for (int jj = (int) floor(q[1] - d[1] - radius); 
				jj <= ceil(q[1] + d[1] + radius); jj++) {

				/*
				 * Determine the boundaries of the region
				 * common to delta frame pixel (ii, jj) and
				 * accumulated image pixel (i, j).
				 */
		
				ale_pos top = (q[0] - d[0] < ii - radius)
					    ? (ii - radius)
					    : (q[0] - d[0]);
				ale_pos bot = (q[0] + d[0] > ii + radius)
					    ? (ii + radius)
					    : (q[0] + d[0]);
				ale_pos lef = (q[1] - d[1] < jj - radius)
					    ? (jj - radius)
					    : (q[1] - d[1]);
				ale_pos rig = (q[1] + d[1] > jj + radius)
					    ? (jj + radius)
					    : (q[1] + d[1]);

				if (top < bot 
				 && lef < rig 
				 && ii >= (int) 0
				 && ii < (int) delta->height()
				 && jj >= (int) 0
				 && jj < (int) delta->width()) {

					pixel value = delta->get_pixel(ii, jj);
					pixel weight = drizzle_weight->get_pixel(i, j);
					pixel thisw  = (bot - top) * (rig - lef) / (d[0] * d[1])
						     * delta->exp().confidence(value);

					drizzle_weight->pix(i, j) = weight + thisw;

					drizzle_image->set_pixel(i, j,
						(weight * drizzle_image->get_pixel(
							i, j) + thisw * value)
						/ (weight + thisw));

				}
			}
		}
	}

public:

	/*
	 * Constructor
	 */
	drizzle(int extend, ale_pos radius, ale_pos scale_factor) {
		step = -1;
		this->extend = extend;
		drizzle_weight = NULL;
		drizzle_image = NULL;
		this->radius = radius;
		this->scale_factor = scale_factor;
	}

	/*
	 * Result of rendering.
	 */

	virtual const image *get_image() {
		assert (step >= 0);
		return drizzle_image;
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() {
		assert (step >= 0);
		return drizzle_weight;
	}

	/*
	 * Perform rendering steps requiring no frames beyond frame N.
	 */

	virtual void sync(int n) {
		assert (step >= -1);
		for (int i = step + 1; i <= n; i++) {
			if (i == 0) {
				const image *im = image_rw::open(0);

				drizzle_image = new image_ale_real(
						(int) floor(im->height() * scale_factor), 
						(int) floor(im->width() * scale_factor), 3);

				drizzle_weight = new image_ale_real(drizzle_image->height(),
						drizzle_image->width(), 3);

				_drizzle(im, transformation::eu_identity(im, scale_factor));

				image_rw::close(0);
			} else if (align::match(i)) {
				transformation t = align::of(i);

				if  (replace)
				for (unsigned int i = 0; i < drizzle_weight->height(); i++)
				for (unsigned int j = 0; j < drizzle_weight->width(); j++)
					drizzle_weight->pix(i, j) *= (1 / (ale_real) 1000);

				if (extend)
					increase_extents(t);

				const image *im = image_rw::open(i);
				_drizzle(im, t);
				image_rw::close(i);
			}
			step = i;
		}
	}

};

#endif
