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

#ifndef __drizzle_h__
#define __drizzle_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "gpt.h"
#include "image.h"
#include "point.h"

class drizzle : public render {
private:
	image *drizzle_image;
	image_weights *drizzle_weight;
	int step;
	int extend;
	double radius;
	double scale_factor;

	/*
	 * Increase extents of image and weight according to a new image to be
	 * drizzled.
	 */

	void increase_extents(transformation t) {

		assert (drizzle_weight != NULL);
		assert (drizzle_image  != NULL);

		double extend_offset_i = drizzle_image->offset()[0];
		double extend_offset_j = drizzle_image->offset()[1];

		assert (drizzle_weight->offset()[0] == extend_offset_i);
		assert (drizzle_weight->offset()[1] == extend_offset_j);

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
	 * specified transformation.
	 */
	void _drizzle(const image *delta, transformation t) {
		int i, j, k;

		assert (drizzle_image != NULL);
		assert (delta != NULL);
		assert (drizzle_weight != NULL);

		for (i = 0; (unsigned int) i < drizzle_image->height(); i++)
		for (j = 0; (unsigned int) j < drizzle_image->width();  j++) {
			int ii, jj;

			/*
			 * Transform
			 */

			point q;
			point offset = drizzle_image->offset();

			q = t.inverse_transform(
				point(i + offset[0], j + offset[1]));

			my_real ti = q[0];
			my_real tj = q[1];

			q = t.inverse_transform(
				point(i + offset[0] + 1, j + offset[1]));

			my_real ui = fabs(q[0] - ti);
			my_real uj = fabs(q[1] - tj);

			q = t.inverse_transform(
				point(i + offset[0], j + offset[1] + 1));

			my_real vi = fabs(q[0] - ti);
			my_real vj = fabs(q[1] - tj);

			/*
			 * We map the area of the drizzle_image pixel onto the delta
			 * frame as a rectangular area oriented on the delta
			 * frame's axes.  Note that this results in an area
			 * that may be the wrong shape or orientation.
			 *
			 * We define two estimates of the rectangle's
			 * dimensions below.  For rotations of 0, 90, 180, or
			 * 270 degrees, max and sum are identical.  For
			 * other orientations, sum is too large and max is
			 * too small.  We use the mean of max and sum, which we
			 * then divide by two to obtain the distance between
			 * the center and the edge.
			 */

			my_real maxi = (ui > vi) ? ui : vi;
			my_real maxj = (uj > vj) ? uj : vj;
			my_real sumi = ui + vi;
			my_real sumj = uj + vj;

			my_real di = (maxi + sumi) / 4;
			my_real dj = (maxj + sumj) / 4;

			ti /= scale_factor;
			tj /= scale_factor;
			di /= scale_factor;
			dj /= scale_factor;

			for (ii = (int) floor(ti - di - radius); 
				ii <= ceil(ti + di + radius); ii++)
			for (jj = (int) floor(tj - dj - radius); 
				jj <= ceil(tj + dj + radius); jj++) {
		
				my_real top = (ti - di < ii - radius)
					    ? (ii - radius)
					    : (ti - di);
				my_real bot = (ti + di > ii + radius)
					    ? (ii + radius)
					    : (ti + di);
				my_real lef = (tj - dj < jj - radius)
					    ? (jj - radius)
					    : (tj - dj);
				my_real rig = (tj + dj > jj + radius)
					    ? (jj + radius)
					    : (tj + dj);

				if (top < bot 
				 && lef < rig 
				 && ii >= 0
				 && ii < (int) delta->height()
				 && jj >= 0
				 && jj < (int) delta->width()) {

					double weight = drizzle_weight->get_pixel_component(i, j, 0);
					double thisw  = (bot - top) * (rig - lef) / (di * dj);

					drizzle_weight->set_pixel_component(i, j, 0, weight + thisw);

					for (k = 0; k < 3; k++)
						drizzle_image->set_pixel_component(i, j, k, 
							(int) ((weight * drizzle_image->get_pixel_component(
								i, j, k) + thisw
								 * delta->get_pixel_component(ii, jj, k))
							/ (weight + thisw)));
				}
			}
		}
	}

public:

	/*
	 * Constructor
	 */
	drizzle(int extend, double radius, double scale_factor) {
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

	virtual const image_weights *get_defined() {
		assert (step >= 0);
		return drizzle_weight;
	}

	/*
	 * Perform rendering steps requiring no frames beyond frame N.
	 */

	virtual void operator()(int n) {
		assert (step >= -1);
		for (int i = step + 1; i <= n; i++) {
			if (i == 0) {
				const image *im = image_rw::open(0);

				drizzle_image = image_rw::copy(0);
				drizzle_weight = new_image_weights(drizzle_image->height(),
						drizzle_image->width(), 1);

				_drizzle(im, transformation::eu_identity(im, scale_factor));

				image_rw::close(0);
			} else if (align::match(i)) {
				transformation t = align::of(i);
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
