// Copyright 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 * d3/et.h: Represent 3D->2D projective transformations.
 */

#ifndef __pt_h__
#define __pt_h__

#include "et.h"

/*
 * Structure to describe a 3D->2D projective transformation.  3D information is
 * preserved by adding a depth element to the result.
 */

struct pt {
private:
	d2::transformation t;
	et euclidean;
	ale_real view_angle;		/* XXX: should be ale_pos */
	ale_pos scale_factor;
	
public:	

	/*
	 * Constructor
	 */

	pt(d2::transformation t, et e, ale_real va, ale_pos sf = 1) {
		this->t = t;
		euclidean = e;
		view_angle = va;
		scale_factor = sf;
	}

	/*
	 * Modify scale factor
	 */
	void scale(ale_pos sf) {
		scale_factor = sf;
	}

	/*
	 * Get the 2D scale factor
	 */
	ale_pos scale_2d() {
		return t.scale();
	}

	/*
	 * Transform point p.
	 */
	struct point generic_transform(struct point p, ale_pos w, ale_pos h) {
		/*
		 * First, apply the euclidean transformation.
		 */

		p = euclidean(p);

		/*
		 * Divide x and y by negative z
		 */

		p[0] /= -p[2];
		p[1] /= -p[2];

		/*
		 * Scale x and y
		 */

		ale_pos scaling_factor = sqrt(w*w + h*h) / (2 * tan(view_angle / 2));
		p[0] *= scaling_factor;
		p[1] *= scaling_factor;

		/*
		 * Add an offset so that the upper-left corner is the origin.
		 */

		p[0] += h / 2;
		p[1] += w / 2;

		return p;
	}

	/*
	 * Width and height
	 */

	ale_pos scaled_width() {
		return t.scaled_width() * scale_factor;
	}

	ale_pos scaled_height() {
		return t.scaled_height() * scale_factor;
	}

	ale_pos unscaled_width() {
		return t.unscaled_width() * scale_factor;
	}

	ale_pos unscaled_height() {
		return t.unscaled_height() * scale_factor;
	}

	/*
	 * Scaled transform
	 */

	point scaled_transform(point p) {
		return generic_transform(p, t.scaled_width(), t.scaled_height());
	}

	/*
	 * Unscaled transform
	 */

	point unscaled_transform(point p) {
		return generic_transform(p, t.unscaled_width(), t.unscaled_height());
	}

	/*
	 * Map point p using the inverse of the transform.
	 */
	point inverse_transform_generic(point p, ale_pos w, ale_pos h) {
		/*
		 * Subtract offset
		 */

		p[0] -= h / 2;
		p[1] -= w / 2;

		/*
		 * Scale x and y
		 */

		ale_pos scaling_factor = sqrt(w*w + h*h) / (2 * tan(view_angle / 2));
		p[0] /= scaling_factor;
		p[1] /= scaling_factor;

		/*
		 * Multiply x and y by negative z
		 */

		p[0] *= -p[2];
		p[1] *= -p[2];

		/*
		 * Apply the inverse euclidean transformation.
		 */

		p = euclidean.inverse_transform(p);

		return p;
	}

	/*
	 * Inverse transform for scaled points.
	 */

	point inverse_transform_scaled(point p) {
		return inverse_transform_generic(p, t.scaled_width(), t.scaled_height());
	}

	/*
	 * Inverse transform for unscaled points.
	 */

	point inverse_transform_unscaled(point p) {
		return inverse_transform_generic(p, t.unscaled_width(), t.unscaled_height());
	}
};

#endif
