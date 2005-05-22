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
 *
 * The following coordinate systems are considered:
 *
 * 	P: projective
 * 	C: local cartesian
 * 	W: world
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
	ale_pos scale_2d() const {
		return t.scale();
	}

	/*
	 * Transform W to C.
	 */
	point wc(point p) const {
		return euclidean(p);
	}

	/*
	 * Transform C to P for given width and height.
	 */
	point cp_generic(point p, ale_pos w, ale_pos h) const {
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
	 * Transform point p.
	 */
	struct point wp_generic(struct point p, ale_pos w, ale_pos h) const {
		return cp_generic(wc(p), w, h);
	}

	/*
	 * Width and height
	 */

	ale_pos scaled_width() const {
		return t.scaled_width() * scale_factor;
	}

	ale_pos scaled_height() const {
		return t.scaled_height() * scale_factor;
	}

	ale_pos unscaled_width() const {
		return t.unscaled_width() * scale_factor;
	}

	ale_pos unscaled_height() const {
		return t.unscaled_height() * scale_factor;
	}

	/*
	 * Scaled transforms
	 */

	point cp_scaled(point p) const {
		return cp_generic(p, scaled_width(), scaled_height());
	}

	point wp_scaled(point p) const {
		return wp_generic(p, scaled_width(), scaled_height());
	}

	/*
	 * Unscaled transforms
	 */

	point cp_unscaled(point p) const {
		return cp_generic(p, unscaled_width(), unscaled_height());
	}

	point wp_unscaled(point p) const {
		return wp_generic(p, unscaled_width(), unscaled_height());
	}

	/*
	 * Transform P to C.
	 */
	point pc_generic(point p, ale_pos w, ale_pos h) const {
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

		return p;
	}

	/*
	 * Transform C to W
	 */
	point cw(point p) const {
		return euclidean.inverse_transform(p);
	}

	/*
	 * Transform P to W
	 */
	point pw_generic(point p, ale_pos w, ale_pos h) const {
		return cw(pc_generic(p, w, h));
	}

	/*
	 * Inverse transforms for scaled points.
	 */

	point pc_scaled(point p) const {
		return pc_generic(p, scaled_width(), scaled_height());
	}

	point pw_scaled(point p) const {
		return pw_generic(p, scaled_width(), scaled_height());
	}

	/*
	 * Inverse transforms for unscaled points.
	 */

	point pc_unscaled(point p) const {
		return pc_generic(p, unscaled_width(), unscaled_height());
	}

	point pw_unscaled(point p) const {
		return pw_generic(p, unscaled_width(), unscaled_height());
	}

};

#endif
