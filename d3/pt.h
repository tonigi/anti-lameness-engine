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

#include "space.h"
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
	ale_real _view_angle;		/* XXX: should be ale_pos */
	ale_pos scale_factor;
	mutable ale_pos diag_per_depth;

public:	

	/*
	 * Constructor
	 */

	pt() {
		_view_angle = M_PI / 4;
		scale_factor = 1;
		diag_per_depth = 0;
	}

	pt(d2::transformation t, et e, ale_real va, ale_pos sf = 1) {
		this->t = t;
		euclidean = e;
		_view_angle = va;
		scale_factor = sf;
		diag_per_depth = 0;
	}

	/*
	 * Output function
	 */

	void debug_output() {
		t.debug_output();
		euclidean.debug_output();
		fprintf(stderr, "[pt.do va=%f sf=%f dpd=%f\n]",
				_view_angle, scale_factor, diag_per_depth);
	}

	/*
	 * Get euclidean transformation reference.
	 */

	et &e() {
		return euclidean;
	}

	/*
	 * Modify scale factor
	 */
	void scale(ale_pos sf) {
		scale_factor = sf;
	}

	/*
	 * Modify or get view angle
	 */
	void view_angle(ale_pos va) {
		diag_per_depth = 0;
		_view_angle = va;
	}

	ale_pos view_angle() {
		return _view_angle;
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

		ale_pos scaling_factor = sqrt(w*w + h*h) / (2 * tan(_view_angle / 2));
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

	int scaled_in_bounds(point p) const {
		return (p[0] >= 0 && p[0] <= scaled_height() - 1
		     && p[1] >= 0 && p[1] <= scaled_width() - 1);
	}

	int unscaled_in_bounds(point p) const {
		return (p[0] >= 0 && p[0] <= unscaled_height() - 1
		     && p[1] >= 0 && p[1] <= unscaled_width() - 1);
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

		ale_pos scaling_factor = sqrt(w*w + h*h) / (2 * tan(_view_angle / 2));
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

	/*
	 * Density calculation
	 */

	ale_pos c_density_scaled(point p) const {
		ale_pos one_density = 1 / (pc_scaled(point(0, 0, -1)).lengthto(pc_scaled(point(0, 1, -1)))
				         * pc_scaled(point(0, 0, -1)).lengthto(pc_scaled(point(1, 0, -1))));

		return one_density / (p[2] * p[2]);
	}

	ale_pos w_density_scaled(point p) const {
		return c_density_scaled(wc(p));
	}

	ale_pos w_density_scaled_max(point w0, point w1, point w2) {
		point c0 = wc(w0);
		point c1 = wc(w1);
		point c2 = wc(w2);

		/*
		 * Select the point closest to the camera.
		 */

		if (c0[2] > c1[2] && c0[2] > c2[2])
			return c_density_scaled(c0);
		else if (c1[2] > c2[2])
			return c_density_scaled(c1);
		else
			return c_density_scaled(c2);
	}

private:
	void calculate_diag_per_depth() const {
		if (diag_per_depth != 0)
			return;
		ale_pos w = unscaled_width();
		ale_pos h = unscaled_height();

		diag_per_depth = sqrt(2) * (2 * tan(_view_angle / 2)) / sqrt(w*w + h*h);
	}

public:

	/*
	 * Get a trilinear coordinate for a given depth.
	 */
	ale_pos trilinear_coordinate(ale_pos depth, ale_pos diagonal) {
		calculate_diag_per_depth();

		return log(diagonal / (fabs(depth) * diag_per_depth)) / log(2);

	}

	/*
	 * Get a trilinear coordinate for a given position in the world and
	 * a given 2D diagonal distance.
	 */
	ale_pos trilinear_coordinate(point w, ale_pos diagonal) {
		return trilinear_coordinate(wc(w)[2], diagonal);
	}

	/*
	 * Get a trilinear coordinate for a given subspace.
	 */
	ale_pos trilinear_coordinate(const space::traverse &st) {
		point min = st.get_min();
		point max = st.get_max();
		point avg = (min + max) / (ale_pos) 2;

		ale_pos diagonal = min.lengthto(max) * sqrt(2) / sqrt(3);

		point diff = (max - min);

		return trilinear_coordinate(avg, diagonal);
	}

	/*
	 * Get a diagonal distance for a given position in the world
	 * and a given trilinear coordinate.
	 */
	ale_pos diagonal_distance(point w, ale_pos coordinate) const {
		calculate_diag_per_depth();

		ale_pos depth = fabs(wc(w)[2]);
		ale_pos diagonal = pow(2, coordinate) * depth * diag_per_depth;

		return diagonal;
	}

	/*
	 * Get the 3D diagonal for a given depth and trilinear coordinate.
	 */
	ale_pos diagonal_distance_3d(ale_pos depth, ale_pos coordinate) const {
		calculate_diag_per_depth();
		return pow(2, coordinate) * fabs(depth) * diag_per_depth * sqrt(3) / sqrt(2);
	}

	/*
	 * Get the 1D distance for a given depth and trilinear coordinate.
	 */
	ale_pos distance_1d(ale_pos depth, ale_pos coordinate) const {
		calculate_diag_per_depth();
		return pow(2, coordinate) * fabs(depth) * diag_per_depth / sqrt(2);
	}

	ale_pos distance_1d(point iw, ale_pos coordinate) const {
		if (wc(iw)[2] >= 0)
			return point::undefined()[0];
		return distance_1d(-wc(iw)[2], coordinate);
	}

	/*
	 * Check for inclusion of a point in the bounding box of projected
	 * vertices.  This function returns non-zero when a point is included,
	 * when one of the vertices is infinite or undefined, or when a vertex
	 * is behind the point of projection.
	 */

	int check_inclusion(point volume_min, point volume_max, d2::point pc, int scaled) const {

		int test[2][2] = {{0, 0}, {0, 0}};

		point wbb[2] = { volume_min, volume_max };

		for (int x = 0; x < 2; x++)
		for (int y = 0; y < 2; y++)
		for (int z = 0; z < 2; z++) {

			point p = scaled ? wp_scaled(point(wbb[x][0], wbb[y][1], wbb[z][2]))
				         : wp_unscaled(point(wbb[x][0], wbb[y][1], wbb[z][2]));

			if (!p.finite())
				return 1;

			if (p[2] > 0)
				return 1;

			for (int d = 0; d < 2; d++) {
				if (p[d] <= pc[d])
					test[d][0] = 1;
				if (p[d] >= pc[d])
					test[d][1] = 1;
			}
		}

		for (int d = 0; d < 2; d++)
		for (int c = 0; c < 2; c++)
			if (test[d][c] == 0)
				return 0;

		return 1;
	}

	int check_inclusion_scaled(point volume_min, point volume_max, d2::point pc) const {
		return check_inclusion(volume_min, volume_max, pc, 1);
	}

	int check_inclusion_unscaled(point volume_min, point volume_max, d2::point pc) const {
		return check_inclusion(volume_min, volume_max, pc, 0);
	}

	int check_inclusion_unscaled(const space::traverse &st, d2::point pc) const {
		return check_inclusion_unscaled(st.get_min(), st.get_max(), pc);
	}

	int check_inclusion_scaled(const space::traverse &st, d2::point pc) const {
		return check_inclusion_scaled(st.get_min(), st.get_max(), pc);
	}



	/*
	 * Get bounding box for projection of a subspace.
	 */

	void get_view_local_bb(point volume_min, point volume_max, point bb[2], int scaled) const {

		point min = point::posinf();
		point max = point::neginf();

		point wbb[2] = { volume_min, volume_max };

		for (int x = 0; x < 2; x++)
		for (int y = 0; y < 2; y++)
		for (int z = 0; z < 2; z++) {
			point p = scaled ? wp_scaled(point(wbb[x][0], wbb[y][1], wbb[z][2]))
				         : wp_unscaled(point(wbb[x][0], wbb[y][1], wbb[z][2]));
			
			for (int d = 0; d < 3; d++) {
				if (p[d] < min[d])
					min[d] = p[d];
				if (p[d] > max[d])
					max[d] = p[d];
			}
		}

		/*
		 * Clip bounding box to image extents.
		 */

		if (min[0] < 0)
			min[0] = 0;
		if (min[1] < 0)
			min[1] = 0;
		if (max[0] > scaled_height() - 1)
			max[0] = scaled_height() - 1;
		if (max[1] > scaled_width() - 1)
			max[1] = scaled_width() - 1;

		bb[0] = min;
		bb[1] = max;
	}

	void get_view_local_bb_unscaled(point volume_min, point volume_max, point bb[2]) const {
		get_view_local_bb(volume_min, volume_max, bb, 0);
	}

	void get_view_local_bb_scaled(point volume_min, point volume_max, point bb[2]) const {
		get_view_local_bb(volume_min, volume_max, bb, 1);
	}

	void get_view_local_bb_scaled(const space::traverse &st, point bb[2]) const {
		get_view_local_bb_scaled(st.get_min(), st.get_max(), bb);
	}

	void get_view_local_bb_unscaled(const space::traverse &t, point bb[2]) const {
		get_view_local_bb_unscaled(t.get_min(), t.get_max(), bb);
	}

	/*
	 * Get the in-bounds centroid for a subspace, if one exists.
	 */

	point centroid(point volume_min, point volume_max) const {

		point bb[2];

		get_view_local_bb_unscaled(volume_min, volume_max, bb);

		point min = bb[0];
		point max = bb[1];

		for (int d = 0; d < 2; d++)
		if (min[d] > max[d])
			return point::undefined();

		if (min[2] >= 0)
			return point::undefined();

		if (max[2] > 0)
			max[2] = 0;

		return (max + min) / 2;
	}

	point centroid(const space::traverse &t) {
		return centroid(t.get_min(), t.get_max());
	}

	/*
	 * Get the local space origin in world space.
	 */

	point origin() {
		return cw(point(0, 0, 0));
	}

};

#endif
