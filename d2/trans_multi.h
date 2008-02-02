// Copyright 2002, 2004, 2007 David Hilvert <dhilvert@auricle.dyndns.org>, 
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
 * trans_multi.h: Represent multiple transformations, affecting different
 * regions of a scene.
 */

#ifndef __trans_multi_h__
#define __trans_multi_h__

#include "trans_abstract.h"
#include "trans_single.h"

struct trans_multi : public trans_abstract {
public:
	struct multi_coordinate {
		int degree;
		int x;
		int y;

	public:
		int operator<(const multi_coordinate &mc) const {
			if (degree < mc.degree
			 || degree == mc.degree && y < mc.y
			 || degree == mc.degree && y == mc.y && x < mc.x)
				return 1;

			return 0;
		}
	};

private:
	static ale_pos _multi_decomp;

	std::vector<trans_single> trans_stack;
	std::vector<multi_coordinate> coord_stack;
	std::map<multi_coordinate, unsigned int> coordinate_map;

	int use_multi;
	unsigned int current_element;

	unsigned int orig_ref_height, orig_ref_width;
	unsigned int cur_ref_height, cur_ref_width;
	point cur_offset;

	unsigned int *spatio_elem_map;
	unsigned int *spatio_elem_map_r;
	
	void push_element() {
		assert (trans_stack.size() > 0);

		if (++current_element == trans_stack.size())
			trans_stack.push_back(trans_stack.back());
	}

	trans_multi() : trans_stack() {
		use_multi = 0;
		current_element = 0;
		orig_ref_height = 0;
		orig_ref_width = 0;
		cur_ref_height = 0;
		cur_ref_width = 0;
		spatio_elem_map = NULL;
		spatio_elem_map_r = NULL;
	}

public:	

	static void set_md(double d) {
		_multi_decomp = d;
	}

	/*
	 * Calculate euclidean identity transform for a given image.
	 */
	static struct trans_multi eu_identity(const image *i = NULL, ale_pos scale_factor = 1) {
		struct trans_multi r;
		multi_coordinate mc;
		
		mc.degree = 0;
		mc.x = 0;
		mc.y = 0;

		r.input_width = i ? i->width() : 2;
		r.input_height = i ? i->height() : 2;
		r.scale_factor = scale_factor;
		r.trans_stack.push_back(trans_single::eu_identity(i, scale_factor));
		r.coord_stack.push_back(mc);
		r.coordinate_map[mc] = r.trans_stack.size() - 1;
		r.current_element = 0;
		return r;
	}

	/*
	 * Generate an array of identity transformations.
	 */
	static trans_multi *new_eu_identity_array(unsigned int size) {
		trans_multi *result = new trans_multi[size];
		for (unsigned int i = 0; i < size; i++)
			result[i] = eu_identity();
		return result;
	}

	/*
	 * Calculate projective transformation parameters from a euclidean
	 * transformation.
	 */
	void eu_to_gpt() {
		for (unsigned int t = 0; t < trans_stack.size(); t++)
			trans_stack[t].eu_to_gpt();
	}

	/*
	 * Calculate projective identity transform for a given image.
	 */
	static trans_multi gpt_identity(const image *i, ale_pos scale_factor) {
		struct trans_multi r = eu_identity(i, scale_factor);
		r.eu_to_gpt();
		return r;
	}

	trans_multi &operator=(const trans_multi &tm) {
		this->trans_abstract::operator=(*((trans_abstract *) &tm));

		trans_stack = tm.trans_stack;
		coord_stack = tm.coord_stack;
		coordinate_map = tm.coordinate_map;

		use_multi = tm.use_multi;
		current_element = tm.current_element;

		orig_ref_height = tm.orig_ref_height;
		orig_ref_width = tm.orig_ref_width;
		cur_ref_height = tm.cur_ref_height;
		cur_ref_width = tm.cur_ref_width;
		cur_offset = tm.cur_offset;

		size_t cur_size = cur_ref_width * cur_ref_height * sizeof(unsigned int);

		if (cur_size > 0) {
			spatio_elem_map = (unsigned int *) malloc(cur_size);
			assert (spatio_elem_map);
			memcpy(spatio_elem_map, tm.spatio_elem_map, cur_size);
		} else {
			spatio_elem_map = NULL;
		}

		cur_size = input_height * input_width * sizeof(unsigned int);
		if (cur_size > 0) {
			spatio_elem_map_r = (unsigned int *) malloc(cur_size);
			assert (spatio_elem_map_r);
			memcpy(spatio_elem_map_r, tm.spatio_elem_map_r, cur_size);
		} else {
			spatio_elem_map_r = NULL;
		}

		return *this;
	}

	trans_multi(const trans_multi &tm) : trans_stack() {
		operator=(tm);
	}

	~trans_multi() {
		free(spatio_elem_map);
		free(spatio_elem_map_r);
	}

	trans_single get_element(unsigned int index) {
		assert (index < trans_stack.size());

		return trans_stack[index];
	}

	trans_single get_current_element() {
		return get_element(current_element);
	}

	void set_element(unsigned int index, trans_single t) {
		assert (index < trans_stack.size());

		trans_stack[index] = t;
	}

	void set_current_element(trans_single t) {
		set_element(current_element, t);
	}

	unsigned int get_current_index() const {
		return current_element;
	}

	void set_current_index(unsigned int i) {
		assert (i < trans_stack.size());

		current_element = i;
	}

	/*
	 * Set the bounds of the reference image after incorporation
	 * of the original frame.
	 */
	void set_original_bounds(const image *i) {
		assert (orig_ref_width == 0);
		assert (orig_ref_height == 0);

		orig_ref_height = i->height();
		orig_ref_width = i->width();

		assert (orig_ref_width != 0);
		assert (orig_ref_height != 0);
	}

	static multi_coordinate parent_mc(multi_coordinate mc) {
		multi_coordinate result;

		assert (mc.degree > 0);

		if (mc.degree == 1) {
			result.degree = 0;
			result.x = 0;
			result.y = 0;

			return result;
		}

		result.degree = mc.degree - 1;
		result.x = (int) floor((double) mc.x / (double) 2);
		result.y = (int) floor((double) mc.y / (double) 2);

		return result;
	}

	/*
	 * Set the bounds of the reference image after incorporation
	 * of the most recent frame.
	 */
	void set_current_bounds(const image *i) {
		use_multi = 0;
		free(spatio_elem_map);
		free(spatio_elem_map_r);
		spatio_elem_map = NULL;
		spatio_elem_map_r = NULL;

		cur_ref_height = i->height();
		cur_ref_width = i->width();
		cur_offset = i->offset();

		for (int d = 1; d < _multi_decomp - 1; d++) {

			ale_pos height_scale = orig_ref_height / pow(2, d);
			ale_pos width_scale  = orig_ref_width  / pow(2, d);

			for (int i = floor(cur_offset[0] / height_scale);
				 i <= ceil((cur_offset[0] + cur_ref_height) / height_scale);
				 i++)
			for (int j = floor(cur_offset[1] / width_scale);
			         j <= ceil((cur_offset[1] + cur_ref_width) / width_scale);
				 j++) {

				 multi_coordinate c;
				 c.degree = d;
				 c.x = j;
				 c.y = i;

				 if (!coordinate_map.count(c)) {
				 	multi_coordinate parent = parent_mc(c);
					assert (coordinate_map.count(parent));
					trans_stack.push_back(trans_stack[coordinate_map[parent]]);
					coord_stack.push_back(c);
					coordinate_map[c] = trans_stack.size() - 1;
				}
			}
		}
	}

	unsigned int stack_depth() const {
		return trans_stack.size();
	}

	struct elem_bounds_t {
		int imin, imax, jmin, jmax;
	};

	elem_bounds_t elem_bounds() {
		elem_bounds_t result;

		result.imin = floor(cur_offset[0]);
		result.imax = ceil(cur_offset[0] + cur_ref_height);
		result.jmin = floor(cur_offset[1]);
		result.jmax = ceil(cur_offset[1] + cur_ref_width);

		if (current_element == 0)
			return result;

		multi_coordinate mc = coord_stack[current_element];
		
		ale_pos height_scale = orig_ref_height / pow(2, mc.degree);
		ale_pos width_scale  = orig_ref_width  / pow(2, mc.degree);

		if (height_scale * mc.y > result.imin)
			result.imin = height_scale * mc.y;
		if ((height_scale + 1) * mc.y < result.imax)
			result.imax = (height_scale + 1) * mc.y;
		if (width_scale * mc.x > result.jmin)
			result.jmin = width_scale * mc.x;
		if ((width_scale + 1) * mc.x < result.jmax)
			result.jmax = (width_scale + 1) * mc.x;

		return result;
	}

	void set_multi() {
		assert(use_multi == 0);
		assert(spatio_elem_map == NULL);
		assert(spatio_elem_map_r == NULL);
		use_multi = 1;
	}

	/*
	 * Returns non-zero if the transformation might be non-Euclidean.
	 */
	int is_projective() const {
		return trans_stack.front().is_projective();
	}

	/*
	 * Projective/Euclidean transformations
	 */
	struct point pe(struct point p) const {
		return trans_stack.front().pe(p);
	}

	/*
	 * Inverse transformations
	 */
	struct point pei(struct point p) const {
		return trans_stack.front().pei(p);
	}


	/*
	 * Modify a euclidean transform in the indicated manner.
	 */
	void eu_modify(int i1, ale_pos diff) {
		trans_stack[current_element].eu_modify(i1, diff);
	}

	/*
	 * Rotate about a given point in the original reference frame.
	 */
	void eu_rotate_about_scaled(point center, ale_pos diff) {
		trans_stack[current_element].eu_rotate_about_scaled(center, diff);
	}

	/*
	 * Modify all euclidean parameters at once.
	 */
	void eu_set(ale_pos eu[3]) {
		trans_stack[current_element].eu_set(eu);
	}

	/*
	 * Get the specified euclidean parameter
	 */
	ale_pos eu_get(int param) const {
		return trans_stack[current_element].eu_get(param);
	}

	/*
	 * Modify a projective transform in the indicated manner.
	 */
	void gpt_modify(int i1, int i2, ale_pos diff) {
		trans_stack[current_element].gpt_modify(i1, i2, diff);
	}

	/*
	 * Modify a projective transform according to the group operation.
	 */
	void gr_modify(int i1, int i2, ale_pos diff) {
		trans_stack[current_element].gr_modify(i1, i2, diff);
	}

	/*
	 * Modify all projective parameters at once.
	 */
	void gpt_set(point x[4]) {
		trans_stack[current_element].gpt_set(x);
	}

	void gpt_set(point x1, point x2, point x3, point x4) {
		trans_stack[current_element].gpt_set(x1, x2, x3, x4);
	}

	void snap(ale_pos interval) {
		trans_stack[current_element].snap(interval);
	}

	/*
	 * Get the specified projective parameter
	 */
	point gpt_get(int point) const {
		return trans_stack[current_element].gpt_get(point);
	}

	/*
	 * Get the specified projective parameter
	 */
	ale_pos gpt_get(int point, int dim) {
		return trans_stack[current_element].gpt_get(point, dim);
	}

	/*
	 * Translate by a given amount
	 */
	void translate(point p) {
		trans_stack[current_element].translate(p);
	}

	/*
	 * Rotate by a given amount about a given point.
	 */
	void rotate(point p, ale_pos degrees) {
		trans_stack[current_element].rotate(p, degrees);
	}

	void reset_memos() {
		for (unsigned int t = 0; t < trans_stack.size(); t++)
			trans_stack[t].reset_memos();
	}

	/*
	 * Rescale a transform with a given factor.
	 */
	void specific_rescale(ale_pos factor) {

		/*
		 * Ensure that no maps exist.
		 */

		assert (use_multi == 0);
		assert (spatio_elem_map == NULL);
		assert (spatio_elem_map_r == NULL);

		for (unsigned int t = 0; t < trans_stack.size(); t++)
			trans_stack[t].rescale(factor);
	}

	/*
	 * Set the dimensions of the image.
	 */
	void specific_set_dimensions(const image *im) {
		for (unsigned int t = 0; t < trans_stack.size(); t++)
			trans_stack[t].set_dimensions(im);
	} 

	/*
	 * Modify all projective parameters at once.  Accommodate bugs in the
	 * version 0 transformation file handler (ALE versions 0.4.0p1 and
	 * earlier).  This code is only called when using a transformation data
	 * file created with an old version of ALE.
	 */
	void gpt_v0_set(point x[4]) {
		trans_stack[current_element].gpt_v0_set(x);
	}

	/*
	 * Modify all euclidean parameters at once.  Accommodate bugs in the
	 * version 0 transformation file handler (ALE versions 0.4.0p1 and 
	 * earlier).  This code is only called when using a transformation data
	 * file created with an old version of ALE.
	 */
	void eu_v0_set(ale_pos eu[3]) {
		trans_stack[current_element].eu_v0_set(eu);
	}

	void debug_output() {
		for (unsigned int t = 0; t < trans_stack.size(); t++)
			trans_stack[t].debug_output();
	}
};

#endif
