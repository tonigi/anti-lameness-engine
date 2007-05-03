// Copyright 2002, 2004, 2007 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                                          <dhilvert@ugcs.caltech.edu>

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
 * trans_multi.h: Represent multiple transformations, affecting different
 * regions of a scene.
 */

#ifndef __trans_multi_h__
#define __trans_multi_h__

#include "trans_abstract.h"
#include "trans_single.h"

struct trans_multi : public trans_abstract {
private:
	std::vector<trans_single> trans_stack;
	int full_support;
	unsigned int current_element;
	
public:	

	trans_multi() : trans_stack() {
		full_support = 0;
		current_element = 0;
	}

	trans_multi(const trans_multi &tm) : trans_abstract(*(trans_abstract *) &tm) {
		trans_stack = tm.trans_stack;
		full_support = tm.full_support;
		current_element = tm.current_element;
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

	void push_element() {
		assert (trans_stack.size() > 0);

		if (++current_element == trans_stack.size())
			trans_stack.push_back(trans_stack.back());
	}

	unsigned int get_current_index() const {
		return current_element;
	}

	void set_current_index(unsigned int i) {
		assert (i < trans_stack.size());

		current_element = i;
	}

	unsigned int stack_depth() {
		return trans_stack.size();
	}

	int supported(int i, int j) {
		if (full_support || current_element == 0)
			return 1;

		return 0;
	}
	
	void use_full_support() {
		full_support = 1;
	}

	void use_restricted_support() {
		full_support = 0;
	}

	int is_nontrivial() {
		return 0;
	}

	void begin_calculate_scaled_region(unsigned int i_max, unsigned int j_max, point offset) {
		assert(0);
	}

	void begin_calculate_unscaled_region(unsigned int i_max, unsigned int j_max, point offset) {
		assert(0);
	}

	void end_calculate_region() {
		assert(0);
	}

	point get_query_point(int i, int j) {
		return point::undefined();
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
	 * Calculate projective transformation parameters from a euclidean
	 * transformation.
	 */
	void eu_to_gpt() {
		for (unsigned int t = 0; t < trans_stack.size(); t++)
			trans_stack[t].eu_to_gpt();
	}

	/*
	 * Calculate euclidean identity transform for a given image.
	 */
	static struct trans_multi eu_identity(const image *i = NULL, ale_pos scale_factor = 1) {
		struct trans_multi r;
		r.input_width = i ? i->width() : 2;
		r.input_height = i ? i->height() : 2;
		r.scale_factor = scale_factor;
		r.trans_stack.push_back(trans_single::eu_identity(i, scale_factor));
		r.current_element = 0;
		return r;
	}

	/*
	 * Calculate projective identity transform for a given image.
	 */
	static trans_multi gpt_identity(const image *i, ale_pos scale_factor) {
		struct trans_multi r = eu_identity(i, scale_factor);
		r.eu_to_gpt();
		return r;
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

		if (trans_stack.size() > 1) {

			/*
			 * Rescale region map.
			 */

#warning d2::trans_multi::specific_rescale() does no region map rescaling.
		}

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
