// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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

#include "scene.h"

/*
 * See scene.h for details on these variables.
 */

scene::lod *scene::cl;

/*
 * Constant for depth tests.
 */

static double depth_quantum = 0.02;

/*
 * Function bodies for recursive cost adjustment.
 */

void scene::local_puts(char *s, unsigned int search_depth, unsigned int f, 
                       unsigned int x, unsigned int y, double value) {

	static int success_bit = 0;
	static int track_bit = 0;
	static unsigned int last_sum = 0;


	// if (success_bit == 0 && x != 201)
	//	return;

	if (search_depth == 0 && (x == 59 || x == 60)) 
		track_bit = 0;
	else if (search_depth == 0) {
		track_bit = 0;
	}

	// if (s != "success" && s != "early_success" && success_bit != 1)
	// 	return;

	if (success_bit != 1 && (s == "success" || s == "early_success"))
		success_bit = 1;

	if (search_depth == 0 && (s == "success" || s == "early_success"))
		success_bit = 0;

	if (!success_bit && !track_bit)
		return;

	if (search_depth + f + x + y != last_sum) {
		last_sum = search_depth + f + x + y;
	} else if (s != "max_struct.cost") {
		return;
	}

	return;

	if (search_depth > 2)
	     return;

	for (unsigned int i = 0; i < search_depth; i++)
		fprintf(stderr, " ");
	
	// fprintf(stderr, "[d=%d f=%d x=%d y=%d color=(%.10f, %.10f, %.10f) depth=%f %s=%f]\n",
	fprintf(stderr, "[d=%d f=%d x=%d y=%d %s=%f]\n",
			search_depth,
			f, x, y, 
//			cl->pz_color[f]->get_pixel(x, y)[0],
//			cl->pz_color[f]->get_pixel(x, y)[1],
//			cl->pz_color[f]->get_pixel(x, y)[2],
//			cl->pz_depth[f]->get_pixel(x, y)[0],
			s, value);

}

/*
 * Macro to abstract similar code chunks in functions try_color_change() and
 * try_depth_change().  See the function bodies for more information on what
 * the variables are.  A mapping between macro variable name and variable name
 * in the function body is given below:
 * 
 * 	bc --> branch_cost	Lowest cost previously found in this branch
 *	sc --> sub_cost		Lowest cost found on a searched sub-branch
 *	sd --> search_depth	Current search depth
 *	f  --> f		Current frame
 *	x  --> x		Current x position
 *	y  --> y		Current y position
 *	vd --> visit_data	Information about visited nodes
 *	
 *	is_color --> is_color	Which function body we are within
 */

#define D3_SCENE_CC_SUCCESS_MACRO(bc,sc,sd,f,x,y,vd,is_color) { \
	bc = sc; \
	local_puts("branch_cost", sd, f, x, y, bc); \
	if (bc < 0) { \
		local_puts("success", sd, f, x, y, (double) 1); \
		leave(vd, f, x, y, is_color); \
		return bc; \
	} }

/*
 * Try a color change at frame F, location (X, Y).  
 *
 * 	DELTA_COST is the current cost represented as an offset from the best
 * 	cost achieved thus far.
 *
 *	The return value is the best cost found on this branch, reported as an
 *	offset from the best cost prior to investigating this branch.  If this
 *	branch contains no visitable nodes, the return value is equal to
 *	DELTA_COST.
 * 
 * 	The best cost found on this branch will be the first cost that beats
 * 	the previous best cost, or, if no such cost is found, the lowest cost
 * 	found.
 */
double scene::try_color_change(unsigned int f, 
		unsigned int x, unsigned int y, 
		unsigned int search_depth, 
		unsigned int max_search_depth, 
		d2::pixel color, double delta_cost,
		void *visit_data) {

	local_puts("try_color_change_begin", search_depth, f, x, y, 1);

	local_puts("delta_cost", search_depth, f, x, y, delta_cost);

	if (!check_visitable(visit_data, f, x, y, 1))
		return delta_cost;

	visit(visit_data, f, x, y, 1);

	assert (x < cl->pz_color[f]->height());
	assert (y < cl->pz_color[f]->width());

	d2::pixel original_color = cl->pz_color[f]->get_pixel(x, y);

	cl->pz_color[f]->pix(x, y) = color;

	/*
	 * Establish the cost difference for the given color change w.r.t. the
	 * reference frame.
	 */

	double color_cost = frc_cost(f, x, y, color)
	                  - frc_cost(f, x, y, original_color);

	local_puts("color_cost", search_depth, f, x, y, color_cost);

	/*
	 * Establish the cost difference for all other frames.
	 */

	unsigned int max_index = 0;
	double max_cost = 0;
	double sum_cost = 0;

	for (unsigned int ff = 0; ff < d2::image_rw::count(); ff++) {
		double diff = 
			ftf_cost(ff, f, x, y, cl->pz_depth[f]->get_pixel(x, y)[0], color)
		      - ftf_cost(ff, f, x, y, cl->pz_depth[f]->get_pixel(x, y)[0], original_color);

		local_puts("color_change_ftf_cost_diff", search_depth, f, x, y, diff);

		if (diff > max_cost) {
			max_index = ff;
			max_cost = diff;
		}

		sum_cost += diff;
	}

	local_puts("sum_cost", search_depth, f, x, y, sum_cost);
	local_puts("max_cost", search_depth, f, x, y, max_cost);
	local_puts("max_index", search_depth, f, x, y, max_index);

	/*
	 * NODE_COST is the cost after making the desired color change.
	 * BRANCH_COST is the lowest cost of all configurations encountered so
	 * far on this branch.  Both variables are expressed as offsets from
	 * the best cost prior to investigating this branch.
	 */

	const double node_cost = delta_cost + sum_cost + color_cost;
	double branch_cost = node_cost;

	local_puts("node_cost", search_depth, f, x, y, node_cost);
	local_puts("branch_cost", search_depth, f, x, y, branch_cost);

	/*
	 * Early success check
	 */
	if (node_cost < 0) {
		local_puts("early_success", search_depth, f, x, y, 1);
		leave(visit_data, f, x, y, 1);
		return node_cost;
	}

	/*
	 * Maximum depth failure check
	 */

	if (search_depth == max_search_depth) {
		local_puts("max_depth_failure", search_depth, f, x, y, 1);
		cl->pz_color[f]->pix(x, y) = original_color;
		leave(visit_data, f, x, y, 1);
		return node_cost;
	}

	/*
	 * Attempt to achieve gains by reducing the cost of the most
	 * expensive frame.
	 */

	unsigned int x2, y2;
	unsigned int in_bound = 0;
	point mapped;

	mapped = frame_to_frame(max_index, f, x, y, cl->pz_depth[f]->get_pixel(x, y)[0]);

	x2 = (int) round(mapped[0]);
	y2 = (int) round(mapped[1]);

	local_puts("x2", search_depth, f, x, y, (double) x2);
	local_puts("y2", search_depth, f, x, y, (double) y2);

	if (x2 < cl->pz_depth[max_index]->height()
	 && y2 < cl->pz_depth[max_index]->width())
	 	in_bound = 1;

	local_puts("in_bound", search_depth, f, x, y, (double) in_bound);

	/*
	 * Variable to temporarily track a sub-branch's cost.
	 */

	double sub_cost = 0;

	/*
	 * Coincidence is the only possible cause of increased cost.  We can
	 * address coincidence either by changing color or by changing depth.
	 */
	
#if 0
	local_puts("try_color_change_for_coincidence", search_depth, f, x, y, (double) 1);

	if (in_bound && 
	    (sub_cost = try_color_change(max_index, x2, y2,
	                     search_depth + 1, max_search_depth,
			     color,
                             node_cost,
			     visit_data)) < branch_cost) 
		D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 1);
#endif

	local_puts("try_local_depth_change_for_coincidence", search_depth, f, x, y, (double) 1);

	if ((sub_cost = try_neighbor_depth_regions(depth_quantum, f, x, y,
			search_depth + 1, max_search_depth,
			// cl->pz_depth[f]->get_pixel(x, y)[0],
			node_cost, visit_data)) < branch_cost)
		D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 1);

	local_puts("try_remote_depth_change_for_coincidence", search_depth, f, x, y, (double) 1);

	if (in_bound &&
	    (sub_cost = try_neighbor_depth_regions(depth_quantum, max_index, x2, y2,
	    			search_depth + 1, max_search_depth,
				// cl->pz_depth[max_index]->get_pixel(x2, y2)[0],
				node_cost,
				visit_data)) < branch_cost) 
		D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 1);



	/*
	 * Failure
	 */

	assert (branch_cost >= 0);

	local_puts("failure_cost", search_depth, f, x, y, branch_cost);

	cl->pz_color[f]->pix(x, y) = original_color;
	leave(visit_data, f, x, y, 1);
	return branch_cost;
}

/*
 * Try a depth change at frame F, location (X, Y).  
 *
 * 	DELTA_COST is the current cost represented as an offset from the best
 * 	cost achieved thus far.
 *
 *	The return value is the best cost found on this branch, reported as an
 *	offset from the best cost prior to investigating this branch.  If this
 *	branch contains no visitable nodes, the return value is equal to
 *	DELTA_COST.
 * 
 * 	The best cost found on this branch will be the first cost that beats
 * 	the previous best cost, or, if no such cost is found, the lowest cost
 * 	found.
 */
double scene::try_depth_change(unsigned int f, 
		unsigned int x, unsigned int y, 
		unsigned int search_depth, 
		unsigned int max_search_depth, 
		double depth, double delta_cost,
		void *visit_data) {

	local_puts("try_depth_change_begin", search_depth, f, x, y, 1);

	if (!check_visitable(visit_data, f, x, y, 0))
		return delta_cost;

	visit(visit_data, f, x, y, 0);

	double original_depth = cl->pz_depth[f]->get_pixel(x, y)[0];

	cl->pz_depth[f]->pix(x, y)[0] = depth;

	local_puts("original_depth", search_depth, f, x, y, (double) original_depth);
	local_puts("depth", search_depth, f, x, y, (double) depth);

	/*
	 * Establish the cost difference for all other frames.
	 */

	unsigned int max_index = f;
	ftf_cost_diff_struct max_struct;
	double sum_cost = 0;

	max_struct.cost = 0;

	for (unsigned int ff = 0; ff < d2::image_rw::count(); ff++) {

		ftf_cost_diff_struct diff_struct =
			ftf_cost_diff(ff, f, x, y, depth, original_depth, 
					cl->pz_color[f]->get_pixel(x, y));

		if (diff_struct.cost > max_struct.cost) {
			max_index = ff;
			max_struct = diff_struct;
		}

		sum_cost += diff_struct.cost;
	}

	local_puts("sum_cost", search_depth, f, x, y, (double) sum_cost);
	local_puts("max_struct.cost", search_depth, f, x, y, (double) max_struct.cost);
	local_puts("max_index", search_depth, f, x, y, (double) max_index);

	/*
	 * NODE_COST is the cost after making the desired color change.
	 * BRANCH_COST is the lowest cost of all configurations encountered so
	 * far on this branch.  Both variables are expressed as offsets from
	 * the best cost prior to investigating this branch.
	 */

	const double node_cost = delta_cost + sum_cost;
	double branch_cost = node_cost;

	local_puts("node_cost", search_depth, f, x, y, node_cost);
	local_puts("branch_cost", search_depth, f, x, y, branch_cost);

	/*
	 * Early success check
	 */
	if (node_cost < 0) {
		local_puts("early_success", search_depth, f, x, y, 1);
		leave(visit_data, f, x, y, 0);
		return node_cost;
	}

	/*
	 * Maximum depth failure check
	 */

	if (search_depth == max_search_depth) {
		local_puts("max_depth_failure", search_depth, f, x, y, 1);
		cl->pz_depth[f]->pix(x, y)[0] = original_depth;
		leave(visit_data, f, x, y, 0);
		return node_cost;
	}

	/*
	 * Attempt to achieve gains by reducing the cost of the most
	 * expensive frame.
	 */

	unsigned int x2, y2;
	int x2o, y2o;
	unsigned int in_bound = 0;
	point mapped;

	mapped = frame_to_frame(max_index, f, x, y, depth);

	x2 = (int) round(mapped[0]);
	y2 = (int) round(mapped[1]);

	x2o = max_struct.xpos;
	y2o = max_struct.ypos;

	// fprintf(stderr, "[x2o=%d y2o=%d]", 
	//		x2o, y2o);

	local_puts("x2", search_depth, f, x, y, (double) x2);
	local_puts("y2", search_depth, f, x, y, (double) y2);

	if (x2 < cl->pz_depth[max_index]->height()
	 && y2 < cl->pz_depth[max_index]->width())
	 	in_bound = 1;

	local_puts("in_bound", search_depth, f, x, y, (double) in_bound);

	double sub_cost = 0;

	/*
	 * First, check for coincidence
	 */
	
	if (fabs(mapped[2] - depth) < obscuration_depth()) {
		local_puts("try_coincidence_local_color_change", search_depth, f, x, y, (double) 1);

#if 0
		if (in_bound
		 && (sub_cost = try_color_change(f, x, y,
				     search_depth + 1, max_search_depth,
				     cl->pz_color[max_index]->get_pixel(x2, y2),
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

		local_puts("try_coincidence_remote_color_change", search_depth, f, x, y, (double) 1);

		if (in_bound
		 && (sub_cost = try_color_change(max_index, x2, y2,
				     search_depth + 1, max_search_depth,
				     cl->pz_color[f]->get_pixel(x, y),
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

#endif
		local_puts("try_coincidence_remote_depth_change", search_depth, f, x, y, (double) 1);

		if (in_bound
		 && (sub_cost = try_neighbor_depth_regions(depth_quantum, max_index, x2, y2,
				     search_depth + 1, max_search_depth,
				     // cl->pz_depth[max_index]->get_pixel(x2, y2)[0], 
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

	}

	/*
	 * Next, assume obscuration
	 */

	if (depth > original_depth) {

		/*
		 * Case where new position is shallower.  Move the point 
		 * corresponding to the new position first.
		 */

		local_puts("try_remote_obscuration_remote_depth_change_to_coincidence_region", search_depth, f, x, y, (double) 1);

		if (in_bound
		 && (sub_cost = try_depth_region(depth_quantum, max_index, x2, y2,
				     search_depth + 1, max_search_depth,
				     mapped[2],
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

		local_puts("try_local_obscuration_remote_depth_change_within_region", search_depth, f, x, y, (double) 1);

		if (x2o >= 0
		 && (sub_cost = try_neighbor_depth_regions(depth_quantum, max_index, x2o, y2o,
				     search_depth + 1, max_search_depth,
				     // cl->pz_depth[max_index]->get_pixel(x2o, y2o)[0],
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

		local_puts("try_local_obscuration_remote_depth_change_to_local_depth_region", search_depth, f, x, y, (double) 1);

		if (x2o >= 0
		 && (sub_cost = try_depth_region(depth_quantum, max_index, x2o, y2o,
				     search_depth + 1, max_search_depth,
				     mapped[2],
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);
	} else {
		/*
		 * Case where new position is deeper.  Move the point corresponding
		 * to the original depth first.
		 */

		local_puts("try_local_obscuration_remote_depth_change_within_region", search_depth, f, x, y, (double) 1);

		if (x2o >= 0
		 && (sub_cost = try_neighbor_depth_regions(depth_quantum, max_index, x2o, y2o,
				     search_depth + 1, max_search_depth,
				     // cl->pz_depth[max_index]->get_pixel(x2o, y2o)[0],
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

		local_puts("try_local_obscuration_remote_depth_change_to_local_depth_region", search_depth, f, x, y, (double) 1);

		if (x2o >= 0
		 && (sub_cost = try_depth_region(depth_quantum, max_index, x2o, y2o,
				     search_depth + 1, max_search_depth,
				     mapped[2],
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

		local_puts("try_remote_obscuration_remote_depth_change_to_coincidence_region", search_depth, f, x, y, (double) 1);

		if (in_bound
		 && (sub_cost = try_depth_region(depth_quantum, max_index, x2, y2,
				     search_depth + 1, max_search_depth,
				     mapped[2],
				     node_cost,
				     visit_data)) < branch_cost)
			D3_SCENE_CC_SUCCESS_MACRO(branch_cost, sub_cost, search_depth, f, x, y, visit_data, 0);

	}

	/*
	 * Failure
	 */

	assert (branch_cost >= 0);

	local_puts("failure_cost", search_depth, f, x, y, branch_cost);

	cl->pz_depth[f]->pix(x, y)[0] = original_depth;
	leave(visit_data, f, x, y, 0);
	return branch_cost;
}

#undef D3_SCENE_CC_SUCCESS_MACRO
