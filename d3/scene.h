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

/*
 * d3/scene.h: Representation of a 3D scene.
 */

#ifndef __scene_h__
#define __scene_h__

#include "point.h"

class scene {

	/*
	 * Constant for depth tests
	 */

	static const double depth_quantum;

	/*
	 * Structure to hold information for a given level of detail.
	 */
	struct lod {

		/*
		 * Color and depth portions of the partial z-buffer for
		 * each frame.
		 */
		d2::image **pz_color;
		d2::image **pz_depth;

		/*
		 * Reference image for each frame.
		 */

		d2::image **reference;

		/*
		 * Scale factor
		 */

		double sf;

		/*
		 * Next element
		 */

		struct lod *next;

	};

	/*
	 * Current level-of-detail
	 */

	static struct lod *cl;


	static void local_puts(char *s, unsigned int search_depth, unsigned int f,
			unsigned int x, unsigned int y, double value);

public:
	/*
	 * Initialize 3D scene from 2D scene, using 2D and 3D alignment
	 * information.
	 */
	static void init_from_d2() {

		/*
		 * Find out how many input frames there are.
		 */

		int N = d2::image_rw::count();

		/*
		 * Initialize the base level of detail
		 */

		cl = new lod;
		cl->pz_color = NULL;
		cl->pz_depth = NULL;
		cl->reference = NULL;
		cl->next = NULL;
		cl->sf = 1;

		/*
		 * Initialize reference images.
		 */

		cl->reference = (d2::image **) malloc(N * sizeof(d2::image *));

		assert(cl->reference);

		for (int n = 0; n < N; n++) {
			cl->reference[n] = d2::image_rw::copy(n, "3D scene reference");
			assert(cl->reference[n]);
		}

		/*
		 * Initialize partial z-buffer array
		 */

		cl->pz_color = (d2::image **) malloc(N * sizeof(d2::image *));
		cl->pz_depth = (d2::image **) malloc(N * sizeof(d2::image *));

		assert(cl->pz_color);
		assert(cl->pz_depth);

		for (int n = 0; n < N; n++) {

			/*
			 * Get d2 and d3 alignment for frame n
			 */

			d2::transformation t2 = d2::align::of(n);
			                et t3 =     align::of(n);

			/*
			 * Get (scaled) height and width for frame n
			 */

			unsigned int height = (int) ceil(t2.scaled_height());
			unsigned int width = (int) ceil(t2.scaled_width());

			/*
			 * Allocate partial z-buffer array for frame n
			 */

			cl->pz_color[n] = cl->reference[n];
			cl->pz_depth[n] = new d2::image_ale_real(height, width, 3);

			assert(cl->pz_color[n]);
			assert(cl->pz_depth[n]);

			/*
			 * For each pixel, assign depths based on the view
			 * transformation data, and assign colors based on the
			 * corresponding reference image.
			 */

			for (unsigned int i = 0; i < height; i++)
			for (unsigned int j = 0; j < width; j++) {

				/*
				 * Obtain the position Q2 and dimensions D of d2
				 * scaled frame pixel (i, j) in the coordinate
				 * system of the d2 rendered output.
				 */

				d2::point p(i, j);
				d2::point q2;
				ale_pos d[2];

				t2.map_area(p, &q2, d);

				/*
				 * Obtain the 3D position Q3 of d2 scaled frame
				 * pixel (i, j) from the perspective of the
				 * camera, and use the third coordinate to set
				 * the corresponding value in the partial
				 * z-buffer depth array.
				 */

				point q3 = t3(point(q2[0], q2[1], 0));
				cl->pz_depth[n]->pix(i, j)[0] = q3[2];
			}

//			fprintf(stderr, "(0, 0) depth (n=%d) %f\n", n, cl->pz_depth[n]->get_pixel(0, 0)[0]);
//			fprintf(stderr, "(%d, 0) depth (n=%d) %f\n", height - 1, n, cl->pz_depth[n]->get_pixel(height - 1, 0)[0]);
//			fprintf(stderr, "(0, %d) depth (n=%d) %f\n", width - 1, n, cl->pz_depth[n]->get_pixel(0, width - 1)[0]);
//			fprintf(stderr, "(%d, %d) depth (n=%d) %f\n", height - 1, width - 1, n, cl->pz_depth[n]->get_pixel(height - 1, width - 1)[0]);
		}
	}

	/*
	 * Initialize visit data
	 */
	static void *init_visit_data(unsigned int max_depth) {

		/*
		 * Data relating to the entry point:
		 *
		 * [ frame x-position y-position ]
		 *
		 * We initialize frame to be negative, since we don't know the
		 * entry point yet.  When we do know the entry point, we will
		 * initialize all three values within their valid ranges.
		 */

		static unsigned int storage[3]; 

		storage[0] = (unsigned int) -1;
		
		// return NULL;
		return (void *) storage;
	}

	/*
	 * Destroy visit data
	 */
	static void destroy_visit_data(void *data) {
	}

	/*
	 * Check whether a parameter is visitable.
	 */
	static int check_visitable(void *visit_data,
			unsigned int frame, unsigned int x, 
			unsigned int y, unsigned int is_color) {

		unsigned int *dp = (unsigned int *) visit_data;

		if (dp[0] == (unsigned int) -1) {
			dp[0] = frame;
			dp[1] = x;
			dp[2] = y;
		} else if (frame == dp[0]
			&& x     == dp[1]
			&& y     == dp[2]
			&& is_color == 1) {

			return 0;
		}

		return 1;
	}

	/*
	 * Record a visit to a parameter.
	 */
	static void visit(void *visit_data,
			unsigned int frame, unsigned int x, 
			unsigned int y, unsigned int is_color) {
		// assert (check_visitable(visit_data, frame, x, y, is_color));
	}

	/*
	 * Record leaving a parameter.
	 */
	static void leave(void *visit_data,
			unsigned int frame, unsigned int x,
			unsigned int y, unsigned int is_color) {

		unsigned int *dp = (unsigned int *) visit_data;

		if (frame == dp[0]
		 && x     == dp[1]
		 && y     == dp[2]
		 && is_color == 1) {

			dp[0] = (unsigned int) -1;
		}
	}

	/*
	 * Evaluate the cost incurred by a color difference
	 * between a frame state and the reference frame, for
	 * a given specified state.
	 */
	static double frc_cost(unsigned int f, 
			unsigned int x, unsigned int y,
			d2::pixel color) {

		assert (x < cl->reference[f]->height());
		assert (y < cl->reference[f]->width());

		d2::pixel reference_color = cl->reference[f]->get_pixel(x, y);

		d2::pixel diff = reference_color - color;

		return (diff * cl->reference[f]->exp().one_sided_confidence(reference_color, diff)).normsq() / 400;
	}

	/*
	 * Evaluate the cost incurred by a color difference
	 * between a frame state and the reference frame.
	 */
	static double frc_cost(unsigned int f, 
			unsigned int x, unsigned int y) {

		assert (x < cl->pz_color[f]->height());
		assert (y < cl->pz_color[f]->width());

		return frc_cost(f, x, y, cl->pz_color[f]->get_pixel(x, y));
	}

	/*
	 * Return an array of colors satisfying a cost bound
	 * for FRC cost.
	 */
	static void frc_samples(unsigned int f,
			unsigned int x, unsigned int y,
			double max_cost, unsigned int num,
			d2::pixel *colors) {

		assert(num == 1);

		assert (x < cl->reference[f]->height());
		assert (y < cl->reference[f]->width());

		colors[0] = cl->reference[f]->get_pixel(x, y);
	}

	/*
	 * Return a depth-specified pixel location in frame f2 for a
	 * depth-specified pixel x in frame f1.
	 */
	static point frame_to_frame(unsigned int f2,
			unsigned int f1, point x) {

		/*
		 * Scale values from current LOD to full LOD.
		 */

		double h1 = cl->pz_color[f1]->height() / cl->sf;
		double w1 = cl->pz_color[f1]->width() / cl->sf;
		double h2 = cl->pz_color[f2]->height() / cl->sf;
		double w2 = cl->pz_color[f2]->width() / cl->sf;

		x /= cl->sf;

		/*
		 * Calculate half of the width and depth values.
		 */

		double half_h1 = (h1 - 1) / (double) 2;
		double half_w1 = (w1 - 1) / (double) 2;
		double half_h2 = (h2 - 1) / (double) 2;
		double half_w2 = (w2 - 1) / (double) 2;

		double canonical_depth_1 = 
			sqrt(half_h1 * half_h1 + half_w1 * half_w1)
		      / tan(align::angle_of(f1) / 2);

		double canonical_depth_2 =
			sqrt(half_h2 * half_h2 + half_w2 * half_w2)
		      / tan(align::angle_of(f2) / 2);


		/*
		 * Adjust 0, 1 coordinates to be zero-centered
		 */

		x[0] -= half_h1;
		x[1] -= half_w1;

		/*
		 * Modify 0, 1 coordinates based on depth
		 */

		x[0] *= fabs(x[2] / canonical_depth_1);
		x[1] *= fabs(x[2] / canonical_depth_1);

		/*
		 * Transform to world coordinates
		 */

		x = align::of(f1).inverse_transform(x);

		/*
		 * Transform to target coordinates
		 */

		x = align::of(f2).transform(x);

		/*
		 * Adjust 0, 1 coordinates to a canonical depth
		 */

		x[0] *= fabs(canonical_depth_2 / x[2]);
		x[1] *= fabs(canonical_depth_2 / x[2]);

		/*
		 * Adjust 0, 1 coordinates to be corner-centered
		 */

		x[0] += half_h2;
		x[1] += half_w2;

		/*
		 * Scale to current LOD.
		 */

		x *= cl->sf;

		return x;
	}

	/*
	 * Return a depth-specified pixel location in frame f2 for a pixel 
	 * (x, y) in frame f1.
	 */
	static point frame_to_frame(unsigned int f2,
			unsigned int f1, unsigned int x,
			unsigned int y, double depth) {

		return frame_to_frame(f2, f1, point(x, y, depth));
	}

	static double obscuration_cost(double distance) {
		return pow(distance, 2) / 100;
	}

	static double obscuration_depth() {
		// return 0.25;
		return 2500;
	}

	/*
	 * Determine whether obscuration of local F1 element (X1, Y1) at a
	 * given DEPTH occurs due to remote frame F2's (X2, Y2) element.
	 * If so, return nonzero cost.
	 */
	static double ftf_local_obscured_cost(unsigned int f2, 
			unsigned int x2, unsigned int y2,
			unsigned int f1,
			unsigned int x1, unsigned int y1, 
			double depth) {

                if (f2 == f1)
                        return 0;
                                                                                                                               
                point mapped = frame_to_frame(f1, f2, x2, y2, cl->pz_depth[f2]->get_pixel(x2, y2)[0]);
                                                                                                                               
                unsigned int x2m = (unsigned int) floor(mapped[0]);
                unsigned int y2m = (unsigned int) floor(mapped[1]);

                if ((mapped[2] - depth) >= obscuration_depth()
		 && x2m == x1
		 && y2m == y1) {
			fprintf(stderr, "obscured.\n");
                        return obscuration_cost(mapped[2] - depth);
		}

		return 0;
	}


	/*
	 * Determine whether obscuration of an F2 element occurs for frame F1's
	 * (X, Y) element at a given depth.
	 */
	static int ftf_remote_obscured(unsigned int f2, unsigned int f1,
			unsigned int x, unsigned int y, 
			double depth) {

                if (f2 == f1)
                        return 0;
                                                                                                                               
                point mapped = frame_to_frame(f2, f1, x, y, depth);
                                                                                                                               
                unsigned int x2 = (unsigned int) floor(mapped[0]);
                unsigned int y2 = (unsigned int) floor(mapped[1]);

                if ((mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0]) >= obscuration_depth()) {
			fprintf(stderr, "obscured.\n");
                        return 1;
		}

		return 0;
	}


	/*
	 * Evaluate the cost incurred between two frames, for a given specified
	 * depth and color in frame F1.
	 *
	 * Note that there are a few components of cost between frames:
	 *
	 * (0) Out-of-bounds
	 * (1) Coincidence
	 * (2) Obscuration of an f2 element
	 * (3) Obscuration of an f1 element
	 *
	 * This method only accounts for (0), (1), and (2).  If (3) is required, use
	 * ftf_cost_diff() instead.
	 */
	static double ftf_cost(unsigned int f2, unsigned int f1, 
			unsigned int x, unsigned int y,
			double depth, d2::pixel color) {

		if (f2 == f1)
			return 0;

		point mapped = frame_to_frame(f2, f1, x, y, depth);

		unsigned int x2 = (unsigned int) floor(mapped[0]);
		unsigned int y2 = (unsigned int) floor(mapped[1]);

		/*
		 * Check for _really_ out of bounds
		 */

		if (depth >= 0 || mapped[2] >= 0)
			return 0;
			// return 1000000;

		/*
		 * Check for out-of-bounds
		 */

		if (mapped[0] < 0 || mapped[0] > cl->pz_color[f2]->height() - 1
		 || mapped[1] < 0 || mapped[1] > cl->pz_color[f2]->width() - 1)
			return 0;
			// return 10;

		/*
		 * Check for coincidence
		 */

//		fprintf(stderr, "color diff %f\n", fabs(cl->pz_color[f2]->get_pixel(x2, y2)[0] - color[0]));
//		fprintf(stderr, "depth diff %f\n", fabs(mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0]));
//		fprintf(stderr, "mapped depth %f\n", mapped[2]);
//		fprintf(stderr, "f2 depth %f\n", cl->pz_depth[f2]->get_pixel(x2, y2)[0]);

//		if (fabs(mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0]) < obscuration_depth()
//		 && (fabs(cl->pz_color[f2]->get_pixel(x2, y2)[0] - color[0]) > 0.01
//		  || fabs(cl->pz_color[f2]->get_pixel(x2, y2)[1] - color[1]) > 0.01 
//		  || fabs(cl->pz_color[f2]->get_pixel(x2, y2)[2] - color[2]) > 0.01)) 
//			return 50;

//		if (mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0] > -obscuration_depth()) {
		if (fabs(mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0]) < obscuration_depth()) {
			d2::pixel c2 = cl->pz_color[f2]->get_bl(d2::point(mapped[0], mapped[1]));
			d2::pixel diff = c2 - color;

			return diff.normsq();
		} else
			fprintf(stderr, "obscured.\n");

		/*
		 * Check for obscuration
		 */

		if ((mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0]) >= obscuration_depth()) {
//			d2::pixel c2 = cl->pz_color[f2]->get_bl(d2::point(mapped[0], mapped[1]));
//			d2::pixel diff = c2 - color;

			fprintf(stderr, "obscured.\n");
			return obscuration_cost(mapped[2] - cl->pz_depth[f2]->get_pixel(x2, y2)[0]);
			
			// return diff.normsq() + 50;
			// return 25;
			// return 10;
			// return 0;
		}

		return 0;
	}

	struct ftf_cost_diff_struct {
		double cost;
		int xpos;
		int ypos;

		ftf_cost_diff_struct() {
			xpos = -1;
			ypos = -1;
		}
	};

	/*
	 * Evaluate the difference in F2-relative cost for a given depth
	 * difference at (X, Y) in frame F1.  Returns a structure that
	 * indicates the most costly obscuring element from f2, since this
	 * is non-trivial to reconstruct.
	 */
	static ftf_cost_diff_struct ftf_cost_diff(unsigned int f2, unsigned int f1,
			unsigned int x, unsigned int y, double d1, double d2,
			d2::pixel color) {

		ftf_cost_diff_struct result;

		/*
		 * Evaluate cost difference due to the f1 element being obscured
		 * by f2 elements.
		 */

		double f1_obs_cost, f1_obs_cost_d1 = 0, f1_obs_cost_d2 = 0;
		d2::point maximum, minimum;

		maximum = minimum = frame_to_frame(f2, f1, x, y, d1).xy();

		for (int dwch = 0; dwch < 2; dwch++)
		for (int xoff = 0; xoff < 2; xoff++)
		for (int yoff = 0; yoff < 2; yoff++) {
			point preimage, image;

			preimage = point(
				x - 0.5 + xoff,
				y - 0.5 + yoff,
				dwch ? d1 : d2);

			image = frame_to_frame(f2, f1, preimage);

			if (image[0] < minimum[0])
				minimum[0] = image[0];
			if (image[1] < minimum[1])
				minimum[1] = image[1];
			if (image[0] > maximum[0])
				maximum[0] = image[0];
			if (image[1] > maximum[1])
				maximum[1] = image[1];
		}

		if (minimum[0] < 0)
			minimum[0] = 0;
		if (maximum[0] > cl->pz_color[f2]->height() - 1)
			maximum[0] = cl->pz_color[f2]->height() - 1;
		if (minimum[1] < 0)
			minimum[1] = 0;
		if (maximum[1] > cl->pz_color[f2]->width() - 1)
			maximum[1] = cl->pz_color[f2]->width() - 1;

		double greatest_cost = 0;

		for (int xx = (int) floor(minimum[0]); xx <= (int) ceil(maximum[0]); xx++)
		for (int yy = (int) floor(minimum[1]); yy <= (int) ceil(maximum[1]); yy++) {
			double temp_d1 = ftf_local_obscured_cost(f2, xx, yy, f1, x, y, d1);
			double temp_d2 = ftf_local_obscured_cost(f2, xx, yy, f1, x, y, d2);

			f1_obs_cost_d1 += temp_d1;
			f1_obs_cost_d2 += temp_d2;

			if ((temp_d2 - temp_d1) > greatest_cost) {
				greatest_cost = temp_d2 - temp_d1;
				result.xpos = xx;
				result.ypos = yy;
			}

//			if (temp_d1 != temp_d2) {
//				fprintf(stderr, "[lo f1=%d x=%d y=%d d1=%f d2=%f f2=%d x2=%d y2=%d]",
//						f1, x, y, d1, d2, f2, xx, yy);
//			}
		}

		f1_obs_cost = f1_obs_cost_d2 - f1_obs_cost_d1;

		/*
		 * Evaluate all other sources of cost.
		 */

		result.cost = f1_obs_cost
		            + ftf_cost(f2, f1, x, y, d2, color)
		            - ftf_cost(f2, f1, x, y, d1, color);

		return result;
	}

	/*
	 * Try a color change at frame F, location (X, Y).  
	 */
	static double try_color_change(unsigned int f, 
			unsigned int x, unsigned int y, 
			unsigned int search_depth, 
			unsigned int max_search_depth, 
			d2::pixel color, double delta_cost,
			void *visit_data);

	/*
	 * Try a depth change at frame F, location (X, Y).  
	 */
	static double try_depth_change(unsigned int f, 
			unsigned int x, unsigned int y, 
			unsigned int search_depth, 
			unsigned int max_search_depth, 
			double depth, double delta_cost,
			void *visit_data);

	static double try_depth_half_region(double depth_change_increment,
			unsigned int f, 
			unsigned int x, unsigned int y, 
			unsigned int search_depth, 
			unsigned int max_search_depth, 
			double depth, double node_cost,
			void *visit_data) {

		assert (search_depth <= max_search_depth);

		double branch_cost = try_depth_change(f, x, y,
				     search_depth, max_search_depth,
				     depth, node_cost,
				     visit_data);

		if (branch_cost < 0)
			return branch_cost;

		depth += depth_change_increment;

		double sub_cost;
		while ((sub_cost = try_depth_change(f, x, y,
				     search_depth, max_search_depth,
				     depth, node_cost,
				     visit_data)) < branch_cost) {
			branch_cost = sub_cost;
			if (branch_cost < 0)
				return branch_cost;
			depth += depth_change_increment;
		}

		return branch_cost;
	}

	static double try_depth_region(double depth_change_increment,
			unsigned int f, 
			unsigned int x, unsigned int y, 
			unsigned int search_depth, 
			unsigned int max_search_depth, 
			double depth, double node_cost,
			void *visit_data) {

		assert (search_depth <= max_search_depth);

		double branch_cost_pos, branch_cost_neg;
		double depth_pos, depth_neg;

		depth_pos = depth_neg = depth;

		branch_cost_pos = branch_cost_neg 
			      = try_depth_change(f, x, y,
				     search_depth, max_search_depth,
				     depth, node_cost,
				     visit_data);

		if (branch_cost_pos < 0)
			return branch_cost_pos;

		depth_pos += depth_change_increment;
		depth_neg -= depth_change_increment;

		double sub_cost;

		while ((sub_cost = try_depth_change(f, x, y,
				     search_depth, max_search_depth,
				     depth_pos, node_cost,
				     visit_data)) < branch_cost_pos) {
			branch_cost_pos = sub_cost;
			if (branch_cost_pos < 0)
				return branch_cost_pos;
			depth_pos += depth_change_increment;
		}

		while ((sub_cost = try_depth_change(f, x, y,
				     search_depth, max_search_depth,
				     depth_neg, node_cost,
				     visit_data)) < branch_cost_neg) {
			branch_cost_neg = sub_cost;
			if (branch_cost_neg < 0)
				return branch_cost_neg;
			depth_neg += depth_change_increment;
		}

		if (branch_cost_neg < branch_cost_pos)
			return branch_cost_neg;

		return branch_cost_pos;
	}

	static double try_neighbor_depth_regions(double depth_change_increment,
			unsigned int f, 
			unsigned int x, unsigned int y, 
			unsigned int search_depth, 
			unsigned int max_search_depth, 
			double node_cost,
			void *visit_data) {
		double branch_cost = 0;
		int uninitialized = 1;

		assert (x >= 0);
		assert (x < cl->pz_color[f]->height());
		assert (y >= 0);
		assert (y < cl->pz_color[f]->width());

		for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++) {
			if (x + i < 0
			 || y + j < 0
			 || x + i >= cl->pz_color[f]->height()
			 || y + j >= cl->pz_color[f]->width())
				continue;

			double sub_cost = try_depth_region(depth_change_increment,
					f, x, y,
					search_depth, max_search_depth,
					cl->pz_depth[f]->get_pixel(x + i, y + j)[0],
					node_cost,
					visit_data);

			if (uninitialized == 1
			 || sub_cost < branch_cost) {
				branch_cost = sub_cost;
				uninitialized = 0;
			}

			if (branch_cost < 0)
				return branch_cost;
		}

		return branch_cost;
	}

	/*
	 * Reduce the level of detail.  Return 0 when no further reduction
	 * is possible.
	 */
	static int reduce_lod() {

		int result = 1;

		/*
		 * Create a new structure for the reduced LOD.
		 */
		struct lod *nl = new lod;

		/*
		 * Find out how many input frames there are.
		 */

		int N = d2::image_rw::count();

		/*
		 * Initialize reference images and partial z-buffer arrays.
		 */

		nl->reference = (d2::image **) malloc(N * sizeof(d2::image *));
		nl->pz_color  = (d2::image **) malloc(N * sizeof(d2::image *));
		nl->pz_depth  = (d2::image **) malloc(N * sizeof(d2::image *));

		assert(nl->reference);
		assert(nl->pz_color);
		assert(nl->pz_depth);

		for (int n = 0; n < N; n++) {
			nl->reference[n] = cl->reference[n]->scale_by_half("3D, reduced LOD");
			nl->pz_color[n]  = nl->reference[n];
			nl->pz_depth[n]  = cl->pz_depth[n]->scale_by_half("3D, reduced LOD");

			for (unsigned int i = 0; i < nl->pz_depth[n]->height(); i++)
			for (unsigned int j = 0; j < nl->pz_depth[n]->width() ; j++) {
				nl->pz_depth[n]->pix(i, j) /= 2;
			}

			assert(nl->reference[n]);
			assert(nl->pz_color[n]);
			assert(nl->pz_depth[n]);

			if (nl->reference[n]->height() < 4
			 || nl->reference[n]->width () < 4
			 || nl->pz_color [n]->height() < 4
			 || nl->pz_color [n]->width () < 4
			 || nl->pz_depth [n]->height() < 4
			 || nl->pz_depth [n]->width () < 4)
				result = 0;

			if (nl->reference[n]->height() < 2
			 || nl->reference[n]->width () < 2
			 || nl->pz_color [n]->height() < 2
			 || nl->pz_color [n]->width () < 2
			 || nl->pz_depth [n]->height() < 2
			 || nl->pz_depth [n]->width () < 2)
				assert(0);
		}

		nl->sf = cl->sf * 0.5;

		nl->next = cl;

		cl = nl;

		return result;
	}

	/*
	 * Increase the level of detail.  
	 */
	static void increase_lod() {
		/*
		 * Pointer to the next higher LOD.
		 */
		struct lod *nl = cl->next;

		assert (nl != NULL);

		/*
		 * Find out how many input frames there are.
		 */

		int N = d2::image_rw::count();

		/*
		 * Update depth information
		 */

		for (int n = 0; n < N; n++) 
		for (unsigned int i = 0; i < nl->pz_depth[n]->height(); i++)
		for (unsigned int j = 0; j < nl->pz_depth[n]->width(); j++)
			nl->pz_depth[n]->pix(i, j)
				= cl->pz_depth[n]->get_scaled_bl(d2::point(i, j), 2) * d2::pixel(2, 2, 2);

		/*
		 * Delete the current LOD.
		 */

		for (int n = 0; n < N; n++) {
			delete cl->pz_depth[n];
			delete cl->reference[n];
			// delete cl->pz_color[n];
		}

		delete cl;

		cl = nl;
	}

	static const d2::image *view(unsigned int n) {

		const double depth_diff_limit = 0.1;

		assert (n < d2::image_rw::count());

		d2::image *result = cl->pz_color[n]->clone("3D view result");
		
		for (unsigned int x = 0; x < cl->pz_color[n]->height(); x++)
		for (unsigned int y = 0; y < cl->pz_color[n]->width();  y++) {
			double weight = 1;
			double d0 = cl->pz_depth[n]->get_pixel(x, y)[0];

			for (unsigned int m = 0; m < d2::image_rw::count(); m++) {
				if (m == n)
					continue;

				point p0 = frame_to_frame(m, n, x, y, d0);

				if (!cl->pz_color[m]->in_bounds(d2::point(p0[0], p0[1])))
					continue;

				d2::pixel pix_0 = cl->pz_color[m]->get_bl(d2::point(p0[0], p0[1]));
				double s0 = cl->pz_depth[m]->get_bl(d2::point(p0[0], p0[1]))[0];

				if (fabs(p0[2] - s0) < depth_diff_limit) {
					result->pix(x, y) += pix_0;
					weight += 1;
				}
			}

			result->pix(x, y) /= weight;
		}

		/*
		 * XXX: this should be freed somewhere, perhaps?
		 */

		return result;
	}

	static const d2::image *depth(unsigned int n) {
		assert (n < d2::image_rw::count());

		return cl->pz_depth[n];
	}

	/*
	 * Evaluate the interframe cost at for a given location.
	 *
	 * m is remote, n is local
	 *
	 * Set est_max to estimate maximum cost.
	 */
	static double evaluate_cost(int m, int n, int x, int y, int est_max, double modified_depth = 0) {
#if 0

		if (m == n)
			return 0;

		if (!cl->pz_depth[n]->in_bounds(d2::point(x, y)))
			return 0;
		
		double depth = (modified_depth == 0)
			     ? cl->pz_depth[n]->get_pixel(x, y)[0]
			     : modified_depth;

		point remote_p = frame_to_frame(m, n, x, y, depth);

		if (!cl->pz_color[m]->in_bounds(d2::point(remote_p[0], remote_p[1])))
			return 0;

		double r_source = cl->pz_color[n]->get_pixel(x, y)[0];
		double r0 = cl->pz_color[m]->get_bl(d2::point(remote_p[0], remote_p[1]))[0];
		double s0 = cl->pz_depth[m]->get_bl(d2::point(remote_p[0], remote_p[1]))[0];
		point pp1 = frame_to_frame(m, n, x, y, dp1);
		point pm1 = frame_to_frame(m, n, x, y, dm1);

		while (sqrt(pow(pm1[0] - remote_p[0], 2) + pow(pm1[1] - remote_p[1], 2)) < lateral_min
		    && sqrt(pow(pp1[0] - remote_p[0], 2) + pow(pp1[1] - remote_p[1], 2)) < lateral_min
		    && dm1 > depth_ceiling / 2
		    && depth + depth_quantum * multiplier * 2 < 0) {
			multiplier *= 2;
			dp1 = depth + depth_quantum * multiplier;
			dm1 = depth - depth_quantum * multiplier;
			pp1 = frame_to_frame(m, n, x, y, dp1);
			pm1 = frame_to_frame(m, n, x, y, dm1);
		}

		if (dp1 >= 0 || dm1 >= 0) {
			// fprintf(stderr, "Neg.\n");
			continue;
		}

		if (!cl->pz_color[m]->in_bounds(d2::point(pp1[0], pp1[1]))
		 || !cl->pz_color[m]->in_bounds(d2::point(pm1[0], pm1[1])))
			continue;

		err += pow(r_source - r0, 2)
		     + ((remote_p[2] < s0) 
		      ? match_weight_deep
		      : match_weight_shallow)
		     * pow(remote_p[2] - s0, 2);

		double rp1 = cl->pz_color[m]->get_bl(d2::point(pp1[0], pp1[1]))[0];
		double sp1 = cl->pz_depth[m]->get_bl(d2::point(pp1[0], pp1[1]))[0];

		p1_err += pow(r_source - rp1, 2) 
			+ ((pp1[2] < sp1)
			 ? match_weight_deep
			 : match_weight_shallow)
			* pow(pp1[2] - sp1, 2);


		double rm1 = cl->pz_color[m]->get_bl(d2::point(pm1[0], pm1[1]))[0];
		double sm1 = cl->pz_depth[m]->get_bl(d2::point(pm1[0], pm1[1]))[0];

		m1_err += pow(r_source - rm1, 2)
			+ ((pm1[2] < sm1)
			 ? match_weight_deep
			 : match_weight_shallow)
			* pow(pm1[2] - sm1, 2);
#endif
		return 0;
	}

	static void evaluate_total_cost() {
#if 0
		for (unsigned int n = 0; n < d2::image_rw::count(); n++)
		for (unsigned int x = 0; x < cl->pz_color[n]->height(); x++)
		for (unsigned int y = 0; y < cl->pz_color[n]->width();  y++) {
#if 1
			double multiplier = 1;
			double lateral_min = 0.5;
			double depth_ceiling = -pow(10, 5);
			double d0 = cl->pz_depth[n]->get_pixel(x, y)[0];
			double dp1 = d0 + depth_quantum;
			double dm1 = d0 - depth_quantum;
			double err = 0;
			double m1_err = 0;
			double p1_err = 0;

			while (dp1 == d0 || dm1 == d0
			    && dp1 < 0
			    && dm1 < 0) {
				// fprintf(stderr, "Multiplying...\n");
				multiplier *= 2;
				dp1 = d0 + depth_quantum * multiplier;
				dm1 = d0 - depth_quantum * multiplier;
			}

			if (d0 >= 0)
				continue;
			if (dp1 >= 0)
				dp1 = 0 - depth_quantum;
			if (dm1 >= 0)
				dm1 = 0 - depth_quantum;

			for (unsigned int m = 0; m < d2::image_rw::count(); m++) {
				if (m == n)
					continue;

				point p0 = frame_to_frame(m, n, x, y, d0);

				if (!cl->pz_color[m]->in_bounds(d2::point(p0[0], p0[1])))
					continue;

				double r_source = cl->pz_color[n]->get_pixel(x, y)[0];
				double r0 = cl->pz_color[m]->get_bl(d2::point(p0[0], p0[1]))[0];
				double s0 = cl->pz_depth[m]->get_bl(d2::point(p0[0], p0[1]))[0];
				point pp1 = frame_to_frame(m, n, x, y, dp1);
				point pm1 = frame_to_frame(m, n, x, y, dm1);

				while (sqrt(pow(pm1[0] - p0[0], 2) + pow(pm1[1] - p0[1], 2)) < lateral_min
				    && sqrt(pow(pp1[0] - p0[0], 2) + pow(pp1[1] - p0[1], 2)) < lateral_min
				    && dm1 > depth_ceiling / 2
				    && d0 + depth_quantum * multiplier * 2 < 0) {
					multiplier *= 2;
					dp1 = d0 + depth_quantum * multiplier;
					dm1 = d0 - depth_quantum * multiplier;
					pp1 = frame_to_frame(m, n, x, y, dp1);
					pm1 = frame_to_frame(m, n, x, y, dm1);
				}

				if (dp1 >= 0 || dm1 >= 0) {
					// fprintf(stderr, "Neg.\n");
					continue;
				}

				if (!cl->pz_color[m]->in_bounds(d2::point(pp1[0], pp1[1]))
				 || !cl->pz_color[m]->in_bounds(d2::point(pm1[0], pm1[1])))
					continue;

				err += pow(r_source - r0, 2)
				     + ((p0[2] < s0) 
				      ? match_weight_deep
				      : match_weight_shallow)
				     * pow(p0[2] - s0, 2);

				double rp1 = cl->pz_color[m]->get_bl(d2::point(pp1[0], pp1[1]))[0];
				double sp1 = cl->pz_depth[m]->get_bl(d2::point(pp1[0], pp1[1]))[0];

				p1_err += pow(r_source - rp1, 2) 
					+ ((pp1[2] < sp1)
					 ? match_weight_deep
					 : match_weight_shallow)
					* pow(pp1[2] - sp1, 2);


				double rm1 = cl->pz_color[m]->get_bl(d2::point(pm1[0], pm1[1]))[0];
				double sm1 = cl->pz_depth[m]->get_bl(d2::point(pm1[0], pm1[1]))[0];

				m1_err += pow(r_source - rm1, 2)
					+ ((pm1[2] < sm1)
					 ? match_weight_deep
					 : match_weight_shallow)
					* pow(pm1[2] - sm1, 2);
			}

			if (m1_err < p1_err && m1_err < err) {
				cl->pz_depth[n]->pix(x, y)[0] = dm1;
				improved = 1;
			} else if (p1_err < err) {
				cl->pz_depth[n]->pix(x, y)[0] = dp1;
				improved = 1;
			}
#else


			d2::pixel options[1];
			double cur_cost = frc_cost(n, x, y);

			// if (cur_cost < 0.0001)
			//	continue;

			frc_samples(n, x, y, cur_cost/2, 1, options);

			if (try_color_change(n, x, y, 0, max_depth, options[0], cur_cost/2, visit_data) < 0) {

				assert (options[0][0] == cl->pz_color[n]->get_pixel(x, y)[0]);
				
//					fprintf(stderr, "[%d %d %d %f (%f %f %f) (%f %f %f)] ", n, x, y, cur_cost, 
//							cl->pz_color[n]->get_pixel(x, y)[0],
//							cl->pz_color[n]->get_pixel(x, y)[1],
//							cl->pz_color[n]->get_pixel(x, y)[2],
//							options[0][0],
//							options[0][1],
//							options[0][2]
//							);
				improved = 1;
			}
#endif
		}
#endif
	}

	static void reduce_cost_to_search_depth(const char *d_out[], const char *v_out[], d2::exposure *exp_out, int inc_bit) {
		int max_depth = 2;
		int improved = 1;
		int count = 0;
		void *visit_data;
		// time_t start_seconds, cur_seconds, max_seconds = 1200;
		double match_weight_shallow = 0.0001;
		double match_weight_deep = 0.0004;

		assert (max_depth > 0);

		visit_data = init_visit_data(max_depth);

		// time(&start_seconds);

		while(reduce_lod());

		while ((improved /*&& count < 40*/) || cl->next) {

			if (inc_bit)
			for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
				if (d_out[i] != NULL) {
					const d2::image *im = depth(i);
					d2::image_rw::write_image(d_out[i], im, exp_out, 1, 1);
				}

				if (v_out[i] != NULL) {
					const d2::image *im = view(i);
					d2::image_rw::write_image(v_out[i], im, exp_out);
					delete im;
				}
			}

			if (!improved/* || count >= 40*/) {
				// fprintf(stderr, "scale factor: %f\n", cl->sf);
				fprintf(stderr, ".");
				assert (cl->next);

				increase_lod();
				count = 0;

#if 0
				if (cl->sf >= 0.125) {
					while (cl->next)
						increase_lod();
					count = 100000;
					improved = 0;
					continue;
				}
#endif
			}

			count++;
			improved = 0;

			// time (&cur_seconds);

			// if (cur_seconds - start_seconds > max_seconds)
			//	continue;

			for (unsigned int n = 0; n < d2::image_rw::count(); n++)
			for (unsigned int x = 0; x < cl->pz_color[n]->height(); x++)
			for (unsigned int y = 0; y < cl->pz_color[n]->width();  y++) {
#if 1
				double multiplier = 1;
				double lateral_min = 0.5;
				double depth_ceiling = -pow(10, 5);
				double d0 = cl->pz_depth[n]->get_pixel(x, y)[0];
				double dp1 = d0 + depth_quantum;
				double dm1 = d0 - depth_quantum;
				double err = 0;
				double m1_err = 0;
				double p1_err = 0;

				while (dp1 == d0 || dm1 == d0
				    && dp1 < 0
				    && dm1 < 0) {
					// fprintf(stderr, "Multiplying...\n");
					multiplier *= 2;
					dp1 = d0 + depth_quantum * multiplier;
					dm1 = d0 - depth_quantum * multiplier;
				}

				if (d0 >= 0)
					continue;
				if (dp1 >= 0)
					dp1 = 0 - depth_quantum;
				if (dm1 >= 0)
					dm1 = 0 - depth_quantum;

				for (unsigned int m = 0; m < d2::image_rw::count(); m++) {
					if (m == n)
						continue;

					point p0 = frame_to_frame(m, n, x, y, d0);

					if (!cl->pz_color[m]->in_bounds(d2::point(p0[0], p0[1])))
						continue;

					double r_source = cl->pz_color[n]->get_pixel(x, y)[0];
					double r0 = cl->pz_color[m]->get_bl(d2::point(p0[0], p0[1]))[0];
					double s0 = cl->pz_depth[m]->get_bl(d2::point(p0[0], p0[1]))[0];
					point pp1 = frame_to_frame(m, n, x, y, dp1);
					point pm1 = frame_to_frame(m, n, x, y, dm1);

					while (sqrt(pow(pm1[0] - p0[0], 2) + pow(pm1[1] - p0[1], 2)) < lateral_min
					    && sqrt(pow(pp1[0] - p0[0], 2) + pow(pp1[1] - p0[1], 2)) < lateral_min
					    && dm1 > depth_ceiling / 2
					    && d0 + depth_quantum * multiplier * 2 < 0) {
						multiplier *= 2;
						dp1 = d0 + depth_quantum * multiplier;
						dm1 = d0 - depth_quantum * multiplier;
						pp1 = frame_to_frame(m, n, x, y, dp1);
						pm1 = frame_to_frame(m, n, x, y, dm1);
					}

					if (dp1 >= 0 || dm1 >= 0) {
						// fprintf(stderr, "Neg.\n");
						continue;
					}

					if (!cl->pz_color[m]->in_bounds(d2::point(pp1[0], pp1[1]))
					 || !cl->pz_color[m]->in_bounds(d2::point(pm1[0], pm1[1])))
						continue;

					err += pow(r_source - r0, 2)
					     + ((p0[2] < s0) 
					      ? match_weight_deep
					      : match_weight_shallow)
					     * pow(p0[2] - s0, 2);

					double rp1 = cl->pz_color[m]->get_bl(d2::point(pp1[0], pp1[1]))[0];
					double sp1 = cl->pz_depth[m]->get_bl(d2::point(pp1[0], pp1[1]))[0];

					p1_err += pow(r_source - rp1, 2) 
						+ ((pp1[2] < sp1)
						 ? match_weight_deep
						 : match_weight_shallow)
						* pow(pp1[2] - sp1, 2);


					double rm1 = cl->pz_color[m]->get_bl(d2::point(pm1[0], pm1[1]))[0];
					double sm1 = cl->pz_depth[m]->get_bl(d2::point(pm1[0], pm1[1]))[0];

					m1_err += pow(r_source - rm1, 2)
						+ ((pm1[2] < sm1)
						 ? match_weight_deep
						 : match_weight_shallow)
						* pow(pm1[2] - sm1, 2);
				}

				if (m1_err < p1_err && m1_err < err) {
					cl->pz_depth[n]->pix(x, y)[0] = dm1;
					improved = 1;
				} else if (p1_err < err) {
					cl->pz_depth[n]->pix(x, y)[0] = dp1;
					improved = 1;
				}
#else


				d2::pixel options[1];
				double cur_cost = frc_cost(n, x, y);

				// if (cur_cost < 0.0001)
				//	continue;

				frc_samples(n, x, y, cur_cost/2, 1, options);

				if (try_color_change(n, x, y, 0, max_depth, options[0], cur_cost/2, visit_data) < 0) {

					assert (options[0][0] == cl->pz_color[n]->get_pixel(x, y)[0]);
					
//					fprintf(stderr, "[%d %d %d %f (%f %f %f) (%f %f %f)] ", n, x, y, cur_cost, 
//							cl->pz_color[n]->get_pixel(x, y)[0],
//							cl->pz_color[n]->get_pixel(x, y)[1],
//							cl->pz_color[n]->get_pixel(x, y)[2],
//							options[0][0],
//							options[0][1],
//							options[0][2]
//							);
					improved = 1;
				}
#endif
			}
		}

#if 0
		
		fprintf(stderr, "\n\n******* DEPTH TRANSLATION DEBUGGING FEATURE ENABLED *******\n\n");

		for (unsigned int x = 0; x < cl->pz_depth[0]->height(); x++)
		for (unsigned int y = 0; y < cl->pz_depth[0]->width(); y++)
			cl->pz_depth[0]->pix(x, y)[0] += 100;
#endif

#if 0

		fprintf(stderr, "\n\n******* FRAME-TO-FRAME DEBUGGING FEATURE ENABLED ********\n\n");

		/*
		 * This is debugging instrumentation that maps pixels from
		 * frame 0 onto areas in frame 1.
		 */

		for (unsigned int x = 0; x < cl->pz_color[0]->height(); x++)
		for (unsigned int y = 0; y < cl->pz_color[0]->width();  y++) {

			point p(x, y, cl->pz_depth[0]->get_pixel(x, y)[0]);

			point q = frame_to_frame(1, 0, p);

			// fprintf(stderr, "ftf [%f %f %f] [%f %f %f]\n", p[0], p[1], p[2], q[0], q[1], q[2]);

			int x2 = (int) floor(q[0]);
			int y2 = (int) floor(q[1]);

			if (x2 >= 0
			 && y2 >= 0
			 && (unsigned) x2 < cl->pz_color[1]->height()
			 && (unsigned) y2 < cl->pz_color[1]->width()) {
#if 1
				cl->pz_color[1]->set_pixel(x2, y2, cl->pz_color[0]->get_pixel(x, y));
#else
				cl->pz_color[1]->set_pixel(x2, y2, d2::pixel(x / (double) cl->pz_color[0]->height(), 
							                     y / (double) cl->pz_color[0]->width(), 0));
				cl->pz_color[0]->set_pixel(x, y, d2::pixel(x2 / (double) cl->pz_color[1]->height(), 
							                   y2 / (double) cl->pz_color[1]->width(), 0));
#endif
			}
		}
#endif

		destroy_visit_data(visit_data);
	}

#if 0
	/*
	 * Describe a scene to a renderer
	 */
	static void describe(render *r) {
	}
#endif
};

#endif
