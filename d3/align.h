// Copyright 2003 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * d3/align.h: Handle alignment of view pyramids.
 *
 * XXX: this class assumes that the view volume is a right rectangular-based
 * pyramid, which holds for most photographic situations.
 *
 * XXX: this class assumes that the horizontal and vertical angles of the view
 * volume pyramid are proportional to the horizontal and vertical pixel
 * resolutions.
 */

#ifndef __d3align_h__
#define __d3align_h__

#include "point.h"
#include "et.h"
#include "pt.h"

class align {
	static ale_pos _init_angle;
	static et *alignment_array;

	/*
	 * Estimate the point at which the pyramidal axis passes through the
	 * plane, based on the given 2-dimensional transformation and the given
	 * point set indicating the corners of the quadrilateral of
	 * intersection.  (Both transformation and point set alone are
	 * sufficient for calculation.)
	 *
	 * Using the point set, the desired point is exactly the point of
	 * intersection between line segments drawn between opposite corners of
	 * the quadrilateral of intersection.  The algebraic details of this
	 * calculation are desribed below.
	 *
	 * If we represent the desired point by (e1, e2), and the four corners
	 * of the quadrilateral of intersection by (a1, a2), ..., (d1, d2) in
	 * clockwise (or counter-clockwise) order, then we get the two-equation
	 * vector system (implicitly four equations)
	 *
	 * (e1, e2) = (a1, a2) + x[(c1, c2) - (a1, a2)]
	 *          = (b1, b2) + y[(d1, d2) - (b1, b2)]
	 *
	 * Solving for x in terms of y, we get the two-equation system
	 *
	 * x = (b1 - yb1 + yd1 - a1) / (c1 - a1)
	 *   = (b2 - yb2 + yd2 - a2) / (c2 - a2)
	 *
	 * And
	 *
	 * y = (c1b2 - c1a2 - a1b2 - c2b1 + c2a1 - a2b1)
	 *   /(-c2b1 + c2d1 + a2b1 - a2d1 + c1b2 - c1d2 - a1b2 + a1d2)
	 *
	 * However, it's much easier just to project the center point of the
	 * source image onto the quadrilateral of intersection using the
	 * given 2-dimensional transformation.
	 */
	static d2::point axis_intersection(d2::point a, d2::point b, d2::point c, d2::point d, d2::transformation t) {
#if 0
		ale_pos y = (c[0]*b[1] - c[0]*a[1] - a[0]*b[1] - c[1]*b[0] 
			   + c[1]*a[0] - a[1]*b[0])
			  /(-c[1]*b[0] + c[1]*d[0] + a[1]*b[0] - a[1]*d[0] 
			   + c[0]*b[1] - c[0]*d[1] - a[0]*b[1] + a[0]*d[1]);

		return b + y * (d - b);
#else
		return t.transform_scaled(d2::point((t.scaled_height() - 1) / 2, (t.scaled_width() - 1) / 2));
#endif
	}

	/*
	 * gd_position is a gradient-descent automatic positioning algorithm
	 * used to determine the location of the view pyramid apex.  
	 *
	 * Starting from an initial camera position ESTIMATE, this method
	 * determines the angles between each pair of legs in the view volume
	 * (there are six pairs).  The sum of squares of the differences
	 * between the thusly determined angles and the known desired angles is
	 * then used as an error metric, and a kind of gradient descent is used
	 * to find a position for which this metric is minimized.
	 */
	static point gd_position(point estimate, d2::point a, d2::point b, d2::point c, d2::point d, d2::transformation t) {
		ale_pos w = t.scaled_width();
		ale_pos h = t.scaled_height();

		/*
		 * The desired diagonal angle is given (as init_angle).
		 * Calculate the desired side angles as follows:
		 *
		 * The distance to the apex is 
		 *
		 * 	D = sqrt(h*h + w*w) / (2 * tan(init_angle/2))
		 *
		 * Then
		 *
		 * 	desired_h_angle = 2 * arctan(h / (2 * sqrt(D*D + w*w/4)))
		 * 	desired_w_angle = 2 * arctan(w / (2 * sqrt(D*D + h*h/4)))
		 */

		ale_pos D = sqrt(h * h + w * w)
			         / (2 * tan(_init_angle/2));
		ale_pos desired_h_angle = 2 * atan(h / (2 * sqrt(D*D + w*w/4)));
		ale_pos desired_w_angle = 2 * atan(w / (2 * sqrt(D*D + h*h/4)));

		ale_pos estimate_h1_angle = estimate.anglebetw(a, b);
		ale_pos estimate_h2_angle = estimate.anglebetw(c, d);
		ale_pos estimate_w1_angle = estimate.anglebetw(a, d);
		ale_pos estimate_w2_angle = estimate.anglebetw(b, c);
		ale_pos estimate_d1_angle = estimate.anglebetw(a, c);
		ale_pos estimate_d2_angle = estimate.anglebetw(b, d);

		ale_pos error = sqrt(pow(estimate_h1_angle - desired_h_angle, 2)
			           + pow(estimate_h2_angle - desired_h_angle, 2)
			           + pow(estimate_w1_angle - desired_w_angle, 2)
			           + pow(estimate_w2_angle - desired_w_angle, 2)
			           + pow(estimate_d1_angle - _init_angle    , 2)
			           + pow(estimate_d2_angle - _init_angle    , 2));

		/*
		 * Vary the magnitude by which each coordinate of the position
		 * can be changed at each step.  
		 */

		double view_angle = _init_angle;

		for (ale_pos magnitude = estimate[2] / 2; 
			magnitude >= 1;
			magnitude /= 2) {

			/*
			 * Continue searching at this magnitude while error <
			 * old_error.  (Initialize old_error accordingly.)
			 */

			ale_pos old_error = error * 2;

			while(old_error > error) {

//				ale_pos D = sqrt(h * h + w * w)
//						 / (2 * tan(view_angle/2));
//				ale_pos desired_h_angle = 2 * atan(h / (2 * sqrt(D*D + w*w/4)));
//				ale_pos desired_w_angle = 2 * atan(w / (2 * sqrt(D*D + h*h/4)));

//				fprintf(stderr, ".");


//				fprintf(stderr, "estimate: [%f %f %f %f %f %f]\n", 
//						estimate_h1_angle,
//						estimate_h2_angle,
//						estimate_w1_angle,
//						estimate_w2_angle,
//						estimate_d1_angle,
//						estimate_d2_angle);
//
//				fprintf(stderr, "desired : [%f %f %f]\n",
//						desired_h_angle,
//						desired_w_angle,
//						view_angle);
					
				old_error = error;

				for (double c0 = -magnitude; c0 <= magnitude; c0 += magnitude)
				for (double c1 = -magnitude; c1 <= magnitude; c1 += magnitude) 
				for (double c2 = -magnitude; c2 <= magnitude; c2 += magnitude) 
				for (double c3 = -magnitude; c3 <= magnitude; c3 += magnitude) {

					if (c3 > 10)
						c3 = 10;

					// fprintf(stderr, "[%f %f %f]\n", c0, c1, c2);

					estimate[0] += c0;
					estimate[1] += c1;
					estimate[2] += c2;
					// view_angle += c3 / 30;

					ale_pos D = sqrt(h * h + w * w)
							 / (2 * tan(view_angle/2));
					ale_pos desired_h_angle = 2 * atan(h / (2 * sqrt(D*D + w*w/4)));
					ale_pos desired_w_angle = 2 * atan(w / (2 * sqrt(D*D + h*h/4)));

					estimate_h1_angle = estimate.anglebetw(a, b);
					estimate_h2_angle = estimate.anglebetw(c, d);
					estimate_w1_angle = estimate.anglebetw(a, d);
					estimate_w2_angle = estimate.anglebetw(b, c);
					estimate_d1_angle = estimate.anglebetw(a, c);
					estimate_d2_angle = estimate.anglebetw(b, d);

					ale_pos perturbed_error = 
							sqrt(pow(estimate_h1_angle - desired_h_angle, 2)
							   + pow(estimate_h2_angle - desired_h_angle, 2)
							   + pow(estimate_w1_angle - desired_w_angle, 2)
							   + pow(estimate_w2_angle - desired_w_angle, 2)
							   + pow(estimate_d1_angle - view_angle    , 2)
							   + pow(estimate_d2_angle - view_angle    , 2));

					if (perturbed_error < error) {
						error = perturbed_error;
					} else {
						estimate[0] -= c0;
						estimate[1] -= c1;
						estimate[2] -= c2;
						// view_angle -= c3 / 30;
					}
				}
			}
		}

//		fprintf(stderr, "error %f\n", error);


		return estimate;
	}

public:

	/* 
	 * Set the initial estimated diagonal viewing angle of the view
	 * pyramid.  This is the angle formed, e.g., between the top-left and
	 * bottom-right legs of the view pyramid.
	 */
	static void init_angle(ale_pos a) {
		_init_angle = a;
	}

	/*
	 * Initialize d3 transformations from d2 transformations.
	 *
	 * Here are three possible approaches for calculating camera position
	 * based on a known diagonal viewing angle and known area of
	 * intersection of the view volume with a fixed plane:
	 *
	 * (1) divide the rectangular-based pyramidal view volume into two
	 * tetrahedra.  A coordinate system is selected so that the triangle of
	 * intersection between one of the tetrahedra and the fixed plane has
	 * coordinates (0, 0, 0), (1, 0, 0), and (a, b, 0), for some (a, b).
	 * The law of cosines is then used to derive three equations
	 * associating the three known angles at the apex with the lengths of
	 * the edges of the tetrahedron.  The solution of this system of
	 * equations gives the coordinates of the apex in terms of a, b, and
	 * the known angles.  These coordinates are then transformed back to
	 * the original coordinate system.
	 *
	 * (2) The gradient descent approach taken by the method gd_position().
	 *
	 * (3) Assume that the camera is aimed (roughly) perpendicular to the
	 * plane.  In this case, the pyramidal axis forms with each adjacent
	 * pyramidal edge, in combination with a segment in the plane, a right
	 * triangle.  The distance of the camera from the plane is d/(2 *
	 * tan(theta/2)), where d is the distance between opposite corners of
	 * the quadrilateral of intersection and theta is the angle between
	 * opposite pairs of edges of the view volume.
	 *
	 * As an initial approach to the problem, we use (3) followed by (2).
	 *
	 * After position is estimated, we determine orientation from position.
	 * In order to do this, we determine the point at which the axis of
	 * the view pyramid passes through the plane, as described in the
	 * method axis_intersection().
	 */

	static void init_from_d2() { 
		assert (alignment_array == NULL);
		
		/*
		 * Initialize the alignment array.
		 */

		alignment_array = (et *) malloc (d2::image_rw::count() * sizeof(et));

		assert (alignment_array);

		if (!alignment_array) {
			fprintf(stderr, "\n\n*** Unable to allocate memory for 3D alignment array. ***\n\n");
			exit(1);
		}

		/*
		 * Iterate over input frames
		 */

		for (unsigned int i = 0; i < d2::image_rw::count(); i++) {

			/*
			 * Obtain the four mapped corners of the quadrilateral
			 * of intersection.
			 */

			d2::transformation t = d2::align::of(i);

			d2::point a = t.transform_scaled(d2::point(0, 0));
			d2::point b = t.transform_scaled(d2::point(t.scaled_height() - 1, 0));
			d2::point c = t.transform_scaled(d2::point(t.scaled_height() - 1, t.scaled_width() - 1));
			d2::point d = t.transform_scaled(d2::point(0, t.scaled_width() - 1));

			/*
			 * Estimate the point at which the pyramidal axis
			 * passes through the plane.  We denote the point as
			 * 'e'.
			 */

			d2::point e = axis_intersection(a, b, c, d, t);

			/*
			 * Determine the average distance between opposite
			 * corners, and calculate the distance to the camera
			 * based on method (3) described above.
			 */

			ale_pos dist1 = sqrt((a[0] - c[0])
					   * (a[0] - c[0])
					  +  (a[1] - c[1])
					   * (a[1] - c[1]));
			
			ale_pos dist2 = sqrt((b[0] - d[0])
					   * (b[0] - d[0])
					  +  (b[1] - d[1])
					   * (b[1] - d[1]));

			ale_pos avg_dist = (dist1 + dist2) / 2;


			ale_pos tangent = tan(_init_angle / 2);

			ale_pos distance = avg_dist / (2 * tangent);

			/*
			 * Set the position of the camera based on the estimate
			 * from method (3).  This is used as the initial
			 * position for method (2).  We assume that the camera 
			 * is placed so that its 0 and 1 coordinates match those
			 * of e.
			 */

			point estimate;

			estimate[0] = e[0];
			estimate[1] = e[1];
			estimate[2] = distance;

			// fprintf(stderr, "position (n=%d) %f %f %f\n", i, estimate[0], estimate[1], estimate[2]);

			/*
			 * Perform method (2).
			 */

			estimate = gd_position(estimate, a, b, c, d, t);

			// fprintf(stderr, ".");

			/*
			 * Assign transformation values based on the output of
			 * method (2), by modifying transformation parameters
			 * from the identity transformation.
			 */

			alignment_array[i] = et::identity();

			// fprintf(stderr, "position (n=%d) %f %f %f\n", i, estimate[0], estimate[1], estimate[2]);

			/*
			 * Modify translation values
			 */

			alignment_array[i].modify_translation(0, -estimate[0]);
			alignment_array[i].modify_translation(1, -estimate[1]);
			alignment_array[i].modify_translation(2, -estimate[2]);

			/*
			 * Assert that the z-axis translation is nonzero.  This
			 * is important for determining rotation values.
			 */

			assert (estimate[2] != 0);

			/*
			 * Modify rotation values
			 *
			 * The axis of the view pyramid should be mapped onto
			 * the z-axis.  Hence, the point of intersection
			 * between the world-coordinates xy-plane and the axis
			 * of the pyramid should occur at (0, 0) in local x-y
			 * coordinates.
			 *
			 * If we temporarily assume no rotation about the
			 * z-axis (e2 == 0), then the Euler angles (e0, e1, e2)
			 * used to determine rotation give us:
			 *
			 *	   x' = x cos e1 - (z cos e0 - y sin e0) sin e1
			 *	   y' = y cos e0 + z sin e0
			 *	   z' = (z cos e0 - y sin e0) cos e1 + x sin e1
			 *
			 * Setting x' = y' = 0, we get:
			 *
			 *	   e0 = arctan (-y/z)
			 *	   e1 = arctan (x/(z cos e0 - y sin e0))
			 *	        [optionally add pi to achieve z' < 0]
			 *
			 * To set e2, taking T to be the 3D transformation as a 
			 * function of e2, we first ensure that 
			 *
			 * 	tan(T(e2, a)[1] / T(e2, a)[0])
			 *    ==
			 *      tan(width / height)
			 *
			 * by assigning
			 *
			 * 	e2 += tan(width / height)
			 * 	    - tan(T(e2, a)[1] / T(e2, a)[0])
			 *
			 * Once this condition is satisfied, we conditionally
			 * assign
			 *
			 * 	e2 += pi 
			 *
			 * to ensure that
			 *
			 * 	T(e2, a)[0] < 0
			 * 	T(e2, a)[1] < 0
			 *
			 * which places the transformed point in the lower-left
			 * quadrant.  This is correct, since point A is mapped
			 * from (0, 0), the most negative point in the image.
			 */

			point e_translated = alignment_array[i](e);

//			fprintf(stderr, "axis-intersection (n=%d) e1 %f e0 %f e_t1 %f e_t0 %f e_t2 %f\n", 
//				i,
//				e[1],
//				e[0],
//				alignment_array[i](e)[1],
//				alignment_array[i](e)[0],
//				alignment_array[i](e)[2]);
//
//			fprintf(stderr, "camera origin (n=%d) o1 %f o0 %f o2 %f o_t1 %f o_t0 %f o_t2 %f\n", 
//				i,
//				estimate[1],
//				estimate[0],
//				estimate[2],
//				alignment_array[i](estimate)[1],
//				alignment_array[i](estimate)[0],
//				alignment_array[i](estimate)[2]);

			ale_pos e0 = atan(-e_translated[1] / e_translated[2]);
			ale_pos e1 = atan(e_translated[0]
					/ (e_translated[2] * cos(e0)
					 - e_translated[1] * sin(e0)));

			alignment_array[i].modify_rotation(0, e0);
			alignment_array[i].modify_rotation(1, e1);
			
			if (alignment_array[i](e)[2] > 0)
				alignment_array[i].modify_rotation(1, M_PI);

//			fprintf(stderr, "axis-intersection (n=%d) e1 %f e0 %f e_t1 %f e_t0 %f e_t2 %f\n", 
//				i,
//				e[1],
//				e[0],
//				alignment_array[i](e)[1],
//				alignment_array[i](e)[0],
//				alignment_array[i](e)[2]);
//
//			fprintf(stderr, "camera origin (n=%d) o_t1 %f o_t0 %f o_t2 %f\n", 
//				i,
//				alignment_array[i](estimate)[1],
//				alignment_array[i](estimate)[0],
//				alignment_array[i](estimate)[2]);

			point a_transformed = alignment_array[i](a);

			ale_pos e2 = atan((t.scaled_width() - 1) / (t.scaled_height() - 1)) 
				   - atan(a_transformed[1] / a_transformed[0]);

//			fprintf(stderr, "a components (n=%d) a1 %f a0 %f\n", i, a[1], a[0]);
//			fprintf(stderr, "e2 components (n=%d) tw %f th %f a_t1 %f a_t0 %f\n", i, 
//					t.scaled_width(), t.scaled_height(), a_transformed[1], a_transformed[0]);

			alignment_array[i].modify_rotation(2, e2);

//			fprintf(stderr, "rotation (n=%d) e0 %f e1 %f e2 %f\n", i, e0, e1, e2);

			if (alignment_array[i](a)[0] > 0) {
				e2 += M_PI;
				alignment_array[i].modify_rotation(2, M_PI);

				// fprintf(stderr, "adding 2pi to e2.");
			}

			// fprintf(stderr, "rotation (n=%d) e0 %f e1 %f e2 %f\n", i, e0, e1, e2);

//			fprintf(stderr, "(0, 0) output (n=%d) a1 %f a0 %f a_t1 %f a_t0 %f\n", 
//				i,
//				a[1],
//				a[0],
//				alignment_array[i](a)[1],
//				alignment_array[i](a)[0]);

			assert (alignment_array[i](a)[0] < 0);
			assert (alignment_array[i](a)[1] < 0);
		}
	}

	/*
	 * Get alignment for frame N.
	 */
	static et of(unsigned int n) {
		assert (n < d2::image_rw::count());
		assert (alignment_array);

		return alignment_array[n];
	}

	static double angle_of(unsigned int n) {
		return _init_angle;
	}

	static pt projective(unsigned int n) {
		return pt(d2::align::of(n), of(n), angle_of(n));
	}

	static void adjust_translation(unsigned int n, int d, ale_pos x) {
		alignment_array[n].modify_translation(d, x);
	}

	static void adjust_rotation(unsigned int n, int d, ale_pos x) {
		alignment_array[n].modify_rotation(d, x);
	}

	static void adjust_view_angle(ale_pos x) {
		_init_angle += x;
	}
};

#endif
