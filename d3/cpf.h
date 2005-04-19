// Copyright 2003, 2004, 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * d3/cpf.h: Control point file interface.
 */

#ifndef __cpf_h__
#define __cpf_h__

#include "point.h"

#define CPF_VERSION 0
#define CPF_VERSION_MAX 0

class cpf {
	static FILE *load_f;
	static FILE *save_f;
	static int load_version;

	static const char *load_n;
	static const char *save_n;
	static int save_version;

	static ale_pos cpp_upper;
	static ale_pos cpp_lower;

	/*
	 * TYPE is:
	 * 	0 type A
	 * 	1 type B
	 * 	2 type C
	 */
	struct control_point {
		int type;
		point *d2;
		point d3;
	};

	static struct control_point *cp_array;
	static unsigned int cp_array_max;
	static unsigned int cp_index;

	static void error(const char *message) {
		fprintf(stderr, "cpf: Error: %s", message);
		exit(1);
	}

	static void get_integer(int *i) {
		if(fscanf(load_f, " %d", i) != 1)
			error("Could not get integer.");
	}

	static void get_double(double *d) {
		if (fscanf(load_f, " %lf", d) != 1)
			error("Could not get double.");
	}

	static void get_new_line() {
		int next_character = 0;

		while (next_character != EOF
		    && next_character != '\n')
			next_character = fgetc(load_f);
	}

	static void check_version(int v) {
		if (v > CPF_VERSION_MAX)
			error("Incompatible version number.");
	}

	/*
	 * Calculate the centroid of a set of points.
	 */
	static point calculate_centroid(point *points, int n) {
		point sum(0, 0, 0);
		int count = 0;

		for (int i = 0; i < n; i++) {
			if (!points[i].defined())
				continue;
			sum += points[i];
			count++;
		}

		return sum / count;
	}

	/*
	 * Measure the error between a projected system and a solved coordinate.
	 */
	static ale_accum measure_projected_error(point solved, const point coords[], int n) {
		ale_accum error = 0;
		ale_accum divisor = 0;

		for (int i = 0; i < n; i++) {
			if (!coords[i].defined())
				continue;

			pt t = align::projective(i);

			point sp = t.wp_unscaled(solved);

			sp[2] = coords[i][2];

			error += (sp - coords[i]).normsq();
			divisor += 1;
		}

		return sqrt(error / divisor);
	}

	static ale_accum measure_total_error() {
		ale_accum result = 0;
		ale_accum divisor = 0;

		for (unsigned int i = 0; i < cp_array_max; i++) {
			ale_accum e = measure_projected_error(cp_array[i].d3, cp_array[i].d2, d2::image_rw::count());
			if (!finite(e)) 
				continue;

			result += e * e;
			divisor += 1;
		}

		return sqrt(result / divisor);
	}

	/*
	 * Solve for a 3D point from a system of projected coordinates.
	 *
	 * The algorithm is this:
	 *
	 * First, convert all 2D points to 3D points by projecting them onto a
	 * plane perpendicular to the view axis of the corresponding camera,
	 * and passing through the origin of the world coordinate system.
	 * Then, for each point, move it along the view ray to the point
	 * closest to the centroid of all of the points.  Repeat this last loop
	 * until the largest adjustment is smaller than some specified lower
	 * bound.
	 */
	static point solve_projected_system(const point points[], int n) {

		/*
		 * Copy the passed array.
		 */
		point *points_copy = new point[n];
		for (int i = 0; i < n; i++)
			points_copy[i] = points[i];

		/*
		 * Set an initial depth for each point, and convert it to world
		 * coordinates.
		 */

		for (int i = 0; i < n; i++) {
			pt t = align::projective(i);

			point p = t.wc(point(0, 0, 0));
			ale_pos plane_depth = p[2];

			points_copy[i][2] = plane_depth;

			points_copy[i] = t.pw_unscaled(points_copy[i]);
		}

		/*
		 * For each point, adjust the depth along the view ray to
		 * minimize the distance from the centroid of the points_copy.
		 */

		ale_pos max_diff = 0;
		ale_pos prev_max_diff = 0;
		ale_pos diff_bound = 0.99999;

		while (!(max_diff / prev_max_diff > diff_bound) || !finite(max_diff / prev_max_diff)) {

			/*
			 * Calculate the centroid of all points_copy.
			 */

			point centroid = calculate_centroid(points_copy, n);

//			fprintf(stderr, "md %f pmd %f ratio %f ", max_diff, prev_max_diff, max_diff / prev_max_diff);
//			fprintf(stderr, "centroid %f %f %f ", centroid[0], centroid[1], centroid[2]);
//			fprintf(stderr, "%f ", measure_projected_error(centroid, points, n));

			prev_max_diff = max_diff;
			max_diff = 0;

			for (int i = 0; i < n; i++) {
//				fprintf(stderr, "points_copy[%d] %f %f %f ", i, points_copy[i][0], points_copy[i][1], points_copy[i][2]);

				if (!points_copy[i].defined())
					continue;

				pt t = align::projective(i);
				point camera_origin = t.cw(point(0, 0, 0));

				point v = points_copy[i] - camera_origin;
				ale_pos l = (centroid - camera_origin).norm();
				ale_pos alpha = camera_origin.anglebetw(points_copy[i], centroid);

				point new_point = camera_origin + v / v.norm() * l * cos(alpha);

				ale_pos diff = points_copy[i].lengthto(new_point);

				if (diff > max_diff)
					max_diff = diff;

				points_copy[i] = new_point;
			}

		}
//		fprintf(stderr, "%f\n", measure_projected_error(calculate_centroid(points_copy, n), points, n));
		// fprintf(stderr, "md %f pmd %f ratio %f ", max_diff, prev_max_diff, max_diff / prev_max_diff);
		
		point result = calculate_centroid(points_copy, n);

		delete points_copy;

		return result;
	}

	static void solve_total_system() {
		for (unsigned i = 0; i < cp_array_max; i++) {
			if (cp_array[i].type == 0) 
				cp_array[i].d3 = solve_projected_system(cp_array[i].d2, d2::image_rw::count());
		}
	}

	/*
	 * Type A control point record
	 *
	 * A <frame 0 coord 1> <frame 0 coord 0> ... <frame n coord 0>
	 */
	static void get_type_a() {

		point result;

		/*
		 * Get the number of frames.
		 */
		int n = d2::image_rw::count();

		/*
		 * Allocate storage for N frames.
		 */
		point *coords = new point[n];

		for (int i = 0; i < n; i++)
		for (int j = 0; j < 2; j++) {
			double dp;
			get_double(&dp);
			coords[i][1 - j] = dp;
		}

		result = solve_projected_system(coords, n);

		fprintf(stderr, "Coordinate system error: %f\n", measure_projected_error(result, coords, n));

		cp_array[cp_array_max - 1].d2 = coords;
		cp_array[cp_array_max - 1].d3 = result;
		cp_array[cp_array_max - 1].type = 0;

	}

	/*
	 * Type B control point record
	 *
	 * B <coord 1> <coord 0> <coord 2>
	 */
	static void get_type_b() {
		double d[3];

		get_double(&d[1]);
		get_double(&d[0]);
		get_double(&d[2]);

		cp_array[cp_array_max - 1].d3 = point(d[0], d[1], d[2]);
		cp_array[cp_array_max - 1].type = 1;
	}

	/*
	 * Type C control point record
	 *
	 * C <type A data> <type B data>
	 */
	static void get_type_c() {
		get_type_a();
		get_type_b();
		cp_array[cp_array_max - 1].type = 2;
	}

public:

	static void set_cpp_upper(ale_pos cu) {
		cpp_upper = cu;
	}

	static void set_cpp_lower(ale_pos cl) {
		cpp_lower = cl;
	}

	static void init_loadfile(const char *filename) {
		load_n = filename;
		load_f = fopen(load_n, "r");

		if (!load_f) {
			fprintf(stderr, "cpf: Error: could not open control point file '%s'.", load_n);
			exit(1);
		}
	}

	static void init_savefile(const char *filename) {
		save_n = filename;
		save_f = fopen(save_n, "w");

		if (!save_f) {
			fprintf(stderr, "cpf: Error: could not open control point file '%s'.", save_n);
			exit(1);
		}

		fprintf(save_f, "# created by ALE control point file handler version %d\n", CPF_VERSION);

		fclose(save_f);
	}

	static void read_file() {
		while (load_f && !feof(load_f)) {
			int command_char;

			command_char = fgetc(load_f);

			switch (command_char) {
				case EOF:
					return;
				case '#':
				case ' ':
				case '\t':
					get_new_line();
					break;
				case '\r':
				case '\n':
					break;
				case 'V':
					get_integer(&load_version);
					check_version(load_version);
					get_new_line();
				        break;
				case 'A':
					cp_array = (control_point *) realloc(cp_array, ++cp_array_max * sizeof(control_point));
					assert(cp_array);
					get_type_a();
					get_new_line();
					break;
				case 'B':
					cp_array = (control_point *) realloc(cp_array, ++cp_array_max * sizeof(control_point));
					assert(cp_array);
					get_type_b();
					get_new_line();
					break;
				case 'C':
					cp_array = (control_point *) realloc(cp_array, ++cp_array_max * sizeof(control_point));
					assert(cp_array);
					get_type_c();
					get_new_line();
					break;
				default:
					error("Unrecognized command");
			}
		}
	}

	static void adjust_cameras() {
		ale_accum current_error = measure_total_error();
		unsigned int n = d2::image_rw::count();

		/*
		 * Perturbation measured in pixels or degrees
		 *
		 * XXX: should probably be pixel arclength instead of degrees.
		 */

		ale_pos max_perturbation = pow(2, floor(log(cpp_upper) / log(2)));
		ale_pos min_perturbation = cpp_lower;
		ale_pos perturbation = max_perturbation;

		fprintf(stderr, "Init error %f\n", current_error);

		while (perturbation >= min_perturbation) {

			ale_accum previous_error;
			ale_pos angular_p = perturbation / 180 * M_PI;

			fprintf(stderr, "P %f AP %f ", perturbation, angular_p);

			do {

				/*
				 * Minimum frame to adjust
				 */
				int M = 1;

				previous_error = current_error;

				ale_accum test_error;

				/*
				 * Try adjusting camera positions
				 */

				for (unsigned int i =  M; i <  n;  i++)
				for (unsigned int d =  0; d <  3;  d++) 
				for (         int s = -1; s <= 1; s+=2) {
					align::adjust_translation(i, d, s * perturbation);
					solve_total_system();
					test_error = measure_total_error();
					if (test_error < current_error) {
						current_error = test_error;
					} else {
						align::adjust_translation(i, d, -s * perturbation);
					}
				}

				/*
				 * Try adjusting camera orientations
				 */

				for (unsigned int i =  M; i <  n;  i++)
				for (unsigned int d =  0; d <  3;  d++) 
				for (         int s = -1; s <= 1; s+=2) {
					align::adjust_rotation(i, d, s * angular_p);
					solve_total_system();
					test_error = measure_total_error();
					if (test_error < current_error) {
						current_error = test_error;
					} else {
						align::adjust_rotation(i, d, -s * angular_p);
					}
				}

				/*
				 * Try adjusting view angle
				 */
				for (int s = -1; s <= 1; s+=2) {
					align::adjust_view_angle(s * angular_p);
					solve_total_system();
					test_error = measure_total_error();
					if (test_error < current_error) {
						current_error = test_error;
					} else {
						align::adjust_view_angle(-s * angular_p);
					}
				}
			} while (current_error < previous_error);

			fprintf(stderr, "E %f\n", current_error);

			perturbation /= 2;
		}

		solve_total_system();
	}
	
	static point get() {
		ale_accum error_bound = 4;

		if (cp_array == NULL) {
			read_file();
			adjust_cameras();
		}

		if (cp_array == NULL)
			return point::undefined();

		while (cp_index < cp_array_max 
		    && error_bound < measure_projected_error(cp_array[cp_index].d3, 
		                                             cp_array[cp_index].d2, 
					                     d2::image_rw::count()))
			cp_index++;

		if (cp_index >= cp_array_max)
			return point::undefined();

		fprintf(stderr, "Coordinate system error: %f\n", 
				measure_projected_error(cp_array[cp_index].d3, cp_array[cp_index].d2, d2::image_rw::count()));

		return cp_array[cp_index++].d3;
	}

	static void set(point p) {
		assert(0);
	}

	static void finish_loadfile() {
	}

	static void finish_savefile() {
	}
};

#endif
