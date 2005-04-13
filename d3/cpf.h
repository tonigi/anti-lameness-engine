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
	static point solve_projected_system(point *points, int n) {

		/*
		 * Set an initial depth for each point, and convert it to world
		 * coordinates.
		 */

		for (int i = 0; i < n; i++) {
			pt t = align::projective(i);

			point p = t.wc(point(0, 0, 0));
			ale_pos plane_depth = p[2];

			points[i][2] = plane_depth;

			points[i] = t.pw_unscaled(points[i]);
		}

		/*
		 * For each point, adjust the depth along the view ray to
		 * minimize the distance from the centroid of the points.
		 */

		ale_pos max_diff = 0;
		ale_pos prev_max_diff = 0;
		ale_pos diff_bound = 0.99999;

		fprintf(stderr, "foo!\n");
		
		while (!(max_diff / prev_max_diff > diff_bound) || !finite(max_diff / prev_max_diff)) {

			/*
			 * Calculate the centroid of all points.
			 */

			point centroid = calculate_centroid(points, n);

//			fprintf(stderr, "md %f pmd %f ratio %f ", max_diff, prev_max_diff, max_diff / prev_max_diff);
//			fprintf(stderr, "centroid %f %f %f ", centroid[0], centroid[1], centroid[2]);

			prev_max_diff = max_diff;
			max_diff = 0;

			for (int i = 0; i < n; i++) {
//				fprintf(stderr, "points[%d] %f %f %f ", i, points[i][0], points[i][1], points[i][2]);

				if (!points[i].defined())
					continue;

				pt t = align::projective(i);
				point camera_origin = t.cw(point(0, 0, 0));

				point v = points[i] - camera_origin;
				ale_pos l = (centroid - camera_origin).norm();
				ale_pos alpha = camera_origin.anglebetw(points[i], centroid);

				point new_point = camera_origin + v / v.norm() * l * cos(alpha);

				ale_pos diff = points[i].lengthto(new_point);

				if (diff > max_diff)
					max_diff = diff;

				points[i] = new_point;
			}
		}
		fprintf(stderr, "md %f pmd %f ratio %f ", max_diff, prev_max_diff, max_diff / prev_max_diff);

		return calculate_centroid(points, n);
	}

	/*
	 * Type A control point record
	 *
	 * A <frame 0 coord 1> <frame 0 coord 0> ... <frame n coord 0>
	 */
	static point get_type_a() {

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

		get_new_line();

		result = solve_projected_system(coords, n);

		delete[] coords;

		return result;
	}

	/*
	 * Type B control point record
	 *
	 * B <coord 1> <coord 0> <coord 2>
	 */
	static point get_type_b() {
		double d[3];

		get_double(&d[1]);
		get_double(&d[0]);
		get_double(&d[2]);

		get_new_line();

		return point(d[0], d[1], d[2]);
	}

	/*
	 * Type C control point record
	 *
	 * C <type A data> <type B data>
	 */
	static point get_type_c() {
		get_type_a();
		return get_type_b();
	}

public:
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

	static point get() {
		while (load_f && !feof(load_f)) {
			int command_char;

			command_char = fgetc(load_f);

			switch (command_char) {
				case EOF:
					return point::undefined();
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
				        break;
				case 'A':
					return get_type_a();
				case 'B':
					return get_type_b();
				case 'C':
					return get_type_c();
				default:
					error("Unrecognized command");
			}
		}

		return point::undefined();
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
