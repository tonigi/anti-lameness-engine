// Copyright 2003 David Hilvert <dhilvert@ugcs.caltech.edu>

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

#ifndef __ipc_stdin_h__
#define __ipc_stdin_h__

#include <assert.h>
#include "../point.h"

/*
 * Module for Irani-Peleg image reconstruction.  
 *
 * This response function is configured by input from stdin.  A series of
 * prompts indicates the information required.
 */

class ipc_stdin {
	double _radius;
	int _filter_dim_i;
	int _filter_dim_j;
	double *response_array;
public:
	/*
	 * Result type.
	 */
	class ipc_result {
		friend class ipc_stdin;
		double response[3];
	public:
		/*
		 * Response intensity on channel K2 resulting from
		 * stimulus on channel K1.  Channel R=0, G=1, B=2.
		 */
		double operator()(int k1, int k2) {
			if (k1 == k2)
				return response[k1];
			else
				return 0;
		}
	};

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	float min_i() { return -_radius; }
	float max_i() { return  _radius; }
	float min_j() { return -_radius; }
	float max_j() { return  _radius; }

	/*
	 * Get the response to the quadrilateral (p[0], p[1], p[2], p[3]).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  The position (i, j) of the responding pixel is
	 * provided, in case response is not uniform for all pixels (e.g. some
	 * sensor arrays stagger red, green, and blue sensors).
	 */
	struct ipc_result operator()(const point p[4], int i, int j) {
		assert (response_array != NULL);

		struct ipc_result result;

		result.response[0] = 0;
		result.response[1] = 0;
		result.response[2] = 0;

		/* 
		 * XXX: We make some unsafe assumptions here.
		 */

		float top = p[0][0];
		float bot = p[2][0];
		float lef = p[0][1];
		float rig = p[1][1];

		if (top < min_i())
			top = min_i();
		if (bot > max_i())
			bot = max_i();
		if (lef < min_j())
			lef = min_j();
		if (rig > max_j())
			rig = max_j();

		// fprintf(stderr, "(t, b, l, r) = (%f, %f, %f, %f)\n", top, bot, lef, rig);

		int il = (int) floor((top - min_i()) / (max_i() - min_i()) * _filter_dim_i);
		int ih = (int) floor((bot - min_i()) / (max_i() - min_i()) * (_filter_dim_i - 0.001));
		int jl = (int) floor((lef - min_j()) / (max_j() - min_j()) * _filter_dim_j);
		int jh = (int) floor((rig - min_j()) / (max_j() - min_j()) * (_filter_dim_j - 0.001));

		// fprintf(stderr, "(il, ih, jl, jh) = (%d, %d, %d, %d)\n", il, ih, jl, jh);

		for (int ii = il; ii <= ih; ii++)
		for (int jj = jl; jj <= jh; jj++) {

			float ltop = ((float) ii) / _filter_dim_i * (max_i() - min_i()) + min_i();
			float lbot = ((float) ii + 1) / _filter_dim_i * (max_i() - min_i()) + min_i();
			float llef = ((float) jj) / _filter_dim_j * (max_j() - min_j()) + min_j();
			float lrig = ((float) jj + 1) / _filter_dim_j * (max_j() - min_j()) + min_j();

			if (ltop < top)
				ltop = top;
			if (lbot > bot)
				lbot = bot;
			if (llef < lef)
				llef = lef;
			if (lrig > rig)
				lrig = rig;

			// fprintf(stderr, "(ii, jj, lt, lb, ll, lr) = (%d, %d, %f, %f, %f, %f)\n", 
				//	ii, jj, ltop, lbot, llef, lrig);

			assert (ii >= 0);
			assert (ii < 6);
			assert (jj >= 0);
			assert (jj < 6);

			for (int k = 0; k < 3; k++) {
				result.response[k] += (double) ((lbot - ltop) * (lrig - llef)
					      * response_array[3 * _filter_dim_j * ii + 3 * jj + k]);
			}
		}

		// fprintf(stderr, "(result[0], result[1], result[2]) = (%d, %d, %d)\n",
		//		result.response[0], result.response[1], result.response[2]);

		return result;
	}

	void class_error() {
		fprintf(stderr, "\n\nALE Panic: Error acquiring input.  Exiting.\n");
		exit(1);
	}

	ipc_stdin () {
		printf("Using the Standard Input Irani-Peleg configuration module.\n");

		printf("\nEnter filter diameter, in units of pixels (e.g. 2.5): ");
		if (scanf("%lf", &_radius) != 1) {
			class_error();
		}
		_radius /= 2;

		printf("\nEnter the number of rows in the filter (e.g. 3): ");
		if (scanf("%d", &_filter_dim_i) != 1 || _filter_dim_i < 1) {
			class_error();
		}

		printf("\nEnter the number of columns in the filter (e.g. 3): ");
		if (scanf("%d", &_filter_dim_j) != 1 || _filter_dim_j < 1) {
			class_error();
		}

		response_array = (double *) malloc(_filter_dim_i * _filter_dim_j * 3
				* sizeof(double));

		if (response_array == NULL) {
			fprintf(stderr, "\n\nCould not allocate filter.\n");
			exit(1);
		}

		printf("\nFilter elements are labeled as (row, column, channel).  The red channel of\n");
		printf("the top-left element is (0, 0, 0), and the blue channel of the bottom-right\n");
		printf("element is (%d, %d, 2).\n\n", _filter_dim_i - 1, _filter_dim_j - 1);

		for (int i = 0; i < _filter_dim_i; i++)
		for (int j = 0; j < _filter_dim_j; j++) 
		for (int k = 0; k < 3; k++) {
			printf("Enter value for element (%d, %d, %d) (e.g. 2.5): ",
					i, j, k);
			if (scanf("%lf", &(response_array[i * _filter_dim_j * 3 + j * 3 + k])) != 1)
				class_error();
		}
	}

	virtual ~ipc_stdin() {
		free(response_array);
	}
};

#endif
