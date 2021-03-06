// Copyright 2003 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

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

#ifndef __psf_stdin_vg_h__
#define __psf_stdin_vg_h__

#include "../../point.h"
#include "psf.h"

/*
 * Point-spread function module.
 *
 * This response function is configured by input from stdin.  A series of
 * prompts indicates the information required.
 */

class psf_stdin_vg : public psf {
	ale_real _height;
	ale_real _width;
	ale_real gap_width;
	int _filter_dim_i;
	int _filter_dim_j;
	ale_real *response_array;
public:
	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	ale_real min_i() const { return -_height; }
	ale_real max_i() const { return  _height; }
	ale_real min_j() const { return -_width - fabs(gap_width); }
	ale_real max_j() const { return  _width + fabs(gap_width); }

	/*
	 * Get the number of varieties supported by this PSF.  These usually
	 * correspond to different points in the sensor array.
	 */
	virtual unsigned int varieties() const {
		return 8;
	}

	/*
	 * Select the variety appropriate for a given position in the sensor
	 * array.
	 */
	virtual unsigned int select(unsigned int i, unsigned int j) {
		return (j % 8);
	}

	/*
	 * Response functions
	 *
	 * response_generic() and operator()() 
	 *
	 * Get the response to the rectangle bounded by (top, bot, lef, rig).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  The variety of the responding pixel is provided, in
	 * case response is not uniform for all pixels (e.g. some sensor arrays
	 * stagger red, green, and blue sensors).
	 */
	psf_result response_generic(ale_real *response_array, ale_real top, ale_real bot,
			ale_real lef, ale_real rig, unsigned int variety) const {

		assert (response_array != NULL);
		assert (variety < varieties());

		psf_result result;

		ale_pos offset = (gap_width / 2)
			       - (gap_width / 7) * variety;

		lef -= offset;
		rig -= offset;

		if (top < min_i())
			top = min_i();
		if (bot > max_i())
			bot = max_i();
		if (lef < -_width)
			lef = -_width;
		if (rig >  _width)
			rig =  _width;

		int il = (int) floor((top - min_i()) / (max_i() - min_i()) * _filter_dim_i);
		int ih = (int) floor((bot - min_i()) / (max_i() - min_i()) * _filter_dim_i);
		int jl = (int) floor((lef +  _width) / (_width * 2) * _filter_dim_j);
		int jh = (int) floor((rig +  _width) / (_width * 2) * _filter_dim_j);

		// fprintf(stderr, "(il, ih, jl, jh) = (%d, %d, %d, %d)\n", il, ih, jl, jh);

		for (int ii = il; ii <= ih; ii++)
		for (int jj = jl; jj <= jh; jj++) {

			ale_real ltop = ((ale_real) ii) / (ale_real) _filter_dim_i * (max_i() - min_i()) + min_i();
			ale_real lbot = ((ale_real) ii + 1) / (ale_real) _filter_dim_i * (max_i() - min_i()) + min_i();
			ale_real llef = ((ale_real) jj) / (ale_real) _filter_dim_j * ((ale_real) _width * (ale_real) 2) - _width;
			ale_real lrig = ((ale_real) jj + 1) / (ale_real) _filter_dim_j * (_width * (ale_real) 2) - _width;

			if (ltop < top)
				ltop = top;
			if (lbot > bot)
				lbot = bot;
			if (llef < lef)
				llef = lef;
			if (lrig > rig)
				lrig = rig;

			for (int k = 0; k < 3; k++) {
				result.matrix(k, k) += (ale_real) ((lbot - ltop) * (lrig - llef)
					      * response_array[3 * _filter_dim_j * ii + 3 * jj + k]);
			}
		}

		return result;
	}

	psf_result operator()(ale_real top, ale_real bot, ale_real lef, ale_real
			rig, unsigned int variety) const {
		return response_generic(response_array, top, bot, lef, rig, variety);
	}

#if 0
	psf_result operator()(ale_real top, ale_real bot, ale_real lef, ale_real rig) const {
		return response_generic(response_array, top, bot, lef, rig, 4);
	}
#endif

	void class_error() {
		fprintf(stderr, "\n\nALE Panic: Error acquiring input.  Exiting.\n");
		exit(1);
	}

	psf_stdin_vg () {

		printf("\nEnter vertical gap width, in units of pixels (e.g. 1.0): ");
		fflush(stdout);
		double dgap_width;
		if (scanf("%lf", &dgap_width) != 1) {
			class_error();
		}
		gap_width = dgap_width;

		printf("\nEnter filter support height, in units of pixels (e.g. 2.5): ");
		fflush(stdout);
		double dheight;
		if (scanf("%lf", &dheight) != 1) {
			class_error();
		}
		_height = dheight / 2;

		printf("\nEnter filter support width, in units of pixels (e.g. 2.5): ");
		fflush(stdout);
		double dwidth;
		if (scanf("%lf", &dwidth) != 1) {
			class_error();
		}
		_width = dwidth / 2;

		printf("\nEnter the number of rows in the filter (e.g. 3): ");
		fflush(stdout);
		if (scanf("%d", &_filter_dim_i) != 1 || _filter_dim_i < 1) {
			class_error();
		}

		printf("\nEnter the number of columns in the filter (e.g. 3): ");
		fflush(stdout);
		if (scanf("%d", &_filter_dim_j) != 1 || _filter_dim_j < 1) {
			class_error();
		}

		response_array = (ale_real *) malloc(_filter_dim_i * _filter_dim_j * 3
				* sizeof(ale_real));

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
			fflush(stdout);
			double delem;
			if (scanf("%lf", &delem) != 1)
				class_error();
			response_array[i * _filter_dim_j * 3 + j * 3 + k] = delem;
		}
	}

	virtual ~psf_stdin_vg() {

		/*
		 * Don't free this without creating a copy constructor.
		 */

		free(response_array);
	}
};

#endif
