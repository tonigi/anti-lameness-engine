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

#ifndef __normalizer_h__
#define __normalizer_h__

#include "../../point.h"
#include "rasterizer.h"
#include "raster.h"

/*
 * Normalizer for rasterized PSFs.
 *
 * This class normalizes a rasterized PSF.
 */

class normalizer : public raster {
	raster *input;
public:
	unsigned int varieties() const {
		return input->varieties();
	}

	unsigned int select(unsigned int i, unsigned int j) const {
		return input->select(i, j);
	}

private:
	void initialize_response_array(ale_real *response_array) {
		pixel integral;

		integral = integrate(response_array);

		for (unsigned int i = 0; i < _filter_dim_i; i++)
		for (unsigned int j = 0; j < _filter_dim_j; j++)
		for (unsigned int k = 0; k < 3            ; k++)
			response_array[i * _filter_dim_j * 3 + j * 3 + k] /= integral[k];
	}

public:
	normalizer (raster *input) {
		this->input = input;

		_height = -input->min_i();
		assert (input->max_i() == _height);

		_width = -input->min_j();
		assert (input->max_j() == _width);

		/*
		 * The element structure matches that of the input.
		 */

		_filter_dim_i = input->max_elem_i();
		_filter_dim_j = input->max_elem_j();

		/*
		 * Ensure that the array has an odd number of elements in each
		 * direction.  This allows us to move the center to the right
		 * place when using FFTW.
		 */

		assert (_filter_dim_i % 2 == 1);
		assert (_filter_dim_j % 2 == 1);

		/*
		 * Determine the number of arrays to create.
		 */

		num_arrays = input->varieties();

		/*
		 * Create arrays
		 */

		response_arrays = (ale_real **)malloc(num_arrays * sizeof(ale_real *));

		if (!response_arrays) {
			fprintf(stderr, "Could not allocate in normalizer.\n");
			exit(1);
		}

		for (unsigned int n = 0; n < num_arrays; n++) {
			response_arrays[n] = (ale_real *)malloc(_filter_dim_i * _filter_dim_j * 3
					                    * sizeof(ale_real));
					
			if (!response_arrays[n]) {
				fprintf(stderr, "Could not allocate in normalizer.\n");
				exit(1);
			}

			for (unsigned int i = 0; i < _filter_dim_i; i++)
			for (unsigned int j = 0; j < _filter_dim_j; j++)
			for (unsigned int k = 0; k < 3;             k++) {
				response_arrays[n][i * _filter_dim_j * 3 + j * 3 + k]
					= input->element(n, i, j, k);
			}

			initialize_response_array(response_arrays[n]);
		}

#if 0
		avg_response = (ale_real *)malloc(_filter_dim_i * _filter_dim_j * 3
				              * sizeof(ale_real));
					      
		if (!avg_response) {
			fprintf(stderr, "Could not allocate in normalizer.\n");
			exit(1);
		}

		for (unsigned int i = 0; i < _filter_dim_i; i++)
		for (unsigned int j = 0; j < _filter_dim_j; j++)
		for (unsigned int k = 0; k < 3;             k++) {
			avg_response[i * _filter_dim_j * 3 + j * 3 + k]
				= input->element(i, j, k);
		}

		initialize_response_array(avg_response);
#endif
		compute_integrals();
	}
};

#endif
