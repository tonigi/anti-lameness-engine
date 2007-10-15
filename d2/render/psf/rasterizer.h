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

#ifndef __psf_rasterize_h__
#define __psf_rasterize_h__

#include "../../point.h"
#include "raster.h"
#include "psf.h"

/*
 * Point-spread function rasterizer.
 *
 * These operations rasterize a PSF to a multiple of the resolution of the
 * rendering grid for a given frame.  
 */

class rasterizer : public raster {
	psf *input;

public:
	unsigned int varieties() const {
		return input->varieties();
	}

	unsigned int select(unsigned int i, unsigned int j) const {
		return input->select(i, j);
	}

	rasterizer (psf *input, transformation t) {
		this->input = input;

		_height = -input->min_i();

		if (input->max_i() > _height)
			_height = input->max_i();

		_width = -input->min_j();

		if (input->max_j() > _width)
			_width = input->max_j();

		/*
		 * Approximate the desired resolution.
		 *
		 * Assume that maximum resolution is reached at (at least) one
		 * of the corners of the image.  (This should be true for
		 * projective, Euclidean, and translational transformations,
		 * but it would be worthwhile to check/prove this for the
		 * projective case at some point, since it's a bit less
		 * obvious.)
		 */

		point min_diff;

		/*
		 * XXX: this loop breaks when height <= 1 or width <= 1.
		 */

		for (unsigned int i = 0; i < t.unscaled_height(); i += (t.unscaled_height() - 1))
		for (unsigned int j = 0; j < t.unscaled_width();  j += (t.unscaled_width()  - 1)) {
			point corner = point(i, j);
			point delta1 = corner - t.scaled_inverse_transform(t.transform_scaled(corner) + point(1, 0));
			point delta2 = corner - t.scaled_inverse_transform(t.transform_scaled(corner) + point(0, 1));

			for (int index = 0; index < 2; index++) {
				ale_pos d1 = fabs(delta1[index]);
				ale_pos d2 = fabs(delta2[index]);

				/*
				 * Take the largest change in each direction.
				 */

				ale_pos delta = (d1 > d2) ? d1 : d2;

				if ((i == 0 && j == 0) || delta < min_diff[index])
					min_diff[index] = delta;
			}
		}

		ale_real resolution_multiplier = 20;  /* Arbitrary */

		_filter_dim_i = (int) ceil((ale_real) 2 * _height * resolution_multiplier / ale_pos_to_real(min_diff[0]));
		_filter_dim_j = (int) ceil((ale_real) 2 * _width * resolution_multiplier / ale_pos_to_real(min_diff[1]));

		/*
		 * Ensure that the array has an odd number of elements in each
		 * direction.  This allows us to move the center to the right
		 * place when using FFTW.
		 */

		if (_filter_dim_i % 2 == 0)
			_filter_dim_i++;
		if (_filter_dim_j % 2 == 0)
			_filter_dim_j++;

		/*
		 * Determine the number of arrays to create.
		 */

		num_arrays = input->varieties();

		/*
		 * Create arrays
		 */

		response_arrays = (ale_real **)malloc(num_arrays * sizeof(ale_real *));

		if (!response_arrays) {
			fprintf(stderr, "Could not allocate in rasterizer.\n");
			exit(1);
		}

		ale_real stepsize_i = (2 * _height) / _filter_dim_i;
		ale_real stepsize_j = (2 * _width) / _filter_dim_j;
		ale_real divisor = stepsize_i * stepsize_j;

		for (unsigned int n = 0; n < num_arrays; n++) {
			response_arrays[n] = (ale_real *)malloc(_filter_dim_i * _filter_dim_j * 3
					                    * sizeof(ale_real));
					
			if (!response_arrays[n]) {
				fprintf(stderr, "Could not allocate in rasterizer.\n");
				exit(1);
			}

			for (unsigned int i = 0; i < _filter_dim_i; i++)
			for (unsigned int j = 0; j < _filter_dim_j; j++) {
				psf_result r 
					= (*input)(-_height + stepsize_i * (ale_real) i,
						   -_height + stepsize_i * (ale_real) (i + 1),
						   -_width  + stepsize_j * (ale_real) j,
						   -_width  + stepsize_j * (ale_real) (j + 1), n);

				for (unsigned int k = 0; k < 3; k++) {
					response_arrays[n][i * _filter_dim_j * 3 + j * 3 + k]
						= r.matrix(k, k) / divisor;
				}
			}
		}

#if 0
		avg_response = (ale_real *)malloc(_filter_dim_i * _filter_dim_j * 3
				              * sizeof(ale_real));
					      
		if (!avg_response) {
			fprintf(stderr, "Could not allocate in rasterizer.\n");
			exit(1);
		}

		for (unsigned int i = 0; i < _filter_dim_i; i++)
		for (unsigned int j = 0; j < _filter_dim_j; j++) {
			psf::psf_result r 
				= (*input)(-_height + stepsize_i * i,
					   -_height + stepsize_i * (i + 1),
					   -_width  + stepsize_j * j,
					   -_width  + stepsize_j * (j + 1));

			for (unsigned int k = 0; k < 3; k++)
				avg_response[i * _filter_dim_j * 3 + j * 3 + k]
					= r.matrix(k, k) / divisor;
		}
#endif

		compute_integrals();
		
//		fprintf(stderr, "(w=%f h=%f we=%d he=%d [", _width, _height, _filter_dim_j, _filter_dim_i);
//		for (unsigned int i = 0; i < _filter_dim_i; i++)
//		for (unsigned int j = 0; j < _filter_dim_j; j++)
//			fprintf(stderr, "%f ", response_arrays[0][i * _filter_dim_j * 3 + j * 3 + 0]);
//		fprintf(stderr, "])");
	}
};

#endif
