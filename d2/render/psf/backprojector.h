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

#ifndef __backprojector_h__
#define __backprojector_h__

#include "../../point.h"
#include "rasterizer.h"
#include "raster.h"

/*
 * Backprojector for rasterized PSFs.
 *
 * This class converts a rasterized PSF into a rasterized backprojection array.
 */

class backprojector : public raster {
	raster *input;
public:
	unsigned int varieties() const {
		return input->varieties();
	}

	unsigned int select(unsigned int i, unsigned int j) const {
		return input->select(i, j);
	}

private:
	/*
	 * Backprojection for the Irani-Peleg renderer 
	 *
	 * Applying a special case of theorem 4.1 of the source paper by Irani
	 * and Peleg, convergence can be assured for a single image, with
	 * uniform PSF, no change in sampling rate, and taking the normalizing
	 * divisor c == 1, if H[psf](f)H[aux](f) is real and within the open
	 * interval (0, 2), where H[psf] is the frequency-domain representation
	 * of the point-spread function and H[aux] is the frequency-domain
	 * representation of the backprojection kernel.  We can guarantee that
	 * H[psf](f)H[aux](f) is real by making H[aux](f) == k(f)H[psf](f)*,
	 * where k is a real function and '*' indicates the complex conjugate.
	 * If k(f) is equal to 1 for all f, then this is equivalent to the
	 * condition h[aux](x) == h[psf](-x), where h[] are the time domain
	 * representations of the respective functions.  Since this negation
	 * of position is implicitly performed in ipc.h, we don't perform it
	 * here.
	 *
	 * However, to ensure that the range (0, 2) is satisfied, it may be
	 * necessary for k(f) to assume a value other than 1.  We choose a
	 * constant function k, in accordance with the source paper's
	 * normalizing divisor c, but this is not required.  We use FFTW
	 * when available, but it is likely that common cases will not observe
	 * any speed improvement.
	 */
	void initialize_response_array(ale_real *response_array) {
		int cols = _filter_dim_j;
		int rows = _filter_dim_i;

#ifdef USE_FFTW
		fftw_complex *inout;
		fftw_plan p_forward;
		fftw_plan p_backward;

		inout = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cols * rows);

		p_forward = fftw_plan_dft_2d(rows, cols, inout, inout, FFTW_FORWARD, FFTW_ESTIMATE);
		p_backward = fftw_plan_dft_2d(rows, cols, inout, inout, FFTW_BACKWARD, FFTW_ESTIMATE);

		for (int k = 0; k < 3; k++) {

			for (int i = 0; i < rows * cols; i++) {
				/*
				 * Write the values to the FFTW input array,
				 * shifting by (rows * cols - 1) / 2 in order
				 * to accommodate the implicit translation.
				 */

				inout[i][0] = response_array[((i + (rows * cols - 1)/2) * 3 + k)
					                   % (rows * cols * 3)];
				inout[i][1] = 0;
			}

			fftw_execute(p_forward);

			/*
			 * Find the frequency with maximum magnitude, then
			 * adjust this according to the sampling rate
			 * (filter resolution).
			 */

			ale_real max_magnitude = 0;

			for (int i = 0; i < rows * cols; i++) {
				ale_real input_magnitude;

				input_magnitude = sqrt(pow(inout[i][0], 2) + pow(inout[i][1], 2));
				if (input_magnitude > max_magnitude)
					max_magnitude = input_magnitude;
			}

			max_magnitude *= (4 * _height * _width) / (rows * cols);

			/*
			 * Scale the magnitude of all of the frequencies and perform
			 * conjugation.  
			 */

			for (int i = 0; i < rows * cols; i++) {

				/*
				 * Adjust the magnitude
				 *
				 * Note: since we're currently dividing all frequencies
				 * by the same value, there's no need to divide in the 
				 * frequency domain.  However, we might want to do 
				 * something else in the future, so it might be
				 * good to leave the code like this for now.
				 */

				inout[i][0] = inout[i][0] * pow(0.9 / max_magnitude, 2);
				inout[i][1] = inout[i][1] * pow(0.9 / max_magnitude, 2);

				/*
				 * Perform conjugation
				 *
				 * Note: conjugation is implicit in ipc.h, so we omit the
				 * step here.
				 */

				/* inout[i][1] = -inout[i][1]; */
			}

			fftw_execute(p_backward);

			for (int i = 0; i < rows * cols; i++) {
				/*
				 * Read the values from the FFTW output array,
				 * shifting by (rows * cols - 1) / 2 in order
				 * to accommodate the implicit translation.
				 */

				response_array[((i + (rows * cols - 1)/2) * 3 + k)
					         % (rows * cols * 3)] 
						 = inout[i][0]
						 / (rows * cols);
			}
		}

		fftw_destroy_plan(p_forward);
		fftw_destroy_plan(p_backward);
		fftw_free(inout);
#else

		for (int k = 0; k < 3; k++) {
			ale_real *real1 = (ale_real *) calloc(rows * cols, sizeof(ale_real));
			ale_real *imag1 = (ale_real *) calloc(rows * cols, sizeof(ale_real));
			ale_real *real2 = (ale_real *) calloc(rows * cols, sizeof(ale_real));
			ale_real *imag2 = (ale_real *) calloc(rows * cols, sizeof(ale_real));

			assert (real1 && imag1 && real2 && imag2);

			if (!(real1 && imag1 && real2 && imag2)) {
				fprintf(stderr, "Unable to allocate memory in backprojector.\n");
				exit(1);
			}

			/*
			 * Calculate frequencies.  We implement the equations indicated by
			 * the FFTW3 info page (section "What FFTW Really Computes").
			 */

			for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			for (int jj = 0; jj < cols; jj++) {
				real1[i * cols + j] += response_array[((i * cols + jj +
							     (rows * cols - 1)/2) * 3 + k)
					                   % (rows * cols * 3)]
						     * cos((-2 * M_PI * j * jj) / cols);
				imag1[i * cols + j] += response_array[((i * cols + jj +
							     (rows * cols - 1)/2) * 3 + k)
					                   % (rows * cols * 3)]
						     * sin((-2 * M_PI * j * jj) / cols);
			}

			for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			for (int ii = 0; ii < rows; ii++) {
				real2[i * cols + j] += real1[ii * cols + j]
						     * cos((-2 * M_PI * i * ii) / rows)
						     - imag1[ii * cols + j]
						     * sin((-2 * M_PI * i * ii) / rows);
				imag2[i * cols + j] += real1[ii * cols + j]
						     * sin((-2 * M_PI * i * ii) / rows)
						     + imag1[ii * cols + j]
						     * cos((-2 * M_PI * i * ii) / rows);
			}

			/*
			 * Find the frequency with maximum magnitude, then
			 * adjust this according to the sampling rate
			 * (filter resolution).
			 */

			ale_real max_magnitude = 0;

			for (int i = 0; i < rows * cols; i++) {
				ale_real input_magnitude;

				input_magnitude = sqrt(pow(real2[i], 2) + pow(imag2[i], 2));
				if (input_magnitude > max_magnitude)
					max_magnitude = input_magnitude;
			}

			max_magnitude *= (4 * _height * _width) / (rows * cols);

			for (int i = 0; i < rows * cols; i++)
				response_array[i * 3 + k] *= pow(0.9 / max_magnitude, 2);

			free(real1);
			free(imag1);
			free(real2);
			free(imag2);
		}

#endif
	}

public:
	backprojector (raster *input) {
		this->input = input;

		_height = -input->min_i();
		assert (input->max_i() == _height);

		_width  = -input->min_j();
		assert (input->max_j() == _width);

		/*
		 * The element structure matches that of the input.
		 */

		_filter_dim_i = input->max_elem_i();
		_filter_dim_j = input->max_elem_j();

		/*
		 * Ensure that the array has an odd number of elements in each
		 * direction.  This allows us to move the center to the right
		 * place when using a discrete FT.
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
			fprintf(stderr, "Could not allocate in backprojector.\n");
			exit(1);
		}

		for (unsigned int n = 0; n < num_arrays; n++) {
			response_arrays[n] = (ale_real *)malloc(_filter_dim_i * _filter_dim_j * 3
					                    * sizeof(ale_real));
					
			if (!response_arrays[n]) {
				fprintf(stderr, "Could not allocate in backprojector.\n");
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

		avg_response = (ale_real *)malloc(_filter_dim_i * _filter_dim_j * 3
				              * sizeof(ale_real));
					      
		if (!avg_response) {
			fprintf(stderr, "Could not allocate in backprojector.\n");
			exit(1);
		}

		for (unsigned int i = 0; i < _filter_dim_i; i++)
		for (unsigned int j = 0; j < _filter_dim_j; j++)
		for (unsigned int k = 0; k < 3;             k++) {
			avg_response[i * _filter_dim_j * 3 + j * 3 + k]
				= input->element(i, j, k);
		}

		initialize_response_array(avg_response);

		compute_integrals();
	}
};

#endif
