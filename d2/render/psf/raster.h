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

#ifndef __raster_h__
#define __raster_h__

#include "../../point.h"
#include "psf.h"

/*
 * Raster point-spread function.
 */

class raster : public psf {
protected:
	ale_real _height;
	ale_real _width;
	unsigned int _filter_dim_i;
	unsigned int _filter_dim_j;
	unsigned int num_arrays;
	ale_real **response_arrays;
	ale_real **response_partials;
#if 0
	ale_real *avg_response;
#endif
	pixel *response_integrals;
#if 0
	pixel avg_integral;
#endif
public:

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	float min_i() const { return -_height; }
	float max_i() const { return  _height; }
	float min_j() const { return -_width; }
	float max_j() const { return  _width; }

	/*
	 * Element accessor methods.
	 */
	unsigned int max_elem_i() {
		return _filter_dim_i;
	}
	unsigned int max_elem_j() {
		return _filter_dim_j;
	}
	ale_real element(unsigned int n, unsigned int i, unsigned int j, unsigned int k) {

		assert (n < num_arrays);
		assert (i < _filter_dim_i);
		assert (j < _filter_dim_j);
		assert (k < 3);
		
		return response_arrays[n][i * _filter_dim_j * 3 + j * 3 + k];
	}
#if 0
	ale_real element(unsigned int i, unsigned int j, unsigned int k) {
		assert (i < _filter_dim_i);
		assert (j < _filter_dim_j);
		assert (k < 3);
		
		return avg_response[i * _filter_dim_j * 3 + j * 3 + k];
	}
#endif

	/*
	 * Response function
	 *
	 * Get the response to the rectangle bounded by (top, bot, lef, rig).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  
	 *
	 * generic_response (private):
	 *
	 * A response array for this generic function is provided by the
	 * calling function, and a specific or average response is returned,
	 * based on this array.
	 *
	 * operator():
	 *
	 * The index of the specific response array is provided, from which the
	 * specific response is determined.  Alternatively, if no index is
	 * specified, then the average response is returned.
	 */
private:
	psf_result generic_response(ale_real *response_partial, float top, float
			bot, float lef, float rig, char channels) const {

		assert (response_partial != NULL);

		psf_result result;

		/*
		 * lrintf() can be more efficient than floor() or float-to-int
		 * casts.  For more details, see Erik de Castro Lopo, "Faster
		 * Floating Point to Integer Conversions":
		 *
		 * 	http://mega-nerd.com/FPcast/
		 *
		 * In this case, lrintf() seems to be a bit faster than plain
		 * casting, and much faster than floor(0.5 + ...).  Casting
		 * from round() seems to be an acceptable alternative to
		 * lrintf().
		 *
		 * Early calculation of common floating-point constants in the
		 * following code is based on an initial implementation by HJ
		 * Hornbeck.
		 */

		ale_real i_element_scale = (ale_real) _filter_dim_i / (max_i() - min_i());
		ale_real j_element_scale = (ale_real) _filter_dim_j / (max_j() - min_j());

		int il = (int) lrintf(i_element_scale * (top - min_i()));
		int ih = (int) lrintf(i_element_scale * (bot - min_i()));
		int jl = (int) lrintf(j_element_scale * (lef - min_j()));
		int jh = (int) lrintf(j_element_scale * (rig - min_j()));

		/*
		 * Bounds clamping may be faster when performed in integer 
		 * arithmetic than in floating-point, so we do this after 
		 * float-to-int conversion is complete.
		 */

		if (il < 0)
			il = 0;
		if (jl < 0)
			jl = 0;
		if (ih > (int) _filter_dim_i)
			ih = (int) _filter_dim_i;
		if (jh > (int) _filter_dim_j)
			jh = (int) _filter_dim_j;

		if (!(il < ih) || !(jl < jh))
			return result;

		for (int k = 0; k < 3; k++) {
			if (!((1 << k) & channels))
				continue;
			assert (ih > 0 && jh > 0);
			assert (ih <= (int) _filter_dim_i);
			assert (jh <= (int) _filter_dim_j);

			ale_real result_k = 0;

			if (il > 0 && jl > 0) 
				result_k += response_partial[k + 3 * (jl - 1) + 3 * _filter_dim_j * (il - 1)];
			if (il > 0) 
				result_k -= response_partial[k + 3 * (jh - 1) + 3 * _filter_dim_j * (il - 1)];
			if (jl > 0)
				result_k -= response_partial[k + 3 * (jl - 1) + 3 * _filter_dim_j * (ih - 1)];
			result_k += response_partial[k + 3 * (jh - 1) + 3 * _filter_dim_j * (ih - 1)];
			result.set_matrix(k, k, result_k);
		}

		return result;
	}

public:
	virtual unsigned int varieties() const = 0;

	virtual unsigned int select(unsigned int i, unsigned int j) const = 0;

	/*
	 * Get a specific pixel response.
	 */
	psf_result operator()(float top, float bot, float lef, float rig, unsigned int variety, 
			char channels) const {
		assert (variety < num_arrays);

		ale_real *response_partial = response_partials[variety];
		assert (response_partial != NULL);

		return generic_response(response_partial, top, bot, lef, rig, channels);
	}

	psf_result operator()(float top, float bot, float lef, float rig, 
			unsigned int variety) const {
		return operator()(top, bot, lef, rig, variety, 0x7);
	}

#if 0
	/*
	 * Get the average pixel response.
	 */
	psf_result operator()(float top, float bot, float lef, float rig) const {
		return generic_response(avg_response, top, bot, lef, rig);
	}
#endif

protected:
	/*
	 * Integrate over the whole PSF
	 */
	pixel integrate(ale_real *response_array) {
		pixel result;

		for (unsigned int i = 0; i < _filter_dim_i; i++)
		for (unsigned int j = 0; j < _filter_dim_j; j++)
		for (unsigned int k = 0; k < 3            ; k++)
			result[k] += response_array[i * _filter_dim_j * 3 + j * 3 + k];

		for (unsigned int k = 0; k < 3; k++)
			result[k] *= (4 * _height * _width)
				/ (_filter_dim_i * _filter_dim_j);

		return result;
	}

	void partial_integrate(ale_real *target, ale_real *source) {
		ale_real element_area = (max_i() - min_i())
		                      * (max_j() - min_j())
				      / (ale_real) (_filter_dim_i)
				      / (ale_real) (_filter_dim_j);

		for (unsigned int i = 0; i < _filter_dim_i; i++)
		for (unsigned int j = 0; j < _filter_dim_j; j++)
		for (unsigned int k = 0; k < 3            ; k++) {
			unsigned int index = i * _filter_dim_j * 3 + j * 3 + k;
			target[index] = source[index] * element_area
			              + ((j > 0) ? target[index - 3] : 0)
				      + ((i > 0) ? target[index - _filter_dim_j * 3] : 0)
				      - ((i > 0 
				       && j > 0) ? target[index - _filter_dim_j * 3 - 3] : 0);
		}
	}

	/*
	 * Compute integrals.
	 */
	void compute_integrals() {
		response_integrals = new pixel[num_arrays];

		for (unsigned int n = 0; n < num_arrays; n++)
			response_integrals[n] = integrate(response_arrays[n]);
		
#if 0
		avg_integral = integrate(avg_response);
#endif

		response_partials = (ale_real **) malloc(sizeof(ale_real *) * num_arrays);
		assert(response_partials);
		for (unsigned int n = 0; n < num_arrays; n++) {
			response_partials[n] = (ale_real *) malloc(sizeof(ale_real) * _filter_dim_i 
										    * _filter_dim_j
										    * 3);
			assert(response_partials[n]);
			partial_integrate(response_partials[n], response_arrays[n]);
		}
	}

public:

	/*
	 * Return elements of given integrals
	 */
	pixel integral(unsigned int n) const {
		assert (response_integrals != NULL);
		return response_integrals[n];
	}
#if 0
	pixel integral() const {
		return avg_integral;
	}
#endif

	raster () {
		response_integrals = NULL;
	}

	virtual ~raster() {

		/*
		 * Deallocate data structures.
		 */

		for (unsigned int n = 0; n < num_arrays; n++)
			free(response_arrays[n]);
		free(response_arrays);
#if 0
		free(avg_response);
#endif

		if (response_integrals)
			delete response_integrals;
	}
};

#endif
