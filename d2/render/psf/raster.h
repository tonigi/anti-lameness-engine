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
	psf_result generic_response(ale_real *response_array, float top, float
			bot, float lef, float rig) const {

		assert (response_array != NULL);

		psf_result result;

		if (top < min_i())
			top = min_i();
		if (bot > max_i())
			bot = max_i();
		if (lef < min_j())
			lef = min_j();
		if (rig > max_j())
			rig = max_j();

		int il = (int) floor((top - min_i()) / (max_i() - min_i()) * _filter_dim_i);
		int ih = (int) floor((bot - min_i()) / (max_i() - min_i()) * (_filter_dim_i - 0.001));
		int jl = (int) floor((lef - min_j()) / (max_j() - min_j()) * _filter_dim_j);
		int jh = (int) floor((rig - min_j()) / (max_j() - min_j()) * (_filter_dim_j - 0.001));

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

			for (int k = 0; k < 3; k++) {
				result.matrix(k, k) += (ale_real) ((lbot - ltop) * (lrig - llef)
					      * response_array[3 * _filter_dim_j * ii + 3 * jj + k]);
			}
		}

		return result;
	}

public:
	virtual unsigned int varieties() const = 0;

	virtual unsigned int select(unsigned int i, unsigned int j) const = 0;

	/*
	 * Get a specific pixel response.
	 */
	psf_result operator()(float top, float bot, float lef, float rig, unsigned int variety) const {
		assert (variety < num_arrays);

		ale_real *response_array = response_arrays[variety];
		assert (response_array != NULL);

		return generic_response(response_array, top, bot, lef, rig);
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
