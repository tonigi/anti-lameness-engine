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
       inline ale_real min_i() const { return -_height; }
       inline ale_real max_i() const { return  _height; }
       inline ale_real min_j() const { return -_width; }
       inline ale_real max_j() const { return  _width; }

	/*
	 * Element accessor methods.
	 */
       inline unsigned int max_elem_i() {
		return _filter_dim_i;
	}
       inline unsigned int max_elem_j() {
		return _filter_dim_j;
	}
       inline ale_real element(unsigned int n, unsigned int i, unsigned int j, unsigned int k) {

		assert (n < num_arrays);
		assert (i < _filter_dim_i);
		assert (j < _filter_dim_j);
		assert (k < 3);
		
		return response_arrays[n][i * _filter_dim_j * 3 + j * 3 + k];
	}
#if 0
       inline ale_real element(unsigned int i, unsigned int j, unsigned int k) {
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
		int il, jl, ih, jh;

		/** precompute a few values **/
		ale_real range_i = 2 * _height;
		ale_real range_j = 2 * _width;
		ale_real premult_i = _filter_dim_i / range_i;
		ale_real premult_j = _filter_dim_j / range_j;


		/** handle the outliers in an intelligent manner **/
		if (top < -_height) {

			il = 0;
			top = -_height;
		}
		else
			il = (int) floor( (top + _height) * premult_i );

		if (lef < -_width) {
			jl = 0;
			lef = -_width;
		}
		else
			jl = (int) floor( (lef + _width) * premult_j );

		/**   are the funny-looking subtractions needed?  **/
		/*
		 *   The -0.001 terms set a lower bound for each dimension, and
		 *   were added to avoid loss of precision, and hence image
		 *   artifacts, due to division by small numbers.  An
		 *   alternative approach could be used instead.
		 *   	
		 *       --dhilvert, 03-Feb-2007
		*/

		if (bot > _height) {

			ih = (int) floor( _filter_dim_i - 0.001 );
			bot = _height;
		}
		else
			ih = (int) floor( (bot + _height) / range_i * (_filter_dim_i - 0.001) );
                      
		if (rig > _width) {

			jh = (int) floor( _filter_dim_j - 0.001 );
			rig = _width;
		}
		else
			jh = (int) floor( (rig + _width) / range_j * (_filter_dim_i - 0.001) );


		for (int ii = il; ii <= ih; ii++)
		for (int jj = jl; jj <= jh; jj++) {

			float ltop = ((float) ii) / premult_i - _height;
			float lbot = ((float) ii + 1) / premult_i - _height;
			float llef = ((float) jj) / premult_j - _width;
			float lrig = ((float) jj + 1) / premult_j - _width;

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
