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

#ifndef __psf_h__
#define __psf_h__

#include "../../point.h"

/*
 * Point-spread function module abstract base class.
 */

class psf {
public:
	/*
	 * Result type is a matrix.
	 */
	class psf_result {
		friend class psf;
	protected:
		ale_real _matrix[3][3];

	public:
		psf_result() {

			/*
			 * Simplified version -- diagonal matrix
			 */

			for (int i = 0; i < 3; i++)
				_matrix[i][i] = 0;
		}

		ale_real get_matrix(unsigned int i, unsigned int j) {
			assert (i < 3);
			assert (j < 3);

			/*
			 * Simplified version -- diagonal matrix
			 */

			if (i != j)
				return 0;
			else
				return _matrix[i][j];
		}

		void set_matrix(unsigned int i, unsigned int j, ale_real value) {
			assert (i < 3);
			assert (j < 3);

			/*
			 * Simplified version -- diagonal matrix
			 */

			assert (i == j || value == 0);

			_matrix[i][j] = value;
		}

		ale_real &matrix(unsigned int i, unsigned int j) {
			assert (i < 3);
			assert (j < 3);

			/*
			 * Simplified version -- diagonal matrix
			 */

			assert (i == j);

			return _matrix[i][j];
		}

		pixel operator()(pixel p) {

			/*
			 * Simplified version -- diagonal matrix
			 */

			return pixel(_matrix[0][0] * p[0],
				     _matrix[1][1] * p[1],
				     _matrix[2][2] * p[2]);

		}

		/*
		 * Weights associated with the result
		 */
		pixel weight() {

			/*
			 * Simplified version -- diagonal matrix
			 */

			return pixel(
				_matrix[0][0],
				_matrix[1][1],
				_matrix[2][2]);
		}

		void operator*=(ale_real scale) {
			/*
			 * Simplified version -- diagonal matrix
			 */

			for (int i = 0; i < 3; i++)
				_matrix[i][i] *= scale;
		}
	};

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	virtual ale_real min_i() const = 0;
	virtual ale_real max_i() const = 0;
	virtual ale_real min_j() const = 0;
	virtual ale_real max_j() const = 0;

	/*
	 * Get the number of varieties supported by this PSF.  These usually
	 * correspond to different points in the sensor array.
	 */
	virtual unsigned int varieties() const {
		return 1;
	}

	/*
	 * Select the variety appropriate for a given position in the sensor
	 * array.
	 */
	virtual unsigned int select(unsigned int i, unsigned int j) {
		return 0;
	}

	/*
	 * Get the response to the rectangle bounded by (top, bot, lef, rig).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  One of several varieties can be selected, usually
	 * based on position (e.g. some sensor arrays stagger red, green, and
	 * blue sensors).
	 */
	virtual psf_result operator()(ale_real top, ale_real bot, ale_real lef, ale_real
			rig, unsigned int variety) const = 0;

	virtual psf_result operator()(ale_real top, ale_real bot, ale_real lef, ale_real
			rig, unsigned int variety, char channels) const {
		return operator()(top, bot, lef, rig, variety);
	}


#if 0
	/*
	 * Get the average pixel response.  This function should be overloaded
	 * for PSFs that support multiple varieties.
	 */
	virtual psf_result operator()(ale_real top, ale_real bot, ale_real lef, ale_real rig) const {
		return operator()(top, bot, lef, rig, 0);
	}
#endif

	virtual ~psf() {
	}
};

#endif
