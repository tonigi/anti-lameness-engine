// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

/*
 * exposure_default.h: Default exposure properties.
 */

#ifndef __exposure_default_h__
#define __exposure_default_h__

/*
 * The default exposure is modeled after the simple power transfer function
 * described in
 *
 * http://netpbm.sourceforge.net/doc/pnmgamma.html
 *
 * Note: optimizations in d2/image_rw.h depend on the details of this function.
 */

class exposure_default : public exposure {
public:
	pixel linearize(pixel input) const {
#if 0
		/*
		 * Calling pow() may be expensive on some platforms (e.g.,
		 * those lacking hardware support for floating point).
		 */

		return ppow(input, 1/0.45) * get_multiplier();
#else
		const int table_size = 1024;
		const int table_bits = 10;
		const int interp_bits = 6;
		static int table_is_built = 0;
		static ale_real table[table_size];

		pixel result;

		if (!table_is_built) {
			for (int i = 0; i < table_size; i++) {
				table[i] = pow((float) i / (float) (table_size - 1), 1/0.45);
			}
			table_is_built = 1;
		}

		for (int k = 0; k < 3; k++) {
			/*
			 * Clamp.
			 */
			if (input[k] >= 1) {
				result[k] = 1;
				continue;
			} else if (input[k] <= 0) {
				result[k] = 0;
				continue;
			} else if (isnan(input[k])) {
				result[k] = input[k];
				continue;
			}

			int index1 = ale_real_to_int(input[k], 65535);

			int index2 = index1 >> (16 - table_bits);
			int index3 = (index1 >> (16 - table_bits - interp_bits)) 
			           & ((1 << interp_bits) - 1);

			if (index2 >= table_size - 1) {
				result[k] = 1;
				continue;
			}

			ale_real frac = ale_real_from_int((index3 << (16 - interp_bits)), 65535);

			result[k] = (1 - frac) * table[index2] + frac * table[index2 + 1];
		}

		return result * get_multiplier();
#endif
	}
	pixel unlinearize(pixel input) const {
#if 0
		/*
		 * Calling pow() may be expensive on some platforms (e.g.,
		 * those lacking hardware support for floating point).
		 */

		return ppow(input / get_multiplier(), 0.45);
#else

		input /= get_multiplier();

		const int table_size = 1024;
		const int table_bits = 10;
		const int interp_bits = 6;
		static int table_is_built = 0;
		static ale_real table[table_size];

		pixel result;

		if (!table_is_built) {
			for (int i = 0; i < table_size; i++) {
				table[i] = pow((float) i / (float) (table_size - 1), 0.45);
			}
			table_is_built = 1;
		}

		for (int k = 0; k < 3; k++) {
			/*
			 * Clamp.
			 */
			if (input[k] >= 1) {
				result[k] = 1;
				continue;
			} else if (input[k] <= 0) {
				result[k] = 0;
				continue;
			} else if (isnan(input[k])) {
				result[k] = input[k];
				continue;
			}

			int index1 = ale_real_to_int(input[k], 65535);

			int index2 = index1 >> (16 - table_bits);
			int index3 = (index1 >> (16 - table_bits - interp_bits)) 
			           & ((1 << interp_bits) - 1);

			if (index2 >= table_size - 1) {
				result[k] = 1;
				continue;
			}

			ale_real frac = ale_real_from_int((index3 << (16 - interp_bits)), 65535);

			result[k] = (1 - frac) * table[index2] + frac * table[index2 + 1];
		}

		return result;
#endif
	}
};

#endif
