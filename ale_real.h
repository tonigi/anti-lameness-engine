// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                                    <dhilvert@ugcs.caltech.edu>

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

#ifndef __ale_real_h__
#define __ale_real_h__

#include "ale_fixed.h"

#define SINGLE 1
#define DOUBLE 2
#define HALF 3
#define FIXED 4

/*
 * Real-valued type used to represent the range of an image (colors, weights,
 * etc.).
 *
 * ale_real is used in computation.
 * ale_sreal is used for storage.
 */

#if ALE_COLORS == SINGLE

typedef float ale_real;
typedef float ale_sreal;

#define ALE_REAL_PRECISION_STRING "SINGLE"

#elif ALE_COLORS == DOUBLE

typedef double ale_real;
typedef double ale_sreal;

#define ALE_REAL_PRECISION_STRING "DOUBLE"

#elif ALE_COLORS == HALF

/*
 * What follows is one approach to packing a floating point
 * number into 16 bits.  This implementation is very slow.
 */

#define MANTISSA_BITS (9)
#define EXP_BITS (15 - MANTISSA_BITS)
#define EXP_SPECIAL (1 << (EXP_BITS - 1))
#define EXP_MAX (EXP_SPECIAL - 1)
#define EXP_MIN (-EXP_MAX)

typedef float ale_real;

class ale_sreal {

	union {
		uint16_t bits;
		struct {
			uint16_t sign:1;
			uint16_t mant:MANTISSA_BITS;
			int16_t  exp :EXP_BITS;
		} fpr;
	} u;
public:
	ale_sreal() {
		u.bits = 0;
	}
	ale_sreal operator=(float v) {
		if (v == 0) {
			u.bits = 0;
		} else if (isnan(v)) {
			u.fpr.exp = EXP_SPECIAL;
			u.fpr.mant = 1;
		} else {

			if (v > 0)
				u.fpr.sign = 0;
			else if (v < 0) {
				u.fpr.sign = 1;
				v = -v;
			} else
				assert(0);

			/*
			 * Get the exponent.
			 */

			int log2 = (int) floor (log(v) / log(2));

			/*
			 * Test the exponent against the largest expressible
			 * exponent for ale_sreal.
			 */

			if (log2 > EXP_MAX) {
				/*
				 * Infinity
				 */

				u.fpr.exp = EXP_SPECIAL;
				u.fpr.mant = 0;

				return *this;
			}

			/*
			 * Test the exponent against the smallest expressible
			 * exponent for ale_sreal.
			 */

			if (log2 < EXP_MIN) {
				/*
				 * Zero
				 */

				u.fpr.exp = 0x0;
				u.fpr.mant = 0;

				return *this;
			}

			/*
			 * The exponent is in range, so use it.
			 */

			u.fpr.exp = log2;

			u.fpr.mant = (uint16_t) floor(v / pow(2, log2) * (1 << (MANTISSA_BITS - 1)));
		}

		return *this;
	}

	operator float() const {
		float result = 3.14159;

		if (((uint16_t) u.fpr.exp == EXP_SPECIAL) && (u.fpr.mant == 1)) {

			/*
			 * NaN
			 */

			float a = 0;
			float b = 0;

			result = a / b;

		} else if (((uint16_t) u.fpr.exp == EXP_SPECIAL) && (u.fpr.mant == 0)) {

			/*
			 * Infinity
			 */

			float a = 1;
			float b = 0;

			result = (a / b);

		} else if ((uint16_t) u.fpr.exp != EXP_SPECIAL) {

			/*
			 * Value is finite.
			 */

			result = u.fpr.mant / ((double) (1 << (MANTISSA_BITS - 1)))
			       * pow(2, u.fpr.exp);

		} else
			assert(0);

		if (u.fpr.sign)
			result = -result;

		return result;
	}

	ale_sreal operator-=(float r) {
		*this = (float) *this - (float) r;
		return *this;
	}
	ale_sreal operator/=(float r) {
		*this = (float) *this / (float) r;
		return *this;
	}
	ale_sreal operator*=(float r) {
		*this = (float) *this * (float) r;
		return *this;
	}
	ale_sreal operator+=(float r) {
		*this = (float) *this + (float) r;
		return *this;
	}
};

#undef MANTISSA_BITS
#undef EXP_BITS
#undef EXP_SPECIAL
#undef EXP_MAX
#undef EXP_MIN

#define ALE_REAL_PRECISION_STRING "HALF"

#elif ALE_COLORS == FIXED

typedef ale_fixed<14> ale_real;
typedef ale_fixed<14> ale_sreal;

#define ALE_REAL_PRECISION_STRING "FIXED"

#else

#warning Unknown precision in ale_real.h: Choosing PRECISION=SINGLE.

typedef float ale_real;
typedef float ale_sreal;

#define ALE_REAL_PRECISION_STRING "SINGLE"

#endif

const ale_real ale_real_0 = (ale_real) 0;

#undef SINGLE
#undef DOUBLE
#undef HALF
#undef FIXED

#endif
