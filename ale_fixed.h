// Copyright 2007 David Hilvert <dhilvert@auricle.dyndns.org>, 
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

#ifndef __ale_fixed_h__
#define __ale_fixed_h__

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "ale_math.h"

#define FIXED16 4
#define FIXED32 5

/*
 * Define a fixed point data type.
 */

class ale_fixed_16 {
public:
	typedef short bits_t;
	typedef int mulbits_t;
	static bits_t posinf() {
		return 32767;
	}
	static bits_t neginf() {
		return -32766;
	}
	static bits_t nan() {
		return -32767;
	}
	static bits_t rint(double d) {
		return (bits_t) lrint(d);
	}
};

class ale_fixed_16_accum {
public:
	typedef int bits_t;
	typedef int mulbits_t;
	static bits_t posinf() {
		return 2147483647;
	}
	static bits_t neginf() {
		return -2147483646;
	}
	static bits_t nan() {
		return -2147483647;
	}
	static bits_t rint(double d) {
		return (bits_t) lrint(d);
	}
};

#if ALE_COLORS == FIXED32 || ALE_COORDINATES == FIXED32
class ale_fixed_32 {
public:
	typedef int bits_t;
	typedef long long mulbits_t;
	static bits_t posinf() {
		return 2147483647;
	}
	static bits_t neginf() {
		return -2147483646;
	}
	static bits_t nan() {
		return -2147483647;
	}
	static bits_t rint(double d) {
		return (bits_t) lrint(d);
	}
};
#endif

#if ALE_COLORS == FIXED32
class ale_fixed_32_accum {
public:
	typedef long long bits_t;
	typedef long long mulbits_t;
	static bits_t posinf() {
		return 9223372036854775807LL;
	}
	static bits_t neginf() {
		return -9223372036854775806LL;
	}
	static bits_t nan() {
		return -9223372036854775807LL;
	}
	static bits_t rint(double d) {
		return (bits_t) llrint(d);
	}
};
#endif

#define ALE_FIXED_NAN (fixed_type::nan())
#define ALE_FIXED_POSINF (fixed_type::posinf())
#define ALE_FIXED_NEGINF (fixed_type::neginf())

template<class fixed_type, unsigned int N>
class ale_fixed {
	static int casting_disabled;

public:

	typedef typename fixed_type::bits_t bits_t;
	typedef typename fixed_type::mulbits_t mulbits_t;
	
	bits_t bits;

	/*
	 * Constructors.
	 */
	ale_fixed() {
		bits = 0;
	}

	ale_fixed(const ale_fixed &f) {
		bits = f.bits;
	}

	ale_fixed& operator=(const ale_fixed &f) {
		bits = f.bits;
		
		return (*this);
	}

	/*
	 * Disable casting
	 */

	static void disable_casting() {
		casting_disabled = 1;
	}

	/*
	 * Enable casting
	 */

	static void enable_casting() {
		casting_disabled = 0;
	}

	/*
	 * Casting status.
	 */

	static int casting_status() {
		return !casting_disabled;
	}

	/*
	 * Cast to ordinary numbers
	 */

	operator double() const {
#if 1
		/*
		 * Removed for performance reasons.
		 */

		assert(!casting_disabled);
		if (bits == ALE_FIXED_NAN) {
			double zero = 0;
			double nan = zero / zero;
			assert (isnan(nan));
			return nan;
		} else if (bits == ALE_FIXED_NEGINF) {
			double zero = 0;
			double negone = -1;
			double neginf = negone / zero;
			assert (isinf(neginf));
			assert (neginf < 0);
			return neginf;
		} else if (bits == ALE_FIXED_POSINF) {
			double zero = 0;
			double posone = +1;
			double posinf = posone / zero;
			assert (isinf(posinf));
			assert (posinf > 0);
			return posinf;
		}
#endif

		return (((double) bits) / (1 << N));
	}

	operator float() const {
#if 1
		/*
		 * Removed for performance reasons.
		 */

		assert(!casting_disabled);
		if (bits == ALE_FIXED_NAN) {
			float zero = 0;
			float nan = zero / zero;
			assert (isnan(nan));
			return nan;
		} else if (bits == ALE_FIXED_NEGINF) {
			float zero = 0;
			float negone = -1;
			float neginf = negone / zero;
			assert (isinf(neginf));
			assert (neginf < 0);
			return neginf;
		} else if (bits == ALE_FIXED_POSINF) {
			float zero = 0;
			float posone = +1;
			float posinf = posone / zero;
			assert (isinf(posinf));
			assert (posinf > 0);
			return posinf;
		}
#endif

		return (((float) bits) / (1 << N));
	}

	operator int() const {
#if 1
		/*
		 * Removed for performance reasons.
		 */

		assert (bits != ALE_FIXED_NAN);
		assert (bits != ALE_FIXED_POSINF);
		assert (bits != ALE_FIXED_NEGINF);
#endif

		return bits / (1 << N);
	}

	operator unsigned int() const {
#if 1
		/*
		 * Removed for performance reasons.
		 */

		assert (bits != ALE_FIXED_NAN);
		assert (bits != ALE_FIXED_POSINF);
		assert (bits != ALE_FIXED_NEGINF);
		assert (bits >= 0);
#endif

		return (unsigned int) operator int();
	}

#if 0
	template<class fixed_type_2, unsigned int M>
	operator ale_fixed<fixed_type_2, M>() const {
		ale_fixed<fixed_type_2, M> result;

		if (bits == ALE_FIXED_NAN) {
			result.bits = fixed_type_2::nan();
			return result;
		}

		if (bits == ALE_FIXED_POSINF) {
			result.bits = fixed_type_2::posinf();
			return result;
		}

		if (bits == ALE_FIXED_NEGINF) {
			result.bits = fixed_type_2::neginf();
			return result;
		}

		if (sizeof(ale_fixed<fixed_type_2,M>) > sizeof(ale_fixed<fixed_type,N>)) {
			typedef typename fixed_type_2::bits_t bits_t_calc;

			bits_t_calc type_result;

			if (M >= N)
				type_result = bits << (bits_t_calc) ((int) M - (int) N);
			else
				type_result = bits / ((bits_t_calc) 1 << (bits_t_calc) ((int) N - (int) M));

			result.bits = type_result;

		} else {
			typedef bits_t bits_t_calc;

			bits_t_calc type_result;

			if (M >= N)
				type_result = bits << (bits_t_calc) ((int) M - (int) N);
			else
				type_result = bits / ((bits_t_calc) 1 << (bits_t_calc) ((int) N - (int) M));

			if (type_result > fixed_type_2::posinf())
				result.bits = fixed_type_2::posinf();
			else if (type_result < fixed_type_2::neginf())
				result.bits = fixed_type_2::neginf();
			else
				result.bits = type_result;
		}

		return result;
	}
#endif

	/*
	 * Cast from ordinary numbers
	 */

	template<class fixed_type_2, unsigned int M>
	ale_fixed(const ale_fixed<fixed_type_2,M> &d) {

		/*
		 * XXX: this shouldn't be necessary.
		 */

		bits = 0;

		if (d.bits == fixed_type_2::nan()) {
			bits = ALE_FIXED_NAN;
			return;
		}

		if (bits == fixed_type_2::posinf()) {
			bits = ALE_FIXED_POSINF;
			return;
		}

		if (bits == fixed_type_2::neginf()) {
			bits = ALE_FIXED_NEGINF;
			return;
		}

		if (sizeof(ale_fixed<fixed_type,N>) > sizeof(ale_fixed<fixed_type_2,M>)) {
			if (N >= M)
				bits = d.bits << (bits_t) ((int) N - (int) M);
			else
				bits = d.bits / ((bits_t) 1 << (bits_t) ((int) M - (int) N));
		} else {
			typedef typename ale_fixed<fixed_type_2,M>::bits_t bits_t_calc;

			bits_t_calc type_result;

			if (N >= M)
				type_result = d.bits << (bits_t_calc) ((int) N - (int) M);
			else
				type_result = d.bits / ((bits_t_calc) 1 << (bits_t_calc) ((int) M - (int) N));

			if (type_result > ALE_FIXED_POSINF)
				bits = ALE_FIXED_POSINF;
			else if (type_result < ALE_FIXED_NEGINF)
				bits = ALE_FIXED_NEGINF;
			else
				bits = type_result;
		}
	}

	ale_fixed(double d) {
		assert(!casting_disabled);

		if (isnan(d)) {
			bits = ALE_FIXED_NAN;
		} else if (isinf(d) && d > 0) {
			bits = ALE_FIXED_POSINF;
		} else if (isinf(d) && d < 0) {
			bits = ALE_FIXED_NEGINF;
		} else {
			bits = (bits_t) fixed_type::rint(d * (1 << N));

#if 1
			/*
			 * Removed for performance reasons.
			 */

			assert((double) *this > (d - (double) 1 / (1 << N)));
			assert((double) *this < (d + (double) 1 / (1 << N)));

			assert(bits < ALE_FIXED_POSINF);
			assert(bits > ALE_FIXED_NEGINF);
#endif
		}
	}

	ale_fixed(int d) {
		bits = (bits_t) d << N;
#if 1
		/*
		 * Removed for performance reasons.
		 */

		assert((d >= 0 && bits >> N == d)
		    || (d < 0 && (-bits) >> N == -d));

		assert (bits < ALE_FIXED_POSINF);
		assert (bits > ALE_FIXED_NEGINF);
#endif
	}

	ale_fixed(unsigned int d) {
		bits = (bits_t) d << N;

		assert((unsigned int) (bits >> N) == d);

		assert (bits < ALE_FIXED_POSINF);
		assert (bits > ALE_FIXED_NEGINF);
	}

	/*
	 * Operators.
	 */

	ale_fixed operator-() const {

		ale_fixed result;

		if (bits == ALE_FIXED_NAN || bits == 0)
			return *this;
		else if (bits == ALE_FIXED_POSINF)
			result.bits = ALE_FIXED_NEGINF;
		else if (bits == ALE_FIXED_NEGINF)
			result.bits = ALE_FIXED_POSINF;
		else
			result.bits = -bits;

		return result;
	}

	ale_fixed operator+(ale_fixed f) const {
		ale_fixed result;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN
		 || (bits == ALE_FIXED_POSINF && f.bits == ALE_FIXED_NEGINF)
		 || (bits == ALE_FIXED_NEGINF && f.bits == ALE_FIXED_POSINF)) {
			result.bits = ALE_FIXED_NAN;
			return result;
		}
#endif

		bits_t bits_result = bits + f.bits;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (bits_result >= ALE_FIXED_POSINF
		 || bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_POSINF
		 || bits > 0 && f.bits > 0 && bits_result < 0) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (bits_result <= ALE_FIXED_NEGINF
		      || bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_NEGINF
		      || bits < 0 && f.bits < 0 && bits_result > 0) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		}
#endif

		result.bits = bits_result;

		return result;
	}

	ale_fixed operator+(int i) const {
		return operator+(ale_fixed(i));
	}

	ale_fixed operator+(unsigned int i) const {
		return operator+(ale_fixed(i));
	}

	ale_fixed operator-(ale_fixed f) const {
		ale_fixed result;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN
		 || (bits == ALE_FIXED_POSINF && f.bits == ALE_FIXED_POSINF)
		 || (bits == ALE_FIXED_NEGINF && f.bits == ALE_FIXED_NEGINF)) {
			result.bits = ALE_FIXED_NAN;
			return result;
		}
#endif

		bits_t bits_result = bits - f.bits;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (bits_result >= ALE_FIXED_POSINF
		 || bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_NEGINF
		 || bits > 0 && f.bits < 0 && bits_result < 0) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (bits_result <= ALE_FIXED_NEGINF
		      || bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_POSINF
		      || bits < 0 && f.bits > 0 && bits_result > 0) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		}
#endif

		result.bits = bits_result;

		return result;
	}

	ale_fixed operator-(int i) const {
		return operator-(ale_fixed(i));
	}

	ale_fixed operator-(unsigned int i) const {
		return operator-(ale_fixed(i));
	}

	ale_fixed operator*(ale_fixed f) const {
		ale_fixed result;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN) {
			result.bits = ALE_FIXED_NAN;
			return result;
		}
#endif 
			

		mulbits_t mul_result = ((mulbits_t) bits * (mulbits_t) f.bits) / (1 << N);

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (mul_result > (mulbits_t) ALE_FIXED_POSINF 
		 || mul_result < (mulbits_t) ALE_FIXED_NEGINF
		 || bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_POSINF
		 || bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_NEGINF) {
			if (mul_result > 0)
				result.bits = ALE_FIXED_POSINF;
			else if (mul_result < 0)
				result.bits = ALE_FIXED_NEGINF;
			else if (mul_result == 0)
				result.bits = ALE_FIXED_NAN;
			else
				assert(0);
			return result;
		}
#endif

		result.bits = mul_result;
		return result;
	}

	ale_fixed operator*(int i) const {
		return operator*(ale_fixed(i));
	}

	ale_fixed operator*(unsigned int i) const {
		return operator*(ale_fixed(i));
	}

	ale_fixed operator/(ale_fixed f) const {
		ale_fixed result;

		/*
		 * While this approach may not be suitable for all
		 * applications, it can be a convenient way to detect and
		 * manufacture non-finite values.
		 */
		if ((bits == 0 && f.bits == 0)
#if 1
		/*
		 * Removed for performance reasons.
		 */

		 || bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN
		 || ((bits == ALE_FIXED_NEGINF || bits == ALE_FIXED_POSINF)
		  && (f.bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_POSINF))
#endif
		) {
			result.bits = ALE_FIXED_NAN;
			return result;
		} else if (f.bits == 0 && bits > 0) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (f.bits == 0 && bits < 0) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		} 
		
#if 1
		/*
		 * Removed for performance reasons.
		 */

		else if (f.bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_NEGINF) {
			result.bits = 0;
			return result;
		} 
#endif
			
		mulbits_t div_result = ((mulbits_t) bits << N) / f.bits;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (div_result > (mulbits_t) ALE_FIXED_POSINF) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (div_result < (mulbits_t) ALE_FIXED_NEGINF) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		}
#endif

		result.bits = (bits_t) div_result;
		return result;
	}

	ale_fixed operator/(int i) const {
		return operator/(ale_fixed(i));
	}

	ale_fixed operator/(unsigned int i) const {
		return operator/(ale_fixed(i));
	}

	ale_fixed &operator+=(ale_fixed f) {
		*this = *this + f;
		return *this;
	}

	ale_fixed &operator-=(ale_fixed f) {
		*this = *this - f;
		return *this;
	}

	ale_fixed &operator*=(ale_fixed f) {
		*this = *this * f;
		return *this;
	}

	ale_fixed &operator/=(ale_fixed f) {
		*this = *this / f;
		return *this;
	}

	int operator!=(ale_fixed f) const {
		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN)
			return 1;

		if (bits == f.bits)
			return 0;

		return 1;
	}

	int operator==(ale_fixed f) const {
		return !(operator!=(f));
	}

	int operator<=(ale_fixed f) const {
		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN)
			return 0;

		if (bits <= f.bits)
			return 1;

		return 0;
	}

	int operator>=(ale_fixed f) const {
		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN)
			return 0;

		if (bits >= f.bits)
			return 1;

		return 0;
	}

	int operator>(ale_fixed f) const {
		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN)
			return 0;

		if (bits > f.bits)
			return 1;

		return 0;
	}

	int operator<(ale_fixed f) const {
		if (bits == ALE_FIXED_NAN || f.bits == ALE_FIXED_NAN)
			return 0;

		if (bits < f.bits)
			return 1;

		return 0;
	}

	int operator>=(int d) const {
		return operator>=((ale_fixed) d);
	}

	int operator<=(int d) const {
		return operator<=((ale_fixed) d);
	}

	int operator==(int d) const {
		return operator==((ale_fixed) d);
	}

	int operator!=(int d) const {
		return operator!=((ale_fixed) d);
	}

	int operator>(int d) const {
		return operator>((ale_fixed) d);
	}

	int operator<(int d) const {
		return operator<((ale_fixed) d);
	}

	int operator>=(double d) const {
		return operator>=((ale_fixed) d);
	}

	int operator>=(float d) const {
		return operator>=((ale_fixed) d);
	}

	int operator<=(double d) const {
		return operator<=((ale_fixed) d);
	}

	int operator==(double d) const {
		return operator==((ale_fixed) d);
	}

	int operator!=(double d) const {
		return operator!=((ale_fixed) d);
	}

	int operator>(double d) const {
		return operator>((ale_fixed) d);
	}

	int operator<(double d) const {
		return operator<((ale_fixed) d);
	}

	int operator>=(unsigned int d) const {
		return operator>=((ale_fixed) d);
	}

	int operator<=(unsigned int d) const {
		return operator<=((ale_fixed) d);
	}

	int operator==(unsigned int d) const {
		return operator==((ale_fixed) d);
	}

	int operator!=(unsigned int d) const {
		return operator!=((ale_fixed) d);
	}

	int operator>(unsigned int d) const {
		return operator>((ale_fixed) d);
	}

	int operator<(unsigned int d) const {
		return operator<((ale_fixed) d);
	}

};

#define ALE_FIXED_INCORPORATE_OPERATOR(return_value, op)			\
template<class fixed_type, unsigned int N>					\
return_value operator op(double a, const ale_fixed<fixed_type, N> &f) {		\
	ale_fixed<fixed_type, N> g(a);					\
	return g.operator op(f);				\
}								\
								\
template<class fixed_type, unsigned int N>					\
return_value operator op(int a, const ale_fixed<fixed_type, N> &f) {		\
	return (ale_fixed<fixed_type, N>) a op f;					\
}								\
								\
template<class fixed_type, unsigned int N>					\
return_value operator op(unsigned int a, const ale_fixed<fixed_type, N> &f) {	\
	return (ale_fixed<fixed_type, N>) a op f;					\
}								\

#define STDARGS ale_fixed<fixed_type,N>

ALE_FIXED_INCORPORATE_OPERATOR(STDARGS, +);
ALE_FIXED_INCORPORATE_OPERATOR(STDARGS, -);
ALE_FIXED_INCORPORATE_OPERATOR(STDARGS, *);
ALE_FIXED_INCORPORATE_OPERATOR(STDARGS, /);
ALE_FIXED_INCORPORATE_OPERATOR(int, <=);
ALE_FIXED_INCORPORATE_OPERATOR(int, >=);
ALE_FIXED_INCORPORATE_OPERATOR(int, <);
ALE_FIXED_INCORPORATE_OPERATOR(int, >);
ALE_FIXED_INCORPORATE_OPERATOR(int, !=);
ALE_FIXED_INCORPORATE_OPERATOR(int, ==);

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> fabs(ale_fixed<fixed_type, N> f) {
	if (f < ale_fixed<fixed_type, N>())
		return -f;
	return f;
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> pow(ale_fixed<fixed_type, N> f, double d) {
	return pow((double) f, (double) d);
}

/*
 * sqrt() via the Babylonian method.
 *
 * http://en.wikipedia.org/wiki/Methods_of_computing_square_roots
 */

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> sqrt(ale_fixed<fixed_type, N> f) {
	ale_fixed<fixed_type, N> guess = f;

	for (int i = 0; i < 5; i++) {
		guess.bits >>= 1;

		if (guess.bits <= 0)
			return 0;

		long long sf = (long long) f.bits << (N - 2);
		guess.bits = guess.bits + sf / guess.bits;
	}

	return guess;
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> pow(ale_fixed<fixed_type, N> f, ale_fixed<fixed_type, N> d) {
	if (d == 2) 
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<fixed_type,N>(1);

	return pow((double) f, (double) d);
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> pow(ale_fixed<fixed_type, N> f, int d) {
	if (d == 2)
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<fixed_type, N>(1);

	if (d > 1)
		return pow(f, d / 2) * pow(f, d - d / 2);

	if (d < 0)
		return 1 / pow(f, -d);

	assert(0);
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> pow(ale_fixed<fixed_type, N> f, unsigned int d) {
	if (d == 2)
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<fixed_type, N>(1);

	return pow(f, d / 2) * pow(f, d - d / 2);
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> floor(ale_fixed<fixed_type, N> f) {

#if 1
		/*
		 * Removed for performance reasons.
		 */

	if (N == 0
	 || f.bits == ALE_FIXED_POSINF
	 || f.bits == ALE_FIXED_NEGINF
	 || f.bits == ALE_FIXED_NAN)
		return f;
#endif

	ale_fixed<fixed_type, N> result;

	result.bits = (f.bits & ~((1 << N) - 1));

	/*
	 * XXX: This isn't exactly right.
	 */
	if (f.bits < 0)
		result.bits -= 1;

	return result;
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> lrintf(ale_fixed<fixed_type, N> f) {

#if 1
		/*
		 * Removed for performance reasons.
		 */

	if (N == 0
	 || f.bits == ALE_FIXED_POSINF
	 || f.bits == ALE_FIXED_NEGINF
	 || f.bits == ALE_FIXED_NAN)
		return f;
#endif

	ale_fixed<fixed_type, N> result = floor(f);

	if (f.bits - result.bits >= (1 << N - 1))
		result.bits += (1 << N);

	return result;
}

template<class fixed_type, unsigned int N>
ale_fixed<fixed_type, N> ceil(ale_fixed<fixed_type, N> f) {
	return -floor(-f);
}

template<class fixed_type, unsigned int N>
int ale_isinf(ale_fixed<fixed_type, N> f) {
	return (f.bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_POSINF);
}

template<class fixed_type, unsigned int N>
int ale_isnan(ale_fixed<fixed_type, N> f) {
	return (f.bits == ALE_FIXED_NAN);
}

template<class fixed_type, unsigned int N>
int finite(ale_fixed<fixed_type, N> f) {
	return (f.bits < ALE_FIXED_POSINF && f.bits > ALE_FIXED_NEGINF);
}

template<class fixed_type, unsigned int N, unsigned int M>
ale_fixed<fixed_type, N> convert_precision(ale_fixed<fixed_type, M> m) {

	/*
	 * XXX: Checks should be added that precision is not
	 * lost from most-significant bits.
	 */

	if (N != M) 
		assert (0);
	
	ale_fixed<fixed_type, N> n;

	n.bits = m.bits << (N - M);

	return n;
}

template<class fixed_type, unsigned int N>
int ale_fixed<fixed_type, N>::casting_disabled = 0;

#undef FIXED16
#undef FIXED32

#endif
