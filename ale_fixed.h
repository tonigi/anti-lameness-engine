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

#define FIXED 4

/*
 * Define a fixed point data type.
 */

#define ALE_FIXED_POSINF 2147483647
#define ALE_FIXED_NEGINF -2147483646
#define ALE_FIXED_NAN -2147483647

template <unsigned int N>
class ale_fixed {
	static int casting_disabled;

	static void sanity_check_fail() {
		assert(0);
		fprintf(stderr, "\n\nFailed fixed-point sanity test.\n\n");
		exit(1);
	}

public:
	typedef int i32;
	typedef long long i64;

	i32 bits;

	/*
	 * Sanity checks.
	 */
	static void sanity_check();

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

	/*
	 * Cast from ordinary numbers
	 */

	ale_fixed(double d) {
		assert(!casting_disabled);

		if (isnan(d)) {
			bits = ALE_FIXED_NAN;
		} else if (isinf(d) && d > 0) {
			bits = ALE_FIXED_POSINF;
		} else if (isinf(d) && d < 0) {
			bits = ALE_FIXED_NEGINF;
		} else {
			bits = (int) lrint(d * (1 << N));

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
		bits = d << N;
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
		bits = d << N;

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

		i32 i32_result = bits + f.bits;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (i32_result >= ALE_FIXED_POSINF
		 || bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_POSINF
		 || bits > 0 && f.bits > 0 && i32_result < 0) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (i32_result <= ALE_FIXED_NEGINF
		      || bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_NEGINF
		      || bits < 0 && f.bits < 0 && i32_result > 0) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		}
#endif

		result.bits = i32_result;

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

		i32 i32_result = bits - f.bits;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (i32_result >= ALE_FIXED_POSINF
		 || bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_NEGINF
		 || bits > 0 && f.bits < 0 && i32_result < 0) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (i32_result <= ALE_FIXED_NEGINF
		      || bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_POSINF
		      || bits < 0 && f.bits > 0 && i32_result > 0) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		}
#endif

		result.bits = i32_result;

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
			

		i64 i64_result = ((i64) bits * (i64) f.bits) / (1 << N);

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (i64_result > (i64) ALE_FIXED_POSINF 
		 || i64_result < (i64) ALE_FIXED_NEGINF
		 || bits == ALE_FIXED_POSINF || f.bits == ALE_FIXED_POSINF
		 || bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_NEGINF) {
			if (i64_result > 0)
				result.bits = ALE_FIXED_POSINF;
			else if (i64_result < 0)
				result.bits = ALE_FIXED_NEGINF;
			else if (i64_result == 0)
				result.bits = ALE_FIXED_NAN;
			else
				assert(0);
			return result;
		}
#endif

		result.bits = i64_result;
		return result;
	}

	ale_fixed operator*(int i) const {
		return operator*(ale_fixed(i));
	}

	ale_fixed operator*(unsigned int i) const {
		return operator*(ale_fixed(i));
	}

	ale_fixed operator/(ale_fixed<N> f) const {
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
			
		i64 i64_result = ((i64) bits << N) / f.bits;

#if 1
		/*
		 * Removed for performance reasons.
		 */

		if (i64_result > (i64) ALE_FIXED_POSINF) {
			result.bits = ALE_FIXED_POSINF;
			return result;
		} else if (i64_result < (i64) ALE_FIXED_NEGINF) {
			result.bits = ALE_FIXED_NEGINF;
			return result;
		}
#endif

		result.bits = (i32) i64_result;
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
template<unsigned int N>					\
return_value operator op(double a, const ale_fixed<N> &f) {		\
	ale_fixed<N> g(a);					\
	return g.operator op(f);				\
}								\
								\
template<unsigned int N>					\
return_value operator op(int a, const ale_fixed<N> &f) {		\
	return (ale_fixed<N>) a op f;					\
}								\
								\
template<unsigned int N>					\
return_value operator op(unsigned int a, const ale_fixed<N> &f) {	\
	return (ale_fixed<N>) a op f;					\
}								\

ALE_FIXED_INCORPORATE_OPERATOR(ale_fixed<N>, +);
ALE_FIXED_INCORPORATE_OPERATOR(ale_fixed<N>, -);
ALE_FIXED_INCORPORATE_OPERATOR(ale_fixed<N>, *);
ALE_FIXED_INCORPORATE_OPERATOR(ale_fixed<N>, /);
ALE_FIXED_INCORPORATE_OPERATOR(int, <=);
ALE_FIXED_INCORPORATE_OPERATOR(int, >=);
ALE_FIXED_INCORPORATE_OPERATOR(int, <);
ALE_FIXED_INCORPORATE_OPERATOR(int, >);
ALE_FIXED_INCORPORATE_OPERATOR(int, !=);
ALE_FIXED_INCORPORATE_OPERATOR(int, ==);

template<unsigned int N>
ale_fixed<N> fabs(ale_fixed<N> f) {
	if (f < ale_fixed<N>())
		return -f;
	return f;
}

template<unsigned int N>
ale_fixed<N> pow(ale_fixed<N> f, double d) {
	return pow((double) f, (double) d);
}

/*
 * sqrt() via the Babylonian method.
 *
 * http://en.wikipedia.org/wiki/Methods_of_computing_square_roots
 */

template<unsigned int N>
ale_fixed<N> sqrt(ale_fixed<N> f) {
	ale_fixed<N> guess = f;

	for (int i = 0; i < 5; i++) {
		guess.bits >>= 1;

		if (guess.bits <= 0)
			return 0;

		long long sf = (long long) f.bits << (N - 2);
		guess.bits = guess.bits + sf / guess.bits;
	}

	return guess;
}

template<unsigned int N>
ale_fixed<N> pow(ale_fixed<N> f, ale_fixed<N> d) {
	if (d == 2) 
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<N>(1);

	return pow((double) f, (double) d);
}

template<unsigned int N>
ale_fixed<N> pow(ale_fixed<N> f, int d) {
	if (d == 2)
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<N>(1);

	if (d > 1)
		return pow(f, d / 2) * pow(f, d - d / 2);

	if (d < 0)
		return 1 / pow(f, -d);

	assert(0);
}

template<unsigned int N>
ale_fixed<N> pow(ale_fixed<N> f, unsigned int d) {
	if (d == 2)
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<N>(1);

	return pow(f, d / 2) * pow(f, d - d / 2);
}

template<unsigned int N>
ale_fixed<N> floor(ale_fixed<N> f) {

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

	ale_fixed<N> result;

	result.bits = (f.bits & ~((1 << N) - 1));

	/*
	 * XXX: This isn't exactly right.
	 */
	if (f.bits < 0)
		result.bits -= 1;

	return result;
}

template<unsigned int N>
ale_fixed<N> lrintf(ale_fixed<N> f) {

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

	ale_fixed<N> result = floor(f);

	if (f.bits - result.bits >= (1 << N - 1))
		result.bits += (1 << N);

	return result;
}

template<unsigned int N>
ale_fixed<N> ceil(ale_fixed<N> f) {
	return -floor(-f);
}

#if ALE_COLORS == FIXED || ALE_COORDINATES == FIXED

template<unsigned int N>
int ale_isinf(ale_fixed<N> f) {
	return (f.bits == ALE_FIXED_NEGINF || f.bits == ALE_FIXED_POSINF);
}

template<unsigned int N>
int ale_isnan(ale_fixed<N> f) {
	return (f.bits == ALE_FIXED_NAN);
}

#endif

template<unsigned int N>
int finite(ale_fixed<N> f) {
	return (f.bits < ALE_FIXED_POSINF && f.bits > ALE_FIXED_NEGINF);
}

template<unsigned int N, int M>
ale_fixed<N> convert_precision(ale_fixed<M> m) {

	/*
	 * XXX: Checks should be added that precision is not
	 * lost from most-significant bits.
	 */

	if (N != M) 
		assert (0);
	
	ale_fixed<N> n;

	n.bits = m.bits << (N - M);

	return n;
}

template<unsigned int N>
int ale_fixed<N>::casting_disabled = 0;

template<unsigned int N>
void ale_fixed<N>::sanity_check() {

	/*
	 * i32 should accept 32-bit integers.
	 */

	i32 test_value = ALE_FIXED_POSINF;

	int count = 0;

	while (test_value /= 2)
		count++;

	assert (count == 30);
	if (count != 30)
		sanity_check_fail();

	/*
	 * i32 should be signed.
	 */

	test_value = 0;

	test_value--;

	assert (test_value < 0);
	if (!(test_value < 0))
		sanity_check_fail();

	/*
	 * i64 should accept 64-bit integers.
	 */

	i64 test_value_2 = ALE_FIXED_POSINF;

	test_value_2 *= test_value_2;

	count = 0;

	while (test_value_2 /= 2)
		count++;

	assert (count == 61);
	if (count != 61)
		sanity_check_fail();

	/*
	 * i64 should be signed.
	 */

	test_value_2 = 0;

	test_value_2--;

	assert (test_value_2 < 0);
	if (!(test_value_2 < 0))
		sanity_check_fail();

	/*
	 * Addition should work
	 */

	ale_fixed<N> a(10), b(2.5);
	ale_fixed<N> c = a + b;

	assert ((double) c > 12 && (double) c < 13);
	if (!((double) c > 12 && (double) c < 13))
		sanity_check_fail();

	/*
	 * Multiplication should work
	 */

	a = 11; b = 2.5;
	c = a * b;

	assert ((double) c > 27 && (double) c < 28);
	if (!((double) c > 27 && (double) c < 28))
		sanity_check_fail();

	/*
	 * Division should work.
	 */

	a = 11; b = 2;
	c = a / b;

	assert ((double) c > 5 && (double) c < 6);
	if (!((double) c > 5 && (double) c < 6))
		sanity_check_fail();

	/*
	 * Subtraction should work.
	 */

	a = 11; b = 2.5;
	c = a - b;

	assert ((double) c > 8 && (double) c < 9);
	if (!((double) c > 8 && (double) c < 9))
		sanity_check_fail();

	/*
	 * Infinite values should work.
	 */

	a = +1; b = 0;
	c = a / b;

	if (!(c > 0 && c > a && c > b && b < c && a < c)) {
		assert(0);
		sanity_check_fail();
	}
}

#undef FIXED

#endif
