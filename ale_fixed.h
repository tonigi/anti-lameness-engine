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

/*
 * Define a fixed point data type.
 */

#define ALE_FIXED_POSINF 2147483647
#define ALE_FIXED_NEGINF -2147483646
#define ALE_FIXED_NAN -2147483647

template <int N>
class ale_fixed {
	typedef int i32;
	typedef long long i64;

	static int casting_disabled;

	static void sanity_check_fail() {
		assert(0);
		fprintf(stderr, "\n\nFailed fixed-point sanity test.\n\n");
		exit(1);
	}

public:
	i32 bits;

	/*
	 * Sanity checks.
	 */
	static void sanity_check();

	/*
	 * Constructor.
	 */
	ale_fixed() {
		bits = 0;
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
	 * Cast to ordinary numbers
	 */

	operator double() const {
		assert(!casting_disabled);

		return (((double) bits) * (log(-N) / log(2)));
	}

	/*
	 * Cast from ordinary numbers
	 */

	ale_fixed(double d) {
		assert(!casting_disabled);
		assert(0);
	}

	ale_fixed(int d) {
		assert(!casting_disabled);
		assert(0);
	}

	ale_fixed(unsigned int d) {
		assert(!casting_disabled);
		assert(0);
	}

	/*
	 * Operators.
	 */

	ale_fixed &operator+=(ale_fixed f) {
		assert(0);
		return *this;
	}

	ale_fixed &operator-=(ale_fixed f) {
		assert(0);
		return *this;
	}

	ale_fixed &operator*=(ale_fixed f) {
		assert(0);
		return *this;
	}

	ale_fixed &operator/=(ale_fixed f) {
		assert(0);
		return *this;
	}

	ale_fixed operator-() {
		assert(0);
		return *this;
	}

	ale_fixed operator-(ale_fixed f) {
		assert(0);
		return ale_fixed(0);
	}

	ale_fixed operator-(int i) {
		return operator-(ale_fixed(i));
	}

	ale_fixed operator-(unsigned int i) {
		return operator-(ale_fixed(i));
	}

	ale_fixed operator*(ale_fixed f) {
		assert(0);
		return ale_fixed(0);
	}

	ale_fixed operator*(int i) {
		return operator*(ale_fixed(i));
	}

	ale_fixed operator*(unsigned int i) {
		return operator*(ale_fixed(i));
	}

	ale_fixed operator/(ale_fixed<N> f) {
		assert(0);
		return ale_fixed(0);
	}

	ale_fixed operator/(int i) {
		return operator/(ale_fixed(i));
	}

	ale_fixed operator/(unsigned int i) {
		return operator/(ale_fixed(i));
	}
};

#define ALE_FIXED_INCORPORATE_OPERATOR(op)			\
template<int N>							\
ale_fixed<N> operator op(double a, ale_fixed<N> &f) {		\
	ale_fixed<N> g(a);					\
	return g.operator op(f);				\
}								\
								\
template<int N>							\
ale_fixed<N> operator op(int a, ale_fixed<N> &f) {		\
	return (double) a op f;					\
}								\
								\
template<int N>							\
ale_fixed<N> operator op(unsigned int a, ale_fixed<N> &f) {	\
	return (double) a op f;					\
}								\

ALE_FIXED_INCORPORATE_OPERATOR(-);
ALE_FIXED_INCORPORATE_OPERATOR(*);
ALE_FIXED_INCORPORATE_OPERATOR(/);

template<int N>
ale_fixed<N> fabs(ale_fixed<N> f) {
	if (f < ale_fixed<N>())
		return -f;
	return f;
}

template<int N>
ale_fixed<N> pow(ale_fixed<N> f, double d) {
	if (d == 2) 
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_fixed<N>(1);

	return pow((double) f, d);
}

template<int N>
ale_fixed<N> ceil(ale_fixed<N> f) {
	assert(0);
	return ale_fixed<N>(0);
}

template<int N>
ale_fixed<N> floor(ale_fixed<N> f) {
	assert(0);
	return ale_fixed<N>(0);
}

template<int N, int M>
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

template<int N>
int ale_fixed<N>::casting_disabled = 0;

template<int N>
void ale_fixed<N>::sanity_check() {

	/*
	 * i32 should accept 32-bit integers.
	 */

	i32 test_value = ALE_FIXED_POSINF;

	int count = 0;

	while (test_value /= 2)
		count++;

	if (count != 30)
		sanity_check_fail();

	/*
	 * i32 should be signed.
	 */

	test_value = 0;

	test_value--;

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

	if (count != 61)
		sanity_check_fail();

	/*
	 * i64 should be signed.
	 */

	test_value_2 = 0;

	test_value_2--;

	if (!(test_value_2 < 0))
		sanity_check_fail();

	/*
	 * Addition should work
	 */

	ale_fixed<N> a(10), b(2.5);
	ale_fixed<N> c = a + b;

	if ((double) c <= 12
	 || (double) c >= 13)
		sanity_check_fail();

	/*
	 * Multiplication should work
	 */

	a = 11; b = 2.5;
	c = a * b;

	if ((double) c <= 27
	 || (double) c >= 28)
		sanity_check_fail();

	/*
	 * Division should work.
	 */

	a = 11; b = 2;
	c = a / b;

	if ((double) c <= 5
	 || (double) c >= 6)
		sanity_check_fail();

	/*
	 * Subtraction should work.
	 */

	a = 11; b = 2.5;
	c = a - b;

	if ((double) c <= 8
	 || (double) c >= 9)
		sanity_check_fail();
}


#endif
