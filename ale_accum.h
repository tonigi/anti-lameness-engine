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

#ifndef __ale_accum_h__
#define __ale_accum_h__

#include "ale_fixed.h"

#define FIXED 4

#define ale_accum_disable_casting()
#define ale_accum_enable_casting()

/*
 * Real-valued type used when accumulating over the domain of an image.
 */

#if ALE_COLORS == FIXED

#define ALE_ACCUM_FIXED_POSINF 9223372036854775807LL
#define ALE_ACCUM_FIXED_NEGINF -9223372036854775806LL
#define ALE_ACCUM_FIXED_NAN -9223372036854775807LL

template <unsigned int N>
class ale_accum_fixed {
	static int casting_disabled;

public:
	typedef int i32;
	typedef long long i64;
	typedef unsigned long long ui64;

	i64 bits;

	static const ale_accum_fixed const_0_001;

	/*
	 * Constructor.
	 */
	ale_accum_fixed() {
		bits = 0;
	}

	ale_accum_fixed(const ale_accum_fixed &f) {
		bits = f.bits;
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
		if (bits == ALE_ACCUM_FIXED_NAN) {
			double zero = 0;
			double nan = zero / zero;
			assert (isnan(nan));
			return nan;
		} else if (bits == ALE_ACCUM_FIXED_NEGINF) {
			double zero = 0;
			double negone = -1;
			double neginf = negone / zero;
			assert (isinf(neginf));
			assert (neginf < 0);
			return neginf;
		} else if (bits == ALE_ACCUM_FIXED_POSINF) {
			double zero = 0;
			double posone = +1;
			double posinf = posone / zero;
			assert (isinf(posinf));
			assert (posinf > 0);
			return posinf;
		} else 
			return (((double) bits) / (1 << N));
	}

	template<unsigned int M>
	operator ale_fixed<M>() {
		ale_fixed<M> result;
		i64 i64_result;

		if (bits == ALE_ACCUM_FIXED_NAN)
			result.bits = ALE_FIXED_NAN;
		else if (bits == ALE_ACCUM_FIXED_POSINF)
			result.bits = ALE_FIXED_POSINF;
		else if (bits == ALE_ACCUM_FIXED_NEGINF)
			result.bits = ALE_FIXED_NEGINF;
		else {

			if (bits >= 0 || M >= N)
				i64_result = bits << (i64) (M - N);
			else
				i64_result = -((-bits) >> (i64) (N - M));

			if (i64_result > ALE_FIXED_POSINF)
				result.bits = ALE_FIXED_POSINF;
			else if (i64_result < ALE_FIXED_NEGINF)
				result.bits = ALE_FIXED_NEGINF;
			else
				result.bits = i64_result;
		}

		return result;
	}

	/*
	 * Cast from ordinary numbers
	 */

	ale_accum_fixed(double d) {
		assert(!casting_disabled);

		if (isnan(d)) {
			bits = ALE_ACCUM_FIXED_NAN;
		} else if (isinf(d) && d > 0) {
			bits = ALE_ACCUM_FIXED_POSINF;
		} else if (isinf(d) && d < 0) {
			bits = ALE_ACCUM_FIXED_NEGINF;
		} else {
			bits = (i64) (d * (1 << N));

			assert((double) *this > (d - (double) 1 / (1 << N)));
			assert((double) *this < (d + (double) 1 / (1 << N)));

			assert(bits < ALE_ACCUM_FIXED_POSINF);
			assert(bits > ALE_ACCUM_FIXED_NEGINF);
		}
	}

	ale_accum_fixed(i64 d) {
		assert(!casting_disabled);

		bits = d << N;
		
		assert((d >= 0 && bits >> N == d)
		    || (d < 0 && (-bits) >> N == -d));

		assert (bits < ALE_ACCUM_FIXED_POSINF);
		assert (bits > ALE_ACCUM_FIXED_NEGINF);
	}

	ale_accum_fixed(ui64 d) {
		assert(!casting_disabled);

		bits = d << N;

		assert((ui64) (bits >> N) == d);

		assert (bits < ALE_ACCUM_FIXED_POSINF);
		assert (bits > ALE_ACCUM_FIXED_NEGINF);
	}

	ale_accum_fixed(int i) {
		i64 d = i;

		bits = d << N;
		
#if 0
		assert((d >= 0 && bits >> N == d)
		    || (d < 0 && (-bits) >> N == -d));

		assert (bits < ALE_ACCUM_FIXED_POSINF);
		assert (bits > ALE_ACCUM_FIXED_NEGINF);
#endif
	}
		
	ale_accum_fixed(unsigned int u) {
		ui64 d = u;

		bits = d << N;

		assert((ui64) (bits >> N) == d);

		assert (bits < ALE_ACCUM_FIXED_POSINF);
		assert (bits > ALE_ACCUM_FIXED_NEGINF);
	}


	template<unsigned int M>
	ale_accum_fixed(ale_fixed<M> d) {
		if (d.bits == ALE_FIXED_POSINF)
			bits = ALE_ACCUM_FIXED_POSINF;
		else if (d.bits == ALE_FIXED_NEGINF)
			bits = ALE_ACCUM_FIXED_NEGINF;
		else if (d.bits == ALE_FIXED_NAN)
			bits = ALE_ACCUM_FIXED_NAN;
		else if (d.bits >= 0 || N >= M)
			bits = ((i64) d.bits) << (i64) (N - M);
		else
			bits = -((-(i64) d.bits) >> (i64) (M - N));
	}

	/*
	 * Operators.
	 */

	ale_accum_fixed operator+(ale_accum_fixed f) const {
		ale_accum_fixed result;

		if (bits == ALE_ACCUM_FIXED_NAN || f.bits == ALE_ACCUM_FIXED_NAN
		 || (bits == ALE_ACCUM_FIXED_POSINF && f.bits == ALE_ACCUM_FIXED_NEGINF)
		 || (bits == ALE_ACCUM_FIXED_NEGINF && f.bits == ALE_ACCUM_FIXED_POSINF)) {
			result.bits = ALE_ACCUM_FIXED_NAN;
			return result;
		}

		i64 i64_result = bits + f.bits;

		if (i64_result >= ALE_ACCUM_FIXED_POSINF
		 || bits == ALE_ACCUM_FIXED_POSINF || f.bits == ALE_ACCUM_FIXED_POSINF
		 || bits > 0 && f.bits > 0 && i64_result < 0)
			result.bits = ALE_ACCUM_FIXED_POSINF;
		else if (i64_result <= ALE_ACCUM_FIXED_NEGINF
		      || bits == ALE_ACCUM_FIXED_NEGINF || f.bits == ALE_ACCUM_FIXED_NEGINF
		      || bits < 0 && f.bits < 0 && i64_result > 0)
			result.bits = ALE_ACCUM_FIXED_NEGINF;
		else
			result.bits = i64_result;

		return result;
	}

	ale_accum_fixed operator+(int i) const {
		return operator+(ale_accum_fixed(i));
	}

	ale_accum_fixed operator+(unsigned int i) const {
		return operator+(ale_accum_fixed(i));
	}

	ale_accum_fixed operator-() const {

		ale_accum_fixed result;

		if (bits == ALE_ACCUM_FIXED_NAN || bits == 0)
			return *this;
		else if (bits == ALE_ACCUM_FIXED_POSINF)
			result.bits = ALE_ACCUM_FIXED_NEGINF;
		else if (bits == ALE_ACCUM_FIXED_NEGINF)
			result.bits = ALE_ACCUM_FIXED_POSINF;
		else
			result.bits = -bits;

		return result;
	}

	ale_accum_fixed operator-(ale_accum_fixed f) const {
		ale_accum_fixed result;

		if (bits == ALE_ACCUM_FIXED_NAN || f.bits == ALE_ACCUM_FIXED_NAN
		 || (bits == ALE_ACCUM_FIXED_POSINF && f.bits == ALE_ACCUM_FIXED_POSINF)
		 || (bits == ALE_ACCUM_FIXED_NEGINF && f.bits == ALE_ACCUM_FIXED_NEGINF)) {
			result.bits = ALE_ACCUM_FIXED_NAN;
			return result;
		}

		i64 i64_result = bits - f.bits;

		if (i64_result >= ALE_ACCUM_FIXED_POSINF
		 || bits == ALE_ACCUM_FIXED_POSINF || f.bits == ALE_ACCUM_FIXED_NEGINF
		 || bits > 0 && f.bits < 0 && i64_result < 0)
			result.bits = ALE_ACCUM_FIXED_POSINF;
		else if (i64_result <= ALE_ACCUM_FIXED_NEGINF
		      || bits == ALE_ACCUM_FIXED_NEGINF || f.bits == ALE_ACCUM_FIXED_POSINF
		      || bits < 0 && f.bits > 0 && i64_result > 0)
			result.bits = ALE_ACCUM_FIXED_NEGINF;
		else
			result.bits = i64_result;

		return result;
	}

	ale_accum_fixed operator-(i64 i) const {
		return operator-(ale_accum_fixed(i));
	}

	ale_accum_fixed operator-(ui64 i) const {
		return operator-(ale_accum_fixed(i));
	}

	ale_accum_fixed operator*(ale_accum_fixed f) const {
		ale_accum_fixed result;

		if (bits == ALE_ACCUM_FIXED_NAN || f.bits == ALE_ACCUM_FIXED_NAN) {
			result.bits = ALE_ACCUM_FIXED_NAN;
			return result;
		}
			
		i64 i64_result = ((i64) bits * (i64) f.bits) / (1 << N);

		if (i64_result > (i64) ALE_ACCUM_FIXED_POSINF 
		 || i64_result < (i64) ALE_ACCUM_FIXED_NEGINF
		 || bits == ALE_ACCUM_FIXED_POSINF || f.bits == ALE_ACCUM_FIXED_POSINF
		 || bits == ALE_ACCUM_FIXED_NEGINF || f.bits == ALE_ACCUM_FIXED_NEGINF) {
			if (i64_result > 0)
				result.bits = ALE_ACCUM_FIXED_POSINF;
			else if (i64_result < 0)
				result.bits = ALE_ACCUM_FIXED_NEGINF;
			else if (i64_result == 0)
				result.bits = ALE_ACCUM_FIXED_NAN;
			else
				assert(0);
				
		} else {
			result.bits = i64_result;
		}
		return result;
	}

	ale_accum_fixed operator*(i64 i) const {
		return operator*(ale_accum_fixed(i));
	}

	ale_accum_fixed operator*(ui64 i) const {
		return operator*(ale_accum_fixed(i));
	}

	ale_accum_fixed operator/(ale_accum_fixed<N> f) const {
		ale_accum_fixed result;

		/*
		 * While this approach may not be suitable for all
		 * applications, it can be a convenient way to detect and
		 * manufacture non-finite values.
		 */
		if (bits == ALE_ACCUM_FIXED_NAN || f.bits == ALE_ACCUM_FIXED_NAN
		 || (bits == 0 && f.bits == 0)) {
			result.bits = ALE_ACCUM_FIXED_NAN;
			return result;
		} else if (f.bits == 0 && bits > 0) {
			result.bits = ALE_ACCUM_FIXED_POSINF;
			return result;
		} else if (f.bits == 0 && bits < 0) {
			result.bits = ALE_ACCUM_FIXED_NEGINF;
			return result;
		}
			
		i64 i64_result = ((i64) bits << N) / f.bits;

		if (i64_result > (i64) ALE_ACCUM_FIXED_POSINF)
			result.bits = ALE_ACCUM_FIXED_POSINF;
		else if (i64_result < (i64) ALE_ACCUM_FIXED_NEGINF)
			result.bits = ALE_ACCUM_FIXED_NEGINF;
		else
			result.bits = (i64) i64_result;

		return result;
	}

	ale_accum_fixed operator/(i64 i) const {
		return operator/(ale_accum_fixed(i));
	}

	ale_accum_fixed operator/(ui64 i) const {
		return operator/(ale_accum_fixed(i));
	}

	ale_accum_fixed &operator+=(ale_accum_fixed f) {
		*this = *this + f;
		return *this;
	}

	ale_accum_fixed &operator-=(ale_accum_fixed f) {
		*this = *this - f;
		return *this;
	}

	ale_accum_fixed &operator*=(ale_accum_fixed f) {
		*this = *this * f;
		return *this;
	}

	ale_accum_fixed &operator/=(ale_accum_fixed f) {
		*this = *this / f;
		return *this;
	}

	int operator>(ale_accum_fixed f) {
		if (bits == ALE_ACCUM_FIXED_NAN || f.bits == ALE_ACCUM_FIXED_NAN)
			return 0;

		if (bits > f.bits)
			return 1;

		return 0;
	}

	int operator<(ale_accum_fixed f) const {
		if (bits == ALE_ACCUM_FIXED_NAN || f.bits == ALE_ACCUM_FIXED_NAN)
			return 0;

		if (bits < f.bits)
			return 1;

		return 0;
	}

	int operator>(int d) const {
		return operator>((ale_accum_fixed) d);
	}

	int operator<(int d) const {
		return operator<((ale_accum_fixed) d);
	}

	int operator>(double d) const {
		return operator>((ale_accum_fixed) d);
	}

	int operator<(double d) const {
		return operator<((ale_accum_fixed) d);
	}

	int operator>(unsigned int d) const {
		return operator>((ale_accum_fixed) d);
	}

	int operator<(unsigned int d) const {
		return operator<((ale_accum_fixed) d);
	}
};

#define ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(return_value, op)		\
template<unsigned int N>					\
return_value operator op(double a, const ale_accum_fixed<N> &f) {		\
	ale_accum_fixed<N> g(a);					\
	return g.operator op(f);				\
}								\
								\
template<unsigned int N>					\
return_value operator op(long long a, const ale_accum_fixed<N> &f) {		\
	return (double) a op f;					\
}								\
								\
template<unsigned int N>					\
return_value operator op(unsigned long long a, const ale_accum_fixed<N> &f) {	\
	return (double) a op f;					\
}								\
								\
template<unsigned int N>					\
return_value operator op(int a, const ale_accum_fixed<N> &f) {		\
	return (double) a op f;					\
}								\
								\
template<unsigned int N>					\
return_value operator op(unsigned int a, const ale_accum_fixed<N> &f) {	\
	return (double) a op f;					\
}								\
								\
template<unsigned int N, unsigned int M>					\
return_value operator op(ale_fixed<M> ff, const ale_accum_fixed<N> &f) {	\
	return (ale_accum_fixed<N>) ff op f;					\
}								

ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(ale_accum_fixed<N>, +);
ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(ale_accum_fixed<N>, -);
ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(ale_accum_fixed<N>, *);
ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(ale_accum_fixed<N>, /);
ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(int, <);
ALE_ACCUM_FIXED_INCORPORATE_OPERATOR(int, >);

template<unsigned int N>
ale_accum_fixed<N> fabs(ale_accum_fixed<N> f) {
	if (f < ale_accum_fixed<N>())
		return -f;
	return f;
}

template<unsigned int N>
ale_accum_fixed<N> pow(ale_accum_fixed<N> f, double d) {
	if (d == 2) 
		return f * f;

	if (d == 1) 
		return f;

	if (d == 0)
		return ale_accum_fixed<N>(1);

	return pow((double) f, d);
}

template<unsigned int N>
ale_accum_fixed<N> floor(ale_accum_fixed<N> f) {

	if (N == 0
	 || f.bits == ALE_ACCUM_FIXED_POSINF
	 || f.bits == ALE_ACCUM_FIXED_NEGINF
	 || f.bits == ALE_ACCUM_FIXED_NAN)
		return f;

	ale_accum_fixed<N> result;

	if (f.bits < 0) {
		result.bits = ((f.bits / (1 << N)) - 1) << N;
	} else {
		result.bits = (f.bits / (1 << N)) << N;
	}

	return result;
}

template<unsigned int N>
ale_accum_fixed<N> ceil(ale_accum_fixed<N> f) {
	return -floor(-f);
}

template<unsigned int N>
int ale_accum_fixed<N>::casting_disabled = 0;

template<unsigned int N>
const ale_accum_fixed<N> ale_accum_fixed<N>::const_0_001 = 0.001;

typedef ale_accum_fixed<15> ale_accum;

#undef ale_accum_enable_casting
#undef ale_accum_disable_casting

#define ale_accum_enable_casting() ale_accum::enable_casting()
#define ale_accum_disable_casting() ale_accum::disable_casting()

#else

typedef double ale_accum;

#endif

#undef FIXED

#endif
