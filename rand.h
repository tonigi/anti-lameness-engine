// Copyright 2006 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>
//                              <dhilvert@gmail.com>

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

/*
 * rand.h: random number generator class.
 */

#ifndef __rand_h__
#define __rand_h__

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_PTHREAD
#include <pthread.h>
#endif

class rng_t {
#ifdef USE_PTHREAD
	static pthread_mutex_t rand_mutex;
#endif
	union {
		unsigned int state_ui;
		unsigned short state_us3[3];
	} state_vars;
public:

	void seed(unsigned int s) {
#ifdef USE_PTHREAD
		state_vars.state_us3[0] = (unsigned short) s;
		state_vars.state_ui = s;
#else
		srand(s);
#endif
	}
	int get() {
#ifdef USE_PTHREAD
#if HAVE_NRAND48
		return (int) nrand48(state_vars.state_us3) % RAND_MAX;
#elif HAVE_RAND_R
		return rand_r(&state_vars.state_ui);
#else
		int result;
		pthread_mutex_lock(&rand_mutex);
		result = rand();
		pthread_mutex_unlock(&rand_mutex);
		return result;
#endif
#else
		return rand();
#endif
	}
};

#endif
