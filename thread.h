// Copyright 2006 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>
//                              <dhilvert@gmail.com>

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
 * thread.h: threading details.
 */

#ifndef __thread_h__
#define __thread_h__

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define THREAD_PER_CPU_DEFAULT 1
#define THREAD_COUNT_DEFAULT 4

class thread {
	static unsigned int _count;
	static unsigned int _cpu_count;

	static void try_linux() {
		assert (_cpu_count == 0);

		FILE *cpuinfo;
		char buffer[100];

		cpuinfo = fopen("/proc/cpuinfo", "r");

		if (!cpuinfo)
			return;

		while (!feof(cpuinfo)) {
			fgets(buffer, 100, cpuinfo);
			if (strncmp("processor", buffer, strlen("processor")))
				continue;
			_cpu_count++;
		}
	}

public:
	static void init() {
		if (_cpu_count == 0) {
			try_linux();
		}

		if (_cpu_count > 0) {
			_count = THREAD_PER_CPU_DEFAULT * _cpu_count; 
		} else {
			_count = THREAD_COUNT_DEFAULT;
		}

		assert (_count > 0);
	}

	static void set_per_cpu(unsigned int new_per_cpu) {
		if (_cpu_count == 0) {
			fprintf(stderr, "\n\n");
			fprintf(stderr, "Error: per-cpu thread count specified, but CPU count is unknown.\n");
			fprintf(stderr, "       Try setting the thread count explicitly.\n");

			exit(1);
		}
		if (new_per_cpu == 0) {
			fprintf(stderr, "\n\n");
			fprintf(stderr, "Error: --per-cpu argument must be positive\n");
			fprintf(stderr, "\n");

			exit(1);
		}

		_count = _cpu_count * new_per_cpu;
		assert (_count > 0);
	}

	static void set_count(unsigned int new_count) {
		if (new_count == 0) {
			fprintf(stderr, "\n\n");
			fprintf(stderr, "Error: --thread argument must be positive\n");
			fprintf(stderr, "\n");

			exit(1);
		}

		_count = new_count;
		assert (_count > 0);
	}

	static int count() {
		assert (_count > 0);
		return _count;
	}
};

#endif
