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

#ifdef USE_PTHREAD
#include <pthread.h>
#endif

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

	class rwlock_t {
#ifdef USE_PTHREAD
		pthread_rwlock_t _lock;
#endif
	public:
		rwlock_t() {
#ifdef USE_PTHREAD
			pthread_rwlock_init(&_lock, NULL);
#endif
		}

		void wrlock() {
#ifdef USE_PTHREAD
			pthread_rwlock_wrlock(&_lock);
#endif
		}

		void rdlock() {
#ifdef USE_PTHREAD
			pthread_rwlock_rdlock(&_lock);
#endif
		}

		void unlock() {
#ifdef USE_PTHREAD
			pthread_rwlock_unlock(&_lock);
#endif
		}
	};

	class lock_t {
#ifdef USE_PTHREAD
		pthread_mutex_t _lock;
#endif
	public:
		lock_t() {
#ifdef USE_PTHREAD
			/* _lock = PTHREAD_MUTEX_INITIALIZER; */
			pthread_mutex_init(&_lock, NULL);
#endif
		}

		void lock() {
#ifdef USE_PTHREAD
			pthread_mutex_lock(&_lock);
#endif
		}

		void unlock() {
#ifdef USE_PTHREAD
			pthread_mutex_unlock(&_lock);
#endif
		}
	};

	class decompose_domain {
		lock_t _lock;
		int ilg, ihg, jlg, jhg;

	protected:
		void lock() {
			_lock.lock();
		}

		void unlock() {
			_lock.unlock();
		}

		virtual void prepare_subdomains(unsigned int threads) {
		}
		virtual void subdomain_algorithm(unsigned int thread, 
				int il, int ih, int jl, int jh) {
		}
		virtual void finish_subdomains(unsigned int threads) {
		}
		
	private:
		struct thread_data_t {
			decompose_domain *this_ptr;
			unsigned int thread_index;
			int il, ih, jl, jh;
		};

		static void *run_thread(void *thread_data) {
			thread_data_t *td = (thread_data_t *) thread_data;
			td->this_ptr->subdomain_algorithm(td->thread_index,
				td->il, td->ih, td->jl, td->jh);
			return NULL;
		}

	public:
		decompose_domain(int ilg, int ihg, int jlg, int jhg) {
			this->ilg = ilg;
			this->ihg = ihg;
			this->jlg = jlg;
			this->jhg = jhg;
		}

		void run() {
			int N;
#ifdef USE_PTHREAD
			N = thread::count();

			pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * N);
			pthread_attr_t *thread_attr = (pthread_attr_t *) 
					malloc(sizeof(pthread_attr_t) * N);
#else
			N = 1;
#endif

			prepare_subdomains(N);

			thread_data_t *td = new thread_data_t[N];

			for (int ti = 0; ti < N; ti++) {
				td[ti].this_ptr = this;
				td[ti].thread_index = ti;
				td[ti].il = ilg + ((ihg - ilg) * ti) / N;
				td[ti].ih = ilg + ((ihg - ilg) * (ti + 1)) / N;
				td[ti].jl = jlg;
				td[ti].jh = jhg;
#ifdef USE_PTHREAD
				pthread_attr_init(&thread_attr[ti]);
				pthread_attr_setdetachstate(&thread_attr[ti], 
						PTHREAD_CREATE_JOINABLE);
				pthread_create(&threads[ti], &thread_attr[ti], 
						run_thread, &td[ti]);
#else
				run_thread(&td[ti]);
#endif
			}


#ifdef USE_PTHREAD
			for (int ti = 0; ti < N; ti++) {
				pthread_join(threads[ti], NULL);
			}
			
			free(threads);
			free(thread_attr);
#endif


			delete[] td;

			finish_subdomains(N);
		}

		virtual ~decompose_domain() {
		}
	};
};

#endif
