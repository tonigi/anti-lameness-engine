// Copyright 2003, 2004, 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
 * d3/cpf.h: Control point file interface.
 */

#ifndef __cpf_h__
#define __cpf_h__

#include "point.h"

#define CPF_VERSION 0
#define CPF_VERSION_MAX 0

class cpf {
	static FILE *load_f;
	static FILE *save_f;
	static int load_version;

	static const char *load_n;
	static const char *save_n;
	static int save_version;

	static void error(const char *message) {
		fprintf(stderr, "cpf: Error: %s", message);
		exit(1);
	}

	static void get_integer(int *i) {
		if(fscanf(load_f, " %d", i) != 1)
			error("Could not get integer.");
	}

	static void get_new_line() {
		int next_character = 0;

		while (next_character != EOF
		    && next_character != '\n')
			next_character = fgetc(load_f);
	}

	static void check_version(int v) {
		if (v > CPF_VERSION_MAX)
			error("Incompatible version number.");
	}

	static point get_type_a() {
		return point::undefined();
	}

	static point get_type_b() {
		return point::undefined();
	}

	static point get_type_c() {
		return point::undefined();
	}

public:
	static void init_loadfile(const char *filename) {
		load_n = filename;
		load_f = fopen(load_n, "r");

		if (!load_f) {
			fprintf(stderr, "cpf: Error: could not open control point file '%s'.", load_n);
			exit(1);
		}
	}

	static void init_savefile(const char *filename) {
		save_n = filename;
		save_f = fopen(save_n, "w");

		if (!save_f) {
			fprintf(stderr, "cpf: Error: could not open control point file '%s'.", save_n);
			exit(1);
		}

		fprintf(save_f, "# created by ALE control point file handler version %d\n", CPF_VERSION);

		fclose(save_f);
	}

	static point get() {
		while (load_f && !feof(load_f)) {
			int command_char;

			command_char = fgetc(load_f);

			switch (command_char) {
				case EOF:
					return point::undefined();
				case '#':
				case ' ':
				case '\n':
				case '\r':
				case '\t':
					break;
				case 'V':
					get_integer(&load_version);
					check_version(load_version);
				        break;
				case 'A':
					return get_type_a();
				case 'B':
					return get_type_b();
				case 'C':
					return get_type_c();
				default:
					error("Unrecognized command");
			}

			get_new_line();
		}

		return point::undefined();
	}

	static void set(point p) {
		assert(0);
	}

	static void finish_loadfile() {
	}

	static void finish_savefile() {
	}
};

#endif
