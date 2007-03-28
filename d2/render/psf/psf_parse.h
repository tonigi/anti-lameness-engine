// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

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

#ifndef __psf_parse_h__
#define __psf_parse_h__

#include "psf.h"
#include "box.h"
#include "circle.h"
#include "sgauss.h"
#include "sum.h"
#include "convolution.h"
#include "scalar_mult.h"
#include "stdin.h"
#include "stdin_vg.h"

/*
 * Parse strings describing point-spread functions, and return a psf object
 * satisfying the string.
 */

class psf_parse {
private:
	static int strpfix(const char *a, const char *b) {
		return strncmp(a, b, strlen(a));
	}

	static void nomem() {
		fprintf(stderr, "\n\n*** Error: unable to allocate memory in psf_parse. ***\n\n");
		exit(1);
	}

	static void syntax_error(char *explanation) {
		fprintf(stderr, "\n\n*** Error: PSF syntax: %s ***\n\n", explanation);
		exit(1);
	}

	/*
	 * Evaluate a type string having no remaining binary operators.
	 */
	static psf *get_atomic(int is_linear, const char *orig_type) {
		double param, param2;

		if (!strcmp(orig_type, "stdin")) {

			fprintf(stderr, "\nInitializing ");
			fprintf(stderr, is_linear ? "linear" : "non-linear");
			fprintf(stderr, " PSF.\n");
							
			return new psf_stdin();
		} else if (!strcmp(orig_type, "stdin_vg")) {

			fprintf(stderr, "\nInitializing ");
			fprintf(stderr, is_linear ? "linear" : "non-linear");
			fprintf(stderr, " PSF.\n");
							
			return new psf_stdin_vg();
		} else if (!strpfix("box=", orig_type)) {
			if (sscanf(orig_type + strlen("box="), "%lf", &param) != 1)
				syntax_error("Unable to get box diameter.");
			return new box(param / 2);
		} else if (!strpfix("circle=", orig_type)) {
			if (sscanf(orig_type + strlen("circle="), "%lf", &param) != 1)
				syntax_error("Unable to get circle diameter.");
			return new circle(param / 2);
	       } else if (!strpfix("sgauss=", orig_type)) {
		       if (sscanf(orig_type + strlen("sgauss="), "%lfx%lf", &param, &param2) != 2)
			       syntax_error("Unable to get sgauss diameters.");
		       return new sgauss(param / 2, param2 / 2);
		} else {
			fprintf(stderr, "get_atomic type %s\n", orig_type);
			syntax_error("Unable to get filter.");
		}
		assert(0);
	}

	/*
	 * Get a scalar value
	 */
	static ale_real get_scalar(const char *orig_type) {
		double result;
		if (sscanf(orig_type, "%lf", &result) != 1)
			syntax_error("Unable to get scalar value.");
		return result;
	}

	/*
	 * Split a type string with the binary operator having
	 * third-lowest precedence (i.e., scalar multiplication).
	 */
	static psf *get_scalar_mult(int is_linear, const char *orig_type) {
		char *type = strdup(orig_type);
		char *operator_index = (char *) type;

		assert(type);
		if (!type)
			nomem();

		while (*operator_index != '\0'
		    && *operator_index != '*')
			operator_index++;

		if (*operator_index == '\0') {
			free(type);
			return get_atomic(is_linear, orig_type);
		}
		
		*operator_index = '\0';
		psf *result = new scalar_mult(get_scalar(type), get_scalar_mult(is_linear, operator_index + 1));
		*operator_index = '*';

		free(type);
		return result;
	}

	/*
	 * Split a type string with the binary operator having
	 * second-lowest precedence (i.e., convolution).
	 */
	static psf *get_convolution(int is_linear, const char *orig_type) {
		char *type = strdup(orig_type);
		char *operator_index = (char *) type;

		assert(type);
		if (!type)
			nomem();

		while (*operator_index != '\0'
		    && *operator_index != '^')
			operator_index++;

		if (*operator_index == '\0') {
			free(type);
			return get_scalar_mult(is_linear, orig_type);
		}

		*operator_index = '\0';
		psf *result = new convolution(get_scalar_mult(is_linear, type), get_convolution(is_linear, operator_index + 1));
		*operator_index = '^';

		free(type);
		return result;
	}

	/*
	 * Split the type string using the binary operator with
	 * lowest precedence (addition).
	 */
	static psf *get_summation(int is_linear, const char *orig_type) {
		char *type = strdup(orig_type);
		char *plus_index = (char *) type;

		assert(type);
		if (!type)
			nomem();

		while (*plus_index != '\0'
		    && *plus_index != '+')
			plus_index++;

		if (*plus_index == '\0') {
			free(type);
			return get_convolution(is_linear, orig_type);
		}
		
		*plus_index = '\0';
		psf *result = new sum(get_convolution(is_linear, type), get_summation(is_linear, plus_index + 1));
		*plus_index = '+';

		free(type);
		return result;
	}

public:
	static psf *get(int is_linear, const char *orig_type) {
		return get_summation(is_linear, orig_type);
	}
};

#endif
