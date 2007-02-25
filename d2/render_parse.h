// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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

#ifndef __render_parse_h__
#define __render_parse_h__

#include "render.h"
#include "render/combine.h"
#include "render/invariant.h"
#include "render/incremental.h"
#include "render/zero.h"

/*
 * Parse strings describing renderers, and return a renderer satisfying
 * the string.
 */

class render_parse {
private:
	static int strpfix(const char *a, const char *b) {
		return strncmp(a, b, strlen(a));
	}

	static void nomem() {
		fprintf(stderr, "\n\n*** Error: unable to allocate memory in render_parse. ***\n\n");
		exit(1);
	}

	static void syntax_error(char *explanation) {
		fprintf(stderr, "\n\n*** Error: Render syntax: %s ***\n\n", explanation);
		exit(1);
	}

	static filter::filter *get_SF_atomic(const char *type) {
		double param;

		if (!strcmp("sinc", type)) {
			return new filter::sinc();
		} else if (!strpfix("lanc:", type)) {
			if (sscanf(type + strlen("lanc:"), "%lf", &param) != 1) 
				syntax_error("Unable to get lanczos diameter.");
			return new filter::lanczos(param/2);
		} else if (!strpfix("triangle:", type)) {
			if (sscanf(type + strlen("triangle:"), "%lf", &param) != 1) 
				syntax_error("Unable to get triangle diameter.");
			return new filter::triangle(param/2);
		} else if (!strpfix("box:", type)) {
			if (sscanf(type + strlen("box:"), "%lf", &param) != 1) 
				syntax_error("Unable to get box diameter.");
			return new filter::box(param/2);
        } else if (!strpfix("gauss:", type)) {
            if (sscanf(type + strlen("gauss:"), "%lf", &param) != 1)
                syntax_error("Unable to get gauss deviation.");
            return new filter::gauss(param);
		} else if (!strpfix("zero", type)) {
			return new filter::zero();
		} else {
			fprintf(stderr, "get_SF_atomic type %s\n", type);
			syntax_error("Unable to get filter.");
		}

		assert (0);

		/*
		 * This line should never be reached; it's included to avoid
		 * 'unreachable' messages emitted by some compilers.
		 */

		return NULL;
	}

	static filter::filter *get_SF(const char *orig_type) {
		char *type = strdup(orig_type);
		if (type == NULL)
			nomem();

		char *star_index = (char *) type;
		while (*star_index != '\0'
		    && *star_index != '*')
			star_index++;

		if (*star_index == '\0') {
			free(type);
			return get_SF_atomic(orig_type);
		}
		
		*star_index = '\0';
		filter::filter *result = new filter::mult(
			get_SF_atomic(type),
			get_SF(star_index + 1));
		*star_index = '*';
		free(type);
		return result;
	}	

public:
	static filter::scaled_filter *get_SSF(const char *type) {
		if (!strpfix("coarse:", type)) {
			return new filter::scaled_filter(get_SF(type + strlen("coarse:")), 1);
		} else if (!strpfix("fine:", type)) {
			return new filter::scaled_filter(get_SF(type + strlen("fine:")), 0);
		} else {
			return new filter::scaled_filter(get_SF(type), 1);
		}
	}

	static filter::ssfe *get_SSFE(const char *type) {
		if (!strpfix("ex:", type)) {
			return new filter::ssfe(get_SSF(type + strlen("ex:")), 1);
		} else if (!strpfix("nex:", type)) {
			return new filter::ssfe(get_SSF(type + strlen("nex:")), 0);
		} else {
			return new filter::ssfe(get_SSF(type), 1);
		}
	}

	static render *get_invariant(const char *type) {
		double param;
		int offset;
		invariant *i;
		if (!strpfix("min:", type)) {
			i = new invariant(get_SSFE(type + strlen("min:")));
			i->set_min();
		} else if (!strpfix("max:", type)) {
			i = new invariant(get_SSFE(type + strlen("max:")));
			i->set_max();
		} else if (!strpfix("last:", type)) {
			i = new invariant(get_SSFE(type + strlen("last:")));
			i->set_last();
		} else if (!strpfix("first:", type)) {
			i = new invariant(get_SSFE(type + strlen("first:")));
			i->set_first();
		} else if (!strpfix("avg:", type)) {
			i = new invariant(get_SSFE(type + strlen("avg:")));
			i->set_avg();
		} else if (!strpfix("avgf:", type)) {
			if (sscanf(type + strlen("avgf:"), "%lf%n", &param, &offset) != 1) 
				syntax_error("Unable to get avgf weight criterion.");
			i = new invariant(get_SSFE(type + strlen("avgf:") + offset + strlen(":")));
			i->set_avgf(param);
		} else if (!strpfix("median:", type)) {
			i = new invariant(get_SSFE(type + strlen("median:")));
			i->set_median();
		} else {
			i = new invariant(get_SSFE(type));
			i->set_avg();
		}

		for (int index = 0; index < render::render_count(); index++)
		if  (typeid(*render::render_num(index)) == typeid(incremental)
			 && i->equals(((incremental *)render::render_num(index))->get_invariant())) {
				delete i;
				return (incremental *)render::render_num(index);
		} else if (typeid(*i->ssfe()->get_scaled_filter()->get_filter()) == typeid(filter::zero)) {
			return new zero(i);
		}

		// fprintf(stderr, "  new '%s'.\n", type);

		return new incremental(i);
	}

	static render *get_combination(const char *ptype, const char *dtype) {
		render *partial = get_invariant(ptype);
		render *_default = get(dtype);
		for (int index = 0; index < render::render_count(); index++)
		if  (typeid(*render::render_num(index)) == typeid(combine)
			 && _default == ((combine *)render::render_num(index))->get_default()
			 && partial  == ((combine *)render::render_num(index))->get_partial()) {
				return render::render_num(index);
		}

		// fprintf(stderr, "  new '%s,%s'.\n", ptype, dtype);

		return new combine(_default, partial);
	}

	static render *get(const char *orig_type) {
		char *type = strdup(orig_type);
		if (type == NULL)
			nomem();

		char *comma_index = (char *) type;
		while (*comma_index != '\0'
		    && *comma_index != ',')
			comma_index++;

		if (*comma_index == '\0') {
			free(type);
			return get_invariant(orig_type);
		}
		
		*comma_index = '\0';
		render *result = get_combination(type, comma_index + 1);
		*comma_index = ',';

		free(type);
		return result;
	}

};

#endif
