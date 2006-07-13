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

/*
 * ANSI C and POSIX include files.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <stack>
#include <map>

/*
 * Types
 */

#include "../ale_pos.h"
#include "../ale_real.h"

/*
 * Interface files
 */

#include "ui.h"
#include "unsupported.h"
#include "implication.h"

/*
 * 2D include files
 */

#include "../d2.h"

/*
 * 3D include files
 */

#include "../d3.h"

/*
 * Device configuration files
 */

#include "../device/xvp610_320x240.h"
#include "../device/xvp610_640x480.h"
#include "../device/ov7620_raw_linear.h"
#include "../device/canon_300d_raw_linear.h"
#include "../device/canon_300d_raw_linear_85mm_1_8.h"
#include "../device/canon_300d_raw_linear_50mm_1_8.h"
#include "../device/canon_300d_raw_linear_50mm_1_4.h"
#include "../device/canon_300d_raw_linear_50mm_1_4_1_4.h"

/*
 * Help files
 */

#include "help.h"

/*
 * Argument counter.
 *
 * Counts instances of a given option.
 */
unsigned int arg_count(int argc, const char *argv[], const char *arg) {
	unsigned int count = 0;
	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i], arg))
			count++;
		else if (!strcmp(argv[i], "--"))
			return count;
	}
	return count;
}

/*
 * Argument prefix counter.
 *
 * Counts instances of a given option prefix.
 */
unsigned int arg_prefix_count(int argc, const char *argv[], const char *pfix) {
	unsigned int count = 0;
	for (int i = 0; i < argc; i++) {
		if (!strncmp(argv[i], pfix, strlen(pfix)))
			count++;
		else if (!strcmp(argv[i], "--"))
			return count;
	}
	return count;
}

/*
 * Reallocation function
 */
void *local_realloc(void *ptr, size_t size) {
	void *new_ptr = realloc(ptr, size);

	if (new_ptr == NULL)
		ui::get()->memory_error_location("main()");

	return new_ptr;
}

/*
 * Not enough arguments function.
 */
void not_enough(const char *opt_name) {
	ui::get()->cli_not_enough(opt_name);
}

/*
 * Bad argument function
 */
void bad_arg(const char *opt_name) {
	ui::get()->cli_bad_arg(opt_name);
}

class input {
	/*
	 * Environment structures.
	 *
	 * XXX: It's arguable that these should be public members of the
	 * 'input' class in order to allow passing environment values to other
	 * classes, but, since we're currently using them only to prepare state
	 * for an internal 'input' function, they can stay private for now.  A
	 * more nuanced approach will likely be required later.
	 */

	class environment {
		static std::stack<environment *> environment_stack;

		std::map<const char *, const char *> environment_map;

	public:
		const char *get(const char *name) {
			if (environment_map.count(name) == 0)
				return NULL;

			return environment_map[name];
		}

		void set(const char *name, const char *value) {
			environment_map[name] = value;
		}

		void set_ptr(const char *name, const void *pointer) {
			int chars = sizeof(void *) * 2 + 3;
			char *c = (char *) malloc(sizeof(char) * chars);

			assert(c);

			if (!c)
				ui::get()->memory_error_location("environment::set_ptr");

			int count = snprintf(c, chars, "%p", pointer);

			assert (count >= 0 && count < chars);

			set(name, c);
		}

		/*
		 * Make an environment substructure.  Note that since deep
		 * structures are currently referenced rather than copied when
		 * the stack is pushed, there is no current need for any
		 * chaining mechanism.
		 */
		void make_substructure(const char *name) {
			set_ptr(name, new environment);
		}

		static environment *top() {
			if (environment_stack.empty())
				environment_stack.push(new environment);
			return environment_stack.top();
		}

		static void push() {
			environment *e = new environment;

			e->environment_map = environment_stack.top()->environment_map;

			e->set_ptr("---chain", environment_stack.top());

			environment_stack.push(e);
		}

		static void pop() {

			/*
			 * Execution environments should never be referenced by
			 * structures further up the call chain, so they can
			 * safely be deleted.  (XXX:  In particular, while
			 * lexical scoping may require copying of execution
			 * environments from lower on the call chain, there is
			 * no obvious reason that a reference should be used in
			 * this case; a shallow copy should be used instead.)
			 */

			delete environment_stack.top();

			environment_stack.pop();
		}
	};

	/*
	 * Read tokens from a stream.
	 */
	class token_reader {
		/*
		 * Get the next token
		 */
		virtual char *get() = 0;

		virtual ~token_reader() {
		}
	};

	class cli_token_reader {

		int arg_index;
		int argc;
		const char **argv;

		cli_token_reader(int c, const char *v[]) {
			argc = c;
			argv = v;
			arg_index = 1;
		}

		const char *get() {
			return argv[arg_index++];
		}
	};

public:
	/*
	 * Input handler.
	 *
	 * Does one of two things:
	 *
	 * (1) Output version information if called with '--version'
	 *
	 * (2) Read options and file arguments, and if the arguments are correct, 
	 * write output.  If an error is detected, print the usage statement.
	 *
	 */

	static void handle(int argc, const char *argv[], const char *package, const char *short_version, const char *version) {

		/*
		 * Initialize help object
		 */
		
		help hi(package, argv[0], short_version);

		/*
		 * Output version information if --version appears
		 * on the command line.
		 */

		if (arg_count(argc, argv, "--version")) {
			/*
			 * Output the version
			 */

			fprintf(stdout, "%s", version);

			return;
		}

		/*
		 * Handle help options
		 */

		if  (arg_prefix_count(argc, argv, "--h"))
		for (int i = 1; i < argc; i++) {
			int all = !strcmp(argv[i], "--hA");
			int is_help_option = !strncmp(argv[i], "--h", strlen("--h"));
			int found_help = 0;

			if (!strcmp(argv[i], "--hu") || all)
				hi.usage(), found_help = 1;
			if (!strcmp(argv[i], "--hq") || all)
				hi.defaults(), found_help = 1;
			if (!strcmp(argv[i], "--hf") || all)
				hi.file(), found_help = 1;
			if (!strcmp(argv[i], "--he") || all)
				hi.exclusion(), found_help = 1;
			if (!strcmp(argv[i], "--ha") || all)
				hi.alignment(), found_help = 1;
			if (!strcmp(argv[i], "--hr") || all)
				hi.rendering(), found_help = 1;
			if (!strcmp(argv[i], "--hx") || all)
				hi.exposure(), found_help = 1;
			if (!strcmp(argv[i], "--ht") || all)
				hi.tdf(), found_help = 1;
			if (!strcmp(argv[i], "--hl") || all)
				hi.filtering(), found_help = 1;
			if (!strcmp(argv[i], "--hd") || all)
				hi.device(), found_help = 1;
			if (!strcmp(argv[i], "--hi") || all)
				hi.interface(), found_help = 1;
			if (!strcmp(argv[i], "--hv") || all)
				hi.visp(), found_help = 1;
			if (!strcmp(argv[i], "--hc") || all)
				hi.cp(), found_help = 1;
			if (!strcmp(argv[i], "--h3") || all)
				hi.d3(), found_help = 1;
			if (!strcmp(argv[i], "--hz") || all)
				hi.undocumented(), found_help = 1;

			if (is_help_option && !found_help)
				hi.usage();

			/*
			 * Check for the end-of-options marker, a non-option argument,
			 * or the end of arguments.  In all of these cases, we exit.
			 */

			if (!strcmp(argv[i], "--")
			 || strncmp(argv[i], "--", strlen("--"))
			 || i == argc - 1)
				return;
		}

		/*
		 * Undocumented projective transformation utility
		 */

		if (arg_count(argc, argv, "--ptcalc") > 0) {
			fprintf(stderr, "\n\n*** Warning: this feature is not documented ***\n\n");
			printf("Enter: w h tlx tly blx bly brx bry trx try x y\n\n");
		
			double w, h, tlx, tly, blx, bly, brx, bry, trx, tr_y, x, y;

			printf("> ");

			if (scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
				&w, &h, &tlx, &tly, &blx, &bly, &brx, &bry, &trx, &tr_y, &x, &y) != 12) {

				fprintf(stderr, "Error reading input.\n");
				exit(1);
			}

			d2::image *i = new d2::image_ale_real((int)h, (int)w, 3);
			d2::transformation t = d2::transformation::gpt_identity(i, 1);
			d2::point q[4] = {
				d2::point(tly, tlx),
				d2::point(bly, blx),
				d2::point(bry, brx),
				d2::point(tr_y, trx)
			};
			t.gpt_set(q);

			d2::point a(y, x), b;

			b = t.transform_scaled(a);

			printf("TRANSFORM t(a): (%f, %f)\n", (double) b[1], (double) b[0]);

			b = t.scaled_inverse_transform(a);

			printf("INVERSE t^-1(a): (%f, %f)\n", (double) b[1], (double) b[0]);

			exit(0);
		}

		/*
		 * Flags and variables
		 */

		double scale_factor = 1;
		double vise_scale_factor = 1;
#if 0
		double usm_multiplier = 0.0;
#endif
		int extend = 0;
		struct d2::tload_t *tload = NULL;
		struct d2::tsave_t *tsave = NULL;
		struct d3::tload_t *d3_tload = NULL;
		struct d3::tsave_t *d3_tsave = NULL;
		int ip_iterations = 0;
		int ip_use_median = 0;
		enum { psf_linear, psf_nonlinear, psf_N };
		const char *psf[psf_N] = {NULL, NULL};
		const char *device = NULL;
		int psf_match = 0;
		double psf_match_args[6];
		int inc = 1;
		int exposure_register = 1;
		const char *wm_filename = NULL;
		int wm_offsetx, wm_offsety;
		double cx_parameter = 0;
		double *d3px_parameters = NULL;
		int d3px_count = 0;
		int *ex_parameters = NULL;
		int ex_count = 0;
		int ex_show = 0;
		d2::render *achain;
		const char *achain_type = "triangle:2";
		const char *afilter_type = "internal";
		d2::render **ochain = NULL;
		const char **ochain_names = NULL;
		const char **ochain_types = NULL;
		const char *d3chain_type = NULL;
		int oc_count = 0;
		const char **visp = NULL;
		int vise_count = 0;
		const char **d3_output = NULL;
		const char **d3_depth = NULL;
		unsigned int d3_count = 0;
		double user_view_angle = 0;
		int user_bayer = IMAGE_BAYER_DEFAULT;
		d2::pixel exp_mult = d2::pixel(1, 1, 1);
		std::map<const char *, d3::pt> d3_output_pt;
		std::map<const char *, d3::pt> d3_depth_pt;

		/*
		 * dchain is ochain[0].
		 */

		ochain = (d2::render **) local_realloc(ochain, 
					(oc_count + 1) * sizeof(d2::render *));
		ochain_names = (const char **) local_realloc((void *)ochain_names, 
					(oc_count + 1) * sizeof(const char *));
		ochain_types = (const char **) local_realloc((void *)ochain_types, 
					(oc_count + 1) * sizeof(const char *));

		ochain_types[0] = "sinc*lanc:8";

		oc_count = 1;

		/*
		 * Handle default settings
		 */

		if (arg_prefix_count(argc, argv, "--q") > 1)
			ui::get()->error("more than one default setting option --q* was specified");

		if (arg_count(argc, argv, "--q0")) {
			ochain_types[0] = "fine:box:1,triangle:2";
			achain_type = "triangle:2";
			d2::align::mc(0.3);
			ip_iterations = 0;
			d2::image_rw::exp_noscale();
			cx_parameter = 0;
		} else if (arg_count(argc, argv, "--qn")) {
			ochain_types[0] = "sinc*lanc:6";
			achain_type = "sinc*lanc:6";
			d2::align::mc(0.5);
			ip_iterations = 0;
			d2::image_rw::exp_noscale();
			cx_parameter = 0;
		} else if (arg_count(argc, argv, "--q1")) {
			ochain_types[0] = "median:fine:sinc*lanc:8,triangle:2";
			achain_type = "triangle:2";
			d2::align::mc(0.5);
			ip_iterations = 0;
			d2::image_rw::exp_noscale();
			cx_parameter = 0;
		} else if (arg_count(argc, argv, "--q2")) {
			ochain_types[0] = "sinc*lanc:8";
			achain_type = "sinc*lanc:8";
			d2::align::no_mc();
			ip_iterations = 4;
			d2::image_rw::exp_noscale();
			cx_parameter = 0;
		} else if (arg_count(argc, argv, "--qr")) {
			ochain_types[0] = "sinc*lanc:8";
			achain_type = "sinc*lanc:8";
			d2::align::no_mc();
			ip_iterations = 6;
			d2::image_rw::exp_scale();
			cx_parameter = 0.7;
		} else {
			/*
			 * Same as --q0
			 */
			ochain_types[0] = "fine:box:1,triangle:2";
			achain_type = "triangle:2";
			d2::align::mc(0.3);
			ip_iterations = 0;
			d2::image_rw::exp_noscale();
			cx_parameter = 0;
		}

		/*
		 * Set default d3chain to default ochain[0].
		 */

		d3chain_type = ochain_types[0];

		/* 
		 * Iterate through arguments until we reach the first file 
		 * argument.  After the first file argument, we assume that
		 * all following arguments are files.
		 */

		for (int i = 1; i < argc - 1; i++) {


			if (!strcmp(argv[i], "--q0")
			 || !strcmp(argv[i], "--q1")
			 || !strcmp(argv[i], "--q2")
			 || !strcmp(argv[i], "--qr")
			 || !strcmp(argv[i], "--qn")) {
				/*
				 * Do nothing.  Defaults have already been set.
				 */
			} else if (!strcmp(argv[i], "--8bpc")) {
				d2::image_rw::depth8();
			} else if (!strcmp(argv[i], "--16bpc")) {
				d2::image_rw::depth16();
			} else if (!strcmp(argv[i], "--plain")) {
				d2::image_rw::ppm_plain();
			} else if (!strcmp(argv[i], "--raw")) {
				d2::image_rw::ppm_raw();
			} else if (!strcmp(argv[i], "--auto")) {
				d2::image_rw::ppm_auto();
			} else if (!strcmp(argv[i], "--align-all")) {
				d2::align::all();
			} else if (!strcmp(argv[i], "--align-green")) {
				d2::align::green();
			} else if (!strcmp(argv[i], "--align-sum")) {
				d2::align::sum();
			} else if (!strcmp(argv[i], "--translation")) {
				d2::align::class_translation();
			} else if (!strcmp(argv[i], "--euclidean")) {
				d2::align::class_euclidean();
			} else if (!strcmp(argv[i], "--projective")) {
				d2::align::class_projective();
			} else if (!strcmp(argv[i], "--identity")) {
				d2::align::initial_default_identity();
			} else if (!strcmp(argv[i], "--follow")) {
				d2::align::initial_default_follow();
			} else if (!strcmp(argv[i], "--perturb-output")) {
				d2::align::perturb_output();
			} else if (!strcmp(argv[i], "--perturb-source")) {
				d2::align::perturb_source();
			} else if (!strcmp(argv[i], "--fail-optimal")) {
				d2::align::fail_optimal();
			} else if (!strcmp(argv[i], "--fail-default")) {
				d2::align::fail_default();
			} else if (!strcmp(argv[i], "--no-extend")) {
				extend = 0;
			} else if (!strcmp(argv[i], "--extend")) {
				extend = 1;
			} else if (!strcmp(argv[i], "--no-mc")) {
				d2::align::no_mc();
			} else if (!strcmp(argv[i], "--oc")) {
				d3::scene::oc();
			} else if (!strcmp(argv[i], "--no-oc")) {
				d3::scene::no_oc();
			} else if (!strcmp(argv[i], "--gs")) {
				if (i + 1 >= argc)
					not_enough("--gs");

				d2::align::gs(argv[i+1]);
				i += 1;

			} else if (!strcmp(argv[i], "--gs-mo")) {
				if (i + 1 >= argc)
					not_enough("--gs-mo");

				unsigned int mo_parameter;
				if (sscanf(argv[i+1], "%u", &mo_parameter) != 1)
					bad_arg("--gs-mo");

				d2::align::gs_mo(mo_parameter);
				i += 1;

			} else if (!strcmp(argv[i], "--focus")) {

				/*
				 * Check for argument availability
				 */

				if (i + 1 >= argc) 
					not_enough(argv[i]);

				double one = +1;
				double zero = +0;
				double inf = one / zero;

				assert (isinf(inf) == +1);

				/*
				 * Focus type
				 */

				unsigned int type = 0;
				double distance;
				double px, py;

				if (!strcmp(argv[i+1], "d")) {

					type = 0;

					if (i + 2 > argc) 
						not_enough("--focus d");

					if (sscanf(argv[i+2], "%lf", &distance) != 1)
						bad_arg("--focus d");

					i += 2;

				} else if (!strcmp(argv[i+1], "p")) {

					type = 1;

					if (i + 3 > argc) 
						not_enough("--focus");

					if (sscanf(argv[i+2], "%lf", &px) != 1
					 || sscanf(argv[i+3], "%lf", &py) != 1)
						bad_arg("--focus");

					i += 3;

				} else {
					bad_arg(argv[i]);
				}

				/*
				 * Options
				 */

				unsigned int ci = 0;
				double fr = 0;
				double ht = 0;
				double vt = 0;
				double sd = 0;
				double ed = inf;
				double sx = -inf;
				double ex = inf;
				double sy = -inf;
				double ey = inf;
				double ap = 3;
				unsigned int sc = 3;
				unsigned int fs = 0;
				unsigned int sr = 0;

				int options = 1; 

				i++;

				while (options && i < argc) {
					if (!strncmp(argv[i], "ci=", 3)) {
						if(sscanf(argv[i] + 3, "%u", &ci) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "fr=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &fr) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "ht=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &ht) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "vt=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &vt) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "sy=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &sy) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "ey=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &ey) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "sx=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &sx) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "ex=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &ex) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "sd=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &sd) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "ed=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &ed) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "ap=", 3)) {
						if(sscanf(argv[i] + 3, "%lf", &ap) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "sc=", 3)) {
						if(sscanf(argv[i] + 3, "%u", &sc) != 1)
							bad_arg("--focus");
					} else if (!strncmp(argv[i], "sr=", 3)) {
						if (!strcmp(argv[i], "sr=aperture")) {
							sr = 0;
						} else if (!strcmp(argv[i], "sr=pixel")) {
							sr = 1;
						} else
							bad_arg("--focus");

					} else if (!strncmp(argv[i], "fs=", 3)) {
						if (!strcmp(argv[i], "fs=mean")) {
							fs = 0;
						} else if (!strcmp(argv[i], "fs=median")) {
							fs = 1;
						} else 
							bad_arg("--focus");
					} else
						options = 0;

					if (options)
						i++;
				}

				i--;

				d3::focus::add_region(type, distance, px, py, ci, fr, ht, vt, sd, ed, sx, ex, sy, ey, ap, sc, fs, sr);

			} else if (!strcmp(argv[i], "--3ddp") || !strcmp(argv[i], "--3dvp")) {
				d2::align::keep();

				/*
				 * Unsupported configurations
				 */

				if (ip_iterations)
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (usm_multiplier)
					unsupported::fornow("3D modeling with unsharp mask");
#endif

				/*
				 * Check for argument availability
				 */

				if (i + 10 >= argc) 
					not_enough(argv[i]);

				/*
				 * Initialize if necessary
				 *
				 * Note: because their existence is checked as an
				 * indicator of the presence of 3D arguments, we
				 * initialize these structures here.
				 */

				if (d3_output == NULL) {
					d3_count  = argc - (i + 2) - 1;
					d3_output = (const char **) calloc(d3_count, sizeof(char *));
					d3_depth = (const char **) calloc(d3_count, sizeof(char *));
				}

				unsigned int width, height;
				double view_angle;
				double x, y, z;
				double P, Y, R;

				if (sscanf(argv[i+1], "%u", &width) != 1
				 || sscanf(argv[i+2], "%u", &height) != 1
				 || sscanf(argv[i+3], "%lf", &view_angle) != 1
				 || sscanf(argv[i+4], "%lf", &x) != 1
				 || sscanf(argv[i+5], "%lf", &y) != 1
				 || sscanf(argv[i+6], "%lf", &z) != 1
				 || sscanf(argv[i+7], "%lf", &P) != 1
				 || sscanf(argv[i+8], "%lf", &Y) != 1
				 || sscanf(argv[i+9], "%lf", &R) != 1)
					bad_arg(argv[i]);

				view_angle *= M_PI / 180;
				P *= M_PI / 180;
				Y *= M_PI / 180;
				R *= M_PI / 180;

				d2::transformation t = 
					d2::transformation::eu_identity();
				t.set_domain(height, width);
				d3::pt _pt(t, d3::et(y, x, z, Y, P, R), view_angle);
				
				if (!strcmp(argv[i], "--3dvp")) {
					d3_output_pt[argv[i+10]] = _pt;
				} else if (!strcmp(argv[i], "--3ddp")) {
					d3_depth_pt[argv[i+10]] = _pt;
				} else {
					assert(0);
				}

				i+=10;
			} else if (!strcmp(argv[i], "--3dv")) {
				d2::align::keep();

				unsigned int frame_no;

				/*
				 * Unsupported configurations
				 */

				if (ip_iterations)
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (usm_multiplier)
					unsupported::fornow("3D modeling with unsharp mask");
#endif

				/*
				 * Check for argument availability
				 */

				if (i + 2 >= argc) 
					not_enough("--3dv");

				/*
				 * Initialize if necessary
				 */

				if (d3_output == NULL) {
					d3_count  = argc - (i + 2) - 1;
					d3_output = (const char **) calloc(d3_count, sizeof(char *));
					d3_depth = (const char **) calloc(d3_count, sizeof(char *));
				}

				if (sscanf(argv[i+1], "%d", &frame_no) != 1)
					ui::get()->error("--3dv argument 0 must be an integer");

				if (frame_no >= d3_count)
					ui::get()->error("--3dv argument 0 is too large");

				if (d3_output[frame_no] != NULL) {
					unsupported::fornow ("Writing a single 3D view to more than one output file");
				}

				d3_output[frame_no] = argv[i+2];

				i+=2;
			} else if (!strcmp(argv[i], "--3dd")) {
				d2::align::keep();

				unsigned int frame_no;

				/*
				 * Unsupported configurations
				 */

				if (ip_iterations)
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (usm_multiplier)
					unsupported::fornow("3D modeling with unsharp mask");
#endif

				/*
				 * Check for argument availability
				 */

				if (i + 2 >= argc)
					not_enough("--3dd");

				/*
				 * Initialize if necessary
				 */

				if (d3_output == NULL) {
					d3_count  = argc - (i + 2) - 1;
					d3_output = (const char **) calloc(d3_count, sizeof(char *));
					d3_depth = (const char **) calloc(d3_count, sizeof(char *));
				}

				if (sscanf(argv[i+1], "%d", &frame_no) != 1)
					ui::get()->error("--3dd argument 0 must be an integer");

				if (frame_no >= d3_count)
					ui::get()->error("--3dd argument 0 is too large");

				if (d3_depth[frame_no] != NULL) {
					unsupported::fornow ("Writing a single frame's depth info to more than one output file");
				}

				d3_depth[frame_no] = argv[i+2];

				i+=2;

			} else if (!strcmp(argv[i], "--view-angle")) {
				if (i + 1 >= argc)
					not_enough("--view-angle");

				double va_parameter;
				sscanf(argv[i+1], "%lf", &va_parameter);
				i += 1;
				user_view_angle = va_parameter * M_PI / 180;
			} else if (!strcmp(argv[i], "--cpf-load")) {
				if (i + 1 >= argc)
					not_enough("--cpf-load");

				d3::cpf::init_loadfile(argv[i+1]);
				i += 1;
			} else if (!strcmp(argv[i], "--ui=stream")) {
				ui::set_stream();
			} else if (!strcmp(argv[i], "--ui=tty")) {
				ui::set_tty();
			} else if (!strcmp(argv[i], "--3d-fmr")) {

				if (i + 1 >= argc)
					not_enough("--3d-fmr");

				double fmr_parameter;
				sscanf(argv[i+1], "%lf", &fmr_parameter);
				i += 1;
				d3::scene::fmr(fmr_parameter);

			} else if (!strcmp(argv[i], "--3d-dmr")) {

				if (i + 1 >= argc)
					not_enough("--3d-dmr");

				double dmr_parameter;
				sscanf(argv[i+1], "%lf", &dmr_parameter);
				i += 1;
				d3::scene::dmr(dmr_parameter);

			} else if (!strcmp(argv[i], "--et")) {

				if (i + 1 >= argc)
					not_enough("--et");

				double et_parameter;
				sscanf(argv[i+1], "%lf", &et_parameter);
				i += 1;
				d3::scene::et(et_parameter);

			} else if (!strcmp(argv[i], "--st")) {
				if (i + 1 >= argc)
					not_enough("--st");

				double st_parameter;
				sscanf(argv[i+1], "%lf", &st_parameter);
				i += 1;
				d3::cpf::st(st_parameter);
			} else if (!strcmp(argv[i], "--di-lower")) {
				if (i + 1 >= argc)
					not_enough("--di-lower");

				double di_parameter;
				sscanf(argv[i+1], "%lf", &di_parameter);
				i += 1;
				d3::scene::di_lower(di_parameter);
			} else if (!strcmp(argv[i], "--rc")) {
				if (i + 1 >= argc)
					not_enough("--rc");

				double rc_parameter;
				sscanf(argv[i+1], "%lf", &rc_parameter);
				i += 1;
				d3::scene::rc(rc_parameter);
			} else if (!strcmp(argv[i], "--do-try")) {
				if (i + 1 >= argc)
					not_enough("--do-try");

				double do_parameter;
				sscanf(argv[i+1], "%lf", &do_parameter);
				i += 1;
				d3::scene::do_try(do_parameter);
			} else if (!strcmp(argv[i], "--di-upper")) {
				if (i + 1 >= argc)
					not_enough("--di-upper");

				double dgi_parameter;
				sscanf(argv[i+1], "%lf", &dgi_parameter);
				i += 1;
				d3::scene::di_upper(dgi_parameter);
			} else if (!strcmp(argv[i], "--fc")) {
				if (i + 1 >= argc)
					not_enough("--fc");

				double fc_parameter;
				sscanf(argv[i+1], "%lf", &fc_parameter);
				i += 1;
				d3::scene::fc(fc_parameter);
			} else if (!strcmp(argv[i], "--ecm")) {
				unsupported::discontinued("--ecm <x>");
			} else if (!strcmp(argv[i], "--acm")) {
				unsupported::discontinued("--acm <x>");
			} else if (!strcmp(argv[i], "--def-nn")) {
				if (i + 1 >= argc)
					not_enough("--def-nn");

				double nn_parameter;
				sscanf(argv[i+1], "%lf", &nn_parameter);
				i += 1;
				d2::image_rw::def_nn(nn_parameter);

				if (nn_parameter > 2) {
					fprintf(stderr, "\n\n*** Warning: --def-nn implementation is currently "
							     "inefficient for large radii. ***\n\n");
				}

			} else if (!strcmp(argv[i], "--mc")) {
				if (i + 1 >= argc)
					not_enough("--mc");

				double mc_parameter;
				sscanf(argv[i+1], "%lf", &mc_parameter);
				mc_parameter /= 100;
				i += 1;
				d2::align::mc(mc_parameter);

			} else if (!strcmp(argv[i], "--fx")) {

				if (i + 1 >= argc)
					not_enough("--fx");

				double fx_parameter;
				sscanf(argv[i+1], "%lf", &fx_parameter);
				i += 1;
				d3::scene::fx(fx_parameter);

			} else if (!strcmp(argv[i], "--tcem")) {

				if (i + 1 >= argc)
					not_enough("--tcem");

				double tcem_parameter;
				sscanf(argv[i+1], "%lf", &tcem_parameter);
				i += 1;
				d3::scene::tcem(tcem_parameter);

			} else if (!strcmp(argv[i], "--oui")) {

				if (i + 1 >= argc)
					not_enough("--oui");

				unsigned int oui_parameter;
				sscanf(argv[i+1], "%u", &oui_parameter);
				i += 1;
				d3::scene::oui(oui_parameter);

			} else if (!strcmp(argv[i], "--pa")) {

				if (i + 1 >= argc)
					not_enough("--pa");

				unsigned int pa_parameter;
				sscanf(argv[i+1], "%u", &pa_parameter);
				i += 1;
				d3::scene::pa(pa_parameter);

			} else if (!strcmp(argv[i], "--pc")) {

				if (i + 1 >= argc)
					not_enough("--pc");

				d3::scene::pc(argv[i+1]);
				i += 1;

			} else if (!strcmp(argv[i], "--cw")) {
				d2::align::certainty_weighted(1);
			} else if (!strcmp(argv[i], "--no-cw")) {
				d2::align::certainty_weighted(0);
			} else if (!strcmp(argv[i], "--wm")) {
				if (wm_filename != NULL)
					ui::get()->error("only one weight map can be specified");

				if (i + 3 >= argc)
					not_enough("--wm");
				wm_filename = argv[i+1];

				if (sscanf(argv[i+2], "%d", &wm_offsetx) != 1)
					ui::get()->error("--wm x-argument must be an integer");

				if (sscanf(argv[i+3], "%d", &wm_offsety) != 1)
					ui::get()->error("--wm y-argument must be an integer");

				i += 3;

			} else if (!strcmp(argv[i], "--fl")) {
				if (i + 3 >= argc)
					not_enough("--fl");
				double h, v, a;
				if (sscanf(argv[i+1], "%lf", &h) != 1) 
					ui::get()->error("--fl h-argument must be numerical");
				if (sscanf(argv[i+2], "%lf", &v) != 1)
					ui::get()->error("--fl v-argument must be numerical");
				if (sscanf(argv[i+3], "%lf", &a) != 1)
					ui::get()->error("--fl a-argument must be numerical");
				i += 3;
#ifdef USE_FFTW
				d2::align::set_frequency_cut(h, v, a);
#else
				ui::get()->error_hint("--fl is not supported", "rebuild ALE with FFTW=1");
#endif
			} else if (!strcmp(argv[i], "--wmx")) {
				if (i + 3 >= argc)
					not_enough("--wmx");
#ifdef USE_UNIX
				d2::align::set_wmx(argv[i+1], argv[i+2], argv[i+3]);
#else
				ui::get()->error_hint("--wmx is not supported", "rebuild ALE with POSIX=1");
#endif
				i += 3;
			} else if (!strcmp(argv[i], "--flshow")) {
				if (i + 1 >= argc)
					not_enough("--flshow");
				d2::align::set_fl_show(argv[i+1]);
				i++;
			} else if (!strcmp(argv[i], "--3dpx")) {
				if (i + 6 >= argc)
					not_enough("--3dpx");

				d3px_parameters = (double *) local_realloc(d3px_parameters, (d3px_count + 1) * 6 * sizeof(double));

				for (int param = 0; param < 6; param++)
				if  (sscanf(argv[i + param + 1], "%lf", &(d3px_parameters[6 * d3px_count + param])) != 1)
					bad_arg("--d3px");

				/*
				 * Swap x and y, since their internal meanings differ from their external meanings.
				 */

				for (int param = 0; param < 2; param++) {
					double temp = d3px_parameters[6 * d3px_count + 2 + param];
					d3px_parameters[6 * d3px_count + 2 + param] = d3px_parameters[6 * d3px_count + 0 + param];
					d3px_parameters[6 * d3px_count + 0 + param] = temp;
				}


				/*
				 * Increment counters
				 */

				d3px_count++;
				i += 6;
			} else if (!strcmp(argv[i], "--ex")) {
				if (i + 6 >= argc)
					not_enough("--ex");

				ex_parameters = (int *) local_realloc(ex_parameters, (ex_count + 1) * 6 * sizeof(int));

				for (int param = 0; param < 6; param++)
				if  (sscanf(argv[i + param + 1], "%d", &(ex_parameters[6 * ex_count + param])) != 1)
					bad_arg("--ex");

				/*
				 * Swap x and y, since their internal meanings differ from their external meanings.
				 */

				for (int param = 0; param < 2; param++) {
					int temp = ex_parameters[6 * ex_count + 2 + param];
					ex_parameters[6 * ex_count + 2 + param] = ex_parameters[6 * ex_count + 0 + param];
					ex_parameters[6 * ex_count + 0 + param] = temp;
				}


				/*
				 * Increment counters
				 */

				ex_count++;
				i += 6;
			} else if (!strcmp(argv[i], "--crop")) {
				if (i + 6 >= argc)
					not_enough("--crop");

				ex_parameters = (int *) local_realloc(ex_parameters, (ex_count + 4) * 6 * sizeof(int));
				int crop_args[6];

				for (int param = 0; param < 6; param++)
				if  (sscanf(argv[i + param + 1], "%d", &(crop_args[param])) != 1)
					bad_arg("--crop");

				/*
				 * Construct exclusion regions from the crop area,
				 * swapping x and y, since their internal meanings
				 * differ from their external meanings.
				 */

				/*
				 * Exclusion region 1: low x
				 */

				ex_parameters[6 * ex_count + 0] = INT_MIN;
				ex_parameters[6 * ex_count + 1] = crop_args[2] - 1;
				ex_parameters[6 * ex_count + 2] = INT_MIN;
				ex_parameters[6 * ex_count + 3] = INT_MAX;
				ex_parameters[6 * ex_count + 4] = crop_args[4];
				ex_parameters[6 * ex_count + 5] = crop_args[5];

				/*
				 * Exclusion region 2: low y
				 */

				ex_parameters[6 * ex_count + 6] = INT_MIN;
				ex_parameters[6 * ex_count + 7] = INT_MAX;
				ex_parameters[6 * ex_count + 8] = INT_MIN;
				ex_parameters[6 * ex_count + 9] = crop_args[0] - 1;
				ex_parameters[6 * ex_count + 10] = crop_args[4];
				ex_parameters[6 * ex_count + 11] = crop_args[5];

				/*
				 * Exclusion region 3: high y
				 */

				ex_parameters[6 * ex_count + 12] = INT_MIN;
				ex_parameters[6 * ex_count + 13] = INT_MAX;
				ex_parameters[6 * ex_count + 14] = crop_args[1] + 1;
				ex_parameters[6 * ex_count + 15] = INT_MAX;
				ex_parameters[6 * ex_count + 16] = crop_args[4];
				ex_parameters[6 * ex_count + 17] = crop_args[5];

				/*
				 * Exclusion region 4: high x
				 */

				ex_parameters[6 * ex_count + 18] = crop_args[3] + 1;
				ex_parameters[6 * ex_count + 19] = INT_MAX;
				ex_parameters[6 * ex_count + 20] = INT_MIN;
				ex_parameters[6 * ex_count + 21] = INT_MAX;
				ex_parameters[6 * ex_count + 22] = crop_args[4];
				ex_parameters[6 * ex_count + 23] = crop_args[5];

				/*
				 * Increment counters
				 */

				ex_count += 4;
				i += 6;
			} else if (!strcmp(argv[i], "--exshow")) {
				ex_show = 1;
			} else if (!strcmp(argv[i], "--wt")) {
				if (i + 1 >= argc)
					not_enough("--wt");

				double wt;

				if (sscanf(argv[i + 1], "%lf", &wt) != 1)
					bad_arg("--wt");

				d2::render::set_wt(wt);
				i++;
			} else if (!strcmp(argv[i], "--3d-chain")) {
				if (i + 1 >= argc)
					not_enough("--3d-chain");
				d3chain_type = argv[i+1];
				i++;
			} else if (!strcmp(argv[i], "--dchain")) {
				if (i + 1 >= argc)
					not_enough("--dchain");
				ochain_types[0] = argv[i+1];
				i++;
			} else if (!strcmp(argv[i], "--achain")) {
				if (i + 1 >= argc)
					not_enough("--achain");
				achain_type = argv[i+1];
				i++;
			} else if (!strcmp(argv[i], "--afilter")) {
				if (i + 1 >= argc)
					not_enough("--afilter");
				afilter_type = argv[i+1];
				i++;
			} else if (!strcmp(argv[i], "--ochain")) {
				if (i + 2 >= argc)
					not_enough("--ochain");

				ochain = (d2::render **) local_realloc(ochain, 
							(oc_count + 1) * sizeof(d2::render *));
				ochain_names = (const char **) local_realloc((void *)ochain_names, 
							(oc_count + 1) * sizeof(const char *));
				ochain_types = (const char **) local_realloc((void *)ochain_types, 
							(oc_count + 1) * sizeof(const char *));

				ochain_types[oc_count] = argv[i+1];
				ochain_names[oc_count] = argv[i+2];

				oc_count++;
				i+=2;
			} else if (!strcmp(argv[i], "--visp")) {
				if (i + 5 >= argc)
					not_enough("--visp");

				visp = (const char **) local_realloc((void *)visp, 4 *
							(vise_count + 1) * sizeof(const char *));

				for (int param = 0; param < 4; param++)
					visp[vise_count * 4 + param] = argv[i + 1 + param];

				vise_count++;
				i+=4;
			} else if (!strcmp(argv[i], "--cx")) {
				
				if (i + 1 >= argc)
					not_enough("--cx");

				sscanf(argv[i+1], "%lf", &cx_parameter);
				i += 1;

			} else if (!strcmp(argv[i], "--no-cx")) {
				cx_parameter = 0;
			} else if (!strcmp(argv[i], "--ip")) {
				unsupported::discontinued("--ip <r> <i>", "--lpsf box=<r> --ips <i>");
			} else if (!strcmp(argv[i], "--bayer")) {
				if (i + 1 >= argc)
					not_enough("--bayer");

				/*
				 * External order is clockwise from top-left.  Internal
				 * order is counter-clockwise from top-left.
				 */

				if (!strcmp(argv[i+1], "rgbg")) {
					user_bayer = IMAGE_BAYER_RGBG;
				} else if (!strcmp(argv[i+1], "bgrg")) {
					user_bayer = IMAGE_BAYER_BGRG;
				} else if (!strcmp(argv[i+1], "gbgr")) {
					user_bayer = IMAGE_BAYER_GRGB;
				} else if (!strcmp(argv[i+1], "grgb")) {
					user_bayer = IMAGE_BAYER_GBGR;
				} else if (!strcmp(argv[i+1], "none")) {
					user_bayer = IMAGE_BAYER_NONE;
				} else {
					bad_arg("--bayer");
				}
				i++;
			} else if (!strcmp(argv[i], "--lpsf")) {
				if (i + 1 >= argc)
					not_enough("--lpsf");

				psf[psf_linear] = argv[i+1];
				i++;
			} else if (!strcmp(argv[i], "--nlpsf")) {
				if (i + 1 >= argc)
					not_enough("--nlpsf");

				psf[psf_nonlinear] = argv[i+1];
				i++;

			} else if (!strcmp(argv[i], "--psf-match")) {
				if (i + 6 >= argc)
					not_enough("--psf-match");

				psf_match = 1;

				for (int index = 0; index < 6; index++) {
					if (sscanf(argv[i + 1], "%lf", &psf_match_args[index]) != 1)
						bad_arg("--psf-match");
					i++;
				}

			} else if (!strcmp(argv[i], "--device")) {
				if (i + 1 >= argc) 
					not_enough("--device");

				device = argv[i+1];
				i++;

#if 0
			} else if (!strcmp(argv[i], "--usm")) {

				if (d3_output != NULL)
					unsupported::fornow("3D modeling with unsharp mask");

				if (i + 1 >= argc)
					not_enough("--usm");

				sscanf(argv[i+1], "%lf", &usm_multiplier);
				i++;
#endif

			} else if (!strcmp(argv[i], "--ipr")) {
				
				if (i + 1 >= argc)
					not_enough("--ipr");

				if (sscanf(argv[i+1], "%d", &ip_iterations) != 1)
					ui::get()->error("--ipr requires an integer argument");

				ui::get()->warn("--ipr is deprecated.  Use --ips instead");
				i++;

			} else if (!strcmp(argv[i], "--cpp-err-median")) {
				d3::cpf::err_median();
			} else if (!strcmp(argv[i], "--cpp-err-mean")) {
				d3::cpf::err_mean();
			} else if (!strcmp(argv[i], "--vp-adjust")) {
				d3::align::vp_adjust();
			} else if (!strcmp(argv[i], "--vp-noadjust")) {
				d3::align::vp_noadjust();
			} else if (!strcmp(argv[i], "--vo-adjust")) {
				d3::align::vo_adjust();
			} else if (!strcmp(argv[i], "--vo-noadjust")) {
				d3::align::vo_noadjust();
			} else if (!strcmp(argv[i], "--ip-mean")) {
				ip_use_median = 0;
			} else if (!strcmp(argv[i], "--ip-median")) {
				ip_use_median = 1;
			} else if (!strcmp(argv[i], "--ips")) {

				if (i + 1 >= argc)
					not_enough("--ips");

				if (sscanf(argv[i+1], "%d", &ip_iterations) != 1)
					ui::get()->error("--ips requires an integer argument");
				i++;
			} else if (!strcmp(argv[i], "--ipc")) {
				unsupported::discontinued("--ipc <c> <i>", "--ips <i> --lpsf <c>", "--ips <i> --device <c>");
			} else if (!strcmp(argv[i], "--exp-extend")) {
				d2::image_rw::exp_scale();
			} else if (!strcmp(argv[i], "--exp-noextend")) {
				d2::image_rw::exp_noscale();
			} else if (!strcmp(argv[i], "--exp-register")) {
				exposure_register = 1;
				d2::align::exp_register();
			} else if (!strcmp(argv[i], "--exp-noregister")) {
				exposure_register = 0;
				d2::align::exp_noregister();
			} else if (!strcmp(argv[i], "--exp-meta-only")) {
				exposure_register = 2;
				d2::align::exp_meta_only();
			} else if (!strcmp(argv[i], "--drizzle-only")) {
				unsupported::discontinued("--drizzle-only", "--dchain box:1");
			} else if (!strcmp(argv[i], "--subspace-traverse")) {
				unsupported::undocumented("--subspace-traverse");
				d3::scene::set_subspace_traverse();
			} else if (!strcmp(argv[i], "--3d-nofilter")) {
				d3::scene::nofilter();
			} else if (!strcmp(argv[i], "--3d-filter")) {
				d3::scene::filter();
			} else if (!strcmp(argv[i], "--occ-norm")) {
				d3::scene::nw();
			} else if (!strcmp(argv[i], "--occ-nonorm")) {
				d3::scene::no_nw();
			} else if (!strcmp(argv[i], "--inc")) {
				inc = 1;
			} else if (!strcmp(argv[i], "--no-inc")) {
				inc = 0;
			} else if (!strncmp(argv[i], "--exp-mult=", strlen("--exp-mult="))) {
				double exp_c, exp_r, exp_b;
				sscanf(argv[i] + strlen("--exp-mult="), "%lf,%lf,%lf", &exp_c, &exp_r, &exp_b);
				exp_mult = d2::pixel(1/(exp_r * exp_c), 1/exp_c, 1/(exp_b * exp_c));
			} else if (!strncmp(argv[i], "--visp-scale=", strlen("--visp-scale="))) {

				sscanf(argv[i] + strlen("--visp-scale="), "%lf", &vise_scale_factor);

				if (vise_scale_factor <= 0.0)
					ui::get()->error("VISP scale must be greater than zero");

				if (!finite(vise_scale_factor))
					ui::get()->error("VISP scale must be finite");

			} else if (!strncmp(argv[i], "--scale=", strlen("--scale="))) {

				sscanf(argv[i] + strlen("--scale="), "%lf", &scale_factor);

				if (scale_factor <= 0)
					ui::get()->error("Scale factor must be greater than zero");

				if (!finite(scale_factor))
					ui::get()->error("Scale factor must be finite");

			} else if (!strncmp(argv[i], "--metric=", strlen("--metric="))) {
				double metric;
				sscanf(argv[i] + strlen("--metric="), "%lf", &metric);
				d2::align::set_metric_exponent(metric);
			} else if (!strncmp(argv[i], "--threshold=", strlen("--threshold="))) {
				double match_threshold;
				sscanf(argv[i] + strlen("--threshold="), "%lf", &match_threshold);
				d2::align::set_match_threshold(match_threshold);
			} else if (!strncmp(argv[i], "--drizzle-diam=", strlen("--drizzle-diam="))) {
				unsupported::discontinued("--drizzle-diam=<x>", "--dchain box:1");
				// sscanf(argv[i] + strlen("--drizzle-diam="), "%lf", &drizzle_radius);
				// drizzle_radius /= 2;
			} else if (!strncmp(argv[i], "--perturb-upper=", strlen("--perturb-upper="))) {
				double perturb_upper;
				int characters;
				sscanf(argv[i] + strlen("--perturb-upper="), "%lf%n", &perturb_upper,
										      &characters);
				if (*(argv[i] + strlen("--perturb-upper=") + characters) == '%')
					d2::align::set_perturb_upper(perturb_upper, 1);
				else
					d2::align::set_perturb_upper(perturb_upper, 0);
			} else if (!strncmp(argv[i], "--perturb-lower=", strlen("--perturb-lower="))) {
				double perturb_lower;
				int characters;
				sscanf(argv[i] + strlen("--perturb-lower="), "%lf%n", &perturb_lower,
										      &characters);
				if (perturb_lower <= 0)
					ui::get()->error("--perturb-lower= value is non-positive");

				if (*(argv[i] + strlen("--perturb-lower=") + characters) == '%')
					d2::align::set_perturb_lower(perturb_lower, 1);
				else
					d2::align::set_perturb_lower(perturb_lower, 0);
			} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
				double perturb_lower;
				ui::get()->warn("--stepsize is deprecated.  Use --perturb-lower instead");
				sscanf(argv[i] + strlen("--stepsize="), "%lf", &perturb_lower);
				d2::align::set_perturb_lower(perturb_lower, 0);
			} else if (!strncmp(argv[i], "--va-upper=", strlen("--va-upper="))) {
				double va_upper;
				int characters;
				sscanf(argv[i] + strlen("--va-upper="), "%lf%n", &va_upper,
										      &characters);
				if (*(argv[i] + strlen("--va-upper=") + characters) == '%')
					ui::get()->error("--va-upper= does not accept '%' arguments\n");
				else
					d3::cpf::set_va_upper(va_upper);
			} else if (!strncmp(argv[i], "--cpp-upper=", strlen("--cpp-upper="))) {
				double perturb_upper;
				int characters;
				sscanf(argv[i] + strlen("--cpp-upper="), "%lf%n", &perturb_upper,
										      &characters);
				if (*(argv[i] + strlen("--cpp-upper=") + characters) == '%')
					ui::get()->error("--cpp-upper= does not currently accept '%' arguments\n");
				else
					d3::cpf::set_cpp_upper(perturb_upper);
			} else if (!strncmp(argv[i], "--cpp-lower=", strlen("--cpp-lower="))) {
				double perturb_lower;
				int characters;
				sscanf(argv[i] + strlen("--cpp-lower="), "%lf%n", &perturb_lower,
										      &characters);
				if (*(argv[i] + strlen("--cpp-lower=") + characters) == '%')
					ui::get()->error("--cpp-lower= does not currently accept '%' arguments\n");
				else
					d3::cpf::set_cpp_lower(perturb_lower);
			} else if (!strncmp(argv[i], "--hf-enhance=", strlen("--hf-enhance="))) {
				unsupported::discontinued("--hf-enhance=<x>");
			} else if (!strncmp(argv[i], "--rot-upper=", strlen("--rot-upper="))) {
				double rot_max;
				sscanf(argv[i] + strlen("--rot-upper="), "%lf", &rot_max);
				d2::align::set_rot_max((int) floor(rot_max));
			} else if (!strncmp(argv[i], "--bda-mult=", strlen("--bda-mult="))) {
				double bda_mult;
				sscanf(argv[i] + strlen("--bda-mult="), "%lf", &bda_mult);
				d2::align::set_bda_mult(bda_mult);
			} else if (!strncmp(argv[i], "--bda-rate=", strlen("--bda-rate="))) {
				double bda_rate;
				sscanf(argv[i] + strlen("--bda-rate="), "%lf", &bda_rate);
				d2::align::set_bda_rate(bda_rate);
			} else if (!strncmp(argv[i], "--lod-max=", strlen("--lod-max="))) {
				double lod_max;
				sscanf(argv[i] + strlen("--lod-max="), "%lf", &lod_max);
				d2::align::set_lod_max((int) floor(lod_max));
			} else if (!strncmp(argv[i], "--cpf-load=", strlen("--cpf-load="))) {
				d3::cpf::init_loadfile(argv[i] + strlen("--cpf-load="));
#if 0
			} else if (!strncmp(argv[i], "--model-load=", strlen("--model-load="))) {
				d3::scene::load_model(argv[i] + strlen("--model-load="));
			} else if (!strncmp(argv[i], "--model-save=", strlen("--model-save="))) {
				d3::scene::save_model(argv[i] + strlen("--model-save="));
#endif
			} else if (!strncmp(argv[i], "--trans-load=", strlen("--trans-load="))) {
				d2::tload_delete(tload);
				tload = d2::tload_new(argv[i] + strlen("--trans-load="));
				d2::align::set_tload(tload);
			} else if (!strncmp(argv[i], "--trans-save=", strlen("--trans-save="))) {
				tsave_delete(tsave);
				tsave = d2::tsave_new(argv[i] + strlen("--trans-save="));
				d2::align::set_tsave(tsave);
			} else if (!strncmp(argv[i], "--3d-trans-load=", strlen("--3d-trans-load="))) {
				d3::tload_delete(d3_tload);
				d3_tload = d3::tload_new(argv[i] + strlen("--3d-trans-load="));
				d3::align::set_tload(d3_tload);
			} else if (!strncmp(argv[i], "--3d-trans-save=", strlen("--3d-trans-save="))) {
				d3::tsave_delete(d3_tsave);
				d3_tsave = d3::tsave_new(argv[i] + strlen("--3d-trans-save="));
				d3::align::set_tsave(d3_tsave);
			} else {

				/*
				 * Trap illegal options and end-of-option indicators.
				 */

				if (!strcmp(argv[i], "--"))
					i++;
				else if (!strncmp(argv[i], "--", strlen("--")))
					ui::get()->illegal_option(argv[i]);

				/*
				 * Apply implication logic.
				 */

				if (extend == 0 && vise_count != 0) {
					implication::changed("VISP requires increased image extents.",
							     "Image extension is now enabled.",
							     "--extend");
					extend = 1;
				}

				if (psf_match && ex_count)
					unsupported::fornow("PSF calibration with exclusion regions.");

				
				if (d3_output != NULL && ip_iterations != 0) 
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (extend == 0 && d3_output != NULL) {
					implication::changed("3D modeling requires increased image extents.",
							     "Image extension is now enabled.",
							     "--extend");
					extend = 1;
				}
#endif

				if (cx_parameter != 0 && !exposure_register) {
					implication::changed("Certainty-based rendering requires exposure registration.",
							     "Exposure registration is now enabled.",
							     "--exp-register");
					d2::align::exp_register();
					exposure_register = 1;
				}

				/*
				 * Set alignment class exclusion region static variables
				 */

				d2::align::set_exclusion(ex_parameters, ex_count);

				/*
				 * Initialize renderer class statics.
				 */

				d2::render::render_init(ex_count, ex_parameters, ex_show, extend, scale_factor);

				/*
				 * Set confidence
				 */

				d2::exposure::set_confidence(cx_parameter);

				/*
				 * Keep transformations for Irani-Peleg, psf-match, and
				 * VISE
				 */

				if (ip_iterations > 0 || psf_match || vise_count > 0) {
					d2::align::keep();
				}

				/*
				 * Initialize device-specific variables
				 */

				d2::psf *device_response[psf_N] = { NULL, NULL };
				d2::exposure **input_exposure = NULL;
				ale_pos view_angle = 43.7 * M_PI / 180;  
				// ale_pos view_angle = 90 * M_PI / 180;  
				input_exposure = (d2::exposure **)
					malloc((argc - i - 1) * sizeof(d2::exposure *));

				if (device != NULL) {
					if (!strcmp(device, "xvp610_640x480")) {
						device_response[psf_linear] = new xvp610_640x480::lpsf();
						device_response[psf_nonlinear] = new xvp610_640x480::nlpsf();
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new xvp610_640x480::exposure();
						view_angle = xvp610_640x480::view_angle();
					} else if (!strcmp(device, "xvp610_320x240")) {
						device_response[psf_linear] = new xvp610_320x240::lpsf();
						device_response[psf_nonlinear] = new xvp610_320x240::nlpsf();
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new xvp610_320x240::exposure();
						view_angle = xvp610_320x240::view_angle();
					} else if (!strcmp(device, "ov7620_raw_linear")) {
						device_response[psf_linear] = new ov7620_raw_linear::lpsf();
						device_response[psf_nonlinear] = NULL;
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new ov7620_raw_linear::exposure();
						d2::image_rw::set_default_bayer(IMAGE_BAYER_BGRG);
					} else if (!strcmp(device, "canon_300d_raw_linear")) {
						device_response[psf_linear] = new canon_300d_raw_linear::lpsf();
						device_response[psf_nonlinear] = NULL;
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new canon_300d_raw_linear::exposure();
						d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
					} else if (!strcmp(device, "canon_300d_raw_linear+85mm_1.8")) {
						device_response[psf_linear] = new canon_300d_raw_linear_85mm_1_8::lpsf();
						device_response[psf_nonlinear] = NULL;
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new canon_300d_raw_linear_85mm_1_8::exposure();
						d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
						view_angle = canon_300d_raw_linear_85mm_1_8::view_angle();
					} else if (!strcmp(device, "canon_300d_raw_linear+50mm_1.8")) {
						device_response[psf_linear] = new canon_300d_raw_linear_50mm_1_8::lpsf();
						device_response[psf_nonlinear] = NULL;
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new canon_300d_raw_linear_50mm_1_8::exposure();
						d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
						view_angle = canon_300d_raw_linear_50mm_1_8::view_angle();
					} else if (!strcmp(device, "canon_300d_raw_linear+50mm_1.4")) {
						device_response[psf_linear] = new canon_300d_raw_linear_50mm_1_4::lpsf();
						device_response[psf_nonlinear] = NULL;
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new canon_300d_raw_linear_50mm_1_4::exposure();
						d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
						view_angle = canon_300d_raw_linear_50mm_1_4::view_angle();
					} else if (!strcmp(device, "canon_300d_raw_linear+50mm_1.4@1.4")) {
						device_response[psf_linear] = new canon_300d_raw_linear_50mm_1_4_1_4::lpsf();
						device_response[psf_nonlinear] = NULL;
						for (int ii = 0; ii < argc - i - 1; ii++)
							input_exposure[ii] = new canon_300d_raw_linear_50mm_1_4_1_4::exposure();
						d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
						view_angle = canon_300d_raw_linear_50mm_1_4_1_4::view_angle();
					} else {
						ui::get()->unknown_device(device);
					}
				} else {
					for (int ii = 0; ii < argc - i - 1; ii++)
						input_exposure[ii] = new d2::exposure_default();
				}

				/*
				 * User-specified variables.
				 */

				if (user_view_angle != 0) {
					view_angle = user_view_angle;
				}

				if (user_bayer != IMAGE_BAYER_DEFAULT) {
					d2::image_rw::set_default_bayer(user_bayer);
				}

				/*
				 * PSF-match exposure.
				 */
				if (psf_match) {
					delete input_exposure[argc - i - 2];
					input_exposure[argc - i - 2] = new d2::exposure_default();
				}

				/*
				 * Initialize output exposure
				 */

				d2::exposure *output_exposure = new d2::exposure_default();
				output_exposure->set_multiplier(exp_mult);

				/*
				 * Configure the response function.
				 */

				d2::psf *response[2] = {NULL, NULL};

				for (int n = 0; n < psf_N; n++ ) {
					if (psf[n] != NULL) {

						response[n] = d2::psf_parse::get((n == psf_linear), psf[n]);

					} else if (device_response[n] != NULL) {

						/*
						 * Device-specific response
						 */

						response[n] = device_response[n];

					} else {

						/*
						 * Default point-spread function.
						 */

						if (n == psf_linear) {

							/*
							 * Default lpsf is a box filter
							 * of diameter 1.0 (radius
							 * 0.5).
							 */

							response[n] = new d2::box(0.5);

						} else if (n == psf_nonlinear) {

							/*
							 * nlpsf is disabled by default.
							 */

							 response[n] = NULL;
						}
					}
				}

				/* 
				 * First file argument.  Print general file information as well
				 * as information specific to this argument.  Initialize image
				 * file handler.
				 */

				/*
				 * There should be at least two file arguments.
				 */

				if (i >= argc - 1) {
					hi.usage();
					exit(1);
				}

				d2::image_rw::init(argc - i - 1, argv + i, argv[argc - 1], input_exposure, output_exposure);
				ochain_names[0] = argv[argc - 1];

				/*
				 * Handle control point data for alignment
				 */
				d2::align::set_cp_count(d3::cpf::count());
				for (unsigned int ii = 0; ii < d3::cpf::count(); ii++)
					d2::align::set_cp(ii, d3::cpf::get_2d(ii));

				/*
				 * PSF-match bayer patterns.
				 */

				if (psf_match) {
					d2::image_rw::set_specific_bayer(argc - i - 2, IMAGE_BAYER_NONE);
				}

				/*
				 * Handle alignment weight map, if necessary
				 */

				if (wm_filename != NULL) {
					d2::image *weight_map;
					weight_map = d2::image_rw::read_image(wm_filename, new d2::exposure_linear());
					weight_map->set_offset(wm_offsety, wm_offsetx);
					d2::align::set_weight_map(weight_map);
				}

				/*
				 * Write comment information about original frame and
				 * target image to the transformation save file, if we
				 * have one.
				 */

				const d2::image *im = d2::image_rw::open(0);
				tsave_orig(tsave, argv[i], im->avg_channel_magnitude());
				tsave_target(tsave, argv[argc - 1]);
				d2::image_rw::close(0);

				/*
				 * Initialize alignment interpolant.
				 */

				if (afilter_type != "internal")
					d2::align::set_interpolant(d2::render_parse::get_SSF(afilter_type));

				/*
				 * Initialize achain and ochain.
				 */

				achain = d2::render_parse::get(achain_type);
				
				for (int chain = 0; chain < oc_count; chain++)
					ochain[chain] = d2::render_parse::get(ochain_types[chain]);

				/*
				 * Use merged renderings as reference images in
				 * alignment.
				 */

				d2::align::set_reference(achain);

				/*
				 * Tell the alignment class about the scale factor.
				 */

				d2::align::set_scale(scale_factor);

				/*
				 * Initialize visp.
				 */

				d2::vise_core::set_scale(vise_scale_factor);

				for (int opt = 0; opt < vise_count; opt++) {
					d2::vise_core::add(d2::render_parse::get(visp[opt * 4 + 0]),
							   visp[opt * 4 + 1],
							   visp[opt * 4 + 2],
							   visp[opt * 4 + 3]);
				}

				/*
				 * Initialize non-incremental renderers
				 */

#if 0
				if (usm_multiplier != 0) {

					/*
					 * Unsharp Mask renderer
					 */

					ochain[0] = new d2::usm(ochain[0], scale_factor,
							usm_multiplier, inc, response[psf_linear],
							response[psf_nonlinear], &input_exposure[0]);
				}
#endif

				if (psf_match) {
					
					/*
					 * Point-spread function calibration renderer.
					 * This renderer does not produce image output.
					 * It is reserved for use with the point-spread
					 * function calibration script
					 * ale-psf-calibrate.
					 */

					ochain[0] = new d2::psf_calibrate(ochain[0],
							1, inc, response[psf_linear],
							response[psf_nonlinear],
							psf_match_args);

				} else if (ip_iterations != 0) {

					/*
					 * Irani-Peleg renderer
					 */

					ochain[0] = new d2::ipc( ochain[0], ip_iterations,
							inc, response[psf_linear],
							response[psf_nonlinear],
							(exposure_register == 1), ip_use_median);
				}

				/*
				 * Handle the original frame.
				 */

				ui::get()->original_frame_start(argv[i]);

				for (int opt = 0; opt < oc_count; opt++) {
					ui::get()->set_orender_current(opt);
					ochain[opt]->sync(0);
					if  (inc) {
						ui::get()->writing_output(opt);
						d2::image_rw::write_image(ochain_names[opt], 
							ochain[opt]->get_image(0));
					}
				}

				d2::vise_core::frame_queue_add(0);

				ui::get()->original_frame_done();

				/*
				 * Handle supplemental frames.
				 */

				for (unsigned int j = 1; j < d2::image_rw::count(); j++) {

					const char *name = d2::image_rw::name(j);

					ui::get()->supplemental_frame_start(name);

					/*
					 * Write comment information about the
					 * supplemental frame to the transformation
					 * save file, if we have one.
					 */

					tsave_info (tsave, name);

					const d2::image *im = d2::image_rw::open(j);
					d2::pixel apm = im->avg_channel_magnitude();
					tsave_apm(tsave, apm[0], apm[1], apm[2]);
					d2::image_rw::close(j);

					for (int opt = 0; opt < oc_count; opt++) {
						ui::get()->set_orender_current(opt);
						ochain[opt]->sync(j);
						if (inc) {
							ui::get()->writing_output(opt);
							d2::image_rw::write_image(ochain_names[opt], 
								ochain[opt]->get_image(j));
						}
					}

					d2::vise_core::frame_queue_add(j);

					ui::get()->supplemental_frame_done();
				}

				/*
				 * Do any post-processing and output final image
				 *
				 * XXX: note that non-incremental renderers currently
				 * return zero for ochain[0]->sync(), since they write
				 * output internally when inc != 0.
				 */

				for (int opt = 0; opt < oc_count; opt++) 
				if  ((ochain[opt]->sync() || !inc) && !psf_match)
					d2::image_rw::write_image(ochain_names[opt], ochain[opt]->get_image());

				/*
				 * Output a summary match statistic.
				 */

				ui::get()->ale_2d_done((double) d2::align::match_summary());

				/*
				 * Perform any 3D tasks
				 */

				optimizations::begin_3d_work();

				if (d3_count > 0) {

					ui::get()->d3_start();

					d3::align::init_angle(view_angle);

					ui::get()->d3_init_view_angle(view_angle / M_PI * 180);

					d3::align::init_from_d2();

					if (d3::cpf::count() > 0) {
						ui::get()->d3_control_point_solve();
						d3::cpf::solve_3d();
						ui::get()->d3_control_point_solve_done();
					}

					ui::get()->d3_final_view_angle(d3::align::angle_of(0) / M_PI * 180);

					d3::align::write_alignments();

					d3::scene::set_filter_type(d3chain_type);

					d3::scene::init_from_d2();

					ui::get()->d3_subdividing_space();
					d3::scene::make_space(d3_depth, d3_output, &d3_depth_pt, &d3_output_pt);
					ui::get()->d3_subdividing_space_done();

					ui::get()->d3_updating_occupancy();
					d3::scene::reduce_cost_to_search_depth(output_exposure, inc);
					ui::get()->d3_updating_occupancy_done();

					d3::scene::d3px(d3px_count, d3px_parameters);
					int view_count = 0;
					for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
						assert (i < d3_count);

						if (d3_depth[i] != NULL) {
							ui::get()->d3_writing_output(d3_depth[i]);
							ui::get()->d3_render_status(0, 0, -1, -1, -1, -1, 0);
							const d2::image *im = d3::scene::depth(i);
							d2::image_rw::write_image(d3_depth[i], im, output_exposure, 1, 1);
							delete im;
							ui::get()->d3_writing_output_done();
						}

						if (d3_output[i] != NULL) {
							ui::get()->d3_writing_output(d3_output[i]);
							const d2::image *im = d3::scene::view(i);
							d2::image_rw::write_image(d3_output[i], im, output_exposure);
							delete im;
							d3::focus::set_camera(view_count++);
							ui::get()->d3_writing_output_done();
						}

						for (std::map<const char *, d3::pt>::iterator i = d3_output_pt.begin();
								i != d3_output_pt.end(); i++) {

							ui::get()->d3_writing_output(i->first);
							const d2::image *im = d3::scene::view(i->second);
							d2::image_rw::write_image(i->first, im, output_exposure);
							delete im;
							d3::focus::set_camera(view_count++);
							ui::get()->d3_writing_output_done();
						}

						for (std::map<const char *, d3::pt>::iterator i = d3_depth_pt.begin();
								i != d3_depth_pt.end(); i++) {

							ui::get()->d3_writing_output(i->first);
							ui::get()->d3_render_status(0, 0, -1, -1, -1, -1, 0);
							const d2::image *im = d3::scene::depth(i->second);
							d2::image_rw::write_image(i->first, im, output_exposure, 1, 1);
							delete im;
							ui::get()->d3_writing_output_done();
						}
					}

					for (unsigned int i = d2::image_rw::count(); i < d3_count; i++) {
						if (d3_depth[i] != NULL) {
							fprintf(stderr, "\n\n*** Frame number for --3dd too high. ***\n\n");
						}
						if (d3_output[i] != NULL) {
							fprintf(stderr, "\n\n*** Frame number for --3dv too high. ***\n\n");
						}
					}
				}

				/*
				 * Destroy the image file handler
				 */

				d2::image_rw::destroy();

				/*
				 * Delete the transformation file structures, if any
				 * exist.
				 */

				tsave_delete(tsave);
				tload_delete(tload);

				/*
				 * We're done.
				 */

				exit(0);
			}
		}

		/*
		 * If there was no output, the user might need more information.
		 */

		hi.usage();
		exit(1);

	}
};
