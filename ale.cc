// Copyright 2002, 2003 David Hilvert <dhilvert@auricle.dyndns.org>,
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

/*
 * ale.cc: The main loop of ALE, including the user interface.
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

/*
 * Interface files
 */

#include "unsupported.h"
#include "implication.h"

/*
 * 2D include files
 */

#include "d2.h"

/*
 * Device configuration files
 */

#include "device/xvp610_320x240.h"
#include "device/xvp610_640x480.h"
#include "device/ov7620_raw_linear.h"

/*
 * Help files
 */

#include "help.h"

/*
 * Version Information
 */

char *short_version = "0.7.1-patch1";

char *version = "ALE Version:      0.7.1-patch1\n"
#ifdef USE_MAGICK
		"File handler:     ImageMagick\n"
#else
		"File handler:     PPM\n"
#endif
		"Color data:       " ALE_REAL_PRECISION_STRING "\n"
		"Coordinate data:  " ALE_POS_PRECISION_STRING "\n"
#ifdef USE_FFTW
		"DFT:              FFTW3\n"
#else
		"DFT:              Built-in\n"
#endif
#if defined NDEBUG && !defined DEBUG
                "Assertions:       Disabled\n"
#elif defined DEBUG && !defined NDEBUG
                "Assertions:       Enabled\n"
#elif defined NDEBUG
                "Assertions:       Probably disabled\n"
#else
                "Assertions:       Probably enabled\n"
#endif
;

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
 * Reallocation function
 */
void *local_realloc(void *ptr, size_t size) {
	void *new_ptr = realloc(ptr, size);

	if (new_ptr == NULL) {
		fprintf(stderr, "\n\n*** Unable to allocate memory in main() ***\n\n");
		exit(1);
	}

	return new_ptr;
}

/*
 * Not enough arguments function.
 */
void not_enough(const char *opt_name) {
	fprintf(stderr, "\n\n*** Not enough arguments for %s ***\n\n", opt_name);
	exit(1);
}

/*
 * Bad argument function
 */
void bad_arg(const char *opt_name) {
	fprintf(stderr, "\n\n*** Bad argument to %s ***\n\n", opt_name);
	exit(1);
}

/*
 * Main function.  
 *
 * Does one of two things:
 *
 * (1) Output version information if called with '--version'
 *
 * (2) Read options and file arguments, and if the arguments are correct, 
 * write output.  If an error is detected, print the usage statement.
 *
 */

int main(int argc, const char *argv[]){

	/*
	 * Version information and help.
	 */
	
	help hi(argv[0], short_version);

	if (arg_count(argc, argv, "--version") > 0) {
		/*
		 * Output the version
		 */

		fprintf(stderr, "%s", version);

		return 0;
	} else if (arg_count(argc, argv, "--hu") > 0) {
		hi.usage();
		exit(1);
	} else if (arg_count(argc, argv, "--hq") > 0) {
		hi.defaults();
		exit(1);
	} else if (arg_count(argc, argv, "--hf") > 0) {
		hi.file();
		exit(1);
	} else if (arg_count(argc, argv, "--he") > 0) {
		hi.exclusion();
		exit(1);
	} else if (arg_count(argc, argv, "--ha") > 0) {
		hi.alignment();
		exit(1);
	} else if (arg_count(argc, argv, "--hr") > 0) {
		hi.rendering();
		exit(1);
	} else if (arg_count(argc, argv, "--hx") > 0) {
		hi.exposure();
		exit(1);
	} else if (arg_count(argc, argv, "--ht") > 0) {
		hi.tdf();
		exit(1);
	} else if (arg_count(argc, argv, "--hl") > 0) {
		hi.filtering();
		exit(1);
	} else if (arg_count(argc, argv, "--hd") > 0) {
		hi.device();
		exit(1);
	} else if (arg_count(argc, argv, "--hv") > 0) {
		hi.visp();
		exit(1);
	} else if (arg_count(argc, argv, "--hz") > 0) {
		hi.undocumented();
		exit(1);
	} else if (arg_count(argc, argv, "--hA") > 0) {
                hi.usage();
                hi.defaults();
                hi.file();
                hi.exclusion();
                hi.alignment();
                hi.rendering();
                hi.exposure();
                hi.tdf();
                hi.filtering();
                hi.device();
                hi.visp();
                hi.undocumented();
		exit(1);
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
	 * Flags
	 */

	double scale_factor = 1;
	double vise_scale_factor = 1;
#if 0
	double usm_multiplier = 0.0;
#endif
	int extend = 0;
	struct d2::tload_t *tload = NULL;
	struct d2::tsave_t *tsave = NULL;
	int ip_iterations = 0;
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
	int *ex_parameters = NULL;
	int ex_count = 0;
	int ex_show = 0;
	d2::render *achain;
	const char *achain_type = "triangle:2";
	const char *afilter_type = "internal";
	d2::render **ochain = NULL;
	const char **ochain_names = NULL;
	const char **ochain_types = NULL;
	int oc_count = 0;
	const char **visp = NULL;
	int vise_count = 0;

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
		ochain_types[0] = "fine:sinc*lanc:6,sinc*lanc:6";
		achain_type = "sinc*lanc:6";
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
		} else if (!strcmp(argv[i], "--gs")) {
			if (i + 1 >= argc)
				not_enough("--gs");

			d2::align::gs(argv[i+1]);
			i += 1;

		} else if (!strcmp(argv[i], "--mc")) {
			if (i + 1 >= argc)
				not_enough("--mc");

			double mc_parameter;
			sscanf(argv[i+1], "%lf", &mc_parameter);
			mc_parameter /= 100;
			i += 1;
			d2::align::mc(mc_parameter);

		} else if (!strcmp(argv[i], "--wm")) {
			if (wm_filename != NULL) {
				fprintf(stderr, "\n\n*** Only one weight map can be specified ***\n\n\n");
				exit(1);
			}

			if (i + 3 >= argc)
				not_enough("--wm");
			wm_filename = argv[i+1];

			if (sscanf(argv[i+2], "%d", &wm_offsetx) != 1) {
				fprintf(stderr, "\n\n*** --wm x-argument must be an integer ***\n\n\n");
				exit(1);
			}
			if (sscanf(argv[i+3], "%d", &wm_offsety) != 1) {
				fprintf(stderr, "\n\n*** --wm y-argument must be an integer ***\n\n\n");
				exit(1);
			}
			i += 3;

		} else if (!strcmp(argv[i], "--fl")) {
			if (i + 3 >= argc)
				not_enough("--fl");
			double h, v, a;
			if (sscanf(argv[i+1], "%lf", &h) != 1) {
				fprintf(stderr, "\n\n*** --fl h-argument must be numerical ***\n\n\n");
				exit(1);
			}
			if (sscanf(argv[i+2], "%lf", &v) != 1) {
				fprintf(stderr, "\n\n*** --fl v-argument must be numerical ***\n\n\n");
				exit(1);
			}
			if (sscanf(argv[i+3], "%lf", &a) != 1) {
				fprintf(stderr, "\n\n*** --fl a-argument must be numerical ***\n\n\n");
				exit(1);
			}
			i += 3;
#ifdef USE_FFTW
			d2::align::set_frequency_cut(h, v, a);
#else
			fprintf(stderr, "\n\n*** Error: --fl is not supported                  ***");
			fprintf(stderr,   "\n*** Hint:  rebuild ALE with FFTW=1                ***\n\n\n");
			exit(1);
#endif
		} else if (!strcmp(argv[i], "--wmx")) {
			if (i + 3 >= argc)
				not_enough("--wmx");
#ifdef USE_UNIX
			d2::align::set_wmx(argv[i+1], argv[i+2], argv[i+3]);
#else
			fprintf(stderr, "\n\n*** Error: --wmx is not supported                 ***");
			fprintf(stderr,   "\n*** Hint:  rebuild ALE with POSIX=1               ***");
#endif
			i += 3;
		} else if (!strcmp(argv[i], "--flshow")) {
			if (i + 1 >= argc)
				not_enough("--flshow");
			d2::align::set_fl_show(argv[i+1]);
			i++;
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
			fprintf(stderr, "\n\n*** Error: --ip <r> <i> is no longer supported. ***\n"
					    "*** Use --lpsf box=<r> --ips <i> instead.        ***\n\n");
			exit(1);
		} else if (!strcmp(argv[i], "--bayer")) {
			if (i + 1 >= argc)
				not_enough("--bayer");

			/*
			 * External order is clockwise from top-left.  Internal
			 * order is counter-clockwise from top-left.
			 */

			if (!strcmp(argv[i+1], "rgbg")) {
				d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
			} else if (!strcmp(argv[i+1], "bgrg")) {
				d2::image_rw::set_default_bayer(IMAGE_BAYER_BGRG);
			} else if (!strcmp(argv[i+1], "gbgr")) {
				d2::image_rw::set_default_bayer(IMAGE_BAYER_GRGB);
			} else if (!strcmp(argv[i+1], "grgb")) {
				d2::image_rw::set_default_bayer(IMAGE_BAYER_GBGR);
			} else if (!strcmp(argv[i+1], "none")) {
				d2::image_rw::set_default_bayer(IMAGE_BAYER_NONE);
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

			if (sscanf(argv[i+1], "%d", &ip_iterations) != 1) {
				fprintf(stderr, "\n\n*** --ipr requires an integer argument ***\n\n");
				exit(1);
			}
			fprintf(stderr, "\n\n*** Warning: --ipr is deprecated.  Use --ips instead ***\n\n");
			i++;

		} else if (!strcmp(argv[i], "--ips")) {

			if (i + 1 >= argc)
				not_enough("--ips");

			if (sscanf(argv[i+1], "%d", &ip_iterations) != 1) {
				fprintf(stderr, "\n\n*** --ips requires an integer argument ***\n\n");
				exit(1);
			}
			i++;
		} else if (!strcmp(argv[i], "--ipc")) {
			fprintf(stderr, "\n\n*** Error: --ipc <c> <i> is no longer supported. ***\n"
					    "*** Use either: --ips <i> --lpsf <c>              ***\n"
					    "***         or: --ips <i> --device <c>           ***\n\n");
			exit(1);
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
		} else if (!strcmp(argv[i], "--drizzle-only")) {
			unsupported::discontinued("--drizzle-only", "--dchain box:1");
		} else if (!strcmp(argv[i], "--inc")) {
			inc = 1;
		} else if (!strcmp(argv[i], "--no-inc")) {
			inc = 0;
		} else if (!strncmp(argv[i], "--visp-scale=", strlen("--visp-scale="))) {
			sscanf(argv[i] + strlen("--visp-scale="), "%lf", &vise_scale_factor);
			if (vise_scale_factor <= 0.0) {
				fprintf(stderr, "\n\n*** Error: VISP scale "
						"must be greater than zero. ***\n\n\n");
				exit(1);
			}
		} else if (!strncmp(argv[i], "--scale=", strlen("--scale="))) {
			sscanf(argv[i] + strlen("--scale="), "%lf", &scale_factor);
			if (scale_factor <= 1.0) {
				fprintf(stderr, "\n\n*** Error: Scale "
						"must be greater than one. ***\n\n\n");
				exit(1);
			}
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
			if (*(argv[i] + strlen("--perturb-lower=") + characters) == '%')
				d2::align::set_perturb_lower(perturb_lower, 1);
			else
				d2::align::set_perturb_lower(perturb_lower, 0);
		} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
			double perturb_lower;
			fprintf(stderr, "\n\n*** Warning: --stepsize is deprecated.  "
					"Use --perturb-lower instead. ***\n\n\n");
			sscanf(argv[i] + strlen("--stepsize="), "%lf", &perturb_lower);
			d2::align::set_perturb_lower(perturb_lower, 0);
		} else if (!strncmp(argv[i], "--hf-enhance=", strlen("--hf-enhance="))) {
			fprintf(stderr, "\n\n*** Error: --hf-enhance=<x> is no longer supported.  "
					"Use --usm <x> instead ***\n\n");
			exit(1);
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
		} else if (!strncmp(argv[i], "--trans-load=", strlen("--trans-load="))) {
			d2::tload_delete(tload);
			tload = d2::tload_new(argv[i] + strlen("--trans-load="));
			d2::align::set_tload(tload);
		} else if (!strncmp(argv[i], "--trans-save=", strlen("--trans-save="))) {
			tsave_delete(tsave);
			tsave = d2::tsave_new(argv[i] + strlen("--trans-save="));
			d2::align::set_tsave(tsave);
		} else {

			/*
			 * Trap illegal options and end-of-option indicators.
			 */

			if (!strcmp(argv[i], "--")) {
				i++;
			} else if (!strncmp(argv[i], "--", strlen("--"))) {
				fprintf(stderr, "\n\n*** Error: illegal option %s ***\n", argv[i]);
				exit(1);
			}

			/*
			 * Set exclusion regions
			 */

			d2::align::set_exclusion(ex_parameters, ex_count);
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
			 * Apply implication logic.
			 */
#if 0
			/*
			 * XXX: if this is really desired, it must be moved above render_init().
			 */
			if (extend == 0 && vise_count != 0) {
				implication::changed("VISP requires increased image extents.",
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
			 * Initialize device-specific variables
		         */

			d2::psf *device_response[psf_N] = { NULL, NULL };
			d2::exposure *input_exposure = NULL;
			ale_pos view_angle = 90 * M_PI / 180;

			if (device != NULL) {
				if (!strcmp(device, "xvp610_640x480")) {
					device_response[psf_linear] = new xvp610_640x480::lpsf();
					device_response[psf_nonlinear] = new xvp610_640x480::nlpsf();
					input_exposure = new xvp610_640x480::exposure[argc - i - 1];
					view_angle = xvp610_640x480::view_angle();
				} else if (!strcmp(device, "xvp610_320x240")) {
					device_response[psf_linear] = new xvp610_320x240::lpsf();
					device_response[psf_nonlinear] = new xvp610_320x240::nlpsf();
					input_exposure = new xvp610_320x240::exposure[argc - i - 1];
					view_angle = xvp610_320x240::view_angle();
				} else if (!strcmp(device, "ov7620_raw_linear")) {
					device_response[psf_linear] = new ov7620_raw_linear::lpsf();
					device_response[psf_nonlinear] = NULL;
					input_exposure = new ov7620_raw_linear::exposure[argc - i - 1];
					d2::image_rw::set_default_bayer(IMAGE_BAYER_BGRG);
				} else {
					fprintf(stderr, "\n\n*** Error: Unknown device %s ***\n\n", device);
					exit(1);
				}
			} else {
				input_exposure = new d2::exposure_default[argc - i - 1];
			}

			/*
			 * Initialize output exposure
			 */

			d2::exposure *output_exposure = new d2::exposure_default();

			/*
			 * Configure the response function.
			 */

			d2::psf *response[2] = {NULL, NULL};

			for (int n = 0; n < psf_N; n++ ) {
				if (psf[n] != NULL) while (psf[n] != NULL) {

					/*
					 * Iterate over substrings of the psf command-line
					 * argument separated by '+' characters.  We sum these
					 * response functions.
					 */
#if 0
					/*
					 * Not available everywhere.
					 */
					char *summation_index = index(psf[n], '+');
#else
					/*
					 * This mimics the functionality of index.
					 * Conversion to char * should be OK.  We
					 * use const char * elsewhere because this
					 * is the only place that the string should
					 * be modified.
					 */
					char *summation_index = (char *) psf[n];
					while(*summation_index != '\0'
					   && *summation_index != '+')
						summation_index++;
					if (*summation_index == '\0')
						summation_index = NULL;
					
#endif
					d2::psf *current_response = response[n];

					if (summation_index) {
						summation_index[0] = '\0';
						summation_index++;
					}

					/*
					 * User-specified point-spread function
					 */

					if (!strcmp(psf[n], "stdin")) {

						/*
						 * Standard input filter
						 */

						fprintf(stderr, "\nInitializing ");
						fprintf(stderr, (n == 0) ? "linear" : "non-linear");
						fprintf(stderr, " PSF.\n");
							

						response[n] = new d2::psf_stdin();

					} else if (!strcmp(psf[n], "stdin_vg")) {

						/*
						 * Standard input filter
						 */

						fprintf(stderr, "\nInitializing ");
						fprintf(stderr, (n == 0) ? "linear" : "non-linear");
						fprintf(stderr, " PSF.\n");
							

						response[n] = new d2::psf_stdin_vg();

					} else if (!strncmp(psf[n], "box=", strlen("box="))) {

						/*
						 * Box filter
						 */

						double box_diameter;

						if (sscanf(psf[n] + strlen("box="), "%lf", &box_diameter)
								!= 1) {
							fprintf(stderr, "\n\n*** Error: box= takes a "
									"numerical argument ***\n\n");
							exit(1);
						}

						response[n] = new d2::box(box_diameter / 2);

					} else if (!strncmp(psf[n], "circle=", strlen("circle="))) {

						/*
						 * Circular filter
						 */

						double diameter;

						if (sscanf(psf[n] + strlen("circle="), "%lf", &diameter)
								!= 1) {
							fprintf(stderr, "\n\n*** Error: circle= takes a "
									"numerical argument ***\n\n");
							exit(1);
						}

						response[n] = new d2::circle(diameter / 2);

					} else {
						fprintf(stderr, "Unknown point-spread function %s.\n\n", psf[n]);
						exit(1);
					}

					if (current_response)
						response[n] = new d2::sum(response[n], current_response);

					psf[n] = summation_index;

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

			tsave_orig(tsave, argv[i]);
			tsave_target(tsave, argv[argc - 1]);

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
						exposure_register);
			}

			/*
			 * Handle the original frame.
			 */

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			for (int opt = 0; opt < oc_count; opt++) {
				ochain[opt]->sync(0);
				if  (inc)
					d2::image_rw::write_image(ochain_names[opt], 
						ochain[opt]->get_image(0));
			}

			d2::vise_core::frame_queue_add(0);

			fprintf(stderr, ".\n");

			/*
			 * Handle supplemental frames.
			 */

			for (unsigned int j = 1; j < d2::image_rw::count(); j++) {

				const char *name = d2::image_rw::name(j);

				fprintf(stderr, "Merging supplemental frame '%s'", 
						name);

				/*
				 * Write comment information about the
				 * supplemental frame to the transformation
				 * save file, if we have one.
				 */

				tsave_info (tsave, name);

				for (int opt = 0; opt < oc_count; opt++) {
					ochain[opt]->sync(j);
					if (inc)
						d2::image_rw::write_image(ochain_names[opt], 
							ochain[opt]->get_image(j));
				}

				d2::vise_core::frame_queue_add(j);

				fprintf(stderr, ".\n");
			}

			/*
			 * Do any post-processing and output final image
			 *
			 * XXX: note that non-incremental renderers currently
			 * return zero for ochain[0]->sync(), since they write
			 * output internally when inc != 0.
			 */

			for (int opt = 0; opt < oc_count; opt++) 
			if  (ochain[opt]->sync() || !inc)
				d2::image_rw::write_image(ochain_names[opt], ochain[opt]->get_image());


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
			 * Output a summary match statistic.
			 */

			fprintf(stderr, "Done (%f%% average match).\n", (double) d2::align::match_summary());

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

