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
 * 2D include files
 */

#include "d2.h"

/*
 * Device configuration files
 */

#include "device/xvp610_320x240.h"

/*
 * Version Information
 */

char *version = "ALE Version:      0.6.0\n"
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
 * Describe how to use this program
 */
inline void usage(const char *argv0) {
#define BETWEEN_SECTIONS "\n"
#define HEADER_SPACE ""
	fprintf(stderr, 
		"\n"
		"Usage: %s [<options>] <original-frame> [<supplemental-frame> ...] <output-file>\n"
		"   or: %s --version\n"
		BETWEEN_SECTIONS
		"File output options:\n"
		HEADER_SPACE
		"--8bpc            Write 8 bit per channel output [default]\n"
		"--16bpc           Write 16 bit per channel output\n"
		"\n"
#ifdef USE_MAGICK
		"--auto            Determine file type automatically [default]\n"
		"--raw             Write raw PPM output\n"
		"--plain           Write plain PPM output\n"
#else
		"--raw             Write raw PPM output [default]\n"
		"--plain           Write plain PPM output\n"
#endif
		BETWEEN_SECTIONS
		"Alignment channel options:\n"
		HEADER_SPACE
		"--align-all       Align images using all color channels\n"
		"--align-green     Align images using the green channel\n"
		"--align-sum       Align images using a sum of channels [default]\n"
		BETWEEN_SECTIONS
		"Transformation options:\n"
		HEADER_SPACE
		"--translation     Only adjust the position of images\n"
		"--euclidean       Adjust the position and orientation of images [default]\n"
		"--projective      Use projective transformations.  Best quality, but slow.\n"
		BETWEEN_SECTIONS
		"Image extents:\n"
		HEADER_SPACE
		"--extend          Increase image extents to accommodate all pixel data.\n"
		"--no-extend       Don't increase extents; crop to original frame. [default]\n"
		BETWEEN_SECTIONS
		"Alignment following:\n"
		HEADER_SPACE
                "--identity        Frames align closely with the original frame.  [default]\n"
                "--follow          Frames align closely with their immediate predecessor.\n"
		BETWEEN_SECTIONS
		"Transformation file operations:\n"
		HEADER_SPACE
		"--trans-load=x    Load initial transformation settings from file x\n"
		"--trans-save=x    Save final transformation data in file x\n"
		BETWEEN_SECTIONS
		"Tunable parameters:\n"
		HEADER_SPACE
		"--scale=x         Scale images by the factor x (where x is at least 1.0)\n"
		"--metric=x        Set the alignment error metric exponent.       (2 is default)\n"
		"--threshold=x     Min. match threshold; a perfect match is 100.  (0 is default)\n"
		"--perturb-upper=x Perturbation upper bound in pixels/degrees  (32.0 is default)\n"
		"--perturb-lower=x Perturbation lower bound in pixels/degrees  (.125 is default)\n"
		"--rot-upper=x     Rotation-specific perturbation upper bound  (32.0 is default)\n"
		"--lod-max=x       LOD scale factor is max(1, (2^floor(x))/perturb)  (1 is def.)\n"
		BETWEEN_SECTIONS
		"Drizzling:\n"
		HEADER_SPACE
		"--drizzle-diam=x  Drizzle with input pixel diameter x (where x > 0)\n"
		"--drizzle-only    If drizzling, output black for pixels with no drizzle data.\n"
		BETWEEN_SECTIONS
		"Point-spread functions:\n"
		HEADER_SPACE
		"--lpsf <p>        Set linear colorspace point-spread function to <p>\n"
		"--nlpsf <p>       Set non-linear colorspace point-spread function to <p>\n"
		"                     Available point-spread functions:\n"
		"                        box=<diameter>\n"
		"                        stdin\n"
		"                        <p>+<p> (summation)\n"
		"                     Default lpsf is either 'box=1.0' or device-specific.\n"
		"                     Default nlpsf is either disabled or device-specific.\n"
		"--psf-match       Used by the script ale-psf-calibrate to evaluate PSFs.\n"
		BETWEEN_SECTIONS
		"Device (sets point-spread function defaults):\n"
		HEADER_SPACE
		"--device <d>      Set the capture device to <d>.\n"
		"                     Available devices:\n"
		"                        xvp610_320x240\n"
		BETWEEN_SECTIONS
		"Unsharp Mask (was 'High-frequency Enhancement'):\n"
		HEADER_SPACE
		"--usm <m>         Apply an unsharp mask with multiplier <m>.\n"
		"                     (See also --device, --nlpsf, and --lpsf.)\n"
		BETWEEN_SECTIONS
		"Irani-Peleg iterative solver:\n"
		HEADER_SPACE
		"--ips <i>         Run <i> iterations.  (see also --device, --nlpsf, and --lpsf)\n"
		BETWEEN_SECTIONS
		"Monte Carlo alignment:\n"
		HEADER_SPACE
		"--mc <x>          Align using, on average, x%% of available pixels (0 < x < 100)\n"
		"--no-mc           Align using all pixels.  [default]\n"
		BETWEEN_SECTIONS
		"Certainty-based rendering [Experimental]:\n"
		HEADER_SPACE
		"--cx <x>          Render with certainty exponent <x>\n"
		"--no-cx           Render with uniform certainty.  [default]\n"
		BETWEEN_SECTIONS
		"Exposure options:\n"
		HEADER_SPACE
		"--exp-register    Register exposure between frames.  [default]\n"
		"--exp-noregister  Assume uniform exposure across all frames.\n"
		BETWEEN_SECTIONS
		"File output options:\n"
		HEADER_SPACE
		"--inc             Produce incremental output.  [default]\n"
		"--no-inc          Don't produce incremental output.\n"
		BETWEEN_SECTIONS
		"Pixel replacement:\n"
		HEADER_SPACE
		"--replace         Replace overlapping areas when merging and drizzling.\n"
		"--no-replace      Do not replace.  [default]\n"
		"\n",
		argv0, argv0);
	exit(1);

#undef BETWEEN_SECTIONS
#undef HEADER_SPACE
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
	 * Version information
	 */

	if (argc == 2 && !strcmp(argv[1], "--version")) {

		/*
		 * Output the version
		 */

		fprintf(stderr, "%s", version);

		return 0;
	}

	/*
	 * Undocumented projective transformation utility
	 */

	if (argc == 2 && !strcmp(argv[1], "--ptcalc")) {
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

		b = t(a);

		printf("Result: (%f, %f)\n", (double) b[1], (double) b[0]);

		exit(0);
	}

	/*
	 * Flags
	 */

	double scale_factor = 1;
	double usm_multiplier = 0.0;
	double drizzle_radius = -1;
	int drizzle_only = 0;
	int extend = 0;
	struct d2::tload_t *tload = NULL;
	struct d2::tsave_t *tsave = NULL;
	int ip_iterations = 0;
	enum { psf_linear, psf_nonlinear, psf_N };
	const char *psf[psf_N] = {NULL, NULL};
	const char *device = NULL;
	int psf_match = 0;
	int inc = 1;
	int replace = 0;
	int exposure_register = 1;

	/* 
	 * Iterate through arguments until we reach the first file 
	 * argument.  After the first file argument, we assume that
	 * all following arguments are files.
	 */

	for (int i = 1; i < argc - 1; i++) {


		if (!strcmp(argv[i], "--8bpc")) {
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
		} else if (!strcmp(argv[i], "--no-extend")) {
			d2::align::no_extend();
			extend = 0;
		} else if (!strcmp(argv[i], "--extend")) {
			d2::align::set_extend();
			extend = 1;
		} else if (!strcmp(argv[i], "--no-mc")) {
			d2::align::no_mc();
		} else if (!strcmp(argv[i], "--mc")) {
			double mc_parameter;
			if (i + 1 < argc) {
				sscanf(argv[i+1], "%lf", &mc_parameter);
				mc_parameter /= 100;
				i += 1;
				d2::align::mc(mc_parameter);
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --mc ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--cx")) {
			double cx_parameter;
			if (i + 1 < argc) {
				sscanf(argv[i+1], "%lf", &cx_parameter);
				i += 1;
				d2::exposure::set_confidence(cx_parameter);
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --cx ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--no-cx")) {
			d2::exposure::set_confidence(0);
		} else if (!strcmp(argv[i], "--ip")) {
			fprintf(stderr, "\n\n*** Error: --ip <r> <i> is no longer supported. ***\n"
					    "*** Use --lpsf box=<r> --ips <i> instead.        ***\n\n");
			exit(1);
		} else if (!strcmp(argv[i], "--lpsf")) {
			if (i + 1 < argc) {
				psf[psf_linear] = argv[i+1];
				i++;
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --lpsf ***\n\n\n");
				exit(1);
			}

		} else if (!strcmp(argv[i], "--nlpsf")) {
			if (i + 1 < argc) {
				psf[psf_nonlinear] = argv[i+1];
				i++;
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --nlpsf ***\n\n\n");
				exit(1);
			}

		} else if (!strcmp(argv[i], "--psf-match")) {
			psf_match = 1;
			d2::align::keep();
		} else if (!strcmp(argv[i], "--device")) {
			if (i + 1 < argc) {
				device = argv[i+1];
				i++;
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --device ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--usm")) {
			if (i + 1 < argc) {
				sscanf(argv[i+1], "%lf", &usm_multiplier);
				i++;
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --usm ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--ipr")) {
			if (i + 1 < argc) {
				if (sscanf(argv[i+1], "%d", &ip_iterations) != 1) {
					fprintf(stderr, "\n\n*** --ipr requires an integer argument ***\n\n");
					exit(1);
				}
				fprintf(stderr, "\n\n*** Warning: --ipr is deprecated.  Use --ips instead ***\n\n");
				i++;
				d2::align::keep();
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --ipr ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--ips")) {
			if (i + 1 < argc) {
				if (sscanf(argv[i+1], "%d", &ip_iterations) != 1) {
					fprintf(stderr, "\n\n*** --ips requires an integer argument ***\n\n");
					exit(1);
				}
				i++;
				d2::align::keep();
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --ips ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--ipc")) {
			fprintf(stderr, "\n\n*** Error: --ipc <c> <i> is no longer supported. ***\n"
					    "*** Use either: --ips <i> --lpsf <c>              ***\n"
					    "***         or: --ips <i> --device <c>           ***\n\n");
			exit(1);
		} else if (!strcmp(argv[i], "--exp-register")) {
			exposure_register = 1;
			d2::align::exp_register();
		} else if (!strcmp(argv[i], "--exp-noregister")) {
			exposure_register = 0;
			d2::align::exp_noregister();
		} else if (!strcmp(argv[i], "--drizzle-only")) {
			drizzle_only = 1;
		} else if (!strcmp(argv[i], "--inc")) {
			inc = 1;
		} else if (!strcmp(argv[i], "--no-inc")) {
			inc = 0;
		} else if (!strcmp(argv[i], "--replace")) {
			replace = 1;
		} else if (!strcmp(argv[i], "--no-replace")) {
			replace = 0;
		} else if (!strncmp(argv[i], "--scale=", strlen("--scale="))) {
			sscanf(argv[i] + strlen("--scale="), "%lf", &scale_factor);
			if (scale_factor < 1.0) {
				fprintf(stderr, "\n\n*** Error: Scale "
						"smaller than 1.0 is not supported ***\n\n\n");
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
			sscanf(argv[i] + strlen("--drizzle-diam="), "%lf", &drizzle_radius);
			drizzle_radius /= 2;
		} else if (!strncmp(argv[i], "--perturb-upper=", strlen("--perturb-upper="))) {
			double perturb_upper;
			sscanf(argv[i] + strlen("--perturb-upper="), "%lf", &perturb_upper);
			d2::align::set_perturb_upper(perturb_upper);
		} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
			double perturb_lower;
			fprintf(stderr, "\n\n*** Warning: --stepsize is deprecated.  "
					"Use --perturb-lower instead. ***\n\n\n");
			sscanf(argv[i] + strlen("--stepsize="), "%lf", &perturb_lower);
			d2::align::set_perturb_lower(perturb_lower);
		} else if (!strncmp(argv[i], "--hf-enhance=", strlen("--hf-enhance="))) {
			fprintf(stderr, "\n\n*** Error: --hf-enhance=<x> is no longer supported.  "
					"Use --usm <x> instead ***\n\n");
			exit(1);
		} else if (!strncmp(argv[i], "--perturb-lower=", strlen("--perturb-lower="))) {
			double perturb_lower;
			sscanf(argv[i] + strlen("--perturb-lower="), "%lf", &perturb_lower);
			d2::align::set_perturb_lower(perturb_lower);
		} else if (!strncmp(argv[i], "--rot-upper=", strlen("--rot-upper="))) {
			double rot_max;
			sscanf(argv[i] + strlen("--rot-upper="), "%lf", &rot_max);
			d2::align::set_rot_max((int) floor(rot_max));
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
			 * Initialize device-specific variables
		         */

			d2::psf *device_response[psf_N] = { NULL, NULL };
			d2::exposure *input_exposure = NULL;

			if (device != NULL) {
				if (!strcmp(device, "xvp610_320x240")) {
					device_response[psf_linear] = new xvp610_320x240::lpsf();
					device_response[psf_nonlinear] = new xvp610_320x240::nlpsf();
					input_exposure = new xvp610_320x240::exposure[argc - i - 1];
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

			if (i >= argc - 1)
				usage(argv[0]);

			d2::image_rw::init(argc - i - 1, argv + i, argv[argc - 1], input_exposure, output_exposure);

			/*
			 * Write comment information about original frame and
			 * target image to the transformation save file, if we
			 * have one.
			 */

			tsave_orig(tsave, argv[i]);
			tsave_target(tsave, argv[argc - 1]);

			/*
			 * Create merge renderer.
			 */

			d2::render *_merge;
			
			if (replace == 0)
				_merge = new d2::merge<0>(extend, scale_factor);
			else
				_merge = new d2::merge<1>(extend, scale_factor);

			/*
			 * Use merged renderings as reference images in
			 * alignment.
			 */

			d2::align::set_reference(_merge);

			/*
			 * Tell the alignment class about the scale factor.
			 */

			d2::align::set_scale(scale_factor);

			/*
			 * Configure the output renderer.
			 */

			d2::render *_render = _merge;

			if (drizzle_radius > 0) {
				d2::render *_drizzle;
				
				if (replace == 0)
					_drizzle = new d2::drizzle<0>(extend,
						drizzle_radius, scale_factor);
				else
					_drizzle = new d2::drizzle<1>(extend,
						drizzle_radius, scale_factor);

				if (drizzle_only)
					_render = _drizzle;
				else
					_render = new d2::combine(_render, _drizzle);
			}



			if (usm_multiplier != 0) {

				/*
				 * Unsharp Mask renderer
				 */

				_render = new d2::usm(_render, scale_factor,
						usm_multiplier, inc, response[psf_linear],
						response[psf_nonlinear], &input_exposure[0]);
			}

			if (psf_match) {
				
				/*
				 * Point-spread function calibration renderer.
				 * This renderer does not produce image output.
				 * It is reserved for use with the point-spread
				 * function calibration script
				 * ale-psf-calibrate.
				 */

				_render = new d2::psf_calibrate(_render,
						1, inc, response[psf_linear],
						response[psf_nonlinear]);

			} else if (ip_iterations != 0) {

				/*
				 * Irani-Peleg renderer
				 */

				_render = new d2::ipc( _render, ip_iterations,
						inc, response[psf_linear],
						response[psf_nonlinear],
						exposure_register);
			}

			/*
			 * Handle the original frame.
			 */

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			_render->sync(0);

			if (inc)
				d2::image_rw::output(_render->get_image());

			fprintf(stderr, ".\n");

			/*
			 * Handle supplemental frames.
			 */

			for (int j = 1; j < d2::image_rw::count(); j++) {

				const char *name = d2::image_rw::name(j);

				fprintf(stderr, "Merging supplemental frame '%s'", 
						name);

				/*
				 * Write comment information about the
				 * supplemental frame to the transformation
				 * save file, if we have one.
				 */

				tsave_info (tsave, name);

				_render->sync(j);

				if (inc)
					d2::image_rw::output(_render->get_image());

				fprintf(stderr, ".\n");
			}

			/*
			 * Do any post-processing and output final image
			 *
			 * XXX: note that non-incremental renderers currently
			 * return zero for _render->sync(), since they write
			 * output internally when inc != 0.
			 */

			if (_render->sync() || !inc)
				d2::image_rw::output(_render->get_image());

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

	usage(argv[0]);
}

