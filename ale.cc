// Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

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

#include "image.h"
#include "image_rw.h"
#include "gpt.h"
#include "my_real.h"
#include "tfile.h"
#include "render.h"
#include "combine.h"
#include "align.h"
#include "merge.h"
#include "drizzle.h"
#include "hf_filter.h"
#include "ip.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

/*
 * Version Information
 */

char *version = "0.4.4"
#ifdef USE_MAGICK
		" (File handler: ImageMagick)";
#else
		" (File handler: PPM binary)";
#endif


/*
 * Describe how to use this program
 */
inline void usage(const char *argv0) {
	fprintf(stderr, 
		"\n"
		"Usage: %s [<options>] <original-frame> [<supplemental-frame> ...] <output-file>\n"
		"   or: %s --version\n"
		"\n\n"
		"Alignment channel options:\n\n"
		"--align-all       Align images using all color channels\n"
		"--align-green     Align images using the green channel\n"
		"--align-sum       Align images using a sum of channels [default]\n"
		"\n\n"
		"Transformation options:\n\n"
		"--translation     Only adjust the position of images\n"
		"--euclidean       Adjust the position and orientation of images [default]\n"
		"--projective      Use projective transformations.  Best quality, but slow.\n"
		"\n\n"
		"Image extents:\n\n"
		"--extend          Increase image extents to accommodate all pixel data.\n"
		"--no-extend       Don't increase extents; crop to original frame. [default]\n"
		"\n\n"
		"Transformation defaults:\n\n"
		"--identity        Default alignment begins with identity transform. [default]\n"
		"--follow          Default alignment begins with the previous frame's transform.\n"
		"\n\n"
		"Transformation file operations:\n\n"
		"--trans-load=x    Load initial transformation settings from file x\n"
		"--trans-save=x    Save final transformation data in file x\n"
		"\n\n"
		"Tunable parameters:\n\n"
		"--scale=x         Scale images by the factor x (where x is at least 1.0)\n"
		"--hf-enhance=x    Enhance high frequency details by factor x. (0.0 is default)\n"
		"--metric=x        Set the error metric exponent (2 is default)\n"
		"--threshold=x     Min. match threshold; a perfect match is 100.  (0 is default)\n"
		"--perturb-upper=x Max. correction step, in pixels or degrees (32.0 is default)\n"
		"--perturb-lower=x Min. correction step, in pixels or degrees (0.125 is default)\n"
		"\n\n"
		"Drizzling:\n\n"
		"--drizzle-diam=x  Drizzle with input pixel diameter x (where 0 < x <= 1).\n"
		"--drizzle-only    If drizzling, output black for pixels with no drizzle data.\n"
		"\n\n"
		"Image reconstruction:\n\n"
		"--ip <d> <i>      Irani-Peleg solve with pixel diameter <d> for <i> iterations.\n"
		"\n\n"
		"Monte Carlo alignment:\n\n"
		"--mc <x>          Align using, on average, x%% of available pixels (0 < x < 100)\n"
		"--no-mc           Align using all pixels.  [default]\n"
		"\n\n"
		"File output options:\n\n"
		"--inc             Produce incremental image output.  [default]\n"
		"--no-inc          Don't produce incremental image output.\n"
		"\n",
		argv0, argv0);
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

	if (argc == 2 && !strcmp(argv[1], "--version")) {

		/*
		 * Output the version
		 */

		fprintf(stderr, "%s\n", version);

		return 0;
	}

	/*
	 * Flags
	 */

	double scale_factor = 1;
	double hf_enhance = 0.0;
	double drizzle_radius = -1;
	int drizzle_only = 0;
	int extend = 0;
	struct tload_t *tload = NULL;
	struct tsave_t *tsave = NULL;
	int ip_iterations = 0;
	double ip_radius = 0;
	int inc = 1;

	/* 
	 * Iterate through arguments until we reach the first file 
	 * argument.  After the first file argument, we assume that
	 * all following arguments are files.
	 */

	for (int i = 1; i < argc - 1; i++) {


		if (!strcmp(argv[i], "--scale2")) {
			scale_factor = 2;
			fprintf(stderr, "\n\n*** Warning: --scale2 is deprecated.  "
					"Use --scale=2 instead. ***\n\n\n");
		} else if (!strcmp(argv[i], "--scale4")) {
			scale_factor = 4;
			fprintf(stderr, "\n\n*** Warning: --scale4 is deprecated.  "
					"Use --scale=4 instead. ***\n\n\n");
		} else if (!strcmp(argv[i], "--scale8")) {
			scale_factor = 8;
			fprintf(stderr, "\n\n*** Warning: --scale8 is deprecated.  "
					"Use --scale=8 instead. ***\n\n\n");
		} else if (!strcmp(argv[i], "--align-all")) {
			align::all();
		} else if (!strcmp(argv[i], "--align-green")) {
			align::green();
		} else if (!strcmp(argv[i], "--align-sum")) {
			align::sum();
		} else if (!strcmp(argv[i], "--translation")) {
			align::class_translation();
		} else if (!strcmp(argv[i], "--euclidean")) {
			align::class_euclidean();
		} else if (!strcmp(argv[i], "--projective")) {
			align::class_projective();
		} else if (!strcmp(argv[i], "--identity")) {
			align::initial_default_identity();
		} else if (!strcmp(argv[i], "--follow")) {
			align::initial_default_follow();
		} else if (!strcmp(argv[i], "--no-extend")) {
			align::no_extend();
			extend = 0;
		} else if (!strcmp(argv[i], "--extend")) {
			align::set_extend();
			extend = 1;
		} else if (!strcmp(argv[i], "--no-mc")) {
			align::no_mc();
		} else if (!strcmp(argv[i], "--mc")) {
			double mc_parameter;
			if (i + 1 < argc) {
				sscanf(argv[i+1], "%lf", &mc_parameter);
				mc_parameter /= 100;
				i += 1;
				align::mc(mc_parameter);
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for --mc ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--ip")) {
			if (i + 2 < argc) {
				sscanf(argv[i+1], "%lf", &ip_radius);
				sscanf(argv[i+2], "%d", &ip_iterations);
				ip_radius /= 2;
				i += 2;
				align::keep();
			} else {
				fprintf(stderr, "\n\n*** Not enough arguments for IP-solve ***\n\n\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--drizzle-only")) {
			drizzle_only = 1;
		} else if (!strcmp(argv[i], "--inc")) {
			inc = 1;
		} else if (!strcmp(argv[i], "--no-inc")) {
			inc = 0;
		} else if (!strncmp(argv[i], "--scale=", strlen("--scale="))) {
			sscanf(argv[i] + strlen("--scale="), "%lf", &scale_factor);
			if (scale_factor < 1.0) {
				fprintf(stderr, "\n\n*** Warning: Ignoring scale value "
						"smaller than 1.0 ***\n\n\n");
			}
		} else if (!strncmp(argv[i], "--metric=", strlen("--metric="))) {
			double metric;
			sscanf(argv[i] + strlen("--metric="), "%lf", &metric);
			align::set_metric_exponent(metric);
		} else if (!strncmp(argv[i], "--threshold=", strlen("--threshold="))) {
			double match_threshold;
			sscanf(argv[i] + strlen("--threshold="), "%lf", &match_threshold);
			align::set_match_threshold(match_threshold);
		} else if (!strncmp(argv[i], "--drizzle-diam=", strlen("--drizzle-diam="))) {
			sscanf(argv[i] + strlen("--drizzle-diam="), "%lf", &drizzle_radius);
			drizzle_radius /= 2;
		} else if (!strncmp(argv[i], "--perturb-upper=", strlen("--perturb-upper="))) {
			double perturb_upper;
			sscanf(argv[i] + strlen("--perturb-upper="), "%lf", &perturb_upper);
			align::set_perturb_upper(perturb_upper);
		} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
			double perturb_lower;
			fprintf(stderr, "\n\n*** Warning: --stepsize is deprecated.  "
					"Use --perturb-lower instead. ***\n\n\n");
			sscanf(argv[i] + strlen("--stepsize="), "%lf", &perturb_lower);
			align::set_perturb_lower(perturb_lower);
		} else if (!strncmp(argv[i], "--hf-enhance=", strlen("--hf-enhance="))) {
			sscanf(argv[i] + strlen("--hf-enhance="), "%lf", &hf_enhance);
		} else if (!strncmp(argv[i], "--perturb-lower=", strlen("--perturb-lower="))) {
			double perturb_lower;
			sscanf(argv[i] + strlen("--perturb-lower="), "%lf", &perturb_lower);
			align::set_perturb_lower(perturb_lower);
		} else if (!strncmp(argv[i], "--trans-load=", strlen("--trans-load="))) {
			tload_delete(tload);
			tload = tload_new(argv[i] + strlen("--trans-load="));
			align::set_tload(tload);
		} else if (!strncmp(argv[i], "--trans-save=", strlen("--trans-save="))) {
			tsave_delete(tsave);
			tsave = tsave_new(argv[i] + strlen("--trans-save="));
			align::set_tsave(tsave);
		} else {

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

			image_rw::init(argc - i - 1, argv + i, argv[argc - 1], 
					scale_factor);

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

			merge *_merge = new merge(extend);

			/*
			 * Use merged renderings as reference images in
			 * alignment.
			 */

			align::set_reference(_merge);

			/*
			 * Configure the output renderer.
			 */

			render *_render = _merge;

			if (drizzle_radius > 0) {
				render *_drizzle = new drizzle(extend,
						drizzle_radius, scale_factor);

				if (drizzle_only)
					_render = _drizzle;
				else
					_render = new combine(_render, _drizzle);
			}

			if (hf_enhance != 0)
				_render = new hf_filter(_render, scale_factor, hf_enhance);

			if (ip_iterations != 0)
				_render = new ip(_render, scale_factor,
						ip_radius, ip_iterations);

			/*
			 * Handle the original frame.
			 */

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			_render->operator()(0);

			if (inc)
				image_rw::output(_render->get_image());

			fprintf(stderr, ".\n");

			/*
			 * Handle supplemental frames.
			 */

			for (int j = 1; j < image_rw::count(); j++) {

				const char *name = image_rw::name(j);

				fprintf(stderr, "Merging supplemental frame '%s'", 
						name);

				/*
				 * Write comment information about the
				 * supplemental frame to the transformation
				 * save file, if we have one.
				 */

				tsave_info (tsave, name);

				_render->operator()(j);

				if (inc)
					image_rw::output(_render->get_image());

				fprintf(stderr, ".\n");
			}

			/*
			 * Do any post-processing and output final image
			 */

			if (_render->operator()() || !inc)
				image_rw::output(_render->get_image());

			/*
			 * Destroy the image file handler
			 */

			image_rw::destroy();

			/*
			 * Delete the transformation file structures, if any
			 * exist.
			 */

			tsave_delete(tsave);
			tload_delete(tload);

			/*
			 * Output a summary match statistic.
			 */

			fprintf(stderr, "Done (%f%% average match).\n", align::match_summary());

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

