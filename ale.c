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
#include "hf-filter.h"
#include "image_rw.h"
#include "gpt.h"
#include "my_real.h"
#include "tfile.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

char *version = "0.3.1"
#ifdef USE_MAGICK
		" (File handler: ImageMagick)";
#else
		" (File handler: PPM binary)";
#endif

/*
 * Various character arrays with atavistic names.
 */

image *display_image, *input_image, *weight_image;

/*
 * Global flags
 */

static int transform_code = 1;
static int align_code = 2;
static double metric = 2;
static double match_threshold = 0;
static double perturb_lower = 0.125;
static double perturb_upper = 32;
static double hf_enhance = 0.0;
static double match_sum = 0;
static int match_count = 0;
static struct tload_t *tload = NULL;
static struct tsave_t *tsave = NULL;

/*
 * Not-quite-symmetric difference function.  Determines the difference in areas
 * where the arrays overlap.  Uses the first array's notion of pixel positions.
 */
inline double diff(image *a, image *b, struct gpt t) {
	float result = 0;
	float divisor = 0;
	float apm_a = (float) avg_pixel_magnitude(a);
	float apm_b = (float) avg_pixel_magnitude(b);

	int i, j, k;

	for (i = 0; i < height(a); i++)
		for (j = 0; j < width(a); j++) {

			/*
			 * Transform
			 */

			my_real ti = (t.a * i + t.b * j + t.c)
				   / (t.g * i + t.h * j + 1  );

			my_real tj = (t.d * i + t.e * j + t.f)
				   / (t.g * i + t.h * j + 1  );

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of array b.
			 */

			if (ti >= 0
			 && ti <= height(b) - 1
			 && tj >= 0
			 && tj <= width(b) - 1){ 

				if (align_code == 0) {
					/*
					 * Align based on all channels.
					 */

					for (k = 0; k < 3; k++) {
						float achan = get_pixel_component(a, i, j, k) / apm_a;
						float bchan = get_bl_component(b, ti, tj, k) / apm_b;

						result += pow(fabs(achan - bchan), metric);
						divisor += pow(achan > bchan ? achan : bchan, metric);
					}
				} else if (align_code == 1) {
					/*
					 * Align based on the green channel.  XXX: apm_* shouldn't be used here.
					 */

					float achan = get_pixel_component(a, i, j, 1) / apm_a;
					float bchan = get_bl_component(b, ti, tj, 1) / apm_b;

					result += pow(fabs(achan - bchan), metric);
					divisor += pow(achan > bchan ? achan : bchan, metric);
				} else if (align_code == 2) {
					/*
					 * Align based on the sum of all channels.
					 */

					float asum = 0;
					float bsum = 0;

					for (k = 0; k < 3; k++) {
						asum += get_pixel_component(a, i, j, k) / apm_a;
						bsum += get_bl_component(b, ti, tj, k) / apm_b;
					}

					result += pow(fabs(asum - bsum), metric);
					divisor += pow(asum > bsum ? asum : bsum, metric);
				}
			}
		}

	return pow(result / divisor, 1/metric);
}

/*
 * Merge part of a delta frame with part of a target frame using the specified
 * offset.
 */
inline void
merge(image *target, image *delta, struct gpt t) {
	int i, j, k;

	for (i = 0; i < height(target); i++)
		for (j = 0; j < width(target); j++) {

			/*
			 * Transform
			 */

			my_real ti = (t.a * i + t.b * j + t.c)
				   / (t.g * i + t.h * j + 1  );

			my_real tj = (t.d * i + t.e * j + t.f)
				   / (t.g * i + t.h * j + 1  );

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of the delta frame.
			 */

			if (ti >= 0
			 && ti <= height(delta) - 1
			 && tj >= 0
			 && tj <= width(delta) - 1){ 

				/*
				 * Determine and update merging weight at this pixel
				 */

				int weight = get_pixel_component(weight_image, i, j, 0);
				set_pixel_component(weight_image, i, j, 0, ++weight);

				/*
				 * Update each channel
				 */
				
				for (k = 0; k < 3; k++)
					set_pixel_component(target, i, j, k,
						(weight * get_pixel_component(target, i, j, k)
					     +   get_bl_component(delta, ti, tj, k)) / (weight + 1));
			}
		}
}

void filter(image *target, image *source, double res_scale, double hf_enhance) {
	int i, j, k;

	assert (height(target) == height(source));
	assert (width(target) == width(source));
	assert (depth(target) == depth(source));

	for (i = 0; i < height(target); i++)
		for (j = 0; j < width(target); j++) 
			for (k = 0; k < 3; k++) {
				int result = 0.01*hf_enhance*hf_filter(res_scale*2.5,i,j,k, source)
					   + get_pixel_component(source, i, j, k);

				if (result < 0)
					result = 0;
				if (result > 255)
					result = 255;

				set_pixel_component(target, i, j, k, result);
			}
}


/*
 * Update an accumulated image with a new input image.
 */
inline void update() {
	double perturb = pow(2, floor(log(perturb_upper) / log(2)));
	int w = width(input_image), h = height(input_image);
	struct gpt offset;
	double here;
	int lod;

	/*
	 * Determine how many whole-pixel adjustment steps will occur
	 */

	int steps = (perturb >= 1) ? (int) (log(perturb) / log(2)) + 1 : 1;
	int step;
	image **display_scales = (image **) malloc(steps * sizeof(image *));
	image **input_scales = (image **) malloc(steps * sizeof(image *));

	/*
	 * Level-of-detail difference.  Use a level of detail 2^lod_diff finer
	 * than the adjustment resolution.  Empirically, it's okay to use a
	 * level-of-detail equal to twice the resolution of the perturbation,
	 * so we set lod_diff to 1, as 2^1==2.
	 */

	const int lod_diff = 1;

	/*
	 * Prepare multiple levels of detail.
	 */

	display_scales[0] = clone(display_image);
	input_scales[0] = clone(input_image);
	
	for (step = 1; step < steps; step++) {
		display_scales[step] = clone(display_scales[step - 1]);
		scale(display_scales[step], 0.5);
		input_scales[step] = clone(input_scales[step - 1]);
		scale(input_scales[step], 0.5);
	}

	/*
	 * Initialize variables used in the main loop.
	 */

	lod = (steps - 1) - lod_diff;
	lod = (lod < 0) ? 0 : lod;

	offset = tload_next(tload, w, h, lod, transform_code == 2);

	here = diff(display_scales[lod], input_scales[lod], offset);

	/*
	 * Simulated annealing perturbation adjustment loop.  
	 */

	while (perturb >= perturb_lower) {

		image *di = display_scales[lod];
		image *ii = input_scales[lod];

		double adj_p = (perturb >= pow(2, lod_diff)) 
			     ? pow(2, lod_diff) : perturb;
		double adj_s;

		struct gpt test_t;
		double test_d;
		double old_here = here;

		int i, j;

		if (transform_code < 2 && transform_code >= 0) {

			/* 
			 * Translational or euclidean transformation
			 */

			for (i = 0; i < 2 + transform_code; i++)
			for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

				test_t = eu_resultant(
					eu_modify(offset, i, adj_s));
				test_d = diff(di, ii, test_t);

				if (test_d < here) {
					here = test_d;
					offset = test_t;
					goto done;
				}
			}
		
		} else if (transform_code == 2) {

			/*
			 * Projective transformation
			 */

			for (i = 0; i < 4; i++)
			for (j = 0; j < 2; j++)
			for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

				test_t = gpt_resultant(
					gpt_modify(offset, j, i, adj_s));
				test_d = diff(di, ii, test_t);

				if (test_d < here) {
					here = test_d;
					offset = test_t;
					adj_s += 3 * adj_p;
				}
			}

		} else assert(0);

done:
		if (!(here < old_here)) {
			perturb *= 0.5;

			if (lod > 0) {
				int i, j;

				for (i = 0; i < 2; i++)
					for (j = 0; j < 4; j++)
						offset.x[i][j] *= 2;

				for (i = 0; i < 2; i++)
					offset.eu[i] *= 2;

				offset.input_width *= 2;
				offset.input_height *= 2;

				offset = gpt_resultant(offset);

				lod--;

				if (perturb >= perturb_lower)
					here = diff(display_scales[lod], 
						input_scales[lod], offset);
			}

			/*
			 * Announce that we've dropped a perturbation level.
			 */

			fprintf(stderr, ".");
		}
	}

	/*
	 * Free the level-of-detail structures
	 */
	for (step = 0; step < steps; step++) {
		image_final(display_scales[step]);
		image_final(input_scales[step]);
		free(display_scales[step]);
		free(input_scales[step]);
	}
	free(display_scales);
	free(input_scales);

	/*
	 * Save the transformation information
	 */

	tsave_next(tsave, offset, transform_code == 2);

	/*
	 * Ensure that the match meets the threshold.
	 */

	if ((1 - here) * 100 > match_threshold) {
		merge (display_image, input_image, offset);
		fprintf(stderr, " okay (%f%% match)", (1 - here) * 100);
	} else {
		fprintf(stderr, " no match (%f%% match)", (1 - here) * 100);
	}

	match_sum += (1 - here) * 100;
	match_count++;
}

int main(int argc, char *argv[]){
	double res_scale = 1;
	int i;

	if (argc == 2 && !strcmp(argv[1], "--version")) {

		/*
		 * Output the version
		 */

		fprintf(stderr, "%s\n", version);

		return 0;
	}

	/* 
	 * Iterate through all arguments but the last argument; the last
	 * argument is the output file, which we update incrementally.
	 */

	for (i = 1; i < argc - 1; i++) {


		if (!strcmp(argv[i], "--scale2")) {
			res_scale = 2;
		} else if (!strcmp(argv[i], "--scale4")) {
			res_scale = 4;
		} else if (!strcmp(argv[i], "--scale8")) {
			res_scale = 8;
		} else if (!strcmp(argv[i], "--align-all")) {
			align_code = 0;
		} else if (!strcmp(argv[i], "--align-green")) {
			align_code = 1;
		} else if (!strcmp(argv[i], "--align-sum")) {
			align_code = 2;
		} else if (!strcmp(argv[i], "--translation")) {
			transform_code = 0;
		} else if (!strcmp(argv[i], "--euclidean")) {
			transform_code = 1;
		} else if (!strcmp(argv[i], "--projective")) {
			transform_code = 2;
		} else if (!strncmp(argv[i], "--metric=", strlen("--metric="))) {
			sscanf(argv[i] + strlen("--metric="), "%lf", &metric);
		} else if (!strncmp(argv[i], "--threshold=", strlen("--threshold="))) {
			sscanf(argv[i] + strlen("--threshold="), "%lf", &match_threshold);
		} else if (!strncmp(argv[i], "--perturb-upper=", strlen("--perturb-upper="))) {
			sscanf(argv[i] + strlen("--perturb-upper="), "%lf", &perturb_upper);
		} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
			fprintf(stderr, "\n\n*** Warning: --stepsize is deprecated.  "
					"Use --perturb-lower instead. ***\n\n\n");
			sscanf(argv[i] + strlen("--stepsize="), "%lf", &perturb_lower);
		} else if (!strncmp(argv[i], "--hf-enhance=", strlen("--hf-enhance="))) {
			sscanf(argv[i] + strlen("--hf-enhance="), "%lf", &hf_enhance);
		} else if (!strncmp(argv[i], "--perturb-lower=", strlen("--perturb-lower="))) {
			sscanf(argv[i] + strlen("--perturb-lower="), "%lf", &perturb_lower);
		} else if (!strncmp(argv[i], "--trans-load=", strlen("--trans-load="))) {
			tload_delete(tload);
			tload = tload_new(argv[i] + strlen("--trans-load="));
		} else if (!strncmp(argv[i], "--trans-save=", strlen("--trans-save="))) {
			tsave_delete(tsave);
			tsave = tsave_new(argv[i] + strlen("--trans-save="));
		} else if (display_image == NULL) {

			/* 
			 * First file argument.  Print general file information as well
			 * as information specific to this argument.
			 */

			tsave_orig(tsave, argv[i]);
			tsave_target(tsave, argv[argc - 1]);

			fprintf(stderr, "Output file will be '%s'.\n", 
					argv[argc - 1]);

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			display_image = read_image(argv[i]);

			if (res_scale != 1)
				scale(display_image, res_scale);
			weight_image = new_image(height(display_image), 
						 width(display_image), 1);

			write_image(argv[argc - 1], display_image);

			fprintf(stderr, ".\n");

		} else {

			tsave_info (tsave, argv[i]);

			fprintf(stderr, "Merging supplemental frame '%s'", 
					argv[i]);

			input_image = read_image(argv[i]);
			if (res_scale != 1)
				scale(input_image, res_scale);
			update();
			image_final(input_image);
			free(input_image);

			write_image(argv[argc - 1], display_image);

			fprintf(stderr, ".\n");
			
		}

	}

	if (display_image) {

		if (res_scale != 1) {

			image *pptarget = clone(display_image);
						    
			fprintf(stderr, "Post-enhancing high frequencies");

			filter(pptarget, display_image, res_scale, hf_enhance);

			write_image(argv[argc - 1], pptarget);

			fprintf(stderr, ".\n");
		}

		tsave_delete(tsave);
		tload_delete(tload);
		fprintf(stderr, "Done (%f%% average match).\n", match_sum / match_count);
		
	} else {

		/*
		 * If there was no output, the user might need more information.
		 */

		fprintf(stderr, 
			"\n"
			"Usage: %s [<options>] <input-files> ... <output-file>\n"
			"   or: %s --version\n"
			"\n\n"
			"Scaling options:\n\n"
			"--scale2          Scale input images up by 2\n"
			"--scale4          Scale input images up by 4\n"
			"--scale8          Scale input images up by 8\n"
			"--hf-enhance=n    Post-enhance high freq. details by factor n. (0.0 is default)\n"
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
			"Transformation file operations:\n\n"
			"--trans-load=x    Load initial transformation settings from file x\n"
			"--trans-save=x    Save final transformation data in file x\n"
			"\n\n"
			"Tunable parameters:\n\n"
			"--metric=x        Set the error metric exponent (2 is default)\n"
			"--threshold=x     Min. match threshold; a perfect match is 100.  (0 is default)\n"
			"--perturb-upper=x Max. correction step, in pixels or degrees (32.0 is default)\n"
			"--perturb-lower=x Min. correction step, in pixels or degrees (0.125 is default)\n"
			"\n",
			argv[0], argv[0]);
		return 1;
	}

	return 0;
}

