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
#include "ppm.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

char *version = "0.1.0";

image *display_image, *input_image, *weight_image;

/*
 * Global flags
 */

static int use_rotation = 1;
static int align_code = 2;
static double metric = 2;
static double match_threshold = 0;
static double minimum_stepsize = 0.1;

/*
 * Not-quite-symmetric difference function.  Determines the difference in areas
 * where the arrays overlap.  Uses the first array's notion of pixel positions.
 */
inline double diff(image *a, image *b, double xoff, double yoff, double theta) {
	float result = 0;
	float divisor = 0;
	float apm_a = (float) avg_pixel_magnitude(a);
	float apm_b = (float) avg_pixel_magnitude(b);

	int i, j, k;

	for (i = 0; i < height(a); i++)
		for (j = 0; j < width(a); j++) {

			/*
			 * Translate
			 */

			float ti = i - yoff;
			float tj = j - xoff;

			/*
			 * Rotate
			 */

			float rti = (ti - height(b)/2) * cos(theta)
				   + (tj - width(b)/2)  * sin(theta)
				   + height(b)/2;

			float rtj = (tj - width(b)/2)  * cos(theta)
				   - (ti - height(b)/2) * sin(theta)
				   + width(b)/2;

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of array b.
			 */

			if (rti >= 0
			 && rti <= height(b) - 1
			 && rtj >= 0
			 && rtj <= width(b) - 1){ 

				if (align_code == 0) {
					/*
					 * Align based on all channels.
					 */

					for (k = 0; k < 3; k++) {
						float achan = get_pixel_component(a, i, j, k) / apm_a;
						float bchan = get_bl_component(b, rti, rtj, k) / apm_b;

						result += pow(fabs(achan - bchan), metric);
						divisor += pow(achan > bchan ? achan : bchan, metric);
					}
				} else if (align_code == 1) {
					/*
					 * Align based on the green channel.  XXX: apm_* shouldn't be used here.
					 */

					float achan = get_pixel_component(a, i, j, 1) / apm_a;
					float bchan = get_bl_component(b, rti, rtj, 1) / apm_b;

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
						bsum += get_bl_component(b, rti, rtj, k) / apm_b;
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
merge(image *target, image *delta, double xoff, double yoff, double theta) {
	int i, j, k;

	for (i = 0; i < height(target); i++)
		for (j = 0; j < width(target); j++) {

			/*
			 * Translate
			 */

			double ti = i - yoff;
			double tj = j - xoff;

			/*
			 * Rotate
			 */

			double rti = (ti - height(delta)/2) * cos(theta)
				   + (tj - width(delta)/2)  * sin(theta)
				   + height(delta)/2;

			double rtj = (tj - width(delta)/2)  * cos(theta)
				   - (ti - height(delta)/2) * sin(theta)
				   + width(delta)/2;

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of the delta frame.
			 */

			if (rti >= 0
			 && rti <= height(delta) - 1
			 && rtj >= 0
			 && rtj <= width(delta) - 1){ 

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
					     +   get_bl_component(delta, rti, rtj, k)) / (weight + 1));
			}
		}
}


/*
 * Update an accumulated image with a new input image.
 */
inline void update() {
	double perturb = 16;
	double xoff = 0, yoff = 0, theta = 0;
	double here = diff(display_image, input_image, 0, 0, 0);
	int last = -1;

	while (perturb > minimum_stepsize) {

		/*
		 * Assumes that (perturb < 360) is an invariant.  XXX: this may be
		 * unnecessary.
		 */

		double perturb_radians = (2 * M_PI / 360) * perturb;
		double theta_p_perturb = (theta + perturb_radians >= 2 * M_PI) 
                                       ? (theta + perturb_radians -  2 * M_PI) 
                                       : (theta + perturb_radians); 
		double theta_m_perturb = (theta - perturb_radians < 0) 
                                       ? (theta - perturb_radians + 2 * M_PI) 
                                       : (theta - perturb_radians); 

		double test = 0;

		if (last != 1 && (test = diff(display_image, input_image, xoff, yoff - perturb, theta)) < here) {
			yoff -= perturb;
			here = test;
			last = 0;
		} else if (last != 0 && (test = diff(display_image, input_image, xoff, yoff + perturb, theta)) < here) {
			yoff += perturb;
			here = test;
			last = 1;
		} else if (last != 3 && (test = diff(display_image, input_image, xoff - perturb, yoff, theta)) < here) {
			xoff -= perturb;
			here = test;
			last = 2;
		} else if (last != 2 && (test = diff(display_image, input_image, xoff + perturb, yoff, theta)) < here) {
			xoff += perturb;
			here = test;
			last = 3;
		} else if (last != 5 && (test = diff(display_image, input_image, xoff, yoff, theta_p_perturb)) < here) {
			theta = theta_p_perturb;
			here = test;
			last = 4;
		} else if (last != 4 && (test = diff(display_image, input_image, xoff, yoff, theta_m_perturb)) < here) {
			theta = theta_m_perturb;
			here = test;
			last = 5;
		} else {
			perturb *= 0.5;
			last = -1;

			/*
			 * Announce that we've dropped a perturbation level.
			 */

			fprintf(stderr, ".");

		}

	}

	/*
	 * Ensure that the match meets the threshold.
	 */

	if ((1 - here) * 100 > match_threshold) {
		merge (display_image, input_image, xoff, yoff, theta);
		fprintf(stderr, " okay (%f%% match)", (1 - here) * 100);
	} else {
		fprintf(stderr, " no match (%f%% match)", (1 - here) * 100);
	}
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
		} else if (!strcmp(argv[i], "--rotation")) {
			use_rotation = 1;
		} else if (!strcmp(argv[i], "--no-rotation")) {
			use_rotation = 0;
		} else if (!strncmp(argv[i], "--metric=", strlen("--metric="))) {
			sscanf(argv[i] + strlen("--metric="), "%lf", &metric);
		} else if (!strncmp(argv[i], "--threshold=", strlen("--threshold="))) {
			sscanf(argv[i] + strlen("--threshold="), "%lf", &match_threshold);
		} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
			sscanf(argv[i] + strlen("--stepsize="), "%lf", &minimum_stepsize);
		} else if (display_image == NULL) {

			/* 
			 * First file argument.  Print general file information as well
			 * as information specific to this argument.
			 */

			fprintf(stderr, "Output file will be '%s'.\n", 
					argv[argc - 1]);

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			display_image = read_ppm(argv[i]);

			if (res_scale != 1)
				scale(display_image, res_scale);
			weight_image = new_image(height(display_image), 
						 width(display_image), 1);

			write_ppm(argv[argc - 1], display_image);

			fprintf(stderr, ".\n");

		} else {

			fprintf(stderr, "Merging supplemental frame '%s'", 
					argv[i]);

			input_image = read_ppm(argv[i]);
			if (res_scale != 1)
				scale(input_image, res_scale);
			update();

			write_ppm(argv[argc - 1], display_image);

			fprintf(stderr, ".\n");
			
		}

	}

	/*
	 * If there was no output, the user might need more information.
	 */

	if (display_image == NULL) {
		fprintf(stderr, 
			"Usage: %s [<options>] <input-files> ... <output-file>\n"
			"   or: %s --version\n"
			"\n"
			"Options:\n"
			"--scale2	Scale input images up by 2\n"
			"--scale4	Scale input images up by 4\n"
			"--scale8	Scale input images up by 8\n"
			"--align-all	Align images using all color channels\n"
			"--align-green	Align images using the green channel\n"
			"--align-sum	Align images using a sum of channels [default]\n"
			"--rotation	Make rotational adjustments to images [default]\n"
			"--no-rotation	Don't make rotational adjustments to images\n"
			"--metric=x	Set the error metric exponent (2 is default)\n"
			"--threshold=x  Min. match threshold; a perfect match is 100.  (0 is default)\n"
			"--stepsize=x	Set the min. correction step, in pixels or degrees (0.1 default)\n",
			argv[0], argv[0]);
		return 1;
	} else {
		fprintf(stderr, "Done.\n");
	}

	return 0;
}

