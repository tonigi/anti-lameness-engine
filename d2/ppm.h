// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>, 
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
 * ppm.h: Read and write PPM files.
 */

#ifndef __ppm_h__
#define __ppm_h__

#include "image_ale_real.h"
#include "image_bayer_ale_real.h"
#include "exposure/exposure.h"

/*
 * Extended attributes
 */

struct extended_t {
	int is_extended;
	int black_level;
	ale_real aperture; /* 1 == f/1.0, 1.4 == f/1.4, etc. */
	ale_real shutter;  /* 1 == 1 sec, 0.5 == 1/2 sec, etc. */
	ale_real gain;     /* 1 == ISO 100, 2 == ISO 200, etc. */

	extended_t() {
		is_extended = 0;
		black_level = 0;
		aperture = 1;
		shutter = 1;
		gain = 1;
	}
};

static inline void error_ppm(const char *filename) {
	fprintf(stderr, 
		"\n\n*** '%s' doesn't look like a PPM file.\n"
		"\n*** To handle other file types, compile ALE with\n"
		"*** ImageMagick support ('make IMAGEMAGICK=1').\n"
		"*** (To do this, you must have a source distribution\n"
		"*** of ALE, and you must have ImageMagick installed.)\n\n", 
		filename);
	exit(1);
}

static inline int digest_comment(FILE *f, const char *filename, extended_t *extended) {
	int next = '#';
	int value;
	double fvalue, fvalue2;

	while (next != '\n' && next != '\r' && next != EOF) {
		while (next == ' ' || next == '\t' || next == '#') {
			next = fgetc(f);
			if (feof(f))
				error_ppm(filename);
		}

		if (ungetc(next, f) == EOF) {
			assert(0);
			fprintf(stderr, "Unable to ungetc().");
			exit(1);
		}

		fvalue2 = 1;

		if (extended->is_extended && fscanf(f, "Black-level: %d", &value) == 1)
			extended->black_level = value;
		else if (extended->is_extended && fscanf(f, "ISO: %lf", &fvalue) == 1)
			extended->gain = fvalue / 100;
		else if (extended->is_extended && fscanf(f, "Gain: %lf", &fvalue) == 1)
			extended->gain = fvalue;
		else if (extended->is_extended && fscanf(f, "Aperture: %lf", &fvalue) == 1)
			extended->aperture = fvalue;
		else if (extended->is_extended && fscanf(f, "Shutter: %lf/%lf", &fvalue, &fvalue2) > 0)
			extended->shutter = fvalue / fvalue2;
		else if (next != '\n' && next != '\r' && next != EOF)
			next = fgetc(f);

		next = fgetc(f);
	}

	return next;
}

static inline void eat_comments(FILE *f, const char *filename, extended_t *extended) {
	int next = ' ';

	while (next == ' ' || next == '\n' || next == '\t' || next == '#' || next == '\r') {
		next = fgetc(f);
		if (next == '#')
			next = digest_comment(f, filename, extended);
		if (feof(f))
			error_ppm(filename);
	}

	if (ungetc(next, f) == EOF) {
		assert(0);
		fprintf(stderr, "Unable to ungetc().");
		exit(1);
	}
}

static inline int is_eppm(const char *filename) {
	char m1, m2, m3, m4;
	int n;
	extended_t extended;
	FILE *f = fopen(filename, "rb");

	if (f == NULL)
		return 0;

	/* Magic */

	eat_comments(f, filename, &extended);	/* XXX - should we eat comments here? */
	n = fscanf(f, "%c%c%c%c", &m1, &m2, &m3, &m4);

	fclose(f);

	if (n != 4 || m1 != 'P' || (m2 != '6' && m2 != '3') || m3 != '#' || m4 != 'E') 
		return 0;

	return 1;
}

static inline image *read_ppm(const char *filename, exposure *e, unsigned int bayer, int init_reference_gain = 0) {
	unsigned int i, j, k;
	image *im;
	unsigned char m1, m2, val;
	int m3, m4;
	int ival;
	int w, h, mcv;
	int n;
	struct extended_t extended;
	FILE *f = fopen(filename, "rb");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/* Magic */

	eat_comments(f, filename, &extended);	/* XXX - should we eat comments here? */
	n = fscanf(f, "%c%c", &m1, &m2);
	assert(n == 2 && m1 == 'P' && (m2 == '6' || m2 == '3'));

	if (n != 2 || m1 != 'P' || (m2 != '6' && m2 != '3'))
		error_ppm(filename);

	/* Extended flag */

	m3 = fgetc(f);

	if (m3 == '#') {
		m4 = fgetc(f);
		if (m4 == 'E') 
			extended.is_extended = 1;
		else while (m4 != EOF && m4 != '\n' && m4 != '\r')
			m4 = fgetc(f);
	} else if (ungetc(m3, f) == EOF) {
		assert(0);
		fprintf(stderr, "Unable to ungetc().");
		exit(1);
	}

	/* Width */

	eat_comments(f, filename, &extended);
	n = fscanf(f, " %d", &w);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Height */

	eat_comments(f, filename, &extended);
	n = fscanf(f, "%d", &h);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Maximum component value */

	eat_comments(f, filename, &extended);
	n = fscanf(f, "%d", &mcv);
	assert(n == 1);
	assert(mcv <= 65535 || m2 == '3');

	if (n != 1 || (mcv > 65535 && m2 == '6'))
		error_ppm(filename);

	/* Make a new image */

	if (bayer == IMAGE_BAYER_NONE)
		im = new image_ale_real(h, w, 3, "file", e);
	else 
		im = new image_bayer_ale_real(h, w, 3, bayer, "file", e);

	assert (im);

	/* Trailing whitespace */

	if (fgetc(f) == EOF) {
		assert(0);
		error_ppm(filename);
	}

	/* Pixels */

	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++) {
		pixel p;
		for (k = 0; k < im->depth();  k++) {

			if (m2 == '6') {

				/* Binary data */

				n = fscanf(f, "%c", &val);
				assert (n == 1);

				if (n != 1)
					error_ppm(filename);

				ival = val;

				if (mcv > 255) {
					n = fscanf(f, "%c", &val);
					assert(n == 1);

					if (n != 1)
						error_ppm(filename);

					ival = (ival << 8) | val;
				}

			} else {

				/* ASCII data */

				eat_comments(f, filename, &extended);

				n = fscanf(f, "%d", &ival);

				assert (n == 1);
				if (n != 1)
					error_ppm(filename);
			}

			p[k] = (ale_real) (ival - (extended.is_extended ? extended.black_level : 0)) 
			     / (ale_real) (mcv  - (extended.is_extended ? extended.black_level : 0));
		}

		im->set_pixel(i, j, e->linearize(p));
	}

	/* Handle exposure and gain */

	if (extended.is_extended) {
		double combined_gain = (1 / pow(extended.aperture, 2))
			             * extended.shutter
			             * extended.gain;
		if (init_reference_gain)
			exposure::set_gain_reference(combined_gain);
		else
			e->set_gain_multiplier(exposure::get_gain_reference()
					     / combined_gain);
	}

	/* Done */

	fclose(f);

	return im;
}

static inline void write_ppm(const char *filename, const image *im, exposure *e, 
		unsigned int mcv, int plain, int rezero, int exposure_scale, double nn_defined_radius) {
	unsigned int i, j, k;
	FILE *f = fopen(filename, "wb");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/*
	 * Output a plain (ASCII) or raw (binary) PPM file
	 */

	/* Magic */

	if (plain)
		fprintf(f, "P3 ");
	else
		fprintf(f, "P6 ");

	/* Width */

	fprintf(f, "%d ", im->width());

	/* Height */

	fprintf(f, "%d ", im->height());

	/* Maximum component value */

	fprintf(f, "%d\n", mcv);

	/* Automatic exposure adjustment information */

	ale_real maxval = 1;
	ale_real minval = (rezero ? im->minval() : 0);
	if (minval > 0)
		minval = 0;
	pixel minval_pixel(minval, minval, minval);

	if (exposure_scale) {
		ale_real new_maxval = im->maxval();

		if (new_maxval > maxval)
			maxval = new_maxval;
	}

	/* Pixels */

	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++) {
		pixel value = im->get_pixel(i, j);

		/*
		 * Get nearest-neighbor defined values.
		 */

		for (k = 0; k < 3; k++)
		if (isnan(value[k]))
		for (int radius = 1; radius <= nn_defined_radius; radius++) {
			double nearest_radius_squared = (radius + 1) * (radius + 1);
			for (int ii = -radius; ii <= radius; ii++)
			for (int jj = -radius; jj <= radius; jj++) {
				if (!im->in_bounds(point(i + ii, j + jj)))
					continue;
				if (ii * ii + jj * jj < nearest_radius_squared
				 && finite(im->get_pixel(i + ii, j + jj)[k])) {
					value[k] = im->get_pixel(i + ii, j + jj)[k];
					nearest_radius_squared = ii * ii + jj * jj;
				}
			}
			if (nearest_radius_squared < (radius + 1) * (radius + 1))
				break;
		}

		pixel exposure_adjust = (value - minval_pixel)
			              / (maxval - minval);
		pixel unlinearized = (e->unlinearize(exposure_adjust)).clamp();

		pixel rescaled = unlinearized * (ale_real) mcv;

		for (k = 0; k < im->depth();  k++) {

			uint16_t output_value = (uint16_t) round(rescaled[k]);

			if (plain) {
				fprintf(f, "%d ", output_value);
			} else {
				if (mcv > 255)
					fprintf(f, "%c", output_value >> 8);
				fprintf(f, "%c", 0xff & output_value);
			}
		}

		if (plain)
			fprintf(f, "\n");

	}

	/* Done */

	fclose(f);
}

#endif
