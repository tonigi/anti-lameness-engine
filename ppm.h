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

/*
 * ppm.h: Read and write PPM files.
 */

#ifndef __ppm_h__
#define __ppm_h__

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "channel.h"

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

static inline void eat_comments(FILE *f, const char *filename) {
#if 0

	/*
	 * %[] is apparently not compatible with WINE 20021007.
	 */

	char comment[2];

	while (fscanf(f, " %1[#]", comment)) {
		int a = 0;
		while (a != '\n') {
			if (feof(f))
				error_ppm(filename);
			a = fgetc(f);
		}
	}
#else

	/*
	 * ... so we try a different approach.
	 */

	int next = ' ';

	while (next == ' ' || next == '\n' || next == '\t' || next == '#' || next == '\r') {
		next = fgetc(f);
		if (next == '#')
			while (next != '\n' && next != '\r' && next != EOF)
				next = fgetc(f);
		if (feof(f))
			error_ppm(filename);
	}

	if (ungetc(next, f) == EOF) {
		assert(0);
		fprintf(stderr, "Unable to ungetc().");
		exit(1);
	}
#endif
}

static inline image *read_ppm(const char *filename) {
	unsigned int i, j, k;
	image *im;
	unsigned char m1, m2, val;
	unsigned int ival;
	int w, h, mcv;
	int n;
	FILE *f = fopen(filename, "rb");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/* Magic */

	eat_comments(f, filename);	/* XXX - should we eat comments here? */
	n = fscanf(f, "%c%c", &m1, &m2);
	assert(n == 2 && m1 == 'P' && (m2 == '6' || m2 == '3'));

	if (n != 2 || m1 != 'P' || (m2 != '6' && m2 != '3'))
		error_ppm(filename);

	/* Width */

	eat_comments(f, filename);
	n = fscanf(f, " %d", &w);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Height */

	eat_comments(f, filename);
	n = fscanf(f, "%d", &h);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Maximum component value */

	eat_comments(f, filename);
	n = fscanf(f, "%d", &mcv);
	assert(n == 1);
	assert(mcv <= 65535 || m2 == '3');

	if (n != 1 || (mcv > 65535 && m2 == '6'))
		error_ppm(filename);

	/* Make a new image */

	im = new image(h, w, 3);
	assert (im);

	/* Trailing whitespace */

	if (fgetc(f) == EOF) {
		assert(0);
		error_ppm(filename);
	}

	/* Pixels */

	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++)
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

			eat_comments(f, filename);

			n = fscanf(f, "%d", &ival);

			assert (n == 1);
			if (n != 1)
				error_ppm(filename);
		}

		if (mcv != CHANNEL_MAX)
			ival = (int) round(((double) ival / (double) (mcv)) * (CHANNEL_MAX));

		im->set_pixel_component(i, j, k, ival);
	}

	/* Done */

	fclose(f);

	return im;
}

static inline void write_ppm(const char *filename, const image *im) {
	unsigned int i, j, k;
	FILE *f = fopen(filename, "wb");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

#ifdef PPM_PLAIN

	/*
	 * Output an ASCII PPM of 16 bit depth
	 */

	/* Magic */

	fprintf(f, "P3 ");

	/* Width */

	fprintf(f, "%d ", im->width());

	/* Height */

	fprintf(f, "%d ", im->height());

	/* Maximum component value */

	fprintf(f, "%d\n", CHANNEL_MAX);

	/* Pixels */

	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++) {

		for (k = 0; k < im->depth();  k++)
			fprintf(f, "%d ", im->get_pixel_component(i, j, k));

		fprintf(f, "\n");

	}

#else

	/*
	 * Output a binary PPM
	 */

	/* Magic */

	fprintf(f, "P6 ");

	/* Width */

	fprintf(f, "%d ", im->width());

	/* Height */

	fprintf(f, "%d ", im->height());

	/* Maximum component value */

	fprintf(f, "%d\n", CHANNEL_MAX);

	/* Pixels */

	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++)
	for (k = 0; k < im->depth();  k++) {
		if (CHANNEL_MAX > 255) 
			fprintf(f, "%c", im->get_pixel_component(i, j, k) >> 8);
		fprintf(f, "%c", 0xff & im->get_pixel_component(i, j, k));
	}

#endif

	/* Done */

	fclose(f);
}

#endif
