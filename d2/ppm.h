// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                              <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
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
	ale_real black_level;
	ale_real aperture; /* 1 == f/1.0, 1.4 == f/1.4, etc. */
	ale_real shutter;  /* 1 == 1 sec, 0.5 == 1/2 sec, etc. */
	ale_real gain;     /* 1 == ISO 100, 2 == ISO 200, etc. */

	extended_t() {
		is_extended = 0;
		black_level = 0;
		aperture = 0;
		shutter = 0;
		gain = 0;
	}
};

static inline void error_ppm(const char *filename) {
	fprintf(stderr, 
		"\n\n*** '%s' doesn't look like a PPM file.\n"
		"\n*** (To handle other file types, use a version of ALE with\n"
		"*** ImageMagick support enabled.)\n\n",
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

static int ppm_void_file_close(void *f) {
	return fclose((FILE *) f);
}

static int ppm_short_little_endian_check() {
	/*
	 * Modified from http://unixpapa.com/incnote/byteorder.html
	 */
	
	short one = 1;
	return (*((char *)(&one)));
}

static inline ale_image read_ppm(const char *filename, exposure *e, unsigned int bayer, int init_reference_gain = 0) {
	unsigned int i, j, k;
	ale_image im;
	unsigned char m1, m2, val;
	int m3, m4;
	long ival;
	int w, h;
	unsigned long mcv;  // XXX: could be 32 bit, which could break things.
	int n;
	struct extended_t extended;
	FILE *f = fopen(filename, "rb");

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	assert(f);

	/* Magic */

	eat_comments(f, filename, &extended);	/* XXX - should we eat comments here? */
	n = fscanf(f, "%c%c", &m1, &m2);

	if (n != 2 || m1 != 'P' || (m2 != '6' && m2 != '3'))
		error_ppm(filename);

	assert(n == 2 && m1 == 'P' && (m2 == '6' || m2 == '3'));

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
	n = fscanf(f, "%lu", &mcv);
	assert(n == 1);
	assert(mcv <= 65535 || m2 == '3');

	if (extended.black_level == 0) {
		extended.black_level = e->get_black_level();
	} else {
		extended.black_level /= mcv;
	}

	if (n != 1 || (mcv > 65535 && m2 == '6'))
		error_ppm(filename);

	/* Trailing whitespace */

	if (fgetc(f) == EOF) {
		assert(0);
		error_ppm(filename);
	}

	/* Make a new image */

	im = ale_new_image(accel::context(), 
			(bayer == IMAGE_BAYER_NONE) ? ALE_IMAGE_RGB : ALE_IMAGE_Y,
			(mcv <= 255) ? ALE_TYPE_UINT_8 :
			((mcv <= 65535) ? ALE_TYPE_UINT_16 :
			((mcv <= 4294967295ul) ? ALE_TYPE_UINT_32 : ALE_TYPE_UINT_64)));		// (XXX: mcv type may be too short for this to be meaningful.)

	if (m2 == '6' && bayer == IMAGE_BAYER_NONE 
	 && (mcv <= 255 || ppm_short_little_endian_check())) {

		/*
		 * For 3-channel binary 8-bit, and 16-bit on little-endian
		 * systems, use a file suffix.
		 */

		ale_image_set_file_static(im, w, h, f, ftell(f), ppm_void_file_close, f);

	} else {

		/*
		 * For all others, convert data via a new, temporary file.
		 */

		FILE *converted_f = tmpfile();

		/* Pixels */

		for (i = 0; i < h; i++)
		for (j = 0; j < w;  j++)
		for (k = 0; k < 3;  k++) {

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

				n = fscanf(f, "%ld", &ival);

				assert (n == 1);
				if (n != 1)
					error_ppm(filename);
			}

			if (!ale_has_channel(i, j, k, bayer))
				continue;

			fprintf(converted_f, "%c", ((char *) ival)[0]);
			if (mcv > 255)
				fprintf(converted_f, "%c", ((char *) ival)[1]);
			if (mcv > 65535)
				fprintf(converted_f, "%c%c", ((char *) ival)[2], ((char *) ival)[3]);
			if (mcv > 4294967295ul)
				fprintf(converted_f, "%c%c%c%c", ((char *) ival)[4], ((char *) ival)[5], ((char *) ival)[6], ((char *) ival)[7]);   // XXX: ival might be too short for this
		}

		 ale_image_set_file_static(im, w, h, converted_f, 0, ppm_void_file_close, converted_f);

		 fclose(f);
	}

	/* Handle exposure and gain */

	if (extended.is_extended) {
		if (extended.aperture != 0
		 || extended.shutter != 0
		 || extended.gain != 0) {

			if (extended.aperture == 0)
				extended.aperture = 1;
			if (extended.shutter == 0)
				extended.shutter = 1;
			if (extended.gain == 0)
				extended.gain = 1;

			ale_real combined_gain = (1 / pow(extended.aperture, 2))
					     * extended.shutter
					     * extended.gain;
			
			if (init_reference_gain)
				exposure::set_gain_reference(combined_gain);
			else
				e->set_gain_multiplier(exposure::get_gain_reference()
						     / combined_gain);
		}
	}

	return im;
}

static inline void write_ppm(const char *filename, ale_image im, unsigned int mcv, int plain) {
	unsigned int i, j, k;
	FILE *f = fopen(filename, "wb");

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	assert(f);

	/*
	 * Output a plain (ASCII) or raw (binary) PPM file
	 */

	/* Magic */

	if (plain)
		fprintf(f, "P3 ");
	else
		fprintf(f, "P6 ");

	/* Width */

	fprintf(f, "%d ", ale_image_get_width(im));

	/* Height */

	fprintf(f, "%d ", ale_image_get_height(im));

	/* Maximum component value */

	fprintf(f, "%d\n", mcv);

	/* Pixels */

	FILE *image_data = ale_image_retain_file(im);

	if (plain) {
		while (!feof(image_data)) {
			char c[2];
			c[0] = fgetc(image_data);
			c[1] = fgetc(image_data);

			fprintf(f, "%u", *((unsigned short *) c));
		}
	} else {
		while (!feof(image_data)) {
			char c = fgetc(image_data);
			if (c == EOF)
				break;
		}
	}

	ale_image_release_file(image_data);

	/* Done */

	fclose(f);
}

#endif
