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

#define PPM_CHANNEL_READ(is_binary, transfer, f_source, f_dest, f_name, channel_type, channel_format, channel_mcv, extended_ptr) {\
	channel_type ival;\
\
	if (is_binary) {\
\
		unsigned char val;\
		int n;\
\
		/* Binary data */\
\
		n = fscanf(f_source, "%c", &val);\
		assert (n == 1);\
\
		if (n != 1)\
			error_ppm(f_name);\
\
		ival = val;\
\
		for (channel_type mcv_r = channel_mcv / 256; mcv_r; mcv_r /= 256) {\
\
			n = fscanf(f_source, "%c", &val);\
			assert(n == 1);\
\
			if (n != 1)\
				error_ppm(f_name);\
\
			ival = (ival << 8) | val;\
		}\
\
	} else {\
\
		/* ASCII data */\
\
		eat_comments(f_source, f_name, extended_ptr);\
\
		n = fscanf(f, channel_format, &ival);\
\
		assert (n == 1);\
		if (n != 1)\
			error_ppm(f_name);\
	}\
\
	if (transfer)\
		fwrite(&ival, sizeof(channel_type), 1, f_dest);\
\
}

static inline ale_image read_ppm(const char *filename, exposure *e, unsigned int bayer, int init_reference_gain = 0) {
	unsigned int i, j, k;
	ale_image im;
	unsigned char m1, m2, val;
	int m3, m4;
	int w, h;
	cl_ulong mcv;
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
	n = fscanf(f, "%llu", &mcv);
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
			((mcv <= 4294967295ul) ? ALE_TYPE_UINT_32 : ALE_TYPE_UINT_64)));

	if (m2 == '6' && bayer == IMAGE_BAYER_NONE 
	 && (mcv <= 255 || !ppm_short_little_endian_check())) {

		/*
		 * For 3-channel binary 8-bit, and 16-bit on big-endian
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

			if (mcv <= 255)
				PPM_CHANNEL_READ((m2 == '6'), (ale_has_channel(i, j, k, bayer)), f, converted_f, filename, cl_uchar, "%hhu", ((cl_uchar) mcv), (&extended))
			else if (mcv <= 65535)
				PPM_CHANNEL_READ((m2 == '6'), (ale_has_channel(i, j, k, bayer)), f, converted_f, filename, cl_ushort, "%hu", ((cl_ushort) mcv), (&extended))
			else if (mcv <= 4294967295ul)
				PPM_CHANNEL_READ((m2 == '6'), (ale_has_channel(i, j, k, bayer)), f, converted_f, filename, cl_uint, "%u", ((cl_uint) mcv), (&extended))
			else
				PPM_CHANNEL_READ((m2 == '6'), (ale_has_channel(i, j, k, bayer)), f, converted_f, filename, cl_ulong, "%llu", ((cl_ulong) mcv), (&extended))

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

static inline void write_ppm(const char *filename, ale_image im, unsigned long mcv, int plain) {
	unsigned int i, j, k;
	FILE *f = fopen(filename, "wb");

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	assert(f);

	/*
	 * XXX: For simplicity of implementation, we currently only handle up
	 * to 16-bit output, which should be within the limits of the current
	 * ALE user interface.
	 */

	assert(mcv <= 65535);

	if (mcv > 65535) {
		fprintf(stderr, "error: I don't know how to produce greater than 16-bit output.\n");
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

	fprintf(f, "%d ", ale_image_get_width(im));

	/* Height */

	fprintf(f, "%d ", ale_image_get_height(im));

	/* Maximum component value */

	fprintf(f, "%lu\n", mcv);

	/* Pixels */

	FILE *image_data = ale_image_retain_file(im);

	while (!feof(image_data)) {

		cl_ushort val;

		if (mcv > 255) {
			if (!fread(&val, sizeof(cl_ushort), 1, image_data))
				break;
		} else {
			val = fgetc(image_data);
		}

		if (feof(image_data))
			break;

		if (plain) {
			fprintf(f, "%u ", (unsigned int) val);
		} else if (mcv <= 255) {
			fprintf(f, "%c", (unsigned char) val);
		} else {
			fprintf(f, "%c%c", (unsigned char) (val >> 8), (unsigned char) (val & 0xff));
		}
	}

	ale_image_release_file(im, image_data);

	/* Done */

	fclose(f);
}

#endif
