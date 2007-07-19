// Copyright 2002, 2003, 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
 * d3/tfile.h: Read and write 3D transformation data files.
 */

/*
 * This version of ALE reads and writes version 0 3D transformation data files.
 */

#ifndef __d3tfile_h__
#define __d3tfile_h__

#include "pt.h"

#define D3_TFILE_VERSION     0
#define D3_TFILE_VERSION_MAX 0

extern int tfile_input_version;
extern int tfile_output_version;

/*
 * Structure to describe a transformation data file to load data from.
 */

struct tload_t {
	const char *filename;
	FILE *file;
};

/*
 * Structure to describe a transformation data file to write data to.
 */

struct tsave_t {
	const char *filename;
	const char *target;
	const char *orig;
	FILE *file;
};

/*
 * Create a new tload_t transformation data file structure, used for
 * reading data from transformation data files.
 */

static inline struct tload_t *tload_new(const char *filename) {
	FILE *file = fopen (filename, "r");
	struct tload_t *result = NULL;

	if (!file) {
		fprintf(stderr, "tload: Error: could not open transformation data file '%s'.", filename);
		exit(1);
	}

	result = (struct tload_t *) 
		malloc(sizeof(struct tload_t));
	result->filename = filename;
	result->file = file;

	return result;
}

/*
 * Load the first transformation from a transformation data file associated with
 * transformation data file structure T, or return the default transformation
 * if no transformation is available.
 *
 * 	T is a pointer to the tload_t transformation data file structure.
 *
 * 	DEFAULT_TRANSFORM is the default transformation result.
 *
 *	IS_DEFAULT is used to signal a non-default transformation result.
 */

static inline pt tload_first(struct tload_t *t,
		pt default_transform, int *is_default) {

	pt result = default_transform;

	*is_default = 1;

	/*
	 * If there is no file, return the default.
	 */

	if (t == NULL)
		return result;

	/*
	 * Search through the initial part of the file to determine 
	 * its version.
	 */

	/*
	 * Skip comments
	 */

	int first_character;

	first_character = fgetc(t->file);

	while (first_character == ' '
	    || first_character == 0xa
	    || first_character == 0xd
	    || first_character == '\t'
	    || first_character == '#') {
		ungetc(first_character, t->file);
		char line[1024];
		fgets(line, 1024, t->file);
		if (strlen(line) >= 1023) {
			fprintf(stderr, 
				"\n3d-trans-load: Error: line too long in input file\n");
			exit(1);
		}

		first_character = fgetc(t->file);
	}
		
	if (first_character != EOF)
		ungetc(first_character, t->file);

	if (first_character != 'W') {
		fprintf(stderr, "\n3d-trans-load: First command must be a version number.\n");
		exit(1);
	}

	/*
	 * Obtain version from version command string.
	 */

	char line[1024];
	fgets(line, 1024, t->file);
	if (strlen(line) >= 1023) {
		fprintf(stderr, 
			"\n3d-trans-load: Error: line too long in input file\n");
		exit(1);
	}

	int count = sscanf(line, "W %d", &tfile_input_version); 
	
	if (count < 1) {
		fprintf(stderr, "Error in 3d transformation "
				"file version command.\n");
		exit(1);
	} else if (tfile_input_version > D3_TFILE_VERSION_MAX) {
		fprintf(stderr, "Unsupported 3D transformation "
				"file version %d\n", 
				tfile_input_version);
		exit(1);
	}

	/*
	 * Read each line of the file until we find a transformation
	 * or EOF.
	 */

	while (!feof(t->file)) {
		char line[1024];

		fgets(line, 1024, t->file);

		if (feof(t->file))
			return result;

		if (strlen(line) >= 1023) {
			fprintf(stderr, 
				"\ntrans-load: Error: line too long in input file\n");
			exit(1);
		}

		const double rtod_multiplier = 180 / M_PI;

		switch (line[0]) {
			case ' ':
			case 0xa:
			case 0xd:
			case '\t':
			case '#':
				/* Comment or whitespace */
				break;
			case 'D':
			case 'd':
				/* Default transformation */
				return result;
			case 'V':
			case 'v':
				unsigned int count;
				double view_angle;

				count = sscanf(line, "V %lf", &view_angle);

				if (count < 1) {
					fprintf(stderr, "\n3d-trans-load: Error: "
							"Malformed 'V' command.\n");
					exit(1);
				}

				result.view_angle(view_angle / rtod_multiplier);

				break;
			case 'E':
			case 'e':
				/* Euclidean transformation data */
				*is_default = 0;
				{
					double width, height;
					double values[6] = {0, 0, -1, 0, 0, 0};
					int count;

					count = sscanf(line, "E %lf%lf%lf%lf%lf%lf%lf%lf",
							&width, &height,
							&values[1], &values[0], &values[2],
							&values[4], &values[3], &values[5]);

					if (count < 8) 
						fprintf(stderr, "\n3d-trans-load: warning: "
								"Missing args for 'E'\n");

					if (width != result.scaled_width()
					 || height != result.scaled_height()) {
						fprintf(stderr, "\n3d-trans-load: Error: "
								"Scaled image width and/or height mismatch.");
					}

					for (int i = 3; i < 6; i++) {
						values [i] /= rtod_multiplier;
					}

					result.e().set(values);

					return result;
				}
				break;
			default:
				fprintf(stderr,
					"\ntrans-load: Error in tload_first: unrecognized command '%s'\n",
					line);
				exit(1);
		}
	}

	/*
	 * EOF reached: return default transformation.
	 */

	return result;
}

/*
 * Load the next transformation from a transformation data file associated with
 * transformation data file structure T, or return the default transformation
 * if no transformation is available.
 *
 * 	T is a pointer to the tload_t transformation data file structure.
 *
 * 	IS_P is nonzero if a projective transformation is expected.
 *
 * 	DEFAULT_TRANSFORM is the default transformation result.
 *
 *	IS_DEFAULT is used to signal a non-default transformation result.
 */

static inline pt tload_next(struct tload_t *t, 
		pt default_transform, int *is_default) {

	pt result = default_transform;

	*is_default = 1;

	/*
	 * Read each line of the file until we find a transformation.
	 */

	while (t && !feof(t->file)) {
		char line[1024];

		fgets(line, 1024, t->file);

		if (feof(t->file))
			return result;

		if (strlen(line) >= 1023) {
			fprintf(stderr, 
				"\ntrans-load: warning: line too long in input file\n");
		}

		const double rtod_multiplier = 180 / M_PI;

		switch (line[0]) {
			case ' ':
			case 0xa:
			case 0xd:
			case '\t':
			case '#':
				/* Comment or whitespace */
				break;
			case 'D':
			case 'd':
				/* Default transformation */
				return result;
			case 'V':
			case 'v':
				unsigned int count;
				double view_angle;

				count = sscanf(line, "V %lf", &view_angle);

				if (count < 1) {
					fprintf(stderr, "\n3d-trans-load: Error: "
							"Malformed 'V' command.\n");
					exit(1);
				}

				result.view_angle(view_angle / rtod_multiplier);

				break;
			case 'E':
			case 'e':
				/* Euclidean transformation data */
				*is_default = 0;
				{
					double width, height;
					double values[6] = {0, 0, -1, 0, 0, 0};
					int count;

					count = sscanf(line, "E %lf%lf%lf%lf%lf%lf%lf%lf",
							&width, &height,
							&values[1], &values[0], &values[2],
							&values[4], &values[3], &values[5]);

					if (count < 8) 
						fprintf(stderr, "\n3d-trans-load: warning: "
								"Missing args for 'E'\n");

					if (width != result.scaled_width()
					 || height != result.scaled_height()) {
						fprintf(stderr, "\n3d-trans-load: Error: "
								"Scaled image width and/or height mismatch.");
					}

					for (int i = 3; i < 6; i++) {
						values [i] /= rtod_multiplier;
					}

					result.e().set(values);

					return result;
				}
				break;
			default:
				fprintf(stderr,
					"\ntrans-load: Error in tload_next: unrecognized command '%s'\n",
					line);
				exit(1);
		}
	}

	return result;
}

/*
 * Create a new tsave_t transformation data file structure, used for
 * writing data to transformation data files.
 */

static inline struct tsave_t *tsave_new(const char *filename) {
	FILE *file = fopen (filename, "w");
	struct tsave_t *result = NULL;

	if (!file) {
		fprintf(stderr, "tsave: Error: could not open transformation data file '%s'.", filename);
		exit(1);
	}

	result = (struct tsave_t *) 
		malloc(sizeof(struct tsave_t));
	result->filename = filename;
	result->file = file;
	result->orig = "unknown";
	result->target = "unknown";

	fprintf(file, "# created by ALE 3D transformation file handler version %d\n", 
			D3_TFILE_VERSION);

	fclose(file);

	return result;
}

/*
 * Save the first transformation to a transformation data file associated with
 * transformation data file structure T, or do nothing if T is NULL.  This
 * function also establishes the output file version.
 *
 * 	OFFSET is the transformation to be saved.
 *
 * 	IS_PROJECTIVE indicates whether to write a projective transformation.
 *
 */

static inline void tsave_first(struct tsave_t *t, pt offset) {

	if (t == NULL)
		return;

	t->file = fopen(t->filename, "a");
	
	tfile_output_version = 0;

	fprintf(t->file, "# producing 3D transformation file syntax version %d\n", tfile_output_version);
	fprintf(t->file, "W %d\n", tfile_output_version);

	// fprintf(t->file, "# Comment: Target output file is %s\n", t->target);
	// fprintf(t->file, "# Comment: Original frame is %s\n", t->orig);
	
	const double rtod_multiplier = 180 / M_PI;

	fprintf(t->file, "V %lf\n", offset.view_angle() * rtod_multiplier);

	fprintf(t->file, "E ");
	fprintf(t->file, "%f %f ", (double) offset.scaled_width(), (double) offset.scaled_height());
	fprintf(t->file, "%f ",    (double) offset.e().get_translation(1));
	fprintf(t->file, "%f ",    (double) offset.e().get_translation(0));
	fprintf(t->file, "%f ",    (double) offset.e().get_translation(2));
	fprintf(t->file, "%f ",    (double) offset.e().get_rotation(1) * rtod_multiplier);
	fprintf(t->file, "%f ",    (double) offset.e().get_rotation(0) * rtod_multiplier);
	fprintf(t->file, "%f ",    (double) offset.e().get_rotation(2) * rtod_multiplier);

	fprintf(t->file, "\n");

	fclose(t->file);
}

/*
 * Save the next transformation to a transformation data file associated with
 * transformation data file structure T, or do nothing if T is NULL.
 *
 * 	OFFSET is the transformation to be saved.
 *
 * 	IS_PROJECTIVE indicates whether to write a projective transformation.
 *
 */

static inline void tsave_next(struct tsave_t *t, pt offset) {

	if (t == NULL)
		return;

	t->file = fopen(t->filename, "a");
	
	const double rtod_multiplier = 180 / M_PI;

	fprintf(t->file, "V %lf\n", offset.view_angle() * rtod_multiplier);

	fprintf(t->file, "E ");
	fprintf(t->file, "%f %f ", (double) offset.scaled_width(), (double) offset.scaled_height());
	fprintf(t->file, "%f ",    (double) offset.e().get_translation(1));
	fprintf(t->file, "%f ",    (double) offset.e().get_translation(0));
	fprintf(t->file, "%f ",    (double) offset.e().get_translation(2));
	fprintf(t->file, "%f ",    (double) offset.e().get_rotation(1) * rtod_multiplier);
	fprintf(t->file, "%f ",    (double) offset.e().get_rotation(0) * rtod_multiplier);
	fprintf(t->file, "%f ",    (double) offset.e().get_rotation(2) * rtod_multiplier);

	fprintf(t->file, "\n");

	fclose(t->file);
}

/*
 * Write information to a transformation file indicating the target output
 * file.
 */

static inline void tsave_target(struct tsave_t *t, const char *filename) {
	if (t == NULL)
		return;

	t->target = filename;
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		
		fclose(t->file);
	}
}

/*
 * Write information to a transformation data file indicating the filename
 * of the original frame (i.e. the first frame in the sequence of input
 * frames).
 */

static inline void tsave_orig(struct tsave_t *t, const char *filename) {
	if (t == NULL)
		return;

	t->orig = filename;
}

/*
 * Write information to a transformation data file indicating the filename
 * of a supplemental frame (i.e. a frame in the sequence of input frames 
 * that is not the first frame).
 */

static inline void tsave_info(struct tsave_t *t, const char *filename) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		fprintf(t->file, "# Comment: Supplemental frame %s\n", filename);
		
		fclose(t->file);
	}
}

/*
 * Write information to a transformation data file indicating the tonal
 * registration multiplier.
 */

static inline void tsave_trm(struct tsave_t *t, ale_real r, ale_real g, ale_real b) {
       if (t != NULL) {
               t->file = fopen(t->filename, "a");

               fprintf(t->file, "# Comment: Exposure [r=%f g=%f b=%f]\n", r, g, b);

               fclose(t->file);
       }
}

/*
 * Write information to a transformation data file indicating the average
 * pixel magnitude.
 */

static inline void tsave_apm(struct tsave_t *t, ale_real r, ale_real g, ale_real b) {
       if (t != NULL) {
               t->file = fopen(t->filename, "a");

               fprintf(t->file, "# Comment: Avg magnitude [r=%f g=%f b=%f]\n", r, g, b);

               fclose(t->file);
       }
}


/*
 * Destroy a tload_t transformation data file structure.
 */

static inline void tload_delete(struct tload_t *victim) {
	if (victim)
		fclose(victim->file);
	free(victim);
}

/*
 * Destroy a tsave_t transformation data file structure.
 */

static inline void tsave_delete(struct tsave_t *victim) {
	free(victim);
}

#endif
