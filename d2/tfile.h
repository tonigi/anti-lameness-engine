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
 * tfile.h: Read and write transformation data files.
 */

/*
 * This version of ALE reads transformation data file versions 0, 1, and 2, and
 * writes version 2 transformation data files.  Data file versions 1 and 2 are
 * identified by a version command "V 1" or "V 2", respectively.  Data file
 * version 0 is identified by having no version command.
 */

#ifndef __tfile_h__
#define __tfile_h__

#include "gpt.h"

#define TFILE_VERSION     2
#define TFILE_VERSION_MAX 2

static int tfile_version = 0;

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
	FILE *file;
};

/*
 * Create a new tload_t transformation data file structure, used for
 * reading data from transformation data files.
 */

static inline struct tload_t *tload_new(const char *filename) {
	FILE *file = fopen (filename, "r");
	struct tload_t *result = NULL;

	if (file) {
		result = (struct tload_t *) 
			malloc(sizeof(struct tload_t));
		result->filename = filename;
		result->file = file;
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

	if (file) {
		result = (struct tsave_t *) 
			malloc(sizeof(struct tsave_t));
		result->filename = filename;
		result->file = file;
	}

	fprintf(file, "# created by ALE transformation file handler version %d\n", 
			TFILE_VERSION);
	fprintf(file, "V %d\n", TFILE_VERSION);

	fclose(file);
	
	return result;
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

static inline transformation tload_next(struct tload_t *t, int is_p, 
		transformation default_transform, int *is_default) {

	transformation result = default_transform;

	*is_default = 1;

	/*
	 * Read each line of the file until we find a transformation.
	 */

	while (t && !feof(t->file)) {
		char line[1024];

		fgets(line, 1024, t->file);

		if (strlen(line) >= 1023) {
			fprintf(stderr, 
				"\ntrans-load: warning: line too long in input file\n");
		}

		switch (line[0]) {
			case ' ':
			case 0xa:
			case 0xd:
			case '\t':
			case '#':
				/* Comment or whitespace */
				break;
			case 'V':
			case 'v':
				/* Version number */
				{
					int count = sscanf(line, "V %d", &tfile_version); 
					
					if (count < 1) {
						fprintf(stderr, "Error in transformation "
								"file version command.\n");
						exit(1);
					} else if (tfile_version > TFILE_VERSION_MAX) {
						fprintf(stderr, "Unsupported transformation "
								"file version %d\n", 
								tfile_version);
						exit(1);
					}

					break;
				}
			case 'D':
			case 'd':
				/* Default transformation */
				return result;
			case 'P':
			case 'p':
				/* Projective transformation data */
				*is_default = 0;
				if (is_p == 0) {
					fprintf(stderr, "\ntrans-load: warning: "
						"Ignoring projective data for euclidean "
						"transformation.\n");
					return result;
				} else {
					double width, height, values[8];
					int count, i;
					point x[4];

					count = sscanf(line, "P %lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &width, &height,
							&values[0], &values[1], &values[2], &values[3],
							&values[4], &values[5], &values[6], &values[7]);
				
					int index = 0;
					for (int i = 0; i < 4; i++)
					for (int j = 1; j >= 0; j--) 
						x[i][j] = values[index++];

					if (count < 10) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'P'\n");

					for (i = 0; i < count - 2; i++) {
						ale_pos factor = (i % 2)
							? (result.width() / width)
							: (result.height() / height);
							
						x[i / 2][i % 2] *= factor;
					}

					if (tfile_version < 1) {
						/*
						 * Accommodate older versions
						 * of tfile.
						 */
						for (i = 0; i < 4; i++) {
							ale_pos y = x[i][0];
							  x[i][0] = x[i][1];
							  x[i][1] = y;
						}
						result.gpt_v0_set(x);
					} else {
						result.gpt_set(x);
					}

					return result;
				}
				break;
			case 'E':
			case 'e':
				/* Euclidean transformation data */
				*is_default = 0;
				{
					double width, height;
					double values[3] = {0, 0, 0};
					int count, i;
					ale_pos eu[3];

					count = sscanf(line, "E %lf%lf%lf%lf%lf",
							&width, &height,
							&values[0], &values[1], &values[2]);

					eu[1] = values[0];
					eu[0] = values[1];
					eu[2] = values[2];

					if (tfile_version < 2) {
						ale_pos t = eu[0];
						   eu[0] = eu[1];
						   eu[1] = t;
					}


					if (count < 5) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'E'\n");

					for (i = 0; (i < count - 2) && (i < 2); i++) {
						ale_pos factor = (i % 2)
							? (result.width() / width)
							: (result.height() / height);
							
						eu[i] *= factor;
					}

					if (tfile_version < 1) {
						result.eu_v0_set(eu);
					} else {		
						result.eu_set(eu);
					}

					return result;
				}
				break;
			default:
				fprintf(stderr,
					"\ntrans-load: warning: unrecognized command '%s'\n",
					line);
				break;
		}
	}

	return result;
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

static inline void tsave_next(struct tsave_t *t, transformation offset, int is_projective) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");
		
		if (is_projective) {
			int i, j;

			fprintf(t->file, "P ");
			fprintf(t->file, "%f %f ", (double) offset.width(), (double) offset.height());
			for (i = 0; i < 4; i++)
				for (j = 1; j >= 0; j--)
					fprintf(t->file, "%f ", (double) offset.gpt_get(i, j));
		} else {
			fprintf(t->file, "E ");
			fprintf(t->file, "%f %f ", (double) offset.width(), (double) offset.height());
			fprintf(t->file, "%f ",    (double) offset.eu_get(1));
			fprintf(t->file, "%f ",    (double) offset.eu_get(0));
			fprintf(t->file, "%f ",    (double) offset.eu_get(2));
		}

		fprintf(t->file, "\n");

		fclose(t->file);
	}
}

/*
 * Write information to a transformation file indicating the target output
 * file.
 */

static inline void tsave_target(struct tsave_t *t, const char *filename) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		fprintf(t->file, "# Comment: Target output file is %s\n", filename);
		
		fclose(t->file);
	}
}

/*
 * Write information to a transformation data file indicating the filename
 * of the original frame (i.e. the first frame in the sequence of input
 * frames).
 */

static inline void tsave_orig(struct tsave_t *t, const char *filename) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		fprintf(t->file, "# Comment: Original frame is %s\n", filename);
		
		fclose(t->file);
	}
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

#endif
