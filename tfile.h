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

#include "gpt.h"
#include <stdio.h>

#define TFILE_VERSION     2
#define TFILE_VERSION_MAX 2

static int tfile_version = 0;

struct tload_t {
	char *filename;
	FILE *file;
};

struct tsave_t {
	char *filename;
	FILE *file;
};

static inline struct tload_t *tload_new(char *filename) {
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

static inline struct tsave_t *tsave_new(char *filename) {
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

static inline void tload_delete(struct tload_t *victim) {
	if (victim)
		fclose(victim->file);
	free(victim);
}

static inline void tsave_delete(struct tsave_t *victim) {
	free(victim);
}

static inline transformation tload_next(struct tload_t *t, int lod, int is_p, 
		transformation default_transform) {

	transformation result = default_transform;

	/*
	 * Scale default initial transform for lod
	 */

	result.rescale (1 / pow(2, lod));

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
				if (is_p == 0) {
					fprintf(stderr, "\ntrans-load: warning: "
						"Ignoring projective data for euclidean "
						"transformation.\n");
					return result;
				} else {
					double width, height;
					int count, i;
					point x[4];

					count = sscanf(line, "P %lf %lf " 
							MR MR MR MR MR MR MR MR,
							&width, &height,
							&x[0][1], &x[0][0], 
							&x[1][1], &x[1][0],
							&x[2][1], &x[2][0], 
							&x[3][1], &x[3][0]);

					if (count < 10) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'P'\n");

					for (i = 0; i < count - 2; i++) {
						double factor = (i % 2)
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
							my_real y = x[i][0];
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
				{
					double width, height;
					int count, i;
					my_real eu[3];

					count = sscanf(line, "E %lf %lf " MR MR MR,
							&width, &height,
							&eu[1], &eu[0], 
							&eu[2]);

					if (tfile_version < 2) {
						double t = eu[0];
						   eu[0] = eu[1];
						   eu[1] = t;
					}


					if (count < 5) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'E'\n");

					for (i = 0; (i < count - 2) && (i < 2); i++) {
						double factor = (i % 2)
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

static inline void tsave_next(struct tsave_t *t, transformation offset, int is_projective) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");
		
		if (is_projective) {
			int i, j;

			fprintf(t->file, "P ");
			fprintf(t->file, "%f %f ", offset.width(), offset.height());
			for (i = 0; i < 4; i++)
				for (j = 1; j >= 0; j--)
					fprintf(t->file, "%f ", offset.gpt_get(i, j));
		} else {
			fprintf(t->file, "E ");
			fprintf(t->file, "%f %f ", offset.width(), offset.height());
			fprintf(t->file, "%f ",    offset.eu_get(1));
			fprintf(t->file, "%f ",    offset.eu_get(0));
			fprintf(t->file, "%f ",    offset.eu_get(2));
		}

		fprintf(t->file, "\n");

		fclose(t->file);
	}
}

static inline void tsave_target(struct tsave_t *t, char *filename) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		fprintf(t->file, "# Comment: Target output file is %s\n", filename);
		
		fclose(t->file);
	}
}

static inline void tsave_orig(struct tsave_t *t, char *filename) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		fprintf(t->file, "# Comment: Original frame is %s\n", filename);
		
		fclose(t->file);
	}
}

static inline void tsave_info(struct tsave_t *t, char *filename) {
	if (t != NULL) {
		t->file = fopen(t->filename, "a");

		fprintf(t->file, "# Comment: Supplemental frame %s\n", filename);
		
		fclose(t->file);
	}
}
