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
		result = malloc(sizeof(struct tload_t));
		result->filename = filename;
		result->file = file;
	}

	return result;
}

static inline struct tsave_t *tsave_new(char *filename) {
	FILE *file = fopen (filename, "w");
	struct tsave_t *result = NULL;

	if (file) {
		result = malloc(sizeof(struct tsave_t));
		result->filename = filename;
		result->file = file;
	}

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

static inline struct gpt tload_next(struct tload_t *this, int w, int h, int lod, int is_p) {
	struct gpt result;

	/*
	 * Default results
	 */

	result.eu[0] = 0;
	result.eu[1] = 0;
	result.eu[2] = 0;

//	result.x[0][0] = 0;               result.x[1][0] = 0;
//	result.x[0][1] = w / pow(2, lod); result.x[1][1] = 0;
//	result.x[0][2] = w / pow(2, lod); result.x[1][2] = h / pow(2, lod);
//	result.x[0][3] = 0;               result.x[1][3] = h / pow(2, lod);

	result.input_width = w / pow(2, lod);
	result.input_height = h / pow(2, lod);

	while (this && !feof(this->file)) {
		char line[1024];

		fgets(line, 1024, this->file);

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
			case 'D':
			case 'd':
				/* Default transformation */
				result = eu_resultant(result);
				return result;
			case 'P':
			case 'p':
				/* Projective transformation data */
				if (is_p == 0) {
					fprintf(stderr, "\ntrans-load: warning: "
						"Ignoring projective data for euclidean "
						"transformation.\n");
					result = eu_resultant(result);
					return result;
				} else {
					double width, height;
					int count, i;

					count = sscanf(line, "P %lf %lf " 
							MR MR MR MR MR MR MR MR,
							&width, &height,
							&result.x[0][0], &result.x[1][0], 
							&result.x[0][1], &result.x[1][1],
							&result.x[0][2], &result.x[1][2], 
							&result.x[0][3], &result.x[1][3]);

					if (count < 10) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'P'\n");

					for (i = 0; i < count - 2; i++) {
						double factor = (i % 2)
							? (result.input_width / width)
							: (result.input_height / height);
							
						result.x[i % 2][i / 2] *= factor;
					}

					result = gpt_resultant(result);
					return result;
				}
				break;
			case 'E':
			case 'e':
				/* Euclidean transformation data */
				{
					double width, height;
					int count, i;

					count = sscanf(line, "E %lf %lf " MR MR MR,
							&width, &height,
							&result.eu[0], &result.eu[1], 
							&result.eu[2]);

					if (count < 5) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'E'\n");

					for (i = 0; (i < count - 2) && (i < 2); i++) {
						double factor = (i % 2)
							? (result.input_width / width)
							: (result.input_height / height);
							
						result.eu[i] *= factor;
					}

					result = eu_resultant(result);
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

	result = eu_resultant(result);

	return result;
}

static inline void tsave_next(struct tsave_t *this, struct gpt offset, int is_projective) {
	if (this != NULL) {
		this->file = fopen(this->filename, "a");
		
		if (is_projective) {
			int i, j;

			fprintf(this->file, "P ");
			fprintf(this->file, "%f %f ", offset.input_width, offset.input_height);
			for (i = 0; i < 4; i++)
				for (j = 0; j < 2; j++)
					fprintf(this->file, "%f ", offset.x[j][i]);
		} else {
			int i;
		
			fprintf(this->file, "E ");
			fprintf(this->file, "%f %f ", offset.input_width, offset.input_height);
			for (i = 0; i < 3; i++) 
				fprintf(this->file, "%f ", offset.eu[i]);
		}

		fprintf(this->file, "\n");

		fclose(this->file);
	}
}

static inline void tsave_info(struct tsave_t *this, char *filename) {
	if (this != NULL) {
		this->file = fopen(this->filename, "a");

		fprintf(this->file, "# Comment: Encountered file %s\n", filename);
		
		fclose(this->file);
	}
}
