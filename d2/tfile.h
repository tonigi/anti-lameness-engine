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
 * This version of ALE reads transformation data file versions 0, 1, 2, and 3,
 * and writes version 2 and 3 transformation data files.  Data file versions 1
 * and higher are identified by a version command "V x", where x is the version
 * number, prior to any transformation command.  Data file version 0 is
 * identified by having no version command.
 */

#ifndef __tfile_h__
#define __tfile_h__

#include "transformation.h"

#define TFILE_VERSION     3
#define TFILE_VERSION_MAX 3

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
	pixel orig_apm;
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
 * 	IS_P is nonzero if a projective transformation is expected.
 *
 * 	DEFAULT_TRANSFORM is the default transformation result.
 *
 *	IS_DEFAULT is used to signal a non-default transformation result.
 */

static inline transformation tload_first(struct tload_t *t, int is_p, 
		transformation default_transform, int *is_default) {

	transformation result = default_transform;

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
				"\ntrans-load: Error: line too long in input file\n");
			exit(1);
		}

		first_character = fgetc(t->file);
	}
		
	if (first_character != EOF)
		ungetc(first_character, t->file);

	/*
	 * Check for version 0
	 */

	if (first_character != 'V')

		/*
		 * Must be version 0.
		 */

		return result;

	/*
	 * Obtain version from version command string.
	 */

	char line[1024];
	fgets(line, 1024, t->file);
	if (strlen(line) >= 1023) {
		fprintf(stderr, 
			"\ntrans-load: Error: line too long in input file\n");
		exit(1);
	}

	int count = sscanf(line, "V %d", &tfile_input_version); 
	
	if (count < 1) {
		fprintf(stderr, "Error in transformation "
				"file version command.\n");
		exit(1);
	} else if (tfile_input_version > TFILE_VERSION_MAX) {
		fprintf(stderr, "Unsupported transformation "
				"file version %d\n", 
				tfile_input_version);
		exit(1);
	}

	/*
	 * Handle versions lower than 3.
	 */

	if (tfile_input_version < 3)

		/*
		 * Versions lower than 3 use the default transformation
		 * for the original frame.
		 */

		return result;

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
			case 'B':
			case 'b':
				if (tfile_input_version < 3) {
					fprintf(stderr, "\ntrans-load: Error: "
						"Barrel distortion not supported "
						"for version %d input files.\n"
						"trans-load: Hint:  Use version 3 "
						"file syntax.\n", tfile_input_version);
					exit(1);
				} else {
					unsigned int count;
					unsigned int pos = 0, chars;
					unsigned int bdc;
					double dparameters[BARREL_DEGREE];
					ale_pos parameters[BARREL_DEGREE];

					count = sscanf(line, "B %u%n", &bdc, &chars);
					pos += chars;

					if (count < 1) {
						fprintf(stderr, "\ntrans-load: Error: "
								"Malformed 'B' command.\n");
						exit(1);
					}

					if (bdc > result.bd_max()) {
						fprintf(stderr, "\ntrans-load: Error: "
								"Barrel distortion degree %d "
								"is too large.  (Maximum is %d.)\n"
								"trans-load: Hint:  "
								"Reduce degree or re-compile "
								"with BD_DEGREE=%d\n", bdc, BARREL_DEGREE, bdc);
						exit(1);
					}

					for (unsigned int d = 0; d < bdc; d++) {
						count = sscanf(line + pos, "%lf%n", &dparameters[d], &chars);
						pos += chars;

						if (count < 1) {
							fprintf(stderr, "\ntrans-load: Error: "
									"Malformed 'B' command.\n");
							exit(1);
						}

						parameters[d] = dparameters[d];
					}

					result.bd_set(bdc, parameters);
				}
				break;
			case 'P':
			case 'p':
				/* Projective transformation data */
				*is_default = 0;
				if (is_p == 0) {
					fprintf(stderr, "\ntrans-load: Error: "
						"Projective data for euclidean "
						"transformation.\n"
						"trans-load: Hint:  "
						"Use command-line option --projective.\n");
					exit(1);
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
							? (result.scaled_width() / width)
							: (result.scaled_height() / height);
							
						x[i / 2][i % 2] *= factor;
					}

					result.gpt_set(x);

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

					if (count < 5) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'E'\n");

					for (i = 0; (i < count - 2) && (i < 2); i++) {
						ale_pos factor = (i % 2)
							? (result.scaled_width() / width)
							: (result.scaled_height() / height);
							
						eu[i] *= factor;
					}

					result.eu_set(eu);

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
 *
 *	IS_PRIMARY is used to differentiate primary and non-primary 
 *	transformations
 */

static inline transformation tload_next(struct tload_t *t, int is_p, 
		transformation default_transform, int *is_default, 
		int is_primary) {

	transformation result = default_transform;

	*is_default = 1;

	/*
	 * Read each line of the file until we find a transformation.
	 */

	while (t && !feof(t->file)) {

		char c = fgetc(t->file);
		if (!feof(t->file) && c != EOF)
			ungetc(c, t->file);

		if (feof(t->file)
		 || (!is_primary
		  && c != EOF
		  && c != 'F'
		  && c != 'f'
		  && c != 'Q'
		  && c != 'q')) {
			return result;
		}

		char line[1024];

		fgets(line, 1024, t->file);

		if (feof(t->file))
			return result;

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
				return result;
			case 'B':
			case 'b':
				if (tfile_input_version < 3) {
					fprintf(stderr, "\ntrans-load: Error: "
						"Barrel distortion not supported "
						"for version %d input files.\n"
						"trans-load: Hint:  Use version 3 "
						"file syntax.\n", tfile_input_version);
					exit(1);
				} else {
					unsigned int count;
					unsigned int pos = 0, chars;
					unsigned int bdc;
					ale_pos parameters[BARREL_DEGREE];
					double dparameters[BARREL_DEGREE];

					count = sscanf(line, "B %u%n", &bdc, &chars);
					pos += chars;

					if (count < 1) {
						fprintf(stderr, "\ntrans-load: Error: "
								"Malformed 'B' command.\n");
						exit(1);
					}

					if (bdc > result.bd_max()) {
						fprintf(stderr, "\ntrans-load: Error: "
								"Barrel distortion degree %d "
								"is too large.  (Maximum is %d.)\n"
								"trans-load: Hint:  "
								"Reduce degree or re-compile "
								"with BD_DEGREE=%d\n", bdc, BARREL_DEGREE, bdc);
						exit(1);
					}

					for (unsigned int d = 0; d < bdc; d++) {
						count = sscanf(line + pos, "%lf%n", &dparameters[d], &chars);
						pos += chars;

						if (count < 1) {
							fprintf(stderr, "\ntrans-load: Error: "
									"Malformed 'B' command.\n");
							exit(1);
						}

						parameters[d] = dparameters[d];
					}

					result.bd_set(bdc, parameters);
				}
				break;
			case 'Q':
			case 'q':
				if (is_primary)
					break;
			case 'P':
			case 'p':
				/* Projective transformation data */
				*is_default = 0;
				if (is_p == 0) {
					fprintf(stderr, "\ntrans-load: Error: "
						"Projective data for euclidean "
						"transformation.\n"
						"trans-load: Hint:  "
						"Use command-line option --projective.\n");
					exit(1);
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
							? (result.scaled_width() / width)
							: (result.scaled_height() / height);
							
						x[i / 2][i % 2] *= factor;
					}

					if (tfile_input_version < 1) {
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
			case 'F':
			case 'f':
				if (is_primary)
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

					if (tfile_input_version < 2) {
						ale_pos t = eu[0];
						   eu[0] = eu[1];
						   eu[1] = t;
					}


					if (count < 5) 
						fprintf(stderr, "\ntrans-load: warning:"
								"Missing args for 'E'\n");

					for (i = 0; (i < count - 2) && (i < 2); i++) {
						ale_pos factor = (i % 2)
							? (result.scaled_width() / width)
							: (result.scaled_height() / height);
							
						eu[i] *= factor;
					}

					if (tfile_input_version < 1) {
						result.eu_v0_set(eu);
					} else {		
						result.eu_set(eu);
					}

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

	fprintf(file, "# created by ALE transformation file handler version %d\n", 
			TFILE_VERSION);

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

static inline void tsave_first(struct tsave_t *t, transformation offset, int is_projective) {

	if (t == NULL)
		return;

	t->file = fopen(t->filename, "a");
	
	/*
	 * Determine the output version to use.  We use version 3 output syntax only when 
	 * necessary.  This comprises two cases:
	 *
	 *   (i)  an input file is used, and this file uses version 3 syntax.
	 *   (ii) non-degenerate barrel distortion correction is selected.
	 *
	 * (i) can be directly examined.  When (i) does not hold, (ii) can be
	 * inferred from offset.bd_count(), since this value should be constant
	 * when (i) does not hold.  XXX: This logic should be reviewed.
	 */

	if (tfile_input_version == 3 || offset.bd_count() > 0) 
		tfile_output_version = 3;
	else
		tfile_output_version = 2;
	
	
	fprintf(t->file, "# producing transformation file syntax version %d\n", tfile_output_version);
	fprintf(t->file, "V %d\n", tfile_output_version);

	fprintf(t->file, "# Comment: Target output file is %s\n", t->target);
	fprintf(t->file, "# Comment: Original frame is %s\n", t->orig);
	fprintf(t->file, "# Comment: Avg magnitude [r=%f g=%f b=%f]\n", t->orig_apm[0], t->orig_apm[1], t->orig_apm[2]);

	if (tfile_output_version < 3) {
		fclose(t->file);
		return;
	}

	if (offset.bd_count() > 0) {
		assert (tfile_output_version >= 3);
		unsigned int i;

		fprintf(t->file, "B ");
		fprintf(t->file, "%u ", offset.bd_count());
		for (i = 0; i < offset.bd_count(); i++)
			fprintf(t->file, "%f ", (double) offset.bd_get(i));
		fprintf(t->file, "\n");
	}


	if (is_projective) {
		int i, j;

		fprintf(t->file, "P ");
		fprintf(t->file, "%f %f ", (double) offset.scaled_width(), (double) offset.scaled_height());
		for (i = 0; i < 4; i++)
			for (j = 1; j >= 0; j--)
				fprintf(t->file, "%f ", (double) offset.gpt_get(i, j));
	} else {
		fprintf(t->file, "E ");
		fprintf(t->file, "%f %f ", (double) offset.scaled_width(), (double) offset.scaled_height());
		fprintf(t->file, "%f ",    (double) offset.eu_get(1));
		fprintf(t->file, "%f ",    (double) offset.eu_get(0));
		fprintf(t->file, "%f ",    (double) offset.eu_get(2));
	}

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
 * 	IS_PRIMARY indicates whether to write a primary transformation
 *
 */

static inline void tsave_next(struct tsave_t *t, transformation offset, int is_projective,
		int is_primary) {

	if (t == NULL)
		return;

	t->file = fopen(t->filename, "a");
	
	if (is_primary && offset.bd_count() > 0) {
		assert (tfile_output_version >= 3);
		unsigned int i;

		fprintf(t->file, "B ");
		fprintf(t->file, "%u ", offset.bd_count());
		for (i = 0; i < offset.bd_count(); i++)
			fprintf(t->file, "%f ", (double) offset.bd_get(i));
		fprintf(t->file, "\n");
	}

	if (is_projective) {
		int i, j;

		fprintf(t->file, is_primary ? "P " : "Q ");
		fprintf(t->file, "%f %f ", (double) offset.scaled_width(), (double) offset.scaled_height());
		for (i = 0; i < 4; i++)
			for (j = 1; j >= 0; j--)
				fprintf(t->file, "%f ", (double) offset.gpt_get(i, j));
	} else {
		fprintf(t->file, is_primary ? "E " : "F ");
		fprintf(t->file, "%f %f ", (double) offset.scaled_width(), (double) offset.scaled_height());
		fprintf(t->file, "%f ",    (double) offset.eu_get(1));
		fprintf(t->file, "%f ",    (double) offset.eu_get(0));
		fprintf(t->file, "%f ",    (double) offset.eu_get(2));
	}

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

static inline void tsave_orig(struct tsave_t *t, const char *filename, pixel apm) {
	if (t == NULL)
		return;

	t->orig = filename;
	t->orig_apm = apm;
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
