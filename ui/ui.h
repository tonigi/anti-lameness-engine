// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __ui_h__
#define __ui_h__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "../ale_pos.h"

/*
 * User interface messages.
 */

class ui_wo;

class ui {
private:
	static ui *singleton;

	/*
	 * UI type
	 *
	 * 0. stream
	 * 1. tty
	 */

	static int type;

protected:

	/*
	 * Data
	 */

	FILE *ui_stream;

	struct status_type {
		enum { UNDEFINED, LOAD_FILE, EXPOSURE_PASS_1,
			LODCLUSTER_CREATE, PREMATCH, ALIGN, POSTMATCH,
			EXPOSURE_PASS_2, RENDERA, RENDERD, RENDERO, WRITED,
			WRITEO, FRAME_DONE, SET_DONE, IP_RENDER, 
			IP_STEP_DONE, IP_UPDATE, IP_WRITE } code, orender_current;
		int arender_current;
		double match_value;
		int onum;
		int perturb_steps;
		int perturb_steps_completed;
		double exp_multiplier[3];
		double perturb_size;
		double align_lod;
		unsigned int ip_frame_num;

		status_type() {
			code = UNDEFINED;
			perturb_steps_completed = 0;

			for (int k = 0; k < 3; k++)
				exp_multiplier[k] = 1;
		}
	} status;

	/*
	 * Constructor
	 */

	ui() {
		ui_stream = stderr;
	}

	/*
	 * Print function
	 */

	virtual void printf(char *format, ...) = 0;

	/*
	 * UI update function
	 */

	virtual void update() = 0;

	/*
	 * Match format strings for textual UIs.
	 */

	char *format_string_ok() {
		return " okay (%f%% match)";
	}
	char *format_string_no_match() {
		return " no match (%f%% match)";
	}
	char *format_string_working() {
		return "      (%f%% match)";
	}

public:
	static ui *get();

	static void set_stream() {
		assert(singleton == NULL);
		type = 0;
	}

	static void set_tty() {
		assert(singleton == NULL);
		type = 1;
	}

	/*
	 * Messages from the engine
	 */

	void exp_multiplier(double m0, double m1, double m2) {
		status.exp_multiplier[0] = m0;
		status.exp_multiplier[1] = m1;
		status.exp_multiplier[2] = m2;
	}

	void exp_multiplier(double mult[3]) {
		for (int k = 0; k < 3; k++)
			status.exp_multiplier[k] = mult[k];
	}

	void set_perturb_steps(int count) {
		status.perturb_steps = count;
	}

	void set_match(double match) {
		status.match_value = (1 - match) * 100;
		update();
	}

	void loading_file() {
		status.code = status.LOAD_FILE;
		update();
	}

	void exposure_1() {
		status.code = status.EXPOSURE_PASS_1;
		update();
	}

	void exposure_2() {
		status.code = status.EXPOSURE_PASS_2;
		update();
	}

	void prematching() {
		status.code = status.PREMATCH;
		update();
	}

	void postmatching() {
		status.code = status.POSTMATCH;
		update();
	}

	void constructing_lod_clusters(ale_pos lod) {
		status.code = status.LODCLUSTER_CREATE;
		status.align_lod = lod;
		update();
	}

	void aligning(ale_pos perturb, ale_pos lod) {
		status.perturb_size = perturb;
		status.align_lod = lod;
		status.code = status.ALIGN;
		update();
	}

	void set_orender_current(int num) {
		status.onum = num;
		if (num == 0)
			status.orender_current = status.RENDERD;
		else
			status.orender_current = status.RENDERO;
	}

	void set_arender_current() {
		status.arender_current = 1;
	}

	void clear_arender_current() {
		status.arender_current = 0;
	}

	void rendering() {
		/*
		 * Current alignment rendering tasks must complete
		 * before any current output rendering tasks can
		 * start.
		 */
		if (status.arender_current) {
			status.code = status.RENDERA;
			status.arender_current = 0;
		} else	{
			status.code = status.orender_current;
		}
		update();
	}

	void writing_output(int num) {
		status.onum = num;
		if (num == 0)
			status.code = status.WRITED;
		else
			status.code = status.WRITEO;
		update();
	}

	/*
	 * Informational output
	 */

	void ip_start() {
		printf("Iterating Irani-Peleg");
	}

	void ip_frame_start(unsigned int num) {
		status.code = status.IP_RENDER;
		status.ip_frame_num = num;
		update();
	}

	void ip_update() {
		status.code = status.IP_UPDATE;
		update();
	}

	void ip_write() {
		status.code = status.IP_WRITE;
		update();
	}

	void ip_step_done() {
		status.code = status.IP_STEP_DONE;
		printf(".");
	}

	void ip_done() {
		printf("\n");
	}

	void original_frame_start(const char *name) {
		status.code = status.UNDEFINED;
		printf("Original Frame:\n");
		printf(" '%s'", name);
	}

	void original_frame_done() {
		status.code = status.FRAME_DONE;
		update();
	}

	void supplemental_frame_start(const char *name) {
		static int section_announced = 0;

		if (!section_announced) {
			printf("Supplemental Frames:\n");
			section_announced = 1;
		}

		status.code = status.UNDEFINED;
		status.perturb_steps_completed = 0;
		printf(" '%s'", name);
	}

	void supplemental_frame_done() {
		status.code = status.FRAME_DONE;
		update();
	}

	void alignment_perturbation_level(ale_pos perturb, ale_pos lod) {
		status.perturb_size = perturb;
		status.align_lod = lod;
		status.perturb_steps_completed++;
		printf(".");
	}

	void alignment_match_ok() {
		status.code = status.UNDEFINED;
		printf(format_string_ok(), status.match_value);
	}

	void alignment_no_match() {
		status.code = status.UNDEFINED;
		printf(format_string_no_match(), status.match_value);
	}

	void ale_done(double value) {
		status.code = status.UNDEFINED;
		printf("  average match %f%%", value);
		status.code = status.SET_DONE;
		update();
	}

	/*
	 * Warnings
	 */

	void warn(const char *string) {
		printf("\n\n*** Warning: %s. ***\n\n\n");
	}

	/*
	 * Errors
	 */

	void exec_failure(const char *exec, const char *arg1, const char *arg2) {
		printf("\n\n*** An error occurred while running `%s %s %s`. ***\n\n\n", exec, arg1, arg2);
		exit(1);
	}

	void fork_failure(const char *location) {
		printf("\n\n*** Could not fork in %s.  ***\n\n\n", location);
		exit(1);
	}

	void memory_error(const char *purpose) {
		printf("Unable to allocate memory for %s.\n", purpose);
		exit(1);
	}

	void memory_error_location(const char *location) {
		printf("Unable to allocate memory in %s.\n", location);
		exit(1);
	}

	void cli_not_enough(const char *option) {
		printf("\n\n*** Not enough arguments for %s ***\n\n", option);
		exit(1);
	}

	void cli_bad_arg(const char *option) {
		printf("\n\n*** Bad argument to %s ***\n\n", option);
		exit(1);
	}

	void error(const char *string) {
		printf("\n\n*** Error: %s. ***\n\n\n", string);
		exit(1);
	}

	void illegal_option(const char *string) {
		printf("\n\n*** Error: illegal option %s ***\n\n", string);
		exit(1);
	}

	void unknown_device(const char *string) {
		printf("\n\n*** Error: unknown device %s ***\n\n", string);
		exit(1);
	}

	void error_hint(const char *error, const char *hint) {
		printf("\n\n*** Error: %s", error);
		printf(  "\n*** Hint:  %s\n\n\n", hint);
		exit(1);
	}

	virtual ~ui() {
	}
};

#include "ui_wo.h"

#endif
