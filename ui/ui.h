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
#include "../config.h"
#if HAVE_TIME_H
#include <time.h>
#endif
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include <map>

/*
 * Time structures.
 */

namespace d2 {
	struct trans_multi;
	typedef trans_multi transformation;
}

class ale_timer_t {
#if HAVE_GETTIMEOFDAY
	struct timeval tv;
#else
	time_t tt;
#endif
	double total;

public:

	ale_timer_t() {
		total = 0;
	}

	void start() {
#if HAVE_GETTIMEOFDAY
		gettimeofday(&tv, NULL);
#else
		tt = time(NULL);
#endif
	}

	void stop() {
#if HAVE_GETTIMEOFDAY
		timeval t;

		gettimeofday(&t, NULL);

		t.tv_sec -= tv.tv_sec;
		t.tv_usec -= tv.tv_usec;
		
		total += t.tv_sec + ((double) 1 / (double) 1000000) * t.tv_usec;
#else
		time_t t = time(NULL);
		t -= tt;

		total += t;
#endif
	}

	double get_total() {
		return total;
	}
};

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
	 * 2. log
	 * 3. quiet
	 */

	static int type;

	static int output_performance_data;

protected:

	/*
	 * Data
	 */

	FILE *ui_stream;

	struct status_type {
		enum { 
			/*
			 * Special
			 */
			
			UNDEFINED, 
			
			/*
			 * Incremental rendering.
			 */
			
			LOAD_FILE, EXPOSURE_PASS_1,
			LODCLUSTER_CREATE, PREMATCH, ALIGN, POSTMATCH,
			EXPOSURE_PASS_2, RENDERA, RENDERD, RENDERO, WRITED,
			WRITEO, FRAME_DONE, SET_DONE, 
			
			/*
			 * Irani-Peleg rendering.
			 */
			
			IP_RENDER, IP_STEP_DONE, IP_UPDATE, IP_WRITE,

			/*
			 * 3D.
			 */

			D3_CONTROL_POINT_SOLVE, D3_SUBDIVIDING_SPACE,
			D3_UPDATING_OCCUPANCY, D3_RENDER

		} code, orender_current;

		int arender_current;
		double match_value;
		int onum;
		int steps;
		int steps_completed;
		double exp_multiplier[3];
		double perturb_size;
		double align_lod;
		double mc;
		unsigned int frame_num;
		unsigned int irani_peleg_stage;
		unsigned int secondary_frame_num;
		unsigned int view_num;
		unsigned int x_coordinate, y_coordinate;;
		unsigned int filtering, focusing;
		unsigned int space_num;
		unsigned int total_spaces;
		double cp_max_perturb;
		double cp_min_perturb;
		double cp_cur_perturb;
		double cp_cur_error;

		status_type() {
			code = UNDEFINED;
			steps_completed = 0;

			for (int k = 0; k < 3; k++)
				exp_multiplier[k] = 1;
		}
	} status;

	/*
	 * Performance data
	 */

	ale_timer_t d2_align_sample;
	ale_timer_t d2_align_sim;
	ale_timer_t d2_incremental;
	ale_timer_t d2_irani_peleg;
	std::map<double,ale_timer_t> perturb_timers;

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
		return " okay (%9.6f%% match)";
	}
	char *format_string_no_match() {
		return " no match (%9.6f%% match)";
	}
	char *format_string_working() {
		return "      (%9.6f%% match)";
	}

public:

	/*
	 * Handle options and other user input.
	 */

	static void handle_input(int argc, const char *argv[], const char *package,
			const char *short_version, const char *version);

	static ui *get();

	static void set_stream() {
		assert(singleton == NULL);
		type = 0;
	}

	static void set_tty() {
		assert(singleton == NULL);
		type = 1;
	}

	static void set_log() {
		assert (singleton == NULL);
		type = 2;
	}

	static void set_quiet() {
		assert(singleton == NULL);
		type = 3;
	}

	static void set_profile() {
		output_performance_data = 1;
	}

	/*
	 * Messages from the engine
	 */

	virtual void identify_output(const char *name) {
		printf("Output file will be '%s'.\n", name);
	}

	virtual void d2_align_sim_start() {
		d2_align_sim.start();
	}

	virtual void d2_align_sim_stop() {
		d2_align_sim.stop();
	}

	virtual void d2_align_sample_start() {
		d2_align_sample.start();
	}

	virtual void d2_align_sample_stop() {
		d2_align_sample.stop();
	}

	virtual void d2_incremental_start() {
		d2_incremental.start();
	}

	virtual void d2_incremental_stop() {
		d2_incremental.stop();
	}

	virtual void d2_irani_peleg_start() {
		d2_irani_peleg.start();
	}

	virtual void d2_irani_peleg_stop() {
		d2_irani_peleg.stop();
	}

	virtual void exp_multiplier(double m0, double m1, double m2) {
		status.exp_multiplier[0] = m0;
		status.exp_multiplier[1] = m1;
		status.exp_multiplier[2] = m2;
	}

	void exp_multiplier(double mult[3]) {
		exp_multiplier(mult[0], mult[1], mult[2]);
	}

	virtual void set_steps(int count) {
		status.steps = count;
	}

	virtual void set_steps_completed(int count) {
		status.steps_completed = count;
	}

	virtual void set_match(double match) {
		status.match_value = (1 - match) * 100;
		update();
	}

	virtual void set_offset(d2::transformation offset);

	virtual void loading_file() {
		status.code = status.LOAD_FILE;
		update();
	}

	virtual void exposure_1() {
		status.code = status.EXPOSURE_PASS_1;
		update();
	}

	virtual void exposure_2() {
		status.code = status.EXPOSURE_PASS_2;
		update();
	}

	virtual void prematching() {
		status.code = status.PREMATCH;
		update();
	}

	virtual void postmatching() {
		status.code = status.POSTMATCH;
		update();
	}

	virtual void constructing_lod_clusters(ale_pos lod) {
		status.code = status.LODCLUSTER_CREATE;
		status.align_lod = lod;
		update();
	}

	virtual void aligning(ale_pos perturb, ale_pos lod) {
		perturb_timers[perturb].start();
		status.perturb_size = perturb;
		status.align_lod = lod;
		status.code = status.ALIGN;
		update();
	}

	virtual void set_orender_current(int num) {
		status.onum = num;
		if (num == 0)
			status.orender_current = status.RENDERD;
		else
			status.orender_current = status.RENDERO;
	}

	virtual void set_arender_current() {
		status.arender_current = 1;
	}

	virtual void clear_arender_current() {
		status.arender_current = 0;
	}

	virtual void rendering() {
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

	virtual void writing_output(int num) {
		status.onum = num;
		if (num == 0)
			status.code = status.WRITED;
		else
			status.code = status.WRITEO;
		update();
	}

	virtual void d3_control_point_data(double max_perturbation, double min_perturbation, double cur_perturbation, 
			double current_error) {
		status.cp_max_perturb = max_perturbation;
		status.cp_min_perturb = min_perturbation;
		status.cp_cur_perturb = cur_perturbation;
		status.cp_cur_error   = current_error;
		update();
	}

	virtual void d3_control_point_step() {
		printf(".");
		update();
	}

	virtual void d3_subdivision_status(unsigned int primary_frame, unsigned int secondary_frame,
			unsigned int i, unsigned int j) {
		status.code = status.D3_SUBDIVIDING_SPACE;
		status.frame_num = primary_frame;
		status.secondary_frame_num = secondary_frame;
		status.y_coordinate = i;
		status.x_coordinate = j;

		update();
	}

	virtual void d3_total_spaces(int total_spaces) {
		status.total_spaces = total_spaces;
	}

	virtual void d3_increment_spaces() {
		status.total_spaces++;
	}

	virtual void d3_occupancy_status(int frame) {
		status.code = status.D3_UPDATING_OCCUPANCY;
		status.frame_num = frame;
		status.space_num = 0;
		update();
	}

	virtual void d3_increment_space_num() {
		status.space_num++;
		update();
	}

	virtual void d3_render_status(int filter, int focus, int frame, int view, int i, int j, int space) {
		status.code = status.D3_RENDER;

		status.filtering = filter;
		status.focusing = focus;

		status.frame_num = frame;
		status.view_num = view;
		status.y_coordinate = i;
		status.x_coordinate = j;

		status.space_num = space;

		update();
	}


	/*
	 * Informational output
	 */

	virtual void ip_start() {
		printf("Iterating Irani-Peleg");
	}

	virtual void ip_frame_start(unsigned int num) {
		status.code = status.IP_RENDER;
		status.frame_num = num;
		status.irani_peleg_stage = 0;
	}

	virtual void ip_frame_simulate_start() {
		status.irani_peleg_stage = 1;
		update();
	}

	virtual void ip_frame_correct_start() {
		status.irani_peleg_stage = 2;
		update();
	}

	virtual void ip_update() {
		status.code = status.IP_UPDATE;
		update();
	}

	virtual void ip_write() {
		status.code = status.IP_WRITE;
		update();
	}

	virtual void ip_step_done() {
		status.code = status.IP_STEP_DONE;
		printf(".");
	}

	virtual void ip_done() {
		printf("\n");
	}

	virtual void original_frame_start(const char *name) {
		status.code = status.UNDEFINED;
		printf("Original Frame:\n");
		printf(" '%s'", name);
	}

	virtual void original_frame_done() {
		status.code = status.FRAME_DONE;
		update();
	}

	virtual void supplemental_frame_start(const char *name) {
		static int section_announced = 0;

		if (!section_announced) {
			printf("Supplemental Frames:\n");
			section_announced = 1;
		}

		status.code = status.UNDEFINED;
		status.steps_completed = 0;
		printf(" '%s'", name);
	}

	virtual void supplemental_frame_done() {
		status.code = status.FRAME_DONE;
		update();
	}

	virtual void alignment_monte_carlo_parameter(ale_pos mc) {
		status.mc = (mc > 1) ? 1 : mc;
	}

	virtual void alignment_perturbation_level(ale_pos perturb, ale_pos lod) {
		perturb_timers[status.perturb_size].stop();
		perturb_timers[perturb].start();
		status.perturb_size = perturb;
		status.align_lod = lod;
		status.steps_completed++;
		printf(".");
	}

	virtual void alignment_match_ok() {
		status.code = status.UNDEFINED;
		printf(format_string_ok(), status.match_value);
	}

	virtual void alignment_no_match() {
		status.code = status.UNDEFINED;
		printf(format_string_no_match(), status.match_value);
	}

	virtual void ale_2d_done(double value) {
		status.code = status.UNDEFINED;
		printf("Average match: %f%%", value);
		status.code = status.SET_DONE;
		update();
		if (output_performance_data) {
			printf("\n");
			printf("Real time measurements\n");
			printf("======================\n");
			printf("\n");
			printf("Alignment (sampling)   :  %f s\n", d2_align_sample.get_total());
			printf("Alignment (checking)   :  %f s\n", d2_align_sim.get_total());
			printf("Incremental rendering  :  %f s\n", d2_incremental.get_total());
			printf("Irani-Peleg rendering  :  %f s\n", d2_irani_peleg.get_total());
			printf("\n");

			printf("Details (local alignment)\n");
			printf("-------------------------\n");

			int have_details = 0;
			for (std::map<double,ale_timer_t>::iterator i = perturb_timers.begin(); 
			     i != perturb_timers.end(); i++) {
				if (i->second.get_total() == 0.0
				 && i == perturb_timers.begin())
					continue;

				printf("Alignment (perturb %f): %f s\n", 
						i->first, i->second.get_total());
				
				have_details = 1;
			}

			if (!have_details) {
				printf("No local alignment performed.\n");
			}

			printf("\n");
		}
	}

	virtual void d3_start() {
		status.code = status.UNDEFINED;
		printf("Rendering 3D");
		update();
	}

	virtual void d3_control_point_solve() {
		status.code = status.D3_CONTROL_POINT_SOLVE;
		update();
	}

	virtual void d3_init_view_angle(double angle) {
		status.code = status.UNDEFINED;
		update();
	}

	virtual void d3_final_view_angle(double angle) {
		status.code = status.UNDEFINED;
		update();
	}

	virtual void d3_control_point_solve_done() {
		status.code = status.UNDEFINED;
		update();
	}

	virtual void d3_subdividing_space() {
		status.code = status.D3_SUBDIVIDING_SPACE;
		update();
	}

	virtual void d3_subdividing_space_done() {
		status.code = status.UNDEFINED;
		update();
	}

	virtual void d3_updating_occupancy() {
		status.code = status.D3_UPDATING_OCCUPANCY;
		update();
	}

	virtual void d3_updating_occupancy_done() {
		status.code = status.UNDEFINED;
		update();
	}

	virtual void d3_writing_output(const char *name) {
		static int section_announced = 0;

		if (!section_announced) {
			printf(":\n");
			section_announced = 1;
		}

		printf(" '%s'", name);

		update();
	}

	virtual void d3_writing_output_done() {
		status.code = status.UNDEFINED;
		printf(".\n");
		update();
	}

	/*
	 * Warnings
	 */

	virtual void warn(const char *string) {
		printf("\n\n*** Warning: %s. ***\n\n\n");
	}

	/*
	 * Errors
	 */

	virtual void exec_failure(const char *exec, const char *arg1, const char *arg2) {
		printf("\n\n*** An error occurred while running `%s %s %s`. ***\n\n\n", exec, arg1, arg2);
		exit(1);
	}

	virtual void fork_failure(const char *location) {
		printf("\n\n*** Could not fork in %s.  ***\n\n\n", location);
		exit(1);
	}

	virtual void memory_error(const char *purpose) {
		printf("\n\n*** Unable to allocate memory for %s. ***\n\n\n", purpose);
		exit(1);
	}

	virtual void memory_error_location(const char *location) {
		printf("\n\n*** Unable to allocate memory in %s.\n\n\n", location);
		exit(1);
	}

	virtual void cli_not_enough(const char *option) {
		printf("\n\n*** Not enough arguments for `%s' ***\n\n\n", option);
		exit(1);
	}

	virtual void cli_bad_arg(const char *option) {
		printf("\n\n*** Bad argument to `%s' ***\n\n", option);
		exit(1);
	}

	virtual void error(const char *string) {
		printf("\n\n*** Error: %s. ***\n\n\n", string);
		exit(1);
	}

	virtual void illegal_option(const char *string) {
		printf("\n\n*** Error: illegal option %s ***\n\n", string);
		exit(1);
	}

	virtual void unknown_device(const char *string) {
		printf("\n\n*** Error: unknown device %s ***\n\n", string);
		exit(1);
	}

	virtual void error_hint(const char *error, const char *hint) {
		printf("\n\n*** Error: %s", error);
		printf(  "\n*** Hint:  %s\n\n\n", hint);
		exit(1);
	}

	virtual ~ui() {
	}
};

#include "ui_wo.h"

#endif
