// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __ui_log_h__
#define __ui_log_h__

#include "ui.h"

/*
 * Logging user interface.
 */

class ui_log : public ui {
private:
	void printf(const char *format, ...) {

		fprintf(ui_stream, "ale: %u: ", (unsigned int) time(NULL));

		va_list ap;
		va_start(ap, format);
		vfprintf(ui_stream, format, ap);
		va_end(ap);

		if (format[strlen(format) - 1] != '\n')
			fprintf(ui_stream, "\n");
	}

	void update() {
	}

public:
	ui_log() {
		ui_stream = stdout;
	}

	void exp_multiplier(double m0, double m1, double m2) {
		ui::exp_multiplier(m0, m1, m2);
		printf("Exposure multiplier: %g, %g, %g\n", m0, m1, m2);
	}

	void set_steps(int count) {
		printf("%d steps to complete\n", count);
	}

	void set_steps_completed(int count) {
		printf("%d steps completed\n", count);
	}

	void set_match(double match) {
		printf("match %g / mismatch %g\n", (1 - match), match);
	}

	void set_offset(d2::transformation);
	void set_offset(d2::trans_single);

	void gs_mo(ale_pos gs_mo) {
		printf("Global search minimum overlap is %f pixels\n", (double) gs_mo);
	}

	void loading_file() {
		printf("Loading file.\n");
	}

	void exposure_1() {
		printf("Exposure pass 1.\n");
	}

	void exposure_2() {
		printf("Exposure pass 2.\n");
	}

	void prematching() {
		printf("Prematching.\n");
	}

	void postmatching() {
		printf("Postmatching.\n");
	}

	void constructing_lod_clusters(ale_pos lod) {
		printf("Constructing LOD cluster (%f)\n", (double) lod);
	}

	void global_alignment(ale_pos perturb, ale_pos lod) {
		status.perturb_size = perturb;
		printf("Global alignment (perturb=%f, lod=%f).\n", (double) perturb, (double) lod);
	}

	void aligning(ale_pos perturb, ale_pos lod) {
		perturb_timers[perturb].start();
		status.perturb_size = perturb;
		printf("Aligning (perturb=%f, lod=%f).\n", (double) perturb, (double) lod);
	}

	void following() {
		printf("Applying initial-final following logic.\n");
	}

	void set_orender_current(int num) {
		printf("Preparing to render output (index %d)\n", num);
	}


	void set_arender_current() {
		printf("Preparing to render alignment reference image.\n");
	}

	void rendering() {
		printf("Rendering.\n");
	}

	void writing_output(int num) {
		printf("Writing output (index %d)\n", num);
	}

	void ip_frame_start(unsigned int num) {
		printf("Starting Irani-Peleg frame %d.\n", num);
	}

	void ip_frame_simulate_start() {
		printf("Simulating frame.");
	}

	void ip_frame_correct_start() {
		printf("Correcting frame.\n");
	}

	void ip_write() {
		printf("Writing.\n");
	}

	void ip_step_done() {
		printf("Finished pass.");
	}

	void ip_done() {
		printf("Irani-Peleg done.\n");
	}

	void original_frame_start(const char *name) {
		printf("Starting original frame (%s)\n", name);
	}

	void original_frame_done() {
		printf("Finished original frame\n");
	}

	void supplemental_frame_start(const char *name) {
		printf("Starting supplemental frame (%s)\n", name);
	}

	void supplemental_frame_done() {
		printf("Supplemental frame done.\n");
	}

	void alignment_perturbation_level(ale_pos perturb, ale_pos lod) {
		perturb_timers[status.perturb_size].stop();
		status.perturb_size = perturb;
		perturb_timers[perturb].start();
		printf("Perturbation set to %g; LOD set to %g.\n", (double) perturb, (double) lod);
	}

	void alignment_match_ok() {
		printf("Alignment match OK.\n");
	}

	void alignment_no_match() {
		printf("Alignment failed to match.\n");
	}

	void cache(double usage, double max) {
		printf("Cache usage is %.1f%% of %.0fMB.\n", 100 * usage / max, max);
	}

	void cache_status(unsigned int i) {
		printf("Cache is full.\n");
	}

	void log_message(const char *message) {
		printf(message);
	}

};

#endif
