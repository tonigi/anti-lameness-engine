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

#include "ui_wo.h"
#include "ui_tty.h"
#include "ui_log.h"
#include "ui_quiet.h"
#include "ui_gl.h"
#include "input.h"
#include "ui.h"
#include "d2.h"

/*
 * See ui.h for details on these variables.
 */

ui *ui::singleton = NULL;
int ui::type = 5;   /* auto is default */
int ui::output_performance_data = 0;

ui *ui::get() {
	if (singleton == NULL) {
		switch (type) {
		case 0:
			singleton = new ui_wo();
			break;
		case 1:
			try {
				singleton = new ui_tty();
			} catch (...) {
				singleton = new ui_wo();
			}
			break;
		case 2:
			singleton = new ui_log();
			break;
		case 3:
			singleton = new ui_quiet();
			break;
		case 4:
			try {
				singleton = new ui_gl();
			} catch (const char *error) {
				fprintf(stderr, "Error initializing GL interface: %s.\n", error);
				exit(1);
			} catch (...) {
				fprintf(stderr, "Error initializing GL interface.\n");
				exit(1);
			}
			break;
		case 5:
		{
			const char *gl_default = getenv("ALE_GL_UI_DEFAULT");
			if (gl_default && !strcmp(gl_default, "1")) {
				try {
					singleton = new ui_gl();
				} catch (...) {
					try {
						singleton = new ui_tty();
					} catch (...) {
						singleton = new ui_wo();
					}
				}
			} else {
				try {
					singleton = new ui_tty();
				} catch (...) {
					singleton = new ui_wo();
				}
			}
			break;
		}
		default:
			assert(0);
		}
	}

	return singleton;
}
	
void ui::handle_input(int argc, const char *argv[], const char *package, const char *short_version, const char *version) {
	input::handle(argc, argv, package, short_version, version);
}

void ui::set_offset(d2::trans_single offset) {
}
void ui::set_offset(d2::transformation offset) {
}
