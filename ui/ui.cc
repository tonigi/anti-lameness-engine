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

#include "ui_wo.h"
#include "ui_tty.h"
#include "input.h"
#include "ui.h"

/*
 * See ui.h for details on these variables.
 */

ui *ui::singleton = NULL;
int ui::type = 1;   /* TTY is default */
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
		default:
			assert(0);
		}
	}

	return singleton;
}
	
void ui::handle_input(int argc, const char *argv[], const char *package, const char *short_version, const char *version) {
	input::handle(argc, argv, package, short_version, version);
}
