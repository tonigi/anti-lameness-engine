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

#ifndef __ui_wo_h__
#define __ui_wo_h__

#include "ui.h"

/*
 * Write-once user interface.
 */

class ui_wo : public ui {
private:
	void printf(char *format, ...) {
		va_list ap;
		va_start(ap, format);
		vfprintf(ui_stream, format, ap);
		va_end(ap);
	}

	void update() {
		if (status.code == status.FRAME_DONE
		 || status.code == status.SET_DONE) {
			printf(".\n");
		}
	}

public:
	ui_wo() {
	}
};

#endif
