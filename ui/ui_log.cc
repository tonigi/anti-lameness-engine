// Copyright 2007 David Hilvert <dhilvert@auricle.dyndns.org>, 
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

#include "ui.h"
#include "ui_log.h"
#include "d2.h"

void ui_log::set_offset(d2::trans_single offset, int marker) {
	if (offset.is_projective()) {
		printf("Alignment [%s %f %f %f %f %f %f %f %f %f %f]\n",
				marker ? "Q" : "P",
				(double) offset.scaled_width(),
				(double) offset.scaled_height(),
				(double) offset.gpt_get(0, 1),
				(double) offset.gpt_get(0, 0),
				(double) offset.gpt_get(1, 1),
				(double) offset.gpt_get(1, 0),
				(double) offset.gpt_get(2, 1),
				(double) offset.gpt_get(2, 0),
				(double) offset.gpt_get(3, 1),
				(double) offset.gpt_get(3, 0) );
	} else {
		printf("Alignment [%s %f %f %f %f %f]\n",
				marker ? "F" : "E",
				(double) offset.scaled_width(),
				(double) offset.scaled_height(),
				(double) offset.eu_get(1),
				(double) offset.eu_get(0),
				(double) offset.eu_get(2) );
	}
}

void ui_log::set_offset(d2::trans_single offset) {
	set_offset(offset, 0);
}

void ui_log::set_offset(d2::transformation offset) {
	set_offset(offset.get_element(offset.get_current_index()), 0);
}
