// Copyright 2006 David Hilvert <dhilvert@auricle.dyndns.org>,
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

/*
 * d3/focus.h: Implementation of defocusing logic.
 */

#ifndef __focus_h__
#define __focus_h__

class focus {
private:
	struct entry {
		int type;
		double distance;
		double px, py;
		double focal_range;
		double vertical_gradient;
		double horizontal_gradient;
		double start_depth;
		double end_depth;
		double aperture;
		double sample_count;
	};

	static unsigned int camera_index;
	static std::vector<std::vector<entry> > focus_list;

public:

	struct result {
		double focal_depth;
		double aperture;
		double sample_count;
	};

	static void add_region(unsigned int type, double distance, double px, double py, 
			unsigned int ci, double fr, double ht, double vt, double sd, double ed,
			double ap, double sc) {

		if (focus_list.size() <= ci)
			focus_list.resize(ci + 1);

		entry e = { type, distance, px, py, fr, ht, vt, sd, ed, ap, sc };
		
		focus_list[ci].push_back(e);
	}

	static result get(const d2::image *depth, int i, int j) {

		assert(0);
		
		/*
		 * Trivial implementation
		 */

		result r = { depth->get_pixel(i, j)[0], 0, 1 };

		return r;
	}

	static void set_camera(unsigned int c) {
		camera_index = c;
	}
};

#endif
