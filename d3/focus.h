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
		double start_x;
		double end_x;
		double start_y;
		double end_y;
		double aperture;
		unsigned int sample_count;
		unsigned int focal_statistic;
	};

	static unsigned int _uses_medians;
	static unsigned int _max_samples;
	static unsigned int camera_index;
	static std::vector<std::vector<entry> > focus_list;

public:

	struct result {
		double focal_distance;
		double aperture;
		unsigned int sample_count;
		unsigned int statistic;
	};

	static void add_region(unsigned int type, double distance, double px, double py, 
			unsigned int ci, double fr, double ht, double vt, double sd, double ed,
			double sx, double ex, double sy, double ey, double ap, unsigned int sc, unsigned int fs) {

		if (fs)
			_uses_medians = 1;

		if (sc > _max_samples)
			_max_samples = sc;

		if (focus_list.size() <= ci)
			focus_list.resize(ci + 1);

		entry e = { type, distance, px, py, fr, ht, vt, sd, ed, sx, ex, sy, ey, ap, sc, fs };
		
		focus_list[ci].push_back(e);
	}

	static int is_trivial() {
		return (focus_list.size() == 0);
	}

	static int uses_medians() {
		return _uses_medians;
	}

	static unsigned int max_samples() {
		return _max_samples;
	}

	static result get(const d2::image *depth, int i, int j) {

		ale_pos d = depth->get_pixel(i, j)[0];

		std::vector<entry> *l = &(focus_list[camera_index]);

		/*
		 * Initialize default focus result.
		 */

		result r = { d, 0, 1 };

		/*
		 * Check for relevant focus regions.
		 */

		for (unsigned int n = 0; n < l->size(); n++) {
			entry *e = &((*l)[n]);

			if (i >= e->start_y
			 && i <= e->end_y
			 && j >= e->start_x
			 && j <= e->end_x
			 && ((d >= -e->end_depth
			   && d <= -e->start_depth)
			  || (isnan(d) && (isnan(e->end_depth)
				        || isnan(e->start_depth))))) {
				d2::point focus_origin;
				ale_pos distance_at_focus_origin;

				if (e->type == 0) {
					/*
					 * Distance at frame center.
					 */
					focus_origin = d2::point(depth->height() / 2, depth->width() / 2);
					distance_at_focus_origin = -e->distance;
				} else if (e->type == 1) {
					/*
					 * Distance at a given point.
					 */
					focus_origin = d2::point(e->py, e->px);
					distance_at_focus_origin = depth->get_bl(d2::point(e->py, e->px))[0];
				} else 
					assert(0);

				r.focal_distance = distance_at_focus_origin + (d2::point(i, j) - focus_origin)
					                             .dproduct(d2::point(-e->vertical_gradient,
											 -e->horizontal_gradient));

				/*
				 * Adjust according to focal_range.
				 */

				ale_pos rel_dist = d - r.focal_distance;

				if (fabs(rel_dist) < e->focal_range / 2) {
					r.focal_distance = d;
				} else if (rel_dist > 0) {
					r.focal_distance += e->focal_range / 2;
				} else if (rel_dist < 0) {
					r.focal_distance -= e->focal_range / 2;
				}

				r.aperture = e->aperture;
				r.sample_count = e->sample_count;
				r.statistic = e->focal_statistic;

				break;
			}
		}

		return r;
	}

	static void set_camera(unsigned int c) {
		camera_index = c;
	}
};

#endif
