// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * d3/scene.h: Representation of a 3D scene.
 */

#ifndef __scene_h__
#define __scene_h__

#include "point.h"

class scene {

	/*
	 * Structure to hold a subdivisible triangle.
	 */

	struct triangle {
		point *vertices[3];
		struct triangle *neighbors[3];
		point *division_vertex;
		struct triangle *children[2];
	};

	static struct triangle *triangle_head;

	/*
	 * Structure to hold input frame information for a given level of
	 * detail.
	 */
	struct lod {

		/*
		 * Reference image for each frame.
		 */

		d2::image **reference;

		/*
		 * Scale factor
		 */

		double sf;

		/*
		 * Next element
		 */

		struct lod *next;

	};

	/*
	 * Current level-of-detail
	 */

	static struct lod *cl;

public:
	/*
	 * Initialize 3D scene from 2D scene, using 2D and 3D alignment
	 * information.
	 */
	static void init_from_d2() {

		/*
		 * Find out how many input frames there are.
		 */

		int N = d2::image_rw::count();

		/*
		 * Initialize the base level of detail
		 */

		cl = new lod;
		cl->reference = NULL;
		cl->next = NULL;
		cl->sf = 1;

		/*
		 * Initialize reference images.
		 */

		cl->reference = (d2::image **) malloc(N * sizeof(d2::image *));

		assert(cl->reference);

		for (int n = 0; n < N; n++) {
			cl->reference[n] = d2::image_rw::copy(n, "3D scene reference");
			assert(cl->reference[n]);
		}
	}

	/*
	 * Reduce the level of detail.  Return 0 when no further reduction
	 * is possible.
	 */
	static int reduce_lod() {

		int result = 1;

		/*
		 * Create a new structure for the reduced LOD.
		 */
		struct lod *nl = new lod;

		/*
		 * Find out how many input frames there are.
		 */

		int N = d2::image_rw::count();

		/*
		 * Initialize reference images and partial z-buffer arrays.
		 */

		nl->reference = (d2::image **) malloc(N * sizeof(d2::image *));

		assert(nl->reference);

		for (int n = 0; n < N; n++) {
			nl->reference[n] = cl->reference[n]->scale_by_half("3D, reduced LOD");

			assert(nl->reference[n]);

			if (nl->reference[n]->height() < 4
			 || nl->reference[n]->width () < 4)
				result = 0;

			if (nl->reference[n]->height() < 2
			 || nl->reference[n]->width () < 2)
				assert(0);
		}

		nl->sf = cl->sf * 0.5;

		nl->next = cl;

		cl = nl;

		return result;
	}

	/*
	 * Increase the level of detail.  
	 */
	static void increase_lod() {
		/*
		 * Pointer to the next higher LOD.
		 */
		struct lod *nl = cl->next;

		assert (nl != NULL);

		/*
		 * Find out how many input frames there are.
		 */

		int N = d2::image_rw::count();

		/*
		 * Delete the current LOD.
		 */

		for (int n = 0; n < N; n++)
			delete cl->reference[n];

		delete cl;

		cl = nl;
	}

	static const d2::image *view(unsigned int n) {
		assert (n < d2::image_rw::count());

		assert(0);

		return NULL;
	}

	static const d2::image *depth(unsigned int n) {
		assert (n < d2::image_rw::count());

		assert(0);

		return NULL;
	}

	/*
	 * When using a 3D scene data structure, improvements should occur in 
	 * two passes:
	 *
	 * 	(a) Moving vertices in a direction perpendicular to the surface
	 * 	normal of one of their adjacent triangles
	 *
	 * 	(b) Attempting to subdivide triangles by adding new vertices and
	 * 	moving these in a direction perpendicular to the surface normal
	 * 	of the corresponding triangle's surface normal.
	 *
	 * When neither of these approaches results in improvement, then the
	 * level of detail can be increased.
	 */
	static void reduce_cost_to_search_depth(const char *d_out[], const char *v_out[], d2::exposure *exp_out, int inc_bit) {
		int max_depth = 2;
		int improved = 1;
		int count = 0;
		void *visit_data;
		// time_t start_seconds, cur_seconds, max_seconds = 1200;

		assert (max_depth > 0);

		visit_data = init_visit_data(max_depth);

		// time(&start_seconds);

		while(reduce_lod());

		while ((improved /*&& count < 40*/) || cl->next) {

			if (inc_bit)
			for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
				if (d_out[i] != NULL) {
					const d2::image *im = depth(i);
					d2::image_rw::write_image(d_out[i], im, exp_out, 1, 1);
					delete im;
				}

				if (v_out[i] != NULL) {
					const d2::image *im = view(i);
					d2::image_rw::write_image(v_out[i], im, exp_out);
					delete im;
				}
			}

			if (!improved/* || count >= 40*/) {
				// fprintf(stderr, "scale factor: %f\n", cl->sf);
				fprintf(stderr, ".");
				assert (cl->next);

				increase_lod();
				count = 0;

#if 0
				if (cl->sf >= 0.125) {
					while (cl->next)
						increase_lod();
					count = 100000;
					improved = 0;
					continue;
				}
#endif
			}

			count++;
			improved = 0;

			// time (&cur_seconds);

			// if (cur_seconds - start_seconds > max_seconds)
			//	continue;

		}

		destroy_visit_data(visit_data);
	}

#if 0
	/*
	 * Describe a scene to a renderer
	 */
	static void describe(render *r) {
	}
#endif
};

#endif
