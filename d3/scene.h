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
		d2::pixel color;
		point *vertices[3];
		struct triangle *neighbors[3];
		point *division_vertex;
		struct triangle *children[2];

		triangle() {
			color = d2::pixel(0, 0, 0);

			vertices[0] = NULL;
			vertices[1] = NULL;
			vertices[2] = NULL;

			neighbors[0] = NULL;
			neighbors[1] = NULL;
			neighbors[2] = NULL;

			division_vertex = NULL;

			children[0] = NULL;
			children[1] = NULL;
		}
	};

	static struct triangle *triangle_head[2];

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

		/*
		 * Determine the bounding box of the intersections of the view
		 * pyramids with the estimated scene plane, and construct an
		 * initial (planar) triangular decomposition of the scene.
		 */

		d2::point min = d2::point(0, 0);
		d2::point max = d2::point(0, 0);

		for (int n = 0; n < N; n++) {
			d2::transformation t = d2::align::of(n);

			d2::point a[4];

			a[0] = t.transform_scaled(d2::point(0, 0));
			a[1] = t.transform_scaled(d2::point(t.scaled_height() - 1, 0));
			a[2] = t.transform_scaled(d2::point(t.scaled_height() - 1, t.scaled_width() - 1));
			a[3] = t.transform_scaled(d2::point(0, t.scaled_width() - 1));

			for (int p = 0; p < 4; p++)
			for (int d = 0; d < 2; d++) {
				if (a[p][d] < min[d])
					min[d] = a[p][d];
				if (a[p][d] > max[d])
					max[d] = a[p][d];
			}
		}

		for (int head = 0; head < 2; head++) {
			triangle_head[head] = new triangle;
			assert(triangle_head[head]);
			if (!triangle_head[0] || !triangle_head[1])
				ui::get()->memory_error("triangular approximation of 3D scene");
		}

		triangle_head[0]->vertices[0] = new point(min[0], min[1], 0);
		triangle_head[0]->vertices[1] = new point(max[0], min[1], 0);
		triangle_head[0]->vertices[2] = new point(min[0], max[1], 0);

		triangle_head[0]->neighbors[0] = triangle_head[1];

		triangle_head[0]->color = cl->reference[0]->avg_channel_magnitude();

		triangle_head[1]->vertices[0] = new point(max[0], max[1], 0);
		triangle_head[1]->vertices[1] = new point(min[0], max[1], 0);
		triangle_head[1]->vertices[2] = new point(max[0], min[1], 0);

		triangle_head[1]->neighbors[0] = triangle_head[0];

		triangle_head[1]->color = cl->reference[0]->avg_channel_magnitude();
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

	static void init_zbuf(ale_pos *zbuf, unsigned int size) {
		for (unsigned int i = 0; i < size; i++) {
			zbuf[i] = 0;
			zbuf[i] = 0 / zbuf[i];
			assert(isnan(zbuf[i]));
		}
	}

	static void fill_with_values(pt _pt, d2::image *im, ale_pos *zbuf, triangle *t) {
		if (t->division_vertex) {
			fill_with_values(_pt, im, zbuf, t->children[0]);
			fill_with_values(_pt, im, zbuf, t->children[1]);
			return;
		}

		/*
		 * Map the points of the triangle into the image space.
		 */

		point p[3];

		for (int v = 0; v < 3; v++)
			p[v] = _pt.scaled_transform(*t->vertices[v]);

		/*
		 * Determine the bounding box of the transformed vertices.
		 */

		point max = p[0], min = p[0];

		for (int v = 1; v < 3; v++)
		for (int d = 0; d < 2; d++) {
			if (max[d] < p[v][d])
				max[d] = p[v][d];
			if (min[d] > p[v][d])
				min[d] = p[v][d];
		}

		/*
		 * Intersect the bounding box with the image boundary.
		 */

		if (max[0] >= im->height())
			max[0] = im->height() - 1;
		if (min[0] >= im->height())
			min[0] = im->height() - 1;
		if (max[0] < 0)
			max[0] = 0;
		if (min[0] < 0)
			min[0] = 0;

		/*
		 * Iterate over all points in the bounding box.
		 */

		for (int i = (int) floor(min[0]); i <= (int) ceil(max[0]); i++)
		for (int j = (int) floor(min[1]); j <= (int) ceil(max[1]); j++) {

			ale_real *depth = zbuf + im->width() * i + j;

			/*
			 * Simple test for depth
			 *
			 * XXX: this doesn't work correctly in all cases.
			 */

			if (*depth > p[0][2])
				continue;
			
			/*
			 * Test for interiority
			 */

			int lower[2] = {0, 0};
			int upper[2] = {0, 0};

			for (int v = 0; v < 3; v++) {
				point cv = p[v];
				point nv = p[(v + 1) % 3];
				point test_point = point(i, j, 1);

				for (int d = 0; d < 2; d++)
				if ((test_point[d] - cv[d]) * (i - nv[d]) < 0) {
					int e = (d + 1) % 2;
					ale_pos travel = (i - cv[d]) / (nv[d] - cv[d]);
					ale_pos intersect = cv[e] + travel * (nv[e] - cv[e]);
					if (intersect <= j) 
						lower[e] = 1;
					if (intersect >= j)
						upper[e] = 1;
				}
			}

			if (!lower[0] || !upper[0] || !lower[1] || !upper[1])
				return;

			/*
			 * Assign the color value;
			 */

			im->pix(i, j) = t->color;
		}
	}

	static const d2::image *view(unsigned int n) {
		assert (n < d2::image_rw::count());

		d2::image *im = new d2::image_ale_real((int) floor(d2::align::of(n).scaled_height()),
				               (int) floor(d2::align::of(n).scaled_width()), 3);

		ale_pos *zbuf = (ale_pos *) malloc(im->height() * im->width()
				* sizeof(ale_pos));

		init_zbuf(zbuf, im->height() * im->width());

		fill_with_values(align::projective(n), im, zbuf, triangle_head[0]);
		fill_with_values(align::projective(n), im, zbuf, triangle_head[1]);

		return im;
	}

	static const d2::image *depth(unsigned int n) {
		assert (n < d2::image_rw::count());

		d2::image *im = new d2::image_ale_real((int) floor(d2::align::of(n).scaled_height()),
				               (int) floor(d2::align::of(n).scaled_width()), 3);

		ale_real *zbuf = (ale_real *) malloc(im->height() * im->width()
				* sizeof(ale_real));

		init_zbuf(zbuf, im->height() * im->width());

		fill_with_values(align::projective(n), im, zbuf, triangle_head[0]);
		fill_with_values(align::projective(n), im, zbuf, triangle_head[1]);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++)
			im->pix(i, j) = d2::pixel(1, 1, 1) * zbuf[i * im->width() + 1];

		return im;
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
		// time_t start_seconds, cur_seconds, max_seconds = 1200;

		assert (max_depth > 0);

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
