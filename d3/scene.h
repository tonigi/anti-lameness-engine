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
		d2::pixel weight;
		point *vertices[3];
		struct triangle *neighbors[3];
		point *division_vertex;
		struct triangle *children[2];

		triangle() {
			color = d2::pixel(0, 0, 0);
			weight = d2::pixel(0, 0, 0);

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

	static triangle **init_zbuf(pt _pt) {
		triangle **result = (triangle **) calloc((int) floor(_pt.scaled_width()) *
				(int) floor(_pt.scaled_height()), sizeof(triangle *));

		assert(result);

		if (!result)
			ui::get()->memory_error("Z-buffer");

		return result;
	}

	/*
	 * Z-buffer to determine the closest triangle.
	 */
	static void zbuffer(pt _pt, triangle **zbuf, triangle *t) {
		int height = (int) floor(_pt.scaled_height());
		int width  = (int) floor(_pt.scaled_width());
		
		if (t->division_vertex) {
			zbuffer(_pt, zbuf, t->children[0]);
			zbuffer(_pt, zbuf, t->children[1]);
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
		 *
		 * XXX: this may include additional points
		 *
		 * XXX: this is horribly verbose
		 */

		if (min[0] >= height)
			min[0] = height - 1;
		if (min[1] >= width)
			min[1] = width - 1;
		if (max[0] >= height)
			max[0] = height - 1;
		if (max[1] >= width)
			max[1] = width - 1;
		if (max[0] < 0)
			max[0] = 0;
		if (max[1] < 0)
			max[1] = 0;
		if (min[0] < 0)
			min[0] = 0;
		if (min[1] < 0)
			min[1] = 0;

		/*
		 * Iterate over all points in the bounding box.
		 */

		for (int i = (int) floor(min[0]); i <= (int) ceil(max[0]); i++)
		for (int j = (int) floor(min[1]); j <= (int) ceil(max[1]); j++) {

			triangle **zbuf_tri = zbuf + width * i + j;

			/*
			 * Simple test for depth
			 *
			 * XXX: this doesn't work correctly in all cases.
			 */

			if (*zbuf_tri && _pt.scaled_transform(*(*zbuf_tri)->vertices[0])[2] > p[0][2])
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
				if ((test_point[d] - cv[d]) * (test_point[d] - nv[d]) < 0) {
					int e = (d + 1) % 2;
					ale_pos travel = (test_point[d] - cv[d]) / (nv[d] - cv[d]);
					ale_pos intersect = cv[e] + travel * (nv[e] - cv[e]);
					if (intersect <= test_point[e]) 
						lower[e] = 1;
					if (intersect >= test_point[e])
						upper[e] = 1;
				}
			}

			if (!lower[0] || !upper[0] || !lower[1] || !upper[1])
				continue;

			/*
			 * Assign a new triangle to the zbuffer
			 */

			*zbuf_tri = t;
		}
	}

	/*
	 * Color a 3D scene using value averaging.
	 */
	static void color_average() {
		/*
		 * Iterate over all frames
		 */
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
			/*
			 * Z-buffer to map points to triangles
			 */
			pt _pt = align::projective(n);
			triangle **zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all triangles
			 */
			assert(0);
		}
	}

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
			if (!triangle_head[head])
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

	static const d2::image *view(unsigned int n) {
		assert (n < d2::image_rw::count());

		d2::image *im = new d2::image_ale_real((int) floor(d2::align::of(n).scaled_height()),
				               (int) floor(d2::align::of(n).scaled_width()), 3);

		triangle **zbuf = init_zbuf(align::projective(n));

		zbuffer(align::projective(n), zbuf, triangle_head[0]);
		zbuffer(align::projective(n), zbuf, triangle_head[1]);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++)
		if (zbuf[i * im->width() + j])
			im->pix(i, j) = zbuf[i * im->width() + j]->color;

		free(zbuf);

		return im;
	}

	static const d2::image *depth(unsigned int n) {
		assert (n < d2::image_rw::count());

		d2::image *im = new d2::image_ale_real((int) floor(d2::align::of(n).scaled_height()),
				               (int) floor(d2::align::of(n).scaled_width()), 3);

		pt _pt = align::projective(n);

		triangle **zbuf = init_zbuf(align::projective(n));

		zbuffer(_pt, zbuf, triangle_head[0]);
		zbuffer(_pt, zbuf, triangle_head[1]);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++)
		if (zbuf[i * im->width() + j])
			im->pix(i, j) = d2::pixel(1, 1, 1) * _pt.scaled_transform(*zbuf[i * im->width() + j]->vertices[0])[2];

		free(zbuf);

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
	 * 	of the corresponding triangle.
	 *
	 * When neither of these approaches results in improvement, then the
	 * level of detail can be increased.
	 */
	static void reduce_cost_to_search_depth(const char *d_out[], const char *v_out[], d2::exposure *exp_out, int inc_bit) {
		int max_depth = 2;
		int improved = 1;
		int count = 0;

		assert (max_depth > 0);

		while(reduce_lod());

		while ((improved /*&& count < 40*/) || cl->next) {

			/*
			 * Write output incrementally, if desired.
			 */

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

			/*
			 * Increase LOD if no improvements were achieved in the
			 * most recent pass at the previous LOD.
			 */

			if (!improved) {
				fprintf(stderr, ".");
				assert (cl->next);
				increase_lod();
				count = 0;
			}

			count++;
			improved = 0;

			/*
			 * Try improving the result by moving existing vertices.
			 */

			/*
			 * Try improving the result by adding a vertex.
			 */

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
