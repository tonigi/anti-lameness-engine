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

		/*
		 * Stats
		 */

		d2::pixel color;
		d2::pixel weight;
		d2::pixel aux_stat;

		/*
		 * Connectivity
		 */

		point *vertices[3];
		struct triangle *neighbors[3];
		struct triangle *parent;
		int division_vertex;
		point *division_new_vertex;
		struct triangle *children[2];

		/*
		 * Constructor
		 */

		triangle() {
			color = d2::pixel(0, 0, 0);
			weight = d2::pixel(0, 0, 0);

			vertices[0] = NULL;
			vertices[1] = NULL;
			vertices[2] = NULL;

			neighbors[0] = NULL;
			neighbors[1] = NULL;
			neighbors[2] = NULL;

			parent = NULL;

			division_new_vertex = NULL;

			children[0] = NULL;
			children[1] = NULL;
		}

		void write_tree(int print_vertices = 0, int level = 0) {
			for (int i = 0; i < level; i++)
				fprintf(stderr, " ");
			fprintf(stderr, "%p (%p %p %p)\n", this, neighbors[0], neighbors[1], neighbors[2]);

			if (print_vertices) {
				for (int i = 0; i < level; i++)
					fprintf(stderr, " ");
				fprintf(stderr, "\t[%f %f %f] [%f %f %f] [%f %f %f]\n",
						(*vertices[0])[0],
						(*vertices[0])[1],
						(*vertices[0])[2],
						(*vertices[1])[0],
						(*vertices[1])[1],
						(*vertices[1])[2],
						(*vertices[2])[0],
						(*vertices[2])[1],
						(*vertices[2])[2]);
			}

			if (children[0])
				children[0]->write_tree(print_vertices, level + 1);
			if (children[1])
				children[1]->write_tree(print_vertices, level + 1);
		}

		/*
		 * Statistic initialization
		 */

		void init_color_counters() {
			if (children[0])
				children[0]->init_color_counters();
			if (children[1])
				children[1]->init_color_counters();

			color = d2::pixel(0, 0, 0);
			weight = d2::pixel(0, 0, 0);
		}

		void init_aux_stats(d2::pixel value = d2::pixel(0, 0, 0)) {
			if (children[0])
				children[0]->init_aux_stats(value);
			if (children[1])
				children[1]->init_aux_stats(value);

			aux_stat = value;
		}


		/*
		 * Get the neighbor link from a given neighbor that references
		 * the 'this' object.
		 */

		int self_ref_from_neighbor(int n) {
			triangle *t = neighbors[n];

			assert (t);

			for (int i = 0; i < 3; i++) {
				if (t->neighbors[i] == this)
					return i;
			}

			fprintf(stderr, "error in %p\n", this);
			assert(0);

			return -1;
		}


		/*
		 * Handle the internal data details of splitting.
		 */
		void split_internals(int v, point *nv) {
			assert (children[0] == NULL);
			assert (children[1] == NULL);
			assert (division_new_vertex == NULL);

			division_vertex = v;
			division_new_vertex = nv;

			children[0] = new triangle(*this);
			children[1] = new triangle(*this);

			assert (children[0]);
			assert (children[1]);

			for (int c = 0; c < 2; c++) {
				children[c]->parent = this;
				children[c]->vertices[(division_vertex + 1 + c) % 3] = division_new_vertex;
				children[c]->neighbors[(division_vertex + 2 - c) % 3] = children[(c + 1) % 2];
				children[c]->children[0] = NULL;
				children[c]->children[1] = NULL;
				children[c]->division_new_vertex = NULL;
			}

			for (int i = 0; i < 2; i++) {
				int vv = (v + 1 + i) % 3;
				if (!neighbors[vv])
					continue;

				int self_ref = self_ref_from_neighbor(vv);

				neighbors[vv]->neighbors[self_ref] = children[i];
			}

		}

		void split(int v, point nv) {
			point *nvp = new point(nv);

			assert (nvp);

			split_internals(v, nvp);

			if (!neighbors[v])
				return;

			int self_ref = self_ref_from_neighbor(v);

			neighbors[v]->split_internals(self_ref, nvp);

			children[0]->neighbors[v] = neighbors[v]->children[1];
			children[1]->neighbors[v] = neighbors[v]->children[0];

			neighbors[v]->children[0]->neighbors[self_ref] = children[1];
			neighbors[v]->children[1]->neighbors[self_ref] = children[0];
		}

		/*
		 * Split a triangle between the given vertex and the opposite edge's midpoint.
		 *
		 * XXX: it might be better to use an angle bisector.
		 */

		void split(int v) {
			split(v, (*vertices[(v + 1) % 3] + *vertices[(v + 2) % 3]) / 2);
		}


		/*
		 * Split a triangle at a random vertex among those with largest
		 * angle.
		 */
		void split() {
			ale_pos angles[3] = { vertices[0]->anglebetw(*vertices[1], *vertices[2]),
				              vertices[1]->anglebetw(*vertices[0], *vertices[2]),
					      vertices[2]->anglebetw(*vertices[0], *vertices[1]) };

			for (int v = 0; v < 3; v++)
			if  (angles[v] > angles[(v + 1) % 3]
			  && angles[v] > angles[(v + 2) % 3]) {
				split(v);
				return;
			}

			for (int v = 0; v < 3; v++)
			if  (angles[v] < angles[(v + 1) % 3]
			  && angles[v] < angles[(v + 2) % 3]) {
				split((v + 1 + rand() % 2) % 3);
				return;
			}

			split(rand() % 3);
		}

		void unsplit_internals() {
			assert(children[0]);
			assert(children[1]);

			for (int i = 0; i < 2; i++) {
				int vv = (division_vertex + 1 + i) % 3;

				neighbors[vv] = children[i]->neighbors[vv];

				if (neighbors[vv] == NULL)
					continue;

				int self_ref = children[i]->self_ref_from_neighbor(vv);

				neighbors[vv]->neighbors[self_ref] = this;
			}

			delete children[0];
			delete children[1];

			division_new_vertex = NULL;
			children[0] = NULL;
			children[1] = NULL;
		}

		void unsplit() {
			if (!division_new_vertex) {
				assert (!children[0]);
				assert (!children[1]);
				return;
			}

			assert(children[0]);
			assert(children[1]);

			children[0]->unsplit();
			children[1]->unsplit();

			assert(division_new_vertex);

			delete division_new_vertex;

			unsplit_internals();

			if (!neighbors[division_vertex])
				return;

			assert(neighbors[division_vertex]->children[0]);
			assert(neighbors[division_vertex]->children[1]);

			neighbors[division_vertex]->children[0]->unsplit();
			neighbors[division_vertex]->children[1]->unsplit();

			neighbors[division_vertex]->unsplit_internals();
		}

		int split_on_aux() {

			if (children[0] && children[1]) {

				// int a = rand() % 2;
				int a = 0;
				int b = 1 - a;

				return (children[a]->split_on_aux()
				      | children[b]->split_on_aux());
			}

			assert (!children[0] && !children[1] && !division_new_vertex);

			if (aux_stat != d2::pixel(0, 0, 0)) 
				return 0;

			split();
			return 1;
		}

		int unsplit_on_aux() {
			if (aux_stat != d2::pixel(0, 0, 0)) {
				unsplit();
				return 1;
			} else if (children[0] && children[1]) {
				return children[0]->unsplit_on_aux()
				     | children[1]->unsplit_on_aux();
			}

			return 0;
		}

		ale_pos compute_area(pt _pt) {
			point a = _pt.scaled_transform(*vertices[0]);
			point b = _pt.scaled_transform(*vertices[1]);
			point c = _pt.scaled_transform(*vertices[2]);

			return 0.5 
			     * a.lengthto(b) 
			     * a.lengthto(c)
			     * sin(a.anglebetw(b, c));
		}
	};

	/*
	 * Use a pair of trees to store the triangles.
	 */
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

	/*
	 * Z-buffer initialization function.
	 */
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

		while (t) {

			while (t->division_new_vertex) {
				assert (t->children[0]);
				assert (t->children[1]);
				t = t->children[0];
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
			 * Expand the bounding box to the closest non-interior
			 * integral values.
			 */

			for (int d = 0; d < 2; d++) {
				max[d] = ceil(max[d]);
				min[d] = floor(min[d]);
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

			for (int i = (int) min[0]; i <= (int) max[0]; i++)
			for (int j = (int) min[1]; j <= (int) max[1]; j++) {

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

			while(t->parent && t == t->parent->children[1]) {
				t = t->parent;
			}

			if (t->parent == NULL)
				t = NULL;
			else if (t == t->parent->children[0])
				t = t->parent->children[1];
		}
	}

	/*
	 * Color a 3D scene using value averaging.  This is intended for
	 * producing low level-of-detail output.  Better results for higher
	 * resolution cases could be obtained by using an Irani-Peleg style
	 * backprojection approach.
	 */
	static void color_average() {

		/*
		 * Initialize color counters
		 */

		triangle_head[0]->init_color_counters();
		triangle_head[1]->init_color_counters();

		/*
		 * Iterate over all frames
		 */
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {

			d2::image *im = cl->reference[n];
			ale_pos sf = cl->sf;

			/*
			 * Z-buffer to map points to triangles
			 */
			pt _pt = align::projective(n);
			_pt.scale(sf / _pt.scale_2d());
			assert (im->width() == (unsigned int) floor(_pt.scaled_width()));
			assert (im->height() == (unsigned int) floor(_pt.scaled_height()));
			triangle **zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all points in the frame.
			 */
			for (unsigned int i = 0; i < im->height(); i++)
			for (unsigned int j = 0; j < im->width();  j++) {
				triangle *t = zbuf[i * im->width() + j];

				/*
				 * Check for points without associated triangles.
				 */
				if (!t)
					continue;

				/*
				 * Set new color and weight.
				 */
				t->color = (t->color * t->weight + im->get_pixel(i, j)) / (t->weight + d2::pixel(1, 1, 1));
				t->weight = t->weight + d2::pixel(1, 1, 1);
			}

			free(zbuf);
		}
	}

	/*
	 * Measure the error between reference images and the scene.
	 */
	static ale_accum scene_error() {

		ale_accum error = 0;
		ale_accum max_est = 0;

		/*
		 * Iterate over all frames
		 */
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {

			ale_pos sf = cl->sf;

			/*
			 * Z-buffer to map points to triangles
			 */
			pt _pt = align::projective(n);
			_pt.scale(sf / _pt.scale_2d());
			unsigned int height = (unsigned int) floor(_pt.scaled_height());
			unsigned int width  = (unsigned int) floor(_pt.scaled_width());
			triangle **zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all points in the frame.
			 */
			for (unsigned int i = 0; i < height; i++)
			for (unsigned int j = 0; j < width;  j++) {
				triangle *t = zbuf[i * width + j];

				/*
				 * Check for points without associated triangles.
				 */
				if (!t)
					continue;

				/*
				 * Calculate the difference.
				 */

				d2::pixel ref = cl->reference[n]->pix(i, j);
				d2::pixel scn = t->color;
				d2::pixel diff = ref - scn;

				for (int k = 0; k < 3; k++) {
					error += diff[k] * diff[k];

					if (ref[k] > scn[k])
						max_est += ref[k] * ref[k];
					else
						max_est += scn[k] * scn[k];
				}
			}

			free(zbuf);
		}

		return sqrt(error / max_est);
	}

	/*
	 * Test the density of the mesh for correct sampling in
	 * color_average(), and split (or unsplit) triangles if necessary,
	 * continuing until no more operations can be performed.  
	 */
	static int density_test(int split) {

		// triangle_head[0]->write_tree(1);
		// triangle_head[1]->write_tree(1);

		ale_pos scale = (split ? 1 : 2);
		d2::pixel init_value = d2::pixel(1, 1, 1) * (ale_real) (split ? 1 : 0);

		triangle_head[0]->init_aux_stats(init_value);
		triangle_head[1]->init_aux_stats(init_value);

		/*
		 * Iterate over all frames
		 */
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {

			ale_pos sf = cl->sf;

			/*
			 * Z-buffer to map points to triangles
			 */
			pt _pt = align::projective(n);
			_pt.scale(scale * sf / _pt.scale_2d());
			unsigned int height = (unsigned int) floor(_pt.scaled_height());
			unsigned int width  = (unsigned int) floor(_pt.scaled_width());
			triangle **zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all points in the frame.
			 */
			for (unsigned int i = 0; i < height; i++)
			for (unsigned int j = 0; j < width;  j++) {
				triangle *t = zbuf[i * width + j];

				/*
				 * Check for points without associated triangles.
				 */
				if (!t)
					continue;

				/*
				 * Check for triangles that have already been eliminated.
				 */
				if (t->aux_stat == d2::pixel(-1, -1, -1))
					continue;

				/*
				 * Mark the triangle as under consideration
				 */

				t->aux_stat = d2::pixel(0, 0, 0);

				/*
				 * Check that the triangle area is at least
				 * four.
				 */

				ale_pos area = t->compute_area(_pt);

				if (area < 4 && split)
					t->aux_stat = d2::pixel(-1, -1, -1);
				else if (area < 4 && !split && t->parent)
					t->parent->aux_stat = d2::pixel(-1, -1, -1);
			}

			free(zbuf);
		}

		if (split)
			return (triangle_head[0]->split_on_aux() | triangle_head[1]->split_on_aux());
		if (!split)
			return (triangle_head[0]->unsplit_on_aux() | triangle_head[1]->unsplit_on_aux());

		assert(0);

		return 0;
	}

	static void density_test_split() {
		while(density_test(1));
	}

	static void density_test_unsplit() {
		while(density_test(0));
	}

	/*
	 * Adjust vertices to minimize scene error.  Return 
	 * non-zero if improvements are made.
	 */
	static int adjust_vertices() {
		int improved = 0;

		triangle *t = triangle_head[0];

		while (t) {

			while (t->division_new_vertex) {
				assert (t->children[0]);
				assert (t->children[1]);
				t = t->children[0];
			}

			/*
			 * Determine the cross product of two triangle edges
			 */

			point n = t->vertices[0]->xproduct(*t->vertices[1], *t->vertices[2]);

			/*
			 * Determine the initial error
			 */

			ale_accum err = scene_error();

			/*
			 * Iterate over vertices, checking the effects of
			 * changes on the scene error.
			 */

			for (int v = 0; v < 3; v++) {
				*t->vertices[v] += n;

				if (!(scene_error() < err)) {
					*t->vertices[v] -= 2 * n;

					if (!(scene_error() < err))
						*t->vertices[v] += n;
					else 
						improved = 1;
				} else
					improved = 1;
			}

			while(t->parent && t == t->parent->children[1]) {
				t = t->parent;
			}

			if (t->parent == NULL && t == triangle_head[0])
				t = triangle_head[1];
			else if (t->parent == NULL)
				t = NULL;
			else if (t == t->parent->children[0])
				t = t->parent->children[1];
		}

		return improved;
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

		triangle_head[1]->vertices[0] = new point(max[0], max[1], 0);
		triangle_head[1]->vertices[1] = new point(min[0], max[1], 0);
		triangle_head[1]->vertices[2] = new point(max[0], min[1], 0);

		triangle_head[1]->neighbors[0] = triangle_head[0];

		color_average();
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

		pt _pt = align::projective(n);

		triangle **zbuf = init_zbuf(_pt);

		zbuffer(_pt, zbuf, triangle_head[0]);
		zbuffer(_pt, zbuf, triangle_head[1]);

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

		density_test_split();
		density_test_unsplit();
		color_average();

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
				density_test_split();
				density_test_unsplit();
				color_average();
				count = 0;
			}

			count++;
			improved = 0;

			/*
			 * Try improving the result by moving existing vertices.
			 */

			improved |= adjust_vertices();
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
