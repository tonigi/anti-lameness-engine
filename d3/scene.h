// Copyright 2003, 2004, 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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

/*
 * PLANAR_SUBDIVISION_COUNT must be exactly pow(4, PLANAR_SUBDIVISION_DEPTH)
 */

#define PLANAR_SUBDIVISION_DEPTH 1
#define PLANAR_SUBDIVISION_COUNT 4

class scene {

	/*
	 * Structure to hold a subdivisible triangle.
	 */

	struct triangle {

		/*
		 * Stats and auxiliaries for this element
		 */

		d2::pixel color[PLANAR_SUBDIVISION_COUNT];
		d2::pixel weight[PLANAR_SUBDIVISION_COUNT];
		d2::pixel aux_stat;
		void *aux_var;

		/*
		 * Connectivity
		 */

		point *vertices[3];
		char external_vertices[3];
		struct triangle *neighbors[3];
		struct triangle *parent;
		int division_vertex;
		point *division_new_vertex;
		struct triangle *children[2];

		/*
		 * Constructor
		 */

		triangle() {

			/*
			 * Ensure that the (COUNT, DEPTH) relationship holds for planar elements.
			 */

			assert ((int) PLANAR_SUBDIVISION_COUNT == (int) pow(4, PLANAR_SUBDIVISION_DEPTH));

			for (int sti = 0; sti < PLANAR_SUBDIVISION_COUNT; sti++) {
				color[sti] = d2::pixel(0, 0, 0);
				weight[sti] = d2::pixel(0, 0, 0);
			}

			aux_var = NULL;

			vertices[0] = NULL;
			vertices[1] = NULL;
			vertices[2] = NULL;

			external_vertices[0] = 1;
			external_vertices[1] = 1;
			external_vertices[2] = 1;

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

			for (int sti = 0; sti < PLANAR_SUBDIVISION_COUNT; sti++) {
				color[sti] = d2::pixel(0, 0, 0);
				weight[sti] = d2::pixel(0, 0, 0);
			}
		}

		void init_aux_stats(d2::pixel value = d2::pixel(0, 0, 0)) {
			if (children[0])
				children[0]->init_aux_stats(value);
			if (children[1])
				children[1]->init_aux_stats(value);

			aux_stat = value;
		}

		void free_aux_vars() {
			if (children[0])
				children[0]->free_aux_vars();
			if (children[1])
				children[1]->free_aux_vars();

			free(aux_var);
			aux_var = NULL;
		}
		

		/*
		 * Get the neighbor link from a given neighbor that references
		 * the 'this' object.
		 */

		int self_ref_from_neighbor(int n) const {
			triangle *t = neighbors[n];

			assert (t);

			for (int i = 0; i < 3; i++)
				if (t->neighbors[i] == this)
					return i;

			fprintf(stderr, "error in %p\n", this);
			assert(0);

			return -1;
		}

		/*
		 * Get a reference to a given vertex.
		 */

		int vertex_ref(point *p) const {
			for (int v = 0; v < 3; v++)
				if (vertices[v] == p)
					return v;

			fprintf(stderr, "error in %p\n", this);
			assert(0);

			return -1;
		}

		/*
		 * Get the angle at a given vertex
		 */

		ale_pos vertex_angle(int v) {
			return vertices[v]->anglebetw(*vertices[(v + 1) % 3],
						      *vertices[(v + 2) % 3]);
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
				children[c]->external_vertices[(division_vertex + 1 + c) % 3] 
					= neighbors[division_vertex] ? 0 : 1;
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
		 * Split a number of triangles so as to establish a lower bound
		 * on the size of newly-created angles.
		 */

		void split() {

			int v;
			double angle_lbound = 59 * M_PI / 180;

			for (v = 0; v < 3; v++)
				if (vertex_angle(v) >= angle_lbound)
					break;

			assert (v < 3);
			assert (v >= 0);

			if (!neighbors[v]) {
				split(v);
				return;
			}

			while (neighbors[v]->vertex_angle(self_ref_from_neighbor(v)) 
			     < angle_lbound) {

				neighbors[v]->split();

				/*
				 * Check to see whether this triangle has been split as a
				 * consequence.
				 */

				if (division_new_vertex)
					return;

			}

			split(v);
		}

#if 0
		/*
		 * Split a triangle at a random vertex among those with largest
		 * angle (averaged with any neighbor's opposite angle).
		 */
		void split() {
			
			ale_pos angles[3]; 

			for (int v = 0; v < 3; v++) {
				angles[v] = vertices[v]->anglebetw(*vertices[(v + 1) % 3],
						                   *vertices[(v + 2) % 3]);
				if (neighbors[v]) {
					int self_ref = self_ref_from_neighbor(v);
					point **nvs = neighbors[v]->vertices;
					angles[v] = (angles[v]
					           + nvs[self_ref]->anglebetw(*nvs[(self_ref + 1) % 3],
							                      *nvs[(self_ref + 2) % 3])) / 2;
				}
			}
				

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
#endif

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

			a[2] = 0;
			b[2] = 0;
			c[2] = 0;

			return 0.5 
			     * a.lengthto(b) 
			     * a.lengthto(c)
			     * sin(a.anglebetw(b, c));
		}

		point normal() {
			point unscaled_normal = vertices[0]->xproduct(*vertices[1], *vertices[2]);

			return (unscaled_normal / fabs(unscaled_normal.norm()));
		}

		point centroid() {
			return (*vertices[0] + *vertices[1] + *vertices[2]) / 3;
		}

		/*
		 * Return the maximum angle formed with a neighbor triangle.
		 * Searches two neighbors deep.
		 */
		ale_pos max_neighbor_angle_2(ale_pos displacement = 0) {

			ale_pos max_angle = 0;

			point _normal = normal();

			for (int n = 0; n < 3; n++) {
				if (neighbors[n] == NULL)
					continue;

				int srfn = self_ref_from_neighbor(n);

				point a = *neighbors[n]->vertices[srfn];
				point b = *neighbors[n]->vertices[(srfn + 1) % 3] + displacement * _normal;
				point c = *neighbors[n]->vertices[(srfn + 2) % 3] + displacement * _normal;

				point unscaled_dnn = a.xproduct(b, c);

				ale_pos angle = point(0, 0, 0).anglebetw(unscaled_dnn, _normal);

				if (angle > max_angle)
					max_angle = angle;

				/*
				 * Second level of neighbors
				 */

				for (int m = 0; m < 2; m++) {
					if (neighbors[n]->neighbors[(m + srfn + 1) % 3] == NULL)
						continue;
					int srfn2 = neighbors[n]->self_ref_from_neighbor((m + srfn + 1) % 3);

					point a = *neighbors[n]->neighbors[m]->vertices[srfn2];
					point b = *neighbors[n]->neighbors[m]->vertices[(srfn2 + 1) % 3] + (1 - m) * displacement * _normal;
					point c = *neighbors[n]->neighbors[m]->vertices[(srfn2 + 2) % 3] + m * displacement * _normal;;

					point unscaled_dnn2 = a.xproduct(b, c);

					ale_pos angle = point(0, 0, 0).anglebetw(unscaled_dnn, unscaled_dnn2);

					if (angle > max_angle)
						max_angle = angle;
				}
			}

			return max_angle;
		}

		/*
		 * Adjust vertices according to mapped frames.
		 *
		 * Return non-zero if an adjustment is made.
		 */
		int adjust_vertices() {
			if (children[0] && children[1]) {
				return children[0]->adjust_vertices()
				     | children[1]->adjust_vertices();
			}

			/*
			 * Don't allow triangles to move if they contain
			 * external points.
			 */

			if (external_vertices[0]
			 || external_vertices[1]
			 || external_vertices[2])
				return 0;

			/*
			 * Get the list of frames in which this triangle appears.  If
			 * it doesn't appear in any frames, then we don't adjust any
			 * vertices.
			 */

			char *frame_list = (char *)aux_var;

			if (!frame_list)
				return 0;

			/*
			 * As a basis for determining the step size for
			 * adjustment, determine the average distance between
			 * vertices and divide it by a constant.
			 */

			ale_accum step = (vertices[0]->lengthto(*vertices[1])
				        + vertices[1]->lengthto(*vertices[2])
					+ vertices[2]->lengthto(*vertices[0])) / 21;

			/*
			 * There are three possibilities for each triangle visited:
			 *
			 * Vertices are moved a step in the direction of the normal
			 * Vertices remain in place
			 * Vertices are moved a step against the direction of the normal.
			 *
			 * Try each possibility, and determine the sum-of-squares error
			 * over all frames in which the triangle appears, mapping the triangle
			 * centroid and using bilinear interpolation between reference frame 
			 * pixels to determine error.  The position with least sum-of-squares
			 * error is selected.
			 */

			int best = 0;
			ale_accum lowest_error = +0;
			lowest_error = +1 / lowest_error;
			ale_accum divisors[3];

			assert (lowest_error > 0);
			assert (isinf(lowest_error) == 1);

			for (int dir = 1; dir >= -1; dir--) {

				ale_accum error = 0;
				ale_accum divisor = 0;

				point adjusted_centroid = centroid() + normal() * (step * dir);

				/*
				 * Eliminate from consideration any change that increases to more than
				 * a given amount the angle between the normals of adjacent triangles.
				 */

//				if (max_neighbor_angle_2(step * dir) > 0)
//					continue;
				if (max_neighbor_angle_2(step * dir) > M_PI / 4)
					continue;

				/*
				 * Iterate through all frames.
				 */

				for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
					if (!frame_list[n])
						continue;
					
					pt _pt = align::projective(n);
					_pt.scale(cl->sf / _pt.scale_2d());

					point mapped_centroid = _pt.scaled_transform(adjusted_centroid);

					if (!cl->reference[n]->in_bounds(mapped_centroid.xy()))
						continue;

					int sti = subtriangle_index(_pt, mapped_centroid, this, normal() * step * dir);
					d2::pixel c1 = cl->reference[n]->get_bl(mapped_centroid.xy());
					d2::pixel c2 = color[sti];

					error += (c1 - c2).normsq();

					if (c1.normsq() > c2.normsq()) 
						divisor += c1.normsq();
					else
						divisor += c2.normsq();
				}

				error /= divisor;

				if (error < lowest_error) {
					lowest_error = error;
					best = dir;
				}

				divisors[dir + 1] = divisor;
			}

			/*
			 * Don't move triangles out of the view pyramid of any
			 * frame.
			 */
			if (divisors[best + 1] < divisors[1])
				best = 0;

			for (int v = 0; v < 3; v++)
				(*vertices[v]) += normal() * (step * best);

			if (best != 0)
				return 1;

			return 0;
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
	 * Upper/lower test for interiority.
	 *
	 * P is a set of triangle vertex points mapped into image space.
	 * TEST_POINT is a test point at distance 1 from a camera.
	 */
	 static int is_interior(point p[3], point test_point, int print_on_failure = 0) {

		int lower[2] = {0, 0};
		int upper[2] = {0, 0};

		for (int v = 0; v < 3; v++) {
			point cv = p[v];
			point nv = p[(v + 1) % 3];

			for (int d = 0; d < 2; d++)
			if ((test_point[d] - cv[d]) * (test_point[d] - nv[d]) <= 0
			 && nv[d] - cv[d] != 0) {

				int e = (d + 1) % 2;
				ale_pos travel = (test_point[d] - cv[d]) / (nv[d] - cv[d]);
				ale_pos intersect = cv[e] + travel * (nv[e] - cv[e]);
				if (intersect <= test_point[e]) 
					lower[e] = 1;
				if (intersect >= test_point[e])
					upper[e] = 1;
			}
		}

		if (!lower[0] || !upper[0] || !lower[1] || !upper[1]) {
			if (print_on_failure && 0) {
				 fprintf(stderr, "is_interior({{%f, %f}, {%f, %f}, {%f, %f}}, {%f, %f})",
						 p[0][0], p[0][1],
						 p[1][0], p[1][1],
						 p[2][0], p[2][1],
						 test_point[0], test_point[1]);
				fprintf(stderr, " = 0\n");
			}

			return 0;
		}

		return 1;
	}

	/*
	 * Return a subtriangle index for a given position.
	 */

	static int subtriangle_index(pt _pt, d2::point test_point, point p_w[3], int depth, int prefix) {
		point p_i[3];

		for (int v = 0; v < 3; v++)
			p_i[v] = _pt.scaled_transform(p_w[v]);

		if (!is_interior(p_i, point(test_point[0], test_point[1], 1), 1)) {
			/*
			 * XXX: Error condition ... what should we do here?
			 */
			// assert(0);
			return prefix;
		}

		if (depth == 0)
			return prefix;

		point hp_w[3];
		point hp_i[3];

		for (int v = 0; v < 3; v++)
			hp_w[v] = (p_w[(v + 1) % 3] + p_w[(v + 2) % 3]) / 2;

		for (int v = 0; v < 3; v++)
			hp_i[v] = _pt.scaled_transform(hp_w[v]);

		for (int v = 0; v < 3; v++) {
			point arg_array[3] = {p_i[v], hp_i[(v + 2) % 3], hp_i[(v + 1) % 3]};

			if (is_interior(arg_array, point(test_point[0], test_point[1], 1)))
				return subtriangle_index(_pt, test_point, arg_array, depth - 1, prefix * 4 + v);
		}

		return subtriangle_index(_pt, test_point, hp_w, depth - 1, prefix * 4 + 3);
	}

	static int subtriangle_index(pt _pt, d2::point test_point, triangle *t, point triangle_offset = point(0, 0, 0)) {
		point p[3];

		for (int v = 0; v < 3; v++)
			p[v] = *t->vertices[v] + triangle_offset;

		return subtriangle_index(_pt, test_point, p, PLANAR_SUBDIVISION_DEPTH, 0);
	}

	static int subtriangle_index(pt _pt, point test_point, triangle *t, point triangle_offset = point(0, 0, 0)) {
		return subtriangle_index(_pt, test_point.xy(), t);
	}

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

				point test_point = point(i, j, 1);

				if (!is_interior(p, test_point))
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
				int sti = subtriangle_index(_pt, d2::point(i, j), t);
//				t->color[sti] = (t->color[sti] * t->weight[sti] + im->get_pixel(i, j)) 
//					      / (t->weight[sti] + d2::pixel(1, 1, 1));
//				t->weight[sti] = t->weight[sti] + d2::pixel(1, 1, 1);
				t->color[sti][n % 3] = (t->color[sti][n % 3] * t->weight[sti][n % 3] + im->get_pixel(i, j)[n % 3]) 
					      / (t->weight[sti][n % 3] + 1);
				t->weight[sti][n % 3] = t->weight[sti][n % 3] + 1;
			}

			free(zbuf);
		}
	}

	/*
	 * Test the density of the mesh for correct sampling in
	 * color_average(), and split (or unsplit) triangles if necessary,
	 * continuing until no more operations can be performed.  
	 */
	static int density_test(int split) {

		ale_pos scale = (split ? 1 : 1);
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

				if (area <= 4 && split)
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

	static int density_test_split() {
		return density_test(1);
	}

	static int density_test_unsplit() {
		return density_test(0);
	}

	/*
	 * Determine triangle visibility.
	 */
	static void determine_visibility() {
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
			 * Iterate over all points in the frame, adding this frame
			 * to the list of frames including the associated triangle.
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
				 * Add this frame to the list of frames in which this
				 * triangle is visible, or start a new list if none exists
				 * yet.
				 */


				if (t->aux_var == NULL)
					t->aux_var = calloc(d2::image_rw::count(), sizeof(char));

				char *aux_var = (char *) t->aux_var;

				assert (aux_var);

				aux_var[n] = 1;
			}

			free(zbuf);
		}
	}

	/*
	 * Adjust vertices to minimize scene error.  Return 
	 * non-zero if improvements are made.
	 */
	static int adjust_vertices() {
		/*
		 * Determine the visibility of triangles from frames.
		 */

		determine_visibility();

		int result = triangle_head[0]->adjust_vertices() | triangle_head[1]->adjust_vertices();

		triangle_head[0]->free_aux_vars();
		triangle_head[1]->free_aux_vars();

		return result;
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
		triangle_head[1]->vertices[1] = triangle_head[0]->vertices[2];
		triangle_head[1]->vertices[2] = triangle_head[0]->vertices[1];

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
		if (zbuf[i * im->width() + j]) {
			triangle *t = zbuf[i * im->width() + j];
			int sti = subtriangle_index(_pt, d2::point(i, j), t);
			im->pix(i, j) = t->color[sti];
		}

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

		while ((improved && count < 10) || cl->next) {

			fprintf(stderr, "unsplitting ...\n");

			if (density_test_unsplit())
				continue;

			fprintf(stderr, "splitting/unsplitting ...\n");

			if (density_test_split() && !density_test_unsplit())
				continue;

			fprintf(stderr, "color averaging...\n");

			color_average();

			/*
			 * Write output incrementally, if desired.
			 */

			if (inc_bit && (!improved || !(rand() % 2)))
			for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
				fprintf(stderr, "writing image...\n");
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
			
//			triangle_head[0]->write_tree(1);
//			triangle_head[1]->write_tree(1);
			
			/*
			 * Increase LOD if no improvements were achieved in the
			 * most recent pass at the previous LOD.
			 */

			if (!improved || count > 40 ) {
//				triangle_head[0]->write_tree(1);
//				triangle_head[1]->write_tree(1);
				fprintf(stderr, ".");
				assert (cl->next);
				fprintf(stderr, "increasing LOD ...\n");
				increase_lod();
				count = 0;
				improved = 1;
				continue;
			}

			count++;
			improved = 0;

			/*
			 * Try improving the result by moving existing vertices.
			 */

			fprintf(stderr, "adjusting vertices [disabled] ...\n");

			improved |= adjust_vertices();
		}

		fprintf(stderr, "color averaging...\n");

		color_average();

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
