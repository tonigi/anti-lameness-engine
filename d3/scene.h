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

#define PLANAR_SUBDIVISION_DEPTH 2
#define PLANAR_SUBDIVISION_COUNT 16

class scene {

	/*
	 * Ray-Triangle intersection utility function
	 *
	 * Variables should be specified in cartesian coordinates (not
	 * projective coordinates).
	 *
	 * Return value elements are:
	 *
	 * 	0: ray multiplier
	 * 	1: v[1] - v[0] component
	 * 	2: v[2] - v[0] component
	 */
	static point rt_intersect(point r, point vertices[3]) {
		point a = vertices[0];
		point b = vertices[1];
		point c = vertices[2];
		point d = a - b;
		point e = a - c;

		/*
		ale_pos m[3][3] = {
			{ r[0], r[1], r[2] },
			{ d[0], d[1], d[2] },
			{ e[0], e[1], e[2] }
		};
		*/

		ale_pos m_det = r[0] * d[1] * e[2] 
			      + r[1] * d[2] * e[0]
			      + r[2] * d[0] * e[1]
			      - e[0] * d[1] * r[2]
			      - d[0] * r[1] * e[2]
			      - r[0] * e[1] * d[2];

		ale_pos m_inverse_t[3][3] = {
			{
				(d[1] * e[2] - d[2] * e[1]) / m_det,
				(d[2] * e[0] - d[0] * e[2]) / m_det,
				(d[0] * e[1] - d[1] * e[0]) / m_det
			}, {
				(e[1] * r[2] - e[2] * r[1]) / m_det,
				(e[2] * r[0] - e[0] * r[2]) / m_det,
				(e[0] * r[1] - e[1] * r[0]) / m_det,
			}, {
				(r[1] * d[2] - r[2] * d[1]) / m_det,
				(r[2] * d[0] - r[0] * d[2]) / m_det,
				(r[0] * d[1] - r[1] * d[0]) / m_det,
			}
		};

		point k(a[0] * m_inverse_t[0][0] + a[1] * m_inverse_t[0][1] + a[2] * m_inverse_t[0][2],
		        a[0] * m_inverse_t[1][0] + a[1] * m_inverse_t[1][1] + a[2] * m_inverse_t[1][2],
			a[0] * m_inverse_t[2][0] + a[1] * m_inverse_t[2][1] + a[2] * m_inverse_t[2][2]);

		return k;
	}

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
		 * Color for traversal.
		 */

		char traversal_color;

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

			traversal_color = 0;
		}

		/*
		 * Clear a traversed connected sub-graph to zero.
		 */
		void traversal_clear() {
			if (traversal_color == 0)
				return;

			traversal_color = 0;

			for (int v = 0; v < 3; v++)
			if  (neighbors[v])
				neighbors[v]->traversal_clear();
		}

		/*
		 * Traverse the sub-graph of triangles around a given vertex,
		 * and return an aggregate result for a given function
		 *
		 * TYPE is one of:
		 *
		 * 	0 returns the sum
		 * 	1 returns the maximum
		 *
		 * STARTING_NODE should be left at 1 when this function is
		 * called externally.
		 */
		ale_accum traverse_around_vertex(point *vertex, ale_accum (triangle::*fp)(), int type, 
				int starting_node = 1) {

			/*
			 * Check whether we've already visited this node.
			 */

			if (traversal_color == 1) {
				return 0;
			}

			/*
			 * Mark that we've visited this node.
			 */

			traversal_color = 1;


			/*
			 * Perform the aggregation operation.
			 */

			ale_accum aggregate = (this->*fp)();

			for (int v = 0; v < 3; v++)
			if  (neighbors[v]) {
				ale_accum new_value = neighbors[v]->traverse_sum_around_vertex(vertex, fp, 0);
				switch (type) {
				case 0:
					/*
					 * Sum
					 */
					aggregate += new_value;
					break;
				case 1:
					/*
					 * Maximum
					 */
					if (new_value > aggregate)
						aggregate = new_value;
					break;
				default:
					assert(0);
				}
			}

			/*
			 * Clean up traversal colors.
			 */

			if (starting_node)
				traversal_clear();

			/*
			 * Return the aggregate.
			 */

			return aggregate;
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
		 * Calculate area.
		 */

		ale_accum area() {
			return 0.5 * vertices[0]->xproduct(*vertices[1], *vertices[2]).norm();
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

		int vertex_ref_maybe(point *p) const {
			for (int v = 0; v < 3; v++)
				if (vertices[v] == p)
					return v;

			return -1;
		}

		/*
		 * Get a reference to a given vertex.
		 */

		int vertex_ref(point *p) const {
			int maybe = vertex_ref_maybe(p);

			assert (maybe >= 0);

			return maybe;
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
			point a = _pt.wp_scaled(*vertices[0]);
			point b = _pt.wp_scaled(*vertices[1]);
			point c = _pt.wp_scaled(*vertices[2]);

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

		point element_centroid(int e) {
			point result = centroid();

			/*
			 * Adjust the position, starting with the least
			 * significant digits.
			 */
			ale_pos multiplier = pow(2, -PLANAR_SUBDIVISION_DEPTH);
			for (int d = 0; d < PLANAR_SUBDIVISION_DEPTH; d++) {
				int digit = e % 4;

				for (int v = 0; v < 3; v++)
				if  (digit == v) {
					point offset = *vertices[v] - centroid();
					result += offset * multiplier;
				}

				if (digit == 3)
					result = 2 * centroid() - result;

				multiplier *= 2;
				e /= 4;
			}

			return result;
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
		 * Calculate the error between the model and the reference
		 * images.
		 */
		ale_accum reference_error() {
			assert(!children[0] && !children[1]);

			ale_accum error = 0;
			ale_accum divisor = 0;

			char *frame_list = (char *)aux_var;

			if (!frame_list)
				return 0;

			/*
			 * Iterate over all elements
			 */

			for (unsigned int e = 0; e < PLANAR_SUBDIVISION_COUNT; e++) {

				/*
				 * For all input frames from which this
				 * triangle is visible, determine the color at
				 * the projected element centroid, and accumulate
				 * the difference.
				 */

				point centroid = element_centroid(e);

				for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
					if (!frame_list[n])
						continue;

					pt _pt = align::projective(n);
					_pt.scale(cl->sf / _pt.scale_2d());

					/*
					 * Map the centroid into image space.
					 */

					point mapped_centroid = _pt.wp_scaled(centroid);

					/*
					 * Check the bounds
					 */

					if (!cl->reference[n]->in_bounds(mapped_centroid.xy()))
						continue;

					d2::pixel ca = color[e];
					d2::pixel cb = cl->reference[n]->get_bl(mapped_centroid.xy());

					for (int k = 0; k < 3; k++) {
						if (!finite(ca[k]) || !finite(cb[k]))
							continue;
						error += pow(ca[k] - cb[k], 2);
						divisor += pow(ca[k] > cb[k] ? ca[k] : cb[k], 2);
					}
				}
			}

			error /= divisor;

			return error;

		}

		/*
		 * Recolor, assuming that visibility remains constant.
		 */
		void recolor() {
			// assert(!children[0] && !children[1]);
			if (children[0])
				children[0]->recolor();
			if (children[1])
				children[1]->recolor();

			if (children[0] || children[1])
				return;

			init_color_counters();

			char *frame_list = (char *)aux_var;

			if (!frame_list)
				return;

			/*
			 * Iterate over all elements
			 */

			for (unsigned int e = 0; e < PLANAR_SUBDIVISION_COUNT; e++) {
				d2::pixel color = d2::pixel(0, 0, 0);
				d2::pixel weight = d2::pixel(0, 0, 0);

				/*
				 * For all input frames from which this
				 * triangle is visible, determine the color at
				 * the projected element centroid.  Average all
				 * such colors.
				 */

				point centroid = element_centroid(e);

				for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
					if (!frame_list[n])
						continue;

					pt _pt = align::projective(n);
					_pt.scale(cl->sf / _pt.scale_2d());

					/*
					 * Map the centroid into image space.
					 */

					point mapped_centroid = _pt.wp_scaled(centroid);

					/*
					 * Check the bounds
					 */

					if (!cl->reference[n]->in_bounds(mapped_centroid.xy()))
						continue;

					color += cl->reference[n]->get_bl(mapped_centroid.xy());
					weight += d2::pixel(1, 1, 1);
					// color = weight * (mapped_centroid.xy()[1] / cl->reference[n]->width());

					point centroid_local = _pt.wc(centroid);
					point _vertices[3] = {_pt.wc(*vertices[0]), _pt.wc(*vertices[1]), _pt.wc(*vertices[2])};
					assert (1 || e == subtriangle_index(centroid_local, _vertices, PLANAR_SUBDIVISION_DEPTH, 0));
				}

				color /= weight;

				this->color[e] = color;
				this->weight[e] = weight;
			}
		}

		/*
		 * Color all neighbors (containing at least one of the given
		 * set of vertices) with a negative color.
		 */
		void color_neighbors_negative(point *vertices[3] = NULL, ale_real neg_color = 0) {

			if (vertices == NULL)
				vertices = this->vertices;

			while (neg_color == 0)
				neg_color = -((ale_real) (rand() % 1000));

			if (color[0][0] == neg_color)
				return;

			if (!vertex_ref_maybe(vertices[0])
			 && !vertex_ref_maybe(vertices[1])
			 && !vertex_ref_maybe(vertices[2]))
				return;

			color[0][0] = neg_color;

			for (int n = 0; n < 3; n++) {
				if (neighbors[n])
					neighbors[n]->color_neighbors_negative(vertices, neg_color);
			}
		}

		/*
		 * Accumulate error from neighbors identified with negative coloration
		 */
		ale_accum accumulate_neighbor_error() {

			if (!(color[0][0] < 0))
				return 0;

			recolor();

			assert (!(color[0][0] < 0));


			/*
			 * XXX: This approach to counting could be improved.
			 * Summation of weights is probably the best approach.
			 */

			ale_accum error = reference_error();
			unsigned int error_count = 0;

			if (finite(error))
				error_count = 1;
			else
				error = 0;

			for (int n = 0; n < 3; n++)
			if  (neighbors[n]) {
				ale_accum neighbor_error = neighbors[n]->accumulate_neighbor_error();
				if (finite(neighbor_error)) {
					error += neighbors[n]->accumulate_neighbor_error();
					error_count++;
				}
			}

			return (ale_accum) (error / error_count);
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
//					+ vertices[2]->lengthto(*vertices[0])) / 7;

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

			int improved = 0;
			ale_accum lowest_error = +0;
			lowest_error = +1 / lowest_error;

			assert (lowest_error > 0);
			assert (isinf(lowest_error) == 1);

			point best_vertices[3] = {*vertices[0], *vertices[1], *vertices[2]};

			/*
			 * Evaluate the error at the current position.
			 */

			color_neighbors_negative();
			lowest_error = accumulate_neighbor_error();

			/*
			 * Test modifications to the current position.
			 */

			for (int v = 0; v < 3; v++)
			for (int axis = 0; axis < 3; axis++)
			for (int dir = -1; dir <= 1; dir += 2) {

				/*
				 * Adjust the vertex under consideration
				 */

				*vertices[v] = best_vertices[v] + point::unit(axis) * step * dir;

				/*
				 * Eliminate from consideration any change that increases to more than
				 * a given amount the angle between the normals of adjacent triangles.
				 */

				if (max_neighbor_angle_2() > M_PI / 3) {
					*vertices[v] = best_vertices[v];
					continue;
				}

				/*
				 * Recalculate the color
				 */

				recolor();

				/*
				 * Check the error
				 */

				color_neighbors_negative();
				ale_accum error = accumulate_neighbor_error();
				if (error < lowest_error) {
					lowest_error = error;
					improved = 1;
					best_vertices[v] = *vertices[v];
					break;
				} else {
					*vertices[v] = best_vertices[v];
				}

				if (!finite(error))
					fprintf(stderr, "dir %d, error %e\n", dir, error);
			}

			return improved;
		}
	};

	/*
	 * Use a pair of trees to store the triangles.
	 */
	static struct triangle *triangle_head[2];

	/*
	 * Vector test for interiority in local cartesian space.
	 *
	 * P is a set of triangle vertex points.  R is a ray endpoint.  The
	 * intersection of R and the plane defined by P is the point being
	 * tested for interiority.
	 */
	 static int is_interior_c(point p_c[3], point r_c, int print_on_failure = 0) {

		point multipliers = rt_intersect(r_c, p_c);

		if (multipliers[1] >= 0
		 && multipliers[2] >= 0
		 && (multipliers[1] + multipliers[2] <= 1))
			return 1;

		if (print_on_failure && 0) {
			 fprintf(stderr, "is_interior_c({{%f, %f}, {%f, %f}, {%f, %f}}, {%f, %f})",
					 p_c[0][0], p_c[0][1],
					 p_c[1][0], p_c[1][1],
					 p_c[2][0], p_c[2][1],
					 r_c[0], r_c[1]);
			fprintf(stderr, " = 0\n");
		}

		return 0;
	}

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

			for (int d = 0; d < 2; d++) {
				int e = (d + 1) % 2;

				if ((test_point[d] - cv[d]) * (test_point[d] - nv[d]) <= 0
				 && nv[d] - cv[d] != 0) {
					ale_pos travel = (test_point[d] - cv[d]) / (nv[d] - cv[d]);
					ale_pos intersect = cv[e] + travel * (nv[e] - cv[e]);
					if (intersect <= test_point[e]) 
						lower[e] = 1;
					if (intersect >= test_point[e])
						upper[e] = 1;
				}
				if (nv[d] - cv[d] == 0 && test_point[d] == nv[d]) {
					lower[d] = 1;
					upper[d] = 1;
				}

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
	 * Return a subtriangle index for a given position.  Assumes local
	 * cartesian (not projective) coordinates.
	 *
	 * R_C is the endpoint of a ray whose tail is at zero.
	 * P_C are the vertex coordinates.
	 */
	static unsigned int subtriangle_index(point r_c, point p_c[3], int depth, int prefix) {

		if (depth == 0)
			return prefix;

		if (!is_interior_c(p_c, r_c, 1)) {
			/*
			 * XXX: Error condition ... what should we do here?
			 */
			// assert(0);
			return prefix;
		}

		point hp_c[3];

		for (int v = 0; v < 3; v++)
			hp_c[v] = (p_c[(v + 1) % 3] + p_c[(v + 2) % 3]) / 2;

		for (int v = 0; v < 3; v++) {
			point arg_array_c[3];

			arg_array_c[v] = p_c[v];
			arg_array_c[(v + 1) % 3] = hp_c[(v + 2) % 3];
			arg_array_c[(v + 2) % 3] = hp_c[(v + 1) % 3];

			int interior_t = is_interior_c(arg_array_c, r_c);

			if (interior_t)
				return subtriangle_index(r_c, arg_array_c, depth - 1, prefix * 4 + v);
		}

		return subtriangle_index(r_c, hp_c, depth - 1, prefix * 4 + 3);
	}

	/*
	 * Return a subtriangle index for a given position.
	 */

	static int subtriangle_index(pt _pt, d2::point test_point, point p_w[3], int depth, int prefix) {
		point p_c[3];

		if (depth == 0)
			return prefix;

		for (int v = 0; v < 3; v++)
			p_c[v] = _pt.wc(p_w[v]);

		point r_c = _pt.pc_scaled(point(test_point[0], test_point[1], 1));

		return subtriangle_index(r_c, p_c, depth, prefix);

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
				p[v] = _pt.wp_scaled(*t->vertices[v]);

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

				if (*zbuf_tri && _pt.wp_scaled(*(*zbuf_tri)->vertices[0])[2] > p[0][2])
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

				if (area / PLANAR_SUBDIVISION_COUNT <= 4 && split)
					t->aux_stat = d2::pixel(-1, -1, -1);
				else if (area / PLANAR_SUBDIVISION_COUNT < 4 && !split && t->parent)
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
//			im->pix(i, j) = d2::pixel(sti,sti,sti) / (double) PLANAR_SUBDIVISION_COUNT;
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
			im->pix(i, j) = d2::pixel(1, 1, 1) * _pt.wc(zbuf[i * im->width() + j]->centroid())[2];

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

		for (;;) {

			// fprintf(stderr, "unsplitting ...\n");

//			if (density_test_unsplit())
//				continue;

			// fprintf(stderr, "splitting/unsplitting ...\n");

			if (density_test_split())
				continue;

			// fprintf(stderr, "color averaging...\n");

			color_average();
			
			// fprintf(stderr, "using recolor() to recolor elements...\n");

			determine_visibility();

			triangle_head[0]->recolor();
			triangle_head[1]->recolor();

			triangle_head[0]->free_aux_vars();
			triangle_head[1]->free_aux_vars();

			/*
			 * Write output incrementally, if desired.
			 */

			if (inc_bit && (!improved || !(rand() % 2)))
			for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
				// fprintf(stderr, "writing image...\n");
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

			if (!improved || count > 40) {
//				triangle_head[0]->write_tree(1);
//				triangle_head[1]->write_tree(1);
				fprintf(stderr, ".");

				if (!cl->next)
					break;
				
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

			// fprintf(stderr, "adjusting vertices ...\n");

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
