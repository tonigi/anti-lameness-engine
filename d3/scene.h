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
 * Variables PLANAR_SUBDIVISION_* impose a lower bound on the size of clusters
 * of co-planar equal area triangles.  COUNT is the lower bound on the number
 * of such triangles in a cluster, and DEPTH is the base 4 logarithm of COUNT.
 *
 * NB: PLANAR_SUBDIVISION_COUNT must be exactly pow(4, PLANAR_SUBDIVISION_DEPTH),
 * and both numbers must be integers.
 */

#define PLANAR_SUBDIVISION_COUNT 4
#define PLANAR_SUBDIVISION_DEPTH 1

class scene {

	struct triangle;

	/*
	 * Multiplier used in calculating the edge-length contribution to 
	 * model cost.
	 */
	static ale_pos edge_cost_multiplier;

	/*
	 * Multiplier used in calculating the maximum angle contribution to
	 * model cost.
	 */
	static ale_pos angle_cost_multiplier;

	/*
	 * Vertex Auxiliary Structure
	 *
	 * This can be used for storing miscellaneous information about
	 * vertices.
	 */
	struct vertex_aux {
		point last_position;
		triangle *last_triangle;
		unsigned long update_id;
		ale_pos max_angle;
	};

	/*
	 * Vertex Auxiliary Map.
	 */
	static std::map<point *, vertex_aux> vam;

	/*
	 * Ray-Triangle intersection utility function
	 *
	 * Variables should be specified in cartesian coordinates (not
	 * projective coordinates).
	 *
	 * Return value elements are:
	 *
	 * 	0: v[1] - v[0] component
	 * 	1: v[2] - v[0] component
	 * 	2: ray multiplier
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

		point k(a[0] * m_inverse_t[1][0] + a[1] * m_inverse_t[1][1] + a[2] * m_inverse_t[1][2],
			a[0] * m_inverse_t[2][0] + a[1] * m_inverse_t[2][1] + a[2] * m_inverse_t[2][2],
			a[0] * m_inverse_t[0][0] + a[1] * m_inverse_t[0][1] + a[2] * m_inverse_t[0][2]);

		return k;
	}

	/*
	 * Z-buffer element.  Since we store all triangles encountered (in order
	 * to more easily accommodate scene modifications), we make the element
	 * array dynamically resizable and re-sortable.
	 */
	class zbuf_elem {
		std::set<triangle *> *tset;
		triangle *_nearest;

	public:
		zbuf_elem() {
			tset = new std::set<triangle *>;
			assert (tset);
			_nearest = NULL;
		}

		~zbuf_elem() {
			delete tset;
		}

		void clear_nearest() {
			_nearest = NULL;
		}

		void find_nearest(pt _pt, int i, int j);

		triangle *nearest(pt _pt, int i, int j) {
			find_nearest(_pt, i, j);
			return _nearest;
		}

		void insert(triangle *t) {
			tset->insert(t);
		}

		void insert(std::set<triangle *>::iterator begin, std::set<triangle *>::iterator end) {
			tset->insert(begin, end);
		}
	};

	/*
	 * Frame-to-frame mapping function
	 *
	 * Maps from point P in frame F1 to frame F2.  Returns projective
	 * coordinates.  Returns zero-depth if there is no mutually-visible
	 * model point.
	 */
	static point frame_to_frame(d2::point p, const pt &_pt1, const pt &_pt2, zbuf_elem *z1, zbuf_elem *z2);

	/*
	 * Vertex movement cost
	 *
	 * Evaluate the cost of moving a vertex.  Negative values indicate
	 * improvement; non-negative values indicate lack of improvement.
	 */
	static ale_accum vertex_movement_cost(triangle *t, point *vertex, point new_position, zbuf_elem **z);

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
		char vertex_fixed[3];
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

			vertex_fixed[0] = 1;
			vertex_fixed[1] = 1;
			vertex_fixed[2] = 1;

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

		ale_accum write_tree_wrapper() {
			write_tree();
			return 0;
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
		 * Traverse a tree of triangles and return an aggregate result
		 * for a given function, ignoring non-finite values.
		 *
		 * TYPE is one of:
		 *
		 * 	0 returns the sum
		 * 	1 returns the maximum
		 * 	2 returns the minimum
		 */
		ale_accum traverse_tree(ale_accum (triangle::*fp)(), int type) {

			if (division_new_vertex != NULL) {
				ale_accum c0 = children[0]->traverse_tree(fp, type);
				ale_accum c1 = children[1]->traverse_tree(fp, type);

				if (type == 0)
					return c0 + c1;
				if (type == 1)
					return (c0 > c1) ? c0 : c1;
				if (type == 2)
					return (c0 < c1) ? c0 : c1;

				assert(0);
				return 0;
			}

			return (this->*fp)();

		}

		/*
		 * Traverse the sub-graph of triangles around a given vertex,
		 * and construct a set containing all triangles encountered.
		 *
		 * PARTIAL should be initialized to the empty set by the
		 * non-recursive call, and will reflect the complete subgraph
		 * upon return from this call.
		 */
		void triangles_around_vertex(point *vertex, std::set<triangle *> *partial) {

			assert (partial);

			/*
			 * If this node does not contain the given vertex, or if we've already
			 * visited this node, then return immediately.
			 */

			if (vertex_ref_maybe(vertex) < 0 || partial->count(this) > 0)
				return;

			/*
			 * Add this triangle to the list.
			 */

			partial->insert(this);

			/*
			 * Visit neighbors
			 */

			for (int v = 0; v < 3; v++)
			if  (neighbors[v]) {
				neighbors[v]->triangles_around_vertex(vertex, partial);
			}
		}

		/*
		 * Traverse the sub-graph of triangles around a given vertex,
		 * and return an aggregate result for a given function,
		 * ignoring non-finite values.
		 *
		 * TYPE is one of:
		 *
		 * 	0 returns the sum
		 * 	1 returns the maximum
		 * 	2 returns the minimum
		 *
		 * STARTING_NODE should be left at its default value of 1 when
		 * this function is called externally.
		 */
		ale_accum traverse_around_vertex(point *vertex, ale_accum (triangle::*fp)(), int type, 
				int starting_node = 1) {

			if (starting_node)
				assert(vertex_ref_maybe(vertex) >= 0);

			/*
			 * If we've already visited this node, or if this node
			 * does not contain the given vertex, then return immediately.
			 */

			if (traversal_color == 1 || vertex_ref_maybe(vertex) < 0) {
				double p0 = +0;
				double p1 = +1;
				double n1 = -1;
				double pinf = p1 / p0;
				double ninf = n1 / p0;

				assert (!finite(pinf));
				assert (!finite(ninf));
				assert (pinf > 1);
				assert (ninf < -1);

				if (type == 0)
					return 0;
				if (type == 1)
					return ninf;
				if (type == 2)
					return pinf;

				assert(0);

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

			if (!finite(aggregate) && type == 0)
				aggregate = 0;

			for (int v = 0; v < 3; v++)
			if  (neighbors[v]) {
				ale_accum new_value = neighbors[v]->traverse_around_vertex(vertex, fp, type, 0);
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
				case 2:
					/*
					 * Minimum
					 */
					if (new_value < aggregate)
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
		void split_internals(int v, point *nv, int internal_fixed = 0) {
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
				children[c]->vertex_fixed[(division_vertex + 1 + c) % 3] 
					= neighbors[division_vertex] ? internal_fixed : 1;
			}

			for (int i = 0; i < 2; i++) {
				int vv = (v + 1 + i) % 3;
				if (!neighbors[vv])
					continue;

				int self_ref = self_ref_from_neighbor(vv);

				neighbors[vv]->neighbors[self_ref] = children[i];
			}

		}

		void split(int v, point nv, int internal_fixed = 0) {
			point *nvp = new point(nv);

			assert (nvp);

			split_internals(v, nvp, internal_fixed);

			if (!neighbors[v])
				return;

			int self_ref = self_ref_from_neighbor(v);

			neighbors[v]->split_internals(self_ref, nvp, internal_fixed);

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

		void split_at(int v, int internal_fixed = 0) {
			split(v, (*vertices[(v + 1) % 3] + *vertices[(v + 2) % 3]) / 2, internal_fixed);
		}

		/*
		 * Split a number of triangles so as to establish a lower bound
		 * on the size of newly-created angles.
		 */

		void split(int internal_fixed = 0) {

			int v;
			double angle_lbound = 59 * M_PI / 180;

			for (v = 0; v < 3; v++)
				if (vertex_angle(v) >= angle_lbound)
					break;

			assert (v < 3);
			assert (v >= 0);

			if (!neighbors[v]) {
				split_at(v, internal_fixed);
				return;
			}

			while (neighbors[v]->vertex_angle(self_ref_from_neighbor(v)) 
			     < angle_lbound) {

				neighbors[v]->split();

				/*
				 * Check to see whether this triangle has been split as a
				 * consequence.
				 */

				if (division_new_vertex) {

					if (internal_fixed)
					for (int c = 0; c < 2; c++) {
						int ref = children[c]->vertex_ref(division_new_vertex);
						children[c]->vertex_fixed[ref] = 1;
					}

					return;
				}

			}

			split_at(v, internal_fixed);
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
				split_at(v);
				return;
			}

			for (int v = 0; v < 3; v++)
			if  (angles[v] < angles[(v + 1) % 3]
			  && angles[v] < angles[(v + 2) % 3]) {
				split_at((v + 1 + rand() % 2) % 3);
				return;
			}

			split_at(rand() % 3);
		}
#endif

		void unsplit_internals() {
			assert(children[0]);
			assert(children[1]);
			assert(division_vertex >= 0);
			assert(division_vertex <  3);

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

		/*
		 * Find a triangle to accommodate the new control point NV.
		 * Cost must be lower than *COST; lowest cost triangle pointer
		 * is stored in *T.
		 */
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

		ale_pos compute_projected_area(pt _pt) {
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

		/*
		 * Find the triangle-coplanar point closest to a given point.
		 */
		point nearest_coplanar_point(point p, int offset = 0) {

			point v[3] = { (*vertices[(0 + offset) % 3]) - p, 
				       (*vertices[(1 + offset) % 3]) - p, 
				       (*vertices[(2 + offset) % 3]) - p };

			return rt_intersect(-normal(), v);
		}

		int closest_point_is_internal(point p) {
			point k = nearest_coplanar_point(p);

			if (k[0] >= 0 && k[1] >= 0 && k[0] + k[1] <= 1)
				return 1;

			return 0;
		}

		/*
		 * Distance between a point and a triangle.
		 */
		ale_pos distance_to(point p) {
			point k = nearest_coplanar_point(p);
			point vec1 = (*vertices[1]) - (*vertices[0]);
			point vec2 = (*vertices[2]) - (*vertices[0]);
			point vec3 = (*vertices[2]) - (*vertices[1]);

			/*
			 * If the projected point is internal, then return the
			 * distance to the projected point.
			 */

			if (k[0] >= 0 && k[1] >= 0 && k[0] + k[1] <= 1) {
				return fabs(k[2]);
			}

			/*
			 * Establish a starting length, and then try various
			 * distances to triangle boundaries, searching for a
			 * smaller length.
			 */

			ale_pos min_length = p.lengthto(*vertices[0]);

			/*
			 * Search other vertices for a smaller length
			 */

			for (int v = 1; v < 3; v++)
				if (p.lengthto(*vertices[v]) < min_length)
					min_length = p.lengthto(*vertices[v]);

			/*
			 * Check the edge opposite vertex 2.
			 */

			ale_pos proj1_len = k[0] + vec1.normalize().dproduct(k[1] * vec2) / vec1.norm();
			point proj1 = (*vertices[0]) + vec1 * proj1_len;

			if (proj1_len >= 0 && proj1_len <= 1 && p.lengthto(proj1) < min_length)
				min_length = p.lengthto(proj1);

			/*
			 * Check the edge opposite vertex 1.
			 */

			ale_pos proj2_len = k[1] + vec2.normalize().dproduct(k[0] * vec1) / vec2.norm();
			point proj2 = (*vertices[0]) + vec2 * proj2_len;

			if (proj2_len >= 0 && proj2_len <= 1 && p.lengthto(proj2) < min_length)
				min_length = p.lengthto(proj2);

			/*
			 * Check the edge opposite vertex 0.
			 */

			ale_pos proj3_len = vec3.normalize().dproduct(k[0] * vec1) / vec3.norm()
				        + vec3.normalize().dproduct(k[1] * vec2) / vec3.norm();
			point proj3 = (*vertices[1]) + vec3 * proj3_len;

			if (proj3_len >= 0 && proj3_len <= 1 && p.lengthto(proj3) < min_length)
				min_length = p.lengthto(proj3);

			return min_length;
		}

		/*
		 * Find the triangle for which point addition incurs the lowest cost.
		 * Cost is evaluated based on the distance between the point and the
		 * triangle.
		 */
		void find_least_cost_triangle(point nv, triangle **t, ale_pos *cost) {
			assert(point::defined(nv));

			if (division_new_vertex) {
				children[0]->find_least_cost_triangle(nv, t, cost);
				children[1]->find_least_cost_triangle(nv, t, cost);
				return;
			}

			/*
			 * Evaluate the distance between the triangle and the
			 * new control point.
			 */

			ale_pos this_cost = distance_to(nv);
			assert(this_cost >= 0);
			if (this_cost < *cost) {
				*t = this;
				*cost = this_cost;
			}

			return;
		}

		/*
		 * Add a control point to an existing triangle.
		 */
		void add_control_point(point nv) {
			ale_pos d[3];
			int nearest[3] = {0, 1, 2};

			/*
			 * Establish the distance to each vertex.
			 */
			for (int v = 0; v < 3; v++) {
				d[v] = vertices[v]->lengthto(nv);
				assert (finite(d[v]));
			}

			/*
			 * Bubble sort the vertices according to distance.
			 * (Why not?  There are only three vertices.)
			 */

			while (!(d[nearest[0]] <= d[nearest[1]] && d[nearest[1]] <= d[nearest[2]])) {
				for (int v = 0; v < 2; v++) {
					if (d[nearest[v]] > d[nearest[v + 1]]) {
						int temp = nearest[v + 1];
						nearest[v + 1] = nearest[v];
						nearest[v]     = temp;
					}
				}
			}

			/*
			 * In the case of an internal projected point, we can
			 * use any free vertex, or, if no vertices are free,
			 * then we can split the triangle and use the split
			 * point, unless the split point is on the surface
			 * boundary, in which case we need to recurse.
			 */

			if (closest_point_is_internal(nv)) {

				/*
				 * First, try to use the closest free vertex.
				 */
				for (int v = 0; v < 3; v++) {
					if (!vertex_fixed[nearest[v]]) {
						(*vertices[nearest[v]]) = nv;
						vertex_fixed[nearest[v]] = 1;
						return;
					}
				}

				/*
				 * No free vertices, so we split the triangle.
				 */

				split(1);

				/*
				 * If the split point is shared with a
				 * neighbor, then we can move it.
				 */
				if (neighbors[division_vertex]) {
					(*division_new_vertex) = nv;
					return;
				}

				/*
				 * Otherwise, we need to recurse.
				 */
				ale_pos one = +1;
				ale_pos zero = +0;
				ale_pos cost = one / zero;
				assert (cost > 0);
				assert (!finite(cost));
				assert (isinf(cost));

				triangle *lct = NULL;

				find_least_cost_triangle(nv, &lct, &cost);

				assert (lct);

				if (lct)
					lct->add_control_point(nv);
				
				return;
			}

			/*
			 * Determine the location of the point in edge
			 * coordinates.
			 */

			point k[3];
			for (int v = 0; v < 3; v++)
				k[v] = nearest_coplanar_point(nv, v);

			/*
			 * Check whether we can move an existing free vertex.
			 */
			for (int v = 0; v < 3; v++) {
				if (!vertex_fixed[nearest[v]] && k[nearest[v]][0] + k[nearest[v]][1] < 1) {
					(*vertices[nearest[v]]) = nv;
					vertex_fixed[nearest[v]] = 1;
					return;
				}
			}

			/*
			 * Check whether a vertex must be moved, even if non-free.
			 */
			for (int v = 0; v < 3; v++) {
				if (k[v][0] <= 0 && k[v][1] <= 0) {
					assert (vertex_fixed[v] == 1);
					point old_point = (*vertices[v]);
					(*vertices[v]) = nv;

					if (k[v][0] == 0)
						split_at((v + 1) % 3, 1);
					else
						split_at((v + 2) % 3, 1);

					(*division_new_vertex) = old_point;
					return;
				}
			}

			/*
			 * Now we know that the new point is in an area where
			 * an edge split can accommodate the point.  Find the
			 * correct edge to split.
			 */
			for (int v = 0; v < 3; v++) {
				if (k[v][0] > 0 && k[v][1] > 0) {
					split_at(v, 1);
					(*division_new_vertex) = nv;
					return;
				}
			}

			/*
			 * Ideally, we won't get here.  If we do, we can
			 * probably just throw the point out.
			 */

			assert(0);
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
		 * Return the sum of angles formed with a neighbor triangle.
		 * Searches one neighbor deep.
		 */
		ale_accum sum_neighbor_angle() {

			ale_pos sum = 0;

			point _normal = normal();

			for (int n = 0; n < 3; n++) {
				if (neighbors[n] == NULL)
					continue;

				ale_pos angle = point(0, 0, 0).anglebetw(neighbors[n]->normal(), _normal);

				sum += angle;
			}

			return sum;
		}

		/*
		 * Return the maximum angle formed with a neighbor triangle.
		 * Searches one neighbor deep.
		 */
		ale_accum max_neighbor_angle() {

			ale_pos max_angle = 0;

			point _normal = normal();

			for (int n = 0; n < 3; n++) {
				if (neighbors[n] == NULL)
					continue;

				ale_pos angle = point(0, 0, 0).anglebetw(neighbors[n]->normal(), _normal);

				if (angle > max_angle)
					max_angle = angle;
			}

			return max_angle;
		}

		/*
		 * Return the maximum internal angle of a triangle.
		 */
		ale_accum max_internal_angle() {
			ale_pos max_angle = 0;

			for (int v = 0; v < 3; v++) {
				ale_pos angle = vertices[v]->anglebetw(*vertices[(v + 1) % 3], *vertices[(v + 2) % 3]);

				if (angle > max_angle)
					max_angle = angle;
			}

			return max_angle;
		}

		/*
		 * Return the minimum internal angle of a triangle.
		 */
		ale_accum min_internal_angle() {
			ale_pos min_angle = M_PI;

			for (int v = 0; v < 3; v++) {
				ale_pos angle = vertices[v]->anglebetw(*vertices[(v + 1) % 3], *vertices[(v + 2) % 3]);

				if (angle < min_angle)
					min_angle = angle;
			}

			return min_angle;
		}

		/*
		 * Calculate the error between the model and the reference
		 * images.  Also includes edge length costs.
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

			/*
			 * Add edge length error.
			 */

			error += ((*vertices[0]).lengthto(*vertices[1])
			        + (*vertices[0]).lengthto(*vertices[2])
			        + (*vertices[1]).lengthto(*vertices[2])) * edge_cost_multiplier
				                                         * cl->sf;
			
			/*
			 * Add angle error.
			 */

			error += sum_neighbor_angle() * angle_cost_multiplier;

			return error;

		}

		ale_accum edge_cost() {
			return ((*vertices[0]).lengthto(*vertices[1])
			      + (*vertices[0]).lengthto(*vertices[2])
			      + (*vertices[1]).lengthto(*vertices[2])) * edge_cost_multiplier
				                                       * cl->sf;
		}

		ale_accum angle_cost() {
			return sum_neighbor_angle() * angle_cost_multiplier;
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

//					point centroid_local = _pt.wc(centroid);
//					point _vertices[3] = {_pt.wc(*vertices[0]), _pt.wc(*vertices[1]), _pt.wc(*vertices[2])};
//					assert (1 || e == subtriangle_index(centroid_local, _vertices, PLANAR_SUBDIVISION_DEPTH, 0));
				}

				color /= weight;

				this->color[e] = color;
				this->weight[e] = weight;
			}
		}

		/*
		 * Display information about vertices.
		 */
		ale_accum report_on_vertices() {
			fprintf(stderr, "Triangle %p\n", this);
			for (int v = 0; v < 3; v++) {
				vertex_aux va = vam[vertices[v]];

				fprintf(stderr, " Vertex %d\n", v);
				fprintf(stderr, "  Position       %f %f %f\n", (*vertices[v])[0],
						                         (*vertices[v])[1],
									 (*vertices[v])[2]);
				fprintf(stderr, "  Last position  %f %f %f\n", va.last_position[0],
						                              va.last_position[1],
									      va.last_position[2]);
				fprintf(stderr, "  Last triangle  %p\n", va.last_triangle);
				fprintf(stderr, "  Update ID:     %lu\n", va.update_id);
				fprintf(stderr, "  Maximum angle: %f\n", va.max_angle);
			}

			return 0;
		}

		/*
		 * Step size for adjustment methods.  Determine the average
		 * distance between vertices and divide it by a constant.
		 */
		ale_accum step_size() {
			return (vertices[0]->lengthto(*vertices[1])
			      + vertices[1]->lengthto(*vertices[2])
			      + vertices[2]->lengthto(*vertices[0])) / 9;
		}

		/*
		 * Reduce the sum of maximum angle metrics.
		 */
		int reduce_angle_metrics() {
			if (children[0] && children[1]){
				return children[0]->reduce_angle_metrics()
				     | children[1]->reduce_angle_metrics();
			}

			/*
			 * If all vertices are fixed, then return.
			 */

			if (vertex_fixed[0]
			 && vertex_fixed[1]
			 && vertex_fixed[2])
				return 0;

			/*
			 * Determine the adjustment step size.
			 */

			ale_accum step = step_size();

			/*
			 * Test modifications to the current position.
			 */

			int improved = 0;

			for (int v = 0; v < 3; v++) {

				/*
				 * Check for a fixed vertex.
				 */

				if (vertex_fixed[v])
					continue;

				/*
				 * Evaluate the metrics at the current position
				 */

				ale_accum max_internal = traverse_around_vertex(vertices[v], &triangle::max_internal_angle, 1);
				ale_accum sum_neighbor = traverse_around_vertex(vertices[v], &triangle::sum_neighbor_angle, 0);

				/*
				 * Store the original vertex position
				 */

				point vertex_original_pos = *vertices[v];

				/*
				 * Perturb the position.
				 */

				for (int axis = 0; axis < 3; axis++)
				for (int dir = -1; dir <= 1; dir += 2) {

					/*
					 * Adjust the vertex under consideration
					 */

					*vertices[v] = vertex_original_pos + point::unit(axis) * step * dir;

					/*
					 * Determine the new values of the metrics.
					 */

					ale_accum max_internal_test = traverse_around_vertex(vertices[v],
							&triangle::max_internal_angle, 1);
					ale_accum sum_neighbor_test = traverse_around_vertex(vertices[v],
							&triangle::sum_neighbor_angle, 0);

					/*
					 * Check the new metric values against the old.
					 */

					if (max_internal_test + sum_neighbor_test 
					  < max_internal      + sum_neighbor     ) {
						max_internal = max_internal_test;
						sum_neighbor = sum_neighbor_test;
						vertex_original_pos = *vertices[v];
						improved = 1;
						break;
					} else {
						*vertices[v] = vertex_original_pos;
					}
				}
			}

			return improved;
		}

			

		/*
		 * Adjust vertices according to mapped frames.
		 *
		 * Return non-zero if an adjustment is made.
		 */
		int adjust_vertices(zbuf_elem **z) {

			if (children[0] && children[1]) {
				return children[0]->adjust_vertices(z)
				     | children[1]->adjust_vertices(z);
			}

			/*
			 * If all vertices are fixed, then return.
			 */

			if (vertex_fixed[0]
			 && vertex_fixed[1]
			 && vertex_fixed[2])
				return 0;

			/*
			 * Determine the adjustment step size.
			 */

			ale_accum step = step_size();

			int improved = 0;

			ale_pos allowable_max_neighbor_angle = M_PI / 1.5;
			// ale_pos allowable_max_internal_angle = M_PI / 1.5;
			// ale_pos allowable_min_internal_angle = M_PI / 8;

			/*
			 * Test modifications to the current position of each
			 * vertex in turn.
			 */

			for (int v = 0; v < 3; v++) {

				/*
				 * Check for fixed vertices.
				 */

				if (vertex_fixed[v])
					continue;

				/*
				 * Perturb the position.
				 */

				for (int axis = 0; axis < 3; axis++)
				for (int dir = -1; dir <= 1; dir += 2) {

					point orig = *vertices[v];
					point perturbed = orig + point::unit(axis) * step * dir;

					ale_accum extremum_angle;

					/*
					 * Adjust the vertex under consideration
					 */

					/*
					 * Eliminate from consideration any change that increases to more than
					 * a given amount the angle between the normals of adjacent triangles.
					 */

					*vertices[v] = perturbed;
					extremum_angle = traverse_around_vertex(vertices[v], &triangle::max_neighbor_angle, 1);
					*vertices[v] = orig;
					if (extremum_angle > allowable_max_neighbor_angle)
						continue;

					/*
					 * Check the error
					 */

					ale_accum error = vertex_movement_cost(this, vertices[v], perturbed, z);
					if (error < 0) {
						improved = 1;
						*vertices[v] = perturbed;
						break;
					}
				}
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

		if (multipliers[0] >= 0
		 && multipliers[1] >= 0
		 && (multipliers[0] + multipliers[1] <= 1))
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
	static zbuf_elem *init_zbuf(pt _pt) {
		zbuf_elem *result = new zbuf_elem[(int) floor(_pt.scaled_width()) * (int) floor(_pt.scaled_height())];

		assert(result);

		if (!result)
			ui::get()->memory_error("Z-buffer");

		return result;
	}

	/*
	 * Z-buffer to determine the closest triangle.
	 *
	 * NOTE: we keep track of the depths of all intersected triangles, to
	 * make it easier to determine the closest triangle after the scene has
	 * been changed.
	 */
	static void zbuffer(pt _pt, zbuf_elem *zbuf, triangle *t) {
		int height = (int) floor(_pt.scaled_height());
		int width  = (int) floor(_pt.scaled_width());

		while (t) {

			while (t->division_new_vertex) {
				assert (t->children[0]);
				assert (t->children[1]);
				t = t->children[0];
			}

			/*
			 * Map the points of the triangle into local cartesian and
			 * projective coordinates.
			 */

			point c[3];
			point p[3];

			for (int v = 0; v < 3; v++) {
				c[v] = _pt.wc(*t->vertices[v]);
				p[v] = _pt.cp_scaled(c[v]);
			}

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

				zbuf_elem *element = zbuf + width * i + j;

				/*
				 * Test for interiority
				 *
				 * XXX: It might be better to use the
				 * is_interior_c() test instead...
				 */

				point test_point = point(i, j, 1);

				if (!is_interior(p, test_point))
					continue;

				element->insert(t);
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
			zbuf_elem *zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all points in the frame.
			 */
			for (unsigned int i = 0; i < im->height(); i++)
			for (unsigned int j = 0; j < im->width();  j++) {
				triangle *t = zbuf[i * im->width() + j].nearest(_pt, i, j);

				/*
				 * Check for points without associated triangles.
				 */
				if (!t)
					continue;

				/*
				 * Set new color and weight.
				 */
				int sti = subtriangle_index(_pt, d2::point(i, j), t);
				t->color[sti] = (t->color[sti] * t->weight[sti] + im->get_pixel(i, j)) 
					      / (t->weight[sti] + d2::pixel(1, 1, 1));
				t->weight[sti] = t->weight[sti] + d2::pixel(1, 1, 1);
			}

			delete[] zbuf;
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
			zbuf_elem *zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all points in the frame.
			 */
			for (unsigned int i = 0; i < height; i++)
			for (unsigned int j = 0; j < width;  j++) {
				triangle *t = zbuf[i * width + j].nearest(_pt, i, j);

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

				ale_pos area = t->compute_projected_area(_pt);

				if (area / PLANAR_SUBDIVISION_COUNT <= 4 && split)
					t->aux_stat = d2::pixel(-1, -1, -1);
				else if (area / PLANAR_SUBDIVISION_COUNT < 4 && !split && t->parent)
					t->parent->aux_stat = d2::pixel(-1, -1, -1);
			}

			delete[] zbuf;
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
			zbuf_elem *zbuf = init_zbuf(_pt);
			zbuffer(_pt, zbuf, triangle_head[0]);
			zbuffer(_pt, zbuf, triangle_head[1]);

			/*
			 * Iterate over all points in the frame, adding this frame
			 * to the list of frames including the associated triangle.
			 */
			for (unsigned int i = 0; i < im->height(); i++)
			for (unsigned int j = 0; j < im->width();  j++) {
				triangle *t = zbuf[i * im->width() + j].nearest(_pt, i, j);

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

			delete[] zbuf;
		}
	}

	static zbuf_elem **construct_zbuffers() {
		zbuf_elem **result = (zbuf_elem **) malloc(d2::image_rw::count() * sizeof(zbuf_elem *));

		for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
			pt _pt = align::projective(i);
			_pt.scale(cl->sf / _pt.scale_2d());
			result[i] = init_zbuf(_pt);
			zbuffer(_pt, result[i], triangle_head[0]);
			zbuffer(_pt, result[i], triangle_head[1]);
		}

		return result;
	}

	static void free_zbuffers(zbuf_elem **z) {
		for (unsigned int i = 0; i < d2::image_rw::count(); i++)
			delete[] z[i];
		free(z);
	}

	/*
	 * Adjust vertices to minimize scene error.  Return 
	 * non-zero if improvements are made.
	 */
	static int adjust_vertices() {
		/*
		 * Collect z-buffer data
		 */

		zbuf_elem **z = construct_zbuffers();

		int result = triangle_head[0]->adjust_vertices(z) | triangle_head[1]->adjust_vertices(z);

		free_zbuffers(z);

		return result;
	}

public:

	static void ecm(ale_pos ecm) {
		edge_cost_multiplier = ecm;
	}

	static void acm(ale_pos acm) {
		angle_cost_multiplier = acm;
	}
	
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

		point scene_plane[3] = {point(0, 0, 0), point(1, 0, 0), point(0, 1, 0)};

		for (int n = 0; n < N; n++) {

			pt t = align::projective(n);

			point tsp[3];

			for (int p = 0; p < 3; p++)
				tsp[p] = t.wc(scene_plane[p]);

			point a[4];

			a[0] = point(0, 0, -1);
			a[1] = point(t.scaled_height() - 1, 0, -1);
			a[2] = point(t.scaled_height() - 1, t.scaled_width() - 1, -1);
			a[3] = point(0, t.scaled_width() - 1, -1);

			for (int p = 0; p < 4; p++) {
				a[p] = t.pc_scaled(a[p]);
				a[p] = rt_intersect(a[p], tsp);
			}

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

		zbuf_elem *zbuf = init_zbuf(_pt);

		zbuffer(_pt, zbuf, triangle_head[0]);
		zbuffer(_pt, zbuf, triangle_head[1]);


#if 1
		zbuf_elem **zbs = construct_zbuffers();

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			d2::pixel val;
			ale_real div = 0;

			for (unsigned int f = 0; f < d2::image_rw::count(); f++) {

				pt _ptf = align::projective(f);
				_ptf.scale(cl->sf / _ptf.scale_2d());

				if (f == n) {
					point p = _ptf.wp_scaled(_pt.pw_scaled(point(i, j, -1)));

					if (!cl->reference[f]->in_bounds(p.xy()))
						continue;

					val += cl->reference[f]->get_bl(p.xy());
					div += 1;

					continue;
				}

				point p = frame_to_frame(d2::point(i, j), _pt, _ptf, zbuf, zbs[f]);

				if (!p.defined())
					continue;

				d2::pixel v = cl->reference[f]->get_bl(p.xy());

				val += v;
				div += 1;
			}

			im->pix(i, j) = val / div;
		}

		free_zbuffers(zbs);
#else
		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {

			triangle *t = zbuf[i * im->width() + j].nearest(_pt, i, j);

			if (!t)
				continue;

			int sti = subtriangle_index(_pt, d2::point(i, j), t);
			im->pix(i, j) = t->color[sti];
//			im->pix(i, j) = d2::pixel(sti,sti,sti) / (double) PLANAR_SUBDIVISION_COUNT;
		}
#endif

		delete[] zbuf;

		return im;
	}

	static const d2::image *depth(unsigned int n) {
		assert (n < d2::image_rw::count());

		d2::image *im = new d2::image_ale_real((int) floor(d2::align::of(n).scaled_height()),
				               (int) floor(d2::align::of(n).scaled_width()), 3);

		pt _pt = align::projective(n);

		zbuf_elem *zbuf = init_zbuf(align::projective(n));

		zbuffer(_pt, zbuf, triangle_head[0]);
		zbuffer(_pt, zbuf, triangle_head[1]);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			triangle *t = zbuf[i * im->width() + j].nearest(_pt, i, j);

			if (!t)
				continue;

#if 1
			point r = _pt.pc_scaled(point(i, j, -1));
			point vertices[3];

			for (int v = 0; v < 3; v++)
				vertices[v] = _pt.wc(*t->vertices[v]);

			point intersect = rt_intersect(r, vertices);

			im->pix(i, j) = d2::pixel(1, 1, 1) * -intersect[2];
#else
			im->pix(i, j) = d2::pixel(1, 1, 1) * _pt.wc(t->centroid())[2];
#endif
		}

		delete[] zbuf;

		return im;
	}

	/*
	 * Add specified control points to the model
	 */
	static void add_control_points() {

		for (unsigned int i = 0; i < cpf::count(); i++) {

			point control_point = cpf::get(i);

			if (!control_point.defined())
				continue;

			ale_pos one = +1;
			ale_pos zero = +0;

			ale_pos cost = one / zero;
			assert (cost > 0);
			assert (!finite(cost));
			assert (isinf(cost));

			triangle *lct = NULL;
			
			triangle_head[0]->find_least_cost_triangle(control_point, &lct, &cost);
			triangle_head[1]->find_least_cost_triangle(control_point, &lct, &cost);

			if (lct == NULL)
				continue;

			lct->add_control_point(control_point);
		}

	}

	/*
	 * Relax model to reduce angle sizes.
	 */
	static void relax_triangle_model() {

		if (!triangle_head[0]->reduce_angle_metrics()
		 && !triangle_head[1]->reduce_angle_metrics())
			return;

		fprintf(stderr, "Relaxing triangle model");
		while(triangle_head[0]->reduce_angle_metrics()
		   || triangle_head[1]->reduce_angle_metrics())
			fprintf(stderr, ".");
		fprintf(stderr, "\n");

	}


	/*
	 * When using a 3D scene data structure, improvements should occur in 
	 * two passes:
	 *
	 * 	(a) Moving vertices to reduce error
	 *
	 * 	(b) Attempting to subdivide triangles by adding new vertices 
	 *
	 * When neither of these approaches results in improvement, then the
	 * level of detail can be increased.
	 *
	 * XXX: the name of this function is horribly misleading.  There isn't
	 * even a 'search depth' any longer, since there is no longer any
	 * bounded DFS occurring.
	 */
	static void reduce_cost_to_search_depth(const char *d_out[], const char *v_out[], d2::exposure *exp_out, int inc_bit) {
		int improved = 1;
		int count = 0;

		/*
		 * To start, use the lowest level-of-detail
		 */

		while(reduce_lod());

		/*
		 * Perform cost-reduction.
		 */

		for (;;) {

			if (density_test_split())
				continue;

			color_average();
			
			determine_visibility();

			triangle_head[0]->recolor();
			triangle_head[1]->recolor();

			triangle_head[0]->free_aux_vars();
			triangle_head[1]->free_aux_vars();

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

			if (!improved || count > 40) {
				fprintf(stderr, ".");
				if (!cl->next)
					break;
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

			improved |= adjust_vertices();
		}

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
