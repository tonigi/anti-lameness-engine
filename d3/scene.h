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
 * View angle multiplier.  
 *
 * Setting this to a value larger than one can be useful for debugging.
 */

#define VIEW_ANGLE_MULTIPLIER 1

class scene {

	/*
	 * Geometric model structure.
	 */

	struct triangle;

	/*
	 * Subspace model structure.
	 */

	struct space;


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
	 * Clipping planes
	 */
	static ale_pos front_clip;
	static ale_pos rear_clip;

	/*
	 * Perturb bounds
	 *
	 * Bound types:
	 *
	 * 	0	Absolute bound (projected pixel units [model units])
	 * 	1	Relative bound (model top-level diameter percentage)
	 *
	 */
	static ale_pos mpl_value;
	static int mpl_type;
	static ale_pos mpu_value;
	static int mpu_type;

	/*
	 * Model files
	 */
	static const char *load_model_name;
	static const char *save_model_name;

	/*
	 * Nearness threshold
	 */
	static const ale_real nearness;

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
	 * Change set data types.
	 */
	struct point_2d_lt {
		bool operator()(d2::point a, d2::point b) const {
			for (int d = 0; d < 2; d++) {
				if (a[d] < b[d])
					return true;
				if (a[d] > b[d])
					return false;
			}
			return false;
		}
	};
	struct point_lt {
		bool operator()(point a, point b) const {
			for (int d = 0; d < 2; d++) {
				if (a[d] < b[d])
					return true;
				if (a[d] > b[d])
					return false;
			}
		}
	};
	typedef std::pair<d2::point,d2::point> pp_t;
	struct pp_lt {
		bool operator()(pp_t a, pp_t b) const {
			point_2d_lt lt;
			if (lt(a.first, b.first))
				return true;
			if (lt(b.first, a.first))
				return false;
			if (lt(a.second, b.second))
				return true;
			return false;
		}
	};
	typedef std::set<pp_t,pp_lt> pp_set_t;
	typedef std::set<point,point_lt> change_elem_t;
	struct change_elem_lt {
		bool operator()(change_elem_t a, change_elem_t b) const {
			for (change_elem_t::iterator i = a.begin(), j = b.begin(); 
					i != a.end() && j != b.end(); i++, j++) {
				point_lt lt;
				if (lt(*i, *j))
					return true;
				if (lt(*j, *i))
					return false;
			}
			if (a.size() < b.size())
				return true;

			return false;
		}
	};
	typedef std::set<change_elem_t,change_elem_lt> change_set_t;
	
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
	 * Pyramid bounds check.
	 */
	static int pyramid_bounds_check(point w) {
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
			pt _pt = align::projective(n);
			_pt.scale(1 / _pt.scale_2d());

			point p = _pt.wp_scaled(w);

			/*
			 * Check the four closest points.
			 */

			int p_int[2] = { (int) floor(p[0]), (int) floor(p[1]) };

			if (p_int[0] >= 0
			 && p_int[1] >= 0
			 && p_int[0] <= _pt.scaled_height() - 2
			 && p_int[1] <= _pt.scaled_width() - 2)
				return 1;
		}
		return 0;
	}

	/*
	 * Frame-to-frame mapping function
	 *
	 * Maps from point P in frame F1 to frame F2.  Returns projective
	 * coordinates.  Returns zero-depth if there is no mutually-visible
	 * model point.
	 */
	static point frame_to_frame(d2::point p, const pt &_pt1, const pt &_pt2, zbuf_elem *z1, zbuf_elem *z2);

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
	 * Perturbation amount
	 */

	static ale_pos perturb;

	/*
	 * Current level-of-detail
	 */

	static struct lod *cl;

	/*
	 * Base level-of-detail
	 */

	static struct lod *bl;

	/*
	 * Vertex movement cost
	 *
	 * Evaluate the cost of moving a vertex.  Negative values indicate
	 * improvement; non-negative values indicate lack of improvement.
	 */
	static ale_accum vertex_movement_cost(std::set<triangle *> *t_set, change_elem_t *v_set,
			point *vertex, point new_position, zbuf_elem **z, lod *_lod);

	/*
	 * Structure to hold a subdivisible triangle.
	 */

	struct triangle {

		/*
		 * Color for traversal.
		 */

		char traversal_color;

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
		 * Set of fixed vertices.
		 */

		static std::set<point *> fixed_vertices;

		/*
		 * Constructor
		 */

		triangle() {

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

			traversal_color = 0;
		}

		static void fix_vertex(point *p) {
			fixed_vertices.insert(p);
		}

		static void free_vertex(point *p) {
			fixed_vertices.erase(p);
		}

		static int vertex_is_free(point *p) {
			return !fixed_vertices.count(p);
		}

		/*
		 * Collect information from the triangle tree:
		 *
		 *   1. The set of all non-fixed vertices.
		 *   2. A multimap from non-fixed vertices to triangles.
		 *   3. The set of all fixed vertices (for verification).
		 */
		
		void collect_vertex_data(std::set<point *> *vertex_set,
				         std::multimap<point *, triangle *> *vertex_map,
					 std::set<point *> *fixed_set) {

			if (children[0] && children[1]) {
				children[0]->collect_vertex_data(vertex_set, vertex_map, fixed_set);
				children[1]->collect_vertex_data(vertex_set, vertex_map, fixed_set);
				return;
			}

			assert(!children[0]);
			assert(!children[1]);

			for (int v = 0; v < 3; v++) {
				if (fixed_vertices.count(vertices[v])) {
					fixed_set->insert(vertices[v]);
				} else {
					vertex_set->insert(vertices[v]);
					vertex_map->insert(std::make_pair(vertices[v], this));
				}
			}
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
		void split_internals(int v, point *nv) {
			assert (children[0] == NULL);
			assert (children[1] == NULL);
			assert (division_new_vertex == NULL);

			division_vertex = v;
			division_new_vertex = nv;

			if (!neighbors[division_vertex])
				fix_vertex(nv);

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

		void split(int v, point nv, int internal_fixed = 0) {
			point *nvp = new point(nv);

			assert (nvp);

			if (internal_fixed)
				fix_vertex(nvp);

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
						fix_vertex(division_new_vertex);

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

		int unsplit_on_camera_area_threshold(pt _pt, ale_pos threshold) {

			if (children[0] && children[1]) {
				return (children[0]->unsplit_on_camera_area_threshold(_pt, threshold)
				      | children[1]->unsplit_on_camera_area_threshold(_pt, threshold));
			}

			assert (!children[0] && !children[1] && !division_new_vertex);

			ale_pos _area = area() * _pt.w_density_scaled_max(*vertices[0], *vertices[1], *vertices[2]);

			if (_area > threshold)
				return 0;

			unsplit();
			return 1;
		}

		int split_on_camera_area_threshold(pt _pt, ale_pos threshold) {

			if (children[0] && children[1]) {
				return (children[0]->split_on_camera_area_threshold(_pt, threshold)
				      | children[1]->split_on_camera_area_threshold(_pt, threshold));
			}

			assert (!children[0] && !children[1] && !division_new_vertex);

			ale_pos _area = area() * _pt.w_density_scaled_max(*vertices[0], *vertices[1], *vertices[2]);

			if (_area <= threshold)
				return 0;

			split();
			return 1;
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
					if (vertex_is_free(vertices[nearest[v]])) {
						(*vertices[nearest[v]]) = nv;
						fix_vertex(vertices[nearest[v]]);
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
				if (vertex_is_free(vertices[nearest[v]]) && k[nearest[v]][0] + k[nearest[v]][1] < 1) {
					(*vertices[nearest[v]]) = nv;
					fix_vertex(vertices[nearest[v]]);
					return;
				}
			}

			/*
			 * Check whether a vertex must be moved, even if non-free.
			 */
			for (int v = 0; v < 3; v++) {
				if (k[v][0] <= 0 && k[v][1] <= 0) {
					assert (!vertex_is_free(vertices[v]));
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

		ale_accum edge_cost() {
			return ((*vertices[0]).lengthto(*vertices[1])
			      + (*vertices[0]).lengthto(*vertices[2])
			      + (*vertices[1]).lengthto(*vertices[2])) * edge_cost_multiplier
				                                       * (2 / perturb);
		}

		ale_accum angle_cost() {
			return sum_neighbor_angle() * angle_cost_multiplier;
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

			if (!vertex_is_free(vertices[0])
			 && !vertex_is_free(vertices[1])
			 && !vertex_is_free(vertices[2]))
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

				if (!vertex_is_free(vertices[v]))
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

					point perturbed = vertex_original_pos + point::unit(axis) * step * dir;

					if (perturbed[2] > front_clip
					 || perturbed[2] < rear_clip)
						continue;

					/*
					 * Adjust the vertex under consideration
					 */

					*vertices[v] = perturbed;

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


		int vertex_visible(int v, zbuf_elem **z, lod *_lod) {
			for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
				pt _pt = align::projective(n);
				_pt.scale(_lod->sf / _pt.scale_2d());

				point w = *vertices[v];
				point p = _pt.wp_scaled(w);

				/*
				 * Check the four closest points.
				 */

				int p_int[2] = { (int) floor(p[0]), (int) floor(p[1]) };

				if (p_int[0] < 0
				 || p_int[1] < 0
				 || p_int[0] > _pt.scaled_height() - 2
				 || p_int[1] > _pt.scaled_width() - 2)
					continue;

				for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++) {
					triangle *t = z[n][(p_int[0] + i) * (int) floor(_pt.scaled_width()) + p_int[1] + j].nearest(_pt, p_int[0] + i, p_int[1] + j);
					if (t->vertex_ref_maybe(vertices[v]))
						return 1;
				}
			}
			return 0;
		}
			

#if 0
		/*
		 * Adjust vertices according to mapped frames.
		 *
		 * Return non-zero if an adjustment is made.
		 */
		int adjust_vertices(zbuf_elem **z, lod *_lod) {

			pt _pt0 = d3::align::projective(0);
			// _pt0.scale(1 / _pt0.scale_2d());  // XXX: world coordinates vary with output scale.

			// fprintf(stderr, "%p traversing children [%u] \n", this, time(NULL));
			if (children[0] && children[1]) {
				return children[0]->adjust_vertices(z, _lod)
				     | children[1]->adjust_vertices(z, _lod);
			}
			// fprintf(stderr, "%p leaf node [%u] \n", this, time(NULL));

			/*
			 * If all vertices are fixed, then return.
			 */

			// fprintf(stderr, "%p checking for unfixed vertices [%u] \n", this, time(NULL));

			if (!vertex_is_free(vertices[0])
			 && !vertex_is_free(vertices[1])
			 && !vertex_is_free(vertices[2]))
				return 0;

			int improved = 0;

			ale_pos allowable_max_neighbor_angle = M_PI / 1.5;
			// ale_pos allowable_max_internal_angle = M_PI / 1.5;
			// ale_pos allowable_min_internal_angle = M_PI / 8;

			/*
			 * Test modifications to the current position of each
			 * vertex in turn.
			 */

			// fprintf(stderr, "%p iterating over vertices vertices [%u] \n", this, time(NULL));
			for (int v = 0; v < 3; v++) {

				/*
				 * Check for fixed vertices.
				 */

				// fprintf(stderr, "%p checking that vertex %d is unfixed [%u] \n", this, v, time(NULL));
				if (!vertex_is_free(vertices[v]))
					continue;

				/*
				 * Determine the adjustment step size according to
				 * the perturbation size.
				 */

				ale_pos step = perturb / sqrt(_pt0.w_density_scaled(*vertices[v]));

				/*
				 * Perturb the position.
				 */

				// fprintf(stderr, "%p perturbing vertex %d [%u] \n", this, v, time(NULL));
				for (int axis = 0; axis < 3; axis++)
				for (int dir = -1; dir <= 1; dir += 2) {

					// fprintf(stderr, "%p perturbing vertex %d (%d %d) [%u] \n", this, v, axis, dir, time(NULL));
					point orig = *vertices[v];
					point perturbed = orig + point::unit(axis) * step * dir;

					/*
					 * Check the clipping planes.
					 */

					// fprintf(stderr, "%p checking clipping planes [%u] \n", this, time(NULL));
					if (perturbed[2] > front_clip
					 || perturbed[2] < rear_clip)
						continue;

					/*
					 * Check view pyramid bounds
					 */

					// fprintf(stderr, "%p checking view pyramid bounds [%u] \n", this, time(NULL));
					if (!pyramid_bounds_check(perturbed))
						continue;

					/*
					 * Eliminate from consideration any change that increases to more than
					 * a given amount the angle between the normals of adjacent triangles.
					 */

					ale_accum extremum_angle;
					*vertices[v] = perturbed;
					extremum_angle = traverse_around_vertex(vertices[v], &triangle::max_neighbor_angle, 1);
					*vertices[v] = orig;
					if (extremum_angle > allowable_max_neighbor_angle)
						continue;

					/*
					 * Check the error
					 */

					// fprintf(stderr, "%p calculating error [%u] \n", this, time(NULL));
					ale_accum error = vertex_movement_cost(this, vertices[v], perturbed, z, _lod);
					// fprintf(stderr, "%p done calculating error (%f) [%u] \n", this, error, time(NULL));
					if (error < 0) {

						fprintf(stderr, "%p (%f %f %f) (%f %f %f) [%f]\n",
							vertices[v],
							(*vertices[v])[0],
							(*vertices[v])[1],
							(*vertices[v])[2],
							perturbed[0],
							perturbed[1],
							perturbed[2],
							2 / perturb);

						improved = 1;
						*vertices[v] = perturbed;
						break;
					}
				}
			}
			return improved;
		}
#endif

		void recursive_write_triangle_ascii_ptr(FILE *f, 
				std::map<triangle *, unsigned int> *triangle_map, 
				std::map<point *, unsigned int> *vertex_map, unsigned int *id) {

			if (!triangle_map->count(this))
				(*triangle_map)[this] = (*id)++;

			for (int v = 0; v < 3; v++)
			if (!vertex_map->count(vertices[v]))
				(*vertex_map)[vertices[v]] = (*id)++;

			for (int v = 0; v < 3; v++)
			if (!triangle_map->count(neighbors[v]))
				(*triangle_map)[neighbors[v]] = (*id)++;

			if (!triangle_map->count(parent))
				(*triangle_map)[parent] = (*id)++;

			if (!vertex_map->count(division_new_vertex))
				(*vertex_map)[division_new_vertex] = (*id)++;

			for (int v = 0; v < 2; v++)
			if (!triangle_map->count(children[v]))
				(*triangle_map)[children[v]] = (*id)++;

			fprintf(f, "T %u %u %u %u %u %u %u %u %u %u %u %d %u %u %u\n", 
				(*triangle_map)[this],
				(*vertex_map)[vertices[0]],
				(*vertex_map)[vertices[1]],
				(*vertex_map)[vertices[2]],
				(unsigned int) !vertex_is_free(vertices[0]),
				(unsigned int) !vertex_is_free(vertices[1]),
				(unsigned int) !vertex_is_free(vertices[2]),
				(*triangle_map)[neighbors[0]],
				(*triangle_map)[neighbors[1]],
				(*triangle_map)[neighbors[2]],
				(*triangle_map)[parent],
				division_vertex,
				(*vertex_map)[division_new_vertex],
				(*triangle_map)[children[0]],
				(*triangle_map)[children[1]]);

			if (children[0])
				children[0]->recursive_write_triangle_ascii_ptr(f, triangle_map, vertex_map, id);
			if (children[1])
				children[1]->recursive_write_triangle_ascii_ptr(f, triangle_map, vertex_map, id);
		}

		void read_triangle_ascii_ptr(FILE *f, std::map<unsigned int, triangle *> *triangle_map, 
				std::map<unsigned int, point *> *vertex_map) {
			unsigned int id;

			for (int v = 0; v < 3; v++) {
				if (!fscanf(f, "%u", &id))
					ui::get()->error("Bad model file.");
				if (!vertex_map->count(id))
					(*vertex_map)[id] = new point;
				vertices[v] = (*vertex_map)[id];
			}

			for (int v = 0; v < 3; v++) {
				unsigned int vf;
				if (!fscanf(f, "%u", &vf))
					ui::get()->error("Bad model file.");
				if (vf)
					fix_vertex(vertices[v]);
			}

			for (int v = 0; v < 3; v++) {
				if (!fscanf(f, "%u", &id))
					ui::get()->error("Bad model file.");
				if (!triangle_map->count(id))
					(*triangle_map)[id] = new triangle;
				neighbors[v] = (*triangle_map)[id];
			}

			if (!fscanf(f, "%u", &id))
				ui::get()->error("Bad model file.");
			if (!triangle_map->count(id))
				(*triangle_map)[id] = new triangle;
			parent = (*triangle_map)[id];

			if (!fscanf(f, "%u", &division_vertex))
				ui::get()->error("Bad model file.");

			if (!fscanf(f, "%u", &id))
				ui::get()->error("Bad model file.");
			if (!vertex_map->count(id))
				(*vertex_map)[id] = new point;
			division_new_vertex = (*vertex_map)[id];

			for (int v = 0; v < 2; v++) {
				if (!fscanf(f, "%u", &id))
					ui::get()->error("Bad model file.");
				if (!triangle_map->count(id))
					(*triangle_map)[id] = new triangle;
				children[v] = (*triangle_map)[id];
			}
		}
	};

	/*
	 * Structure to hold a subdivisible region of space.
	 */

	struct space {
		struct space *positive;
		struct space *negative;
	};

	/*
	 * Space root pointer
	 */
	static space *root_space;

	/*
	 * Space traversal and navigation class.
	 */

	class space_traverse {
		space *current;
		point min, max;
		int next_split;

	public:

		static space_traverse root() {

			space_traverse result;

			result.current = root_space;
			result.min = point(-M_PI/2, -M_PI/2, -M_PI/2);
			result.max = point( M_PI/2,  M_PI/2,  M_PI/2);
			result.next_split = 0;

			return result;
		}

		space_traverse positive() {

			if (current->positive == NULL)
				current->positive = new space;
			
			space_traverse result;

			result.current = current->positive;
			result.min = min;
			result.max = max;
			result.next_split = (next_split + 1) % 3;

			result.min[next_split] = (min[next_split] + max[next_split]) / 2;

			return result;
		}

		space_traverse negative() {

			if (current->negative == NULL)
				current->negative = new space;
			
			space_traverse result;

			result.current = current->negative;
			result.min = min;
			result.max = max;
			result.next_split = (next_split + 1) % 3;

			result.max[next_split] = (min[next_split] + max[next_split]) / 2;

			return result;
		}

		point get_min() {
			return point(tan(min[0]), tan(min[1]), tan(min[2]));
		}

		point get_max() {
			return point(tan(max[0]), tan(max[1]), tan(max[2]));
		}

		int includes(point p) {
			point tdp = point(atan(p[0]), atan(p[1]), atan(p[2]));

			for (int d = 0; d < 3; d++) {
				if (tdp[d] > max[d])
					return 0;
				if (tdp[d] < min[d])
					return 0;
				if (isnan(tdp[d]))
					return 0;
			}

			return 1;
		}

		space *get_space() {
			return current;
		}

		int get_next_split() {
			return next_split;
		}
	};

	/*
	 * Class to iterate through subspaces based on proximity to a camera.
	 */

	class space_iterate {
		std::stack<space_traverse> space_stack;
		point camera_origin;
	public:
		space_iterate(pt _pt) {
			camera_origin = _pt.cw(point(0, 0, 0));
			space_stack.push(space_traverse::root());
		}

		int next() {
			if (space_stack.empty())
				return 0;

			space_traverse st = space_stack.top();

			int d = st.get_next_split();

			ale_pos split_point = (st.get_max()[d] + st.get_min()[d]) / 2;

			space *n = st.get_space()->negative;
			space *p = st.get_space()->positive;

			if (camera_origin[d] > split_point) {
				if (n) {
					space_stack.top() = st.negative();
					if (p)
						space_stack.push(st.positive());
				} else {
					if (p)
						space_stack.top() = st.positive();
					else
						space_stack.pop();
				}
			} else {
				if (p) {
					space_stack.top() = st.positive();
					if (n)
						space_stack.push(st.negative());
				} else {
					if (n)
						space_stack.top() = st.negative();
					else
						space_stack.pop();
				}
			}

			return 1;
		}

		space_traverse get() {
			assert (!space_stack.empty());
			return space_stack.top();
		}
	};

	/*
	 * Class for information regarding spatial regions of interest.
	 *
	 * This class is configured for convenience in cases where sampling is
	 * performed using an approximation of the fine:box:1,triangle:2 chain.
	 * In this case, the *_1 variables would store the fine data and the
	 * *_2 variables would store the coarse data.  Other uses are also
	 * possible.
	 */

	class spatial_info {
		/*
		 * Map channel value --> weight.
		 */
		std::map<ale_real, ale_real> color_weights_1[3];
		std::map<ale_real, ale_real> color_weights_2[3];

		/*
		 * Current color.
		 */
		d2::pixel color;

		/*
		 * Map occupancy value --> weight.
		 */
		std::map<ale_real, ale_real> occupancy_weights_1;
		std::map<ale_real, ale_real> occupancy_weights_2;

		/*
		 * Current occupancy value.
		 */
		ale_real occupancy;

		/*
		 * Insert a weight into a list.
		 */
		void insert_weight(std::map<ale_real, ale_real> *m, ale_real v, ale_real w) {
			if (m->count(v))
				(*m)[v] += w;
			else
				(*m)[v] = w;
		}

		/*
		 * Find the median of a weighted list.  Uses NaN for undefined.
		 */
		ale_real find_median(std::map<ale_real, ale_real> *m) {
			std::map<ale_real, ale_real>::iterator a = m->begin(),
				                               b = m->begin();

			ale_real zero1 = 0;
			ale_real zero2 = 0;
			ale_real undefined = zero1 / zero2;

			ale_accum weight_sum = 0;

			while (a != m->end()) {
				weight_sum += a->first;
				a++;
			}

			if (weight_sum == 0)
				return undefined;

			ale_accum midpoint = weight_sum / 2;

			ale_accum weight_sum_2 = 0;

			while (b != m->end() && weight_sum_2 < midpoint) {
				weight_sum_2 += b->first;

				if (weight_sum_2 > midpoint) {
					return b->second;
				} else if (weight_sum_2 == midpoint) {
					ale_accum first_value = b->second;
					b++;
					assert (b != m->end());

					return (first_value + b->second) / 2;
				} 

				b++;
			}

			assert(0);

			return undefined;
		}

	public:
		/*
		 * Constructor.
		 */
		spatial_info() {
			color = d2::pixel::zero();
			occupancy = 0.5;
		}

		/*
		 * Accumulate color; primary data set.
		 */
		void accumulate_color_1(d2::pixel color, d2::pixel weight) {
			for (int k = 0; k < 3; k++)
				insert_weight(&color_weights_1[k], color[k], weight[k]);
		}

		/*
		 * Accumulate color; secondary data set.
		 */
		void accumulate_color_2(d2::pixel color, d2::pixel weight) {
			for (int k = 0; k < 3; k++)
				insert_weight(&color_weights_2[k], color[k], weight[k]);
		}

		/*
		 * Accumulate occupancy; primary data set.
		 */
		void accumulate_occupancy_1(ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_1, occupancy, weight);
		}

		/*
		 * Accumulate occupancy; secondary data set.
		 */
		void accumulate_occupancy_2(ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_2, occupancy, weight);
		}

		/*
		 * Update color (and clear accumulation structures).
		 */
		void update_color() {
			for (int d = 0; d < 3; d++) {
				ale_real c = find_median(&color_weights_1[d]);
				if (isnan(c))
					c = find_median(&color_weights_2[d]);
				if (isnan(c))
					c = 0;

				color[d] = c;

				color_weights_1[d].clear();
				color_weights_2[d].clear();
			}
		}

		/*
		 * Update occupancy (and clear accumulation structures).
		 */
		void update_occupancy() {
			ale_real o = find_median(&occupancy_weights_1);
			if (isnan(o))
				o = find_median(&occupancy_weights_2);
			if (isnan(o))
				o = 0.5;

			occupancy = o;

			occupancy_weights_1.clear();
			occupancy_weights_2.clear();
		}

		/*
		 * Get current color.
		 */
		d2::pixel get_color() {
			return color;
		}

		/*
		 * Get current occupancy.
		 */
		ale_real get_occupancy() {
			return occupancy;
		}
	};

	/*
	 * Map spatial regions of interest to spatial info structures.  XXX:
	 * This may get very poor cache behavior in comparison with, say, an
	 * array.  Unfortunately, there is no immediately obvious array
	 * representation.  If some kind of array representation were adopted,
	 * it would probably cluster regions of similar depth from the
	 * perspective of the typical camera.  In particular, for a
	 * stereoscopic view, depth ordering for two random points tends to be
	 * similar between cameras, I think.  Unfortunately, it is never
	 * identical for all points (unless cameras are co-located).  One
	 * possible approach would be to order based on, say, camera 0's idea
	 * of depth.
	 */

	static std::map<struct space *, spatial_info> spatial_info_map;

	/*
	 * Use a pair of trees to store the triangles.
	 */
	static struct triangle *triangle_head[2];

#if 0
	/*
	 * Map view B and pixel A(i, j) within view A to a set of depths where
	 * A(i, j) matches a point in image B.
	 */
	std::set<ale_pos> depth_match_set(int B, int A, int i, int j) {
	}
#endif

	/*
	 * Read the model file.
	 */
	static void read_model_file() {
		assert(load_model_name);

		FILE *f = fopen(load_model_name, "r");

		if (!f)
			ui::get()->error("Can't open model load file.");

		char command;
		int int_value;
		unsigned int id;
		int root_index;

		/*
		 * Check model file version.
		 */

		if (fscanf(f, "%c %d", &command, &int_value) != 2
		 || command != 'V')
			ui::get()->error("Bad model file.");

		if (int_value != 0)
			ui::get()->error("Unrecognized model file version.");

		/*
		 * Map variables.
		 */

		std::map<unsigned int, triangle *> triangle_map;
		std::map<unsigned int, point *> vertex_map;

		triangle_map[0] = 0;
		vertex_map[0] = 0;

		/*
		 * Read commands
		 */

		while (!feof(f)) {
			if (fscanf(f, "%c", &command) < 1) {
				assert(feof(f));
				break;
			}

			switch(command) {
			case ' ':
			case '\t':
			case '\n':
			case '\r':
				break;
			case 'T':
				/*
				 * Triangle
				 */
				if (!fscanf(f, "%u", &id))
					ui::get()->error("Bad model file.");
				if (!triangle_map.count(id))
					triangle_map[id] = new triangle;
				triangle_map[id]->read_triangle_ascii_ptr(f, &triangle_map, &vertex_map);
				break;
			case 'R':
				/*
				 * Roots
				 */
				if (fscanf(f, "%d %u", &root_index, &id) != 2)
					ui::get()->error("Bad model file.");
				if (root_index < 0 || root_index > 1)
					ui::get()->error("Cannot handle roots other than 0 or 1");
				if (!triangle_map.count(id))
					triangle_map[id] = new triangle;
				triangle_head[root_index] = triangle_map[id];
				break;
			case 'P':
				/*
				 * Vertices
				 */
				if (fscanf(f, "%u", &id) != 1)
					ui::get()->error("Bad model file.");

				if(!vertex_map.count(id))
					vertex_map[id] = new point;

				for (int v = 0; v < 3; v++) {
					double coord;
					if (!fscanf(f, "%lf", &coord))
						ui::get()->error("Bad model file.");
					(*vertex_map[id])[v] = coord;
				}
				break;
			case 'C':
				unsigned int index;
				ale_pos view_angle, x, y, z, P, Y, R;
				if (fscanf(f, "%u %f %f %f %f %f %f %f", &index,
						&view_angle, &y, &x, &z, &Y, &P, &R) != 8)
					ui::get()->error("Bad model file.");

				align::set_view_angle(view_angle * M_PI / 180);
				align::set_translation(index, 0, x);
				align::set_translation(index, 1, y);
				align::set_translation(index, 2, z);
				align::set_rotation(index, 0, P * M_PI / 180);
				align::set_rotation(index, 1, Y * M_PI / 180);
				align::set_rotation(index, 2, R * M_PI / 180);

				break;
			default:
				fprintf(stderr, "Unknown model file command '%c' (%x).\n", command, command);
				ui::get()->error("Bad model file.");
			}
		}

		fclose(f);
	}

	/*
	 * Write the model file.
	 */
	static void write_model_file() {
		if (!save_model_name)
			return;

		FILE *f = fopen(save_model_name, "w");

		if (!f)
			ui::get()->error("Can't open model save file.");

		/*
		 * Print model file version.
		 */

		fprintf(f, "V %d\n", 0);

		/*
		 * Map variables.
		 */

		std::map<triangle *, unsigned int> triangle_map;
		std::map<point *, unsigned int> vertex_map;

		triangle_map[0] = 0;
		vertex_map[0] = 0;

		unsigned int id = 1;

		/*
		 * Print root values.
		 */

		for (int r = 0; r < 2; r++) {
			triangle_map[triangle_head[r]] = id++;
			fprintf(f, "R %d %u\n", r, triangle_map[triangle_head[r]]);
		}

		/*
		 * Print trees
		 */

		for (int r = 0; r < 2; r++)
			triangle_head[r]->recursive_write_triangle_ascii_ptr(f, &triangle_map, &vertex_map, &id);

		/*
		 * Print vertex data
		 */

		for (std::map<point *, unsigned int>::iterator i = 
				vertex_map.begin(); i != vertex_map.end(); i++) {

			if (i->first == NULL)
				continue;

			point p = *i->first;

			fprintf(f, "P %u %f %f %f\n", i->second, p[0], p[1], p[2]);
		}

		/*
		 * Print camera data
		 */

		for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
			et _et = align::of(i);
			pt _pt = align::projective(i);

			fprintf(f, "C %u %f %f %f %f %f %f %f\n", i,
					_pt.view_angle() * 180 / M_PI,
					_et.get_translation(1),
					_et.get_translation(0),
					_et.get_translation(2),
					_et.get_rotation(1) * 180 / M_PI,
					_et.get_rotation(0) * 180 / M_PI,
					_et.get_rotation(2) * 180 / M_PI);
		}

		fclose(f);
	}

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
	 * Test the density of the mesh for correct sampling, and split (or
	 * unsplit) triangles if necessary, continuing until no more operations
	 * can be performed.  
	 */
	static int density_test(int split) {

		int result = 0;

		/*
		 * Iterate over all frames
		 *
		 * XXX: we might not want to do this.  E.g., for splitting, we
		 * might select one frame, since the perturbation magnitude is determined
		 * by a single frame.
		 */
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {

			pt _pt = align::projective(n);
			// _pt.scale(1 / _pt.scale_2d());  // XXX: world coordinates vary with the output scale.
			ale_pos area_threshold = pow(2 * perturb, 2);

			if (split)
				result |= triangle_head[0]->split_on_camera_area_threshold(_pt, area_threshold)
					| triangle_head[1]->split_on_camera_area_threshold(_pt, area_threshold);
			else
				result |= triangle_head[0]->unsplit_on_camera_area_threshold(_pt, area_threshold)
					| triangle_head[1]->unsplit_on_camera_area_threshold(_pt, area_threshold);


		}

		return result;
	}

	static int density_test_split() {
		return density_test(1);
	}

	static int density_test_unsplit() {
		return density_test(0);
	}

	static zbuf_elem **construct_zbuffers(lod *_lod) {
		zbuf_elem **result = (zbuf_elem **) malloc(d2::image_rw::count() * sizeof(zbuf_elem *));

		for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
			pt _pt = align::projective(i);
			_pt.scale(_lod->sf / _pt.scale_2d());
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
	static int adjust_vertices(std::set<point *> *vertex_set, 
			           std::multimap<point *, triangle *> *vertex_map,
				   change_set_t *change_set_prev,
				   change_set_t *change_set_cur) {

		/*
		 * Collect bounding box data.
		 */

		std::vector<pp_set_t> bb_array(d2::image_rw::count());

		for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
			pt _pt = align::projective(f);

			for (change_set_t::iterator i = change_set_prev->begin();
					i != change_set_prev->end(); i++) {
				d2::point max = d2::point::neginf();
				d2::point min = d2::point::posinf();

				for (change_elem_t::iterator j = i->begin();
						j != i->end(); j++) {
					d2::point p = _pt.wp_scaled(*j).xy();
					max.accumulate_max(p);
					min.accumulate_min(p);
				}

				if (max[0] >= min[0]
				 && max[1] >= min[1])
					bb_array[f].insert(std::make_pair(min, max));
			}
		}

		/*
		 * Collect z-buffer data
		 */

		fprintf(stderr, "  begin creating zbuffers[%u] \n", time(NULL));
		zbuf_elem **z = construct_zbuffers(cl);
		fprintf(stderr, "  end creating zbuffers[%u] \n", time(NULL));

#if 0
		fprintf(stderr, "  begin tree traversal[%u] \n", time(NULL));
		int result = triangle_head[0]->adjust_vertices(z, cl) 
			   | triangle_head[1]->adjust_vertices(z, cl);
		fprintf(stderr, "  end tree traversal[%u] \n", time(NULL));
#else
		int result = 0;

		fprintf(stderr, "  begin free vertex set traversal[%u] \n", time(NULL));
		for (std::set<point *>::iterator i = vertex_set->begin(); i != vertex_set->end(); i++) {

			/*
			 * Check for inclusion of the original point in one of
			 * the bounding box regions.
			 */
			int orig_inclusion = 0;
			if (change_set_prev->size() > 0)
			for (unsigned int f = 0; !orig_inclusion && f < d2::image_rw::count(); f++) {
				pt _pt = align::projective(f);

				for (pp_set_t::iterator j = bb_array[f].begin(); 
						!orig_inclusion && j != bb_array[f].end(); j++) {

					point p = (**i);
					point bb_min = j->first;
					point bb_max = j->second;

					if (p[0] >= bb_min[0]
					 && p[1] >= bb_min[1]
					 && p[0] <= bb_max[0]
					 && p[1] <= bb_max[1])
						orig_inclusion = 1;
				}
			}

			/*
			 * Determine the adjustment step size according to
			 * the perturbation size.
			 */

			ale_pos allowable_max_neighbor_angle = M_PI / 1.5;
			ale_pos step = perturb / sqrt(d3::align::projective(0).w_density_scaled(**i));

			/*
			 * Perturb the position.
			 */
			
			for (int axis = 0; axis < 3; axis++)
			for (int dir = -1; dir <= 1; dir += 2) {

				// fprintf(stderr, "%p perturbing vertex %d (%d %d) [%u] \n", this, v, axis, dir, time(NULL));
				point orig = **i;
				point perturbed = orig + point::unit(axis) * step * dir;

				/*
				 * Check the clipping planes.
				 */

				// fprintf(stderr, "%p checking clipping planes [%u] \n", this, time(NULL));
				if (perturbed[2] > front_clip
				 || perturbed[2] < rear_clip)
					continue;

				/*
				 * Check view pyramid bounds
				 */

				// fprintf(stderr, "%p checking view pyramid bounds [%u] \n", this, time(NULL));
				if (!pyramid_bounds_check(perturbed))
					continue;

				/*
				 * Check for inclusion of the perturbed point in one of
				 * the bounding box regions.
				 */
				int perturbed_inclusion = 0;
				if (!orig_inclusion && change_set_prev->size() > 0)
				for (unsigned int f = 0; !perturbed_inclusion && f < d2::image_rw::count(); f++) {
					pt _pt = align::projective(f);

					for (pp_set_t::iterator j = bb_array[f].begin(); 
							!perturbed_inclusion && j != bb_array[f].end(); j++) {

						point p = perturbed;
						point bb_min = j->first;
						point bb_max = j->second;

						if (p[0] >= bb_min[0]
						 && p[1] >= bb_min[1]
						 && p[0] <= bb_max[0]
						 && p[1] <= bb_max[1])
							perturbed_inclusion = 1;
					}
				}

				/*
				 * Bounding box predicate
				 */
				if (change_set_prev->size() > 0 && !orig_inclusion && !perturbed_inclusion)
					continue;

				/*
				 * The set of triangles surrounding the vertex.
				 */

				std::set<triangle *> tset;

				for (std::multimap<point *, triangle *>::iterator t = vertex_map->lower_bound(*i);
						t != vertex_map->upper_bound(*i); t++)
					tset.insert(t->second);

				/*
				 * Eliminate from consideration any change that increases the number of forbidden 
				 * neighbor angles.
				 */

				int orig_forbidden_angles = 0;
				for (std::set<triangle *>::iterator t = tset.begin(); t != tset.end(); t++) {
					if ((*t)->max_neighbor_angle() > allowable_max_neighbor_angle)
						orig_forbidden_angles++;
				}

				**i = perturbed;
				int perturb_forbidden_angles = 0;
				for (std::set<triangle *>::iterator t = tset.begin(); t != tset.end(); t++) {
					if ((*t)->max_neighbor_angle() > allowable_max_neighbor_angle)
						perturb_forbidden_angles++;
				}
				**i = orig;
				if (perturb_forbidden_angles > orig_forbidden_angles)
					continue;

				/*
				 * Determine the set of vertices associated with the triangles
				 * surrounding the vertex being perturbed
				 */

				change_elem_t v_set;
				for (std::set<triangle *>::iterator j = tset.begin(); j != tset.end(); j++)
				for (int v = 0; v < 3; v++)
					v_set.insert(*(*j)->vertices[v]);
				v_set.insert(perturbed);

				/*
				 * Check the error
				 */

				// fprintf(stderr, "%p calculating error [%u] \n", this, time(NULL));
				ale_accum error = vertex_movement_cost(&tset, &v_set, *i, perturbed, z, cl);
				// fprintf(stderr, "%p done calculating error (%f) [%u] \n", this, error, time(NULL));
				if (error < 0) {

					fprintf(stderr, "%p (%f %f %f) (%f %f %f) [%f]\n",
						*i,
						(**i)[0],
						(**i)[1],
						(**i)[2],
						perturbed[0],
						perturbed[1],
						perturbed[2],
						2 / perturb);

					result = 1;
					**i = perturbed;
					change_set_cur->insert(v_set);
					break;
				}
			}

		}
		fprintf(stderr, "  end free vertex set traversal[%u] \n", time(NULL));
#endif

		fprintf(stderr, "  begin freeing zbuffers[%u] \n", time(NULL));
		free_zbuffers(z);
		fprintf(stderr, "  end freeing zbuffers[%u] \n", time(NULL));

		return result;
	}

public:
	static void load_model(const char *name) {
		load_model_name = name;
	}

	static void save_model(const char *name) {
		save_model_name = name;
	}

	static void mpl_absolute(ale_real value) {
		mpl_value = value;
		mpl_type = 0;
	}

	static void mpu_absolute(ale_real value) {
		mpu_value = value;
		mpu_type = 0;
	}

	static void mpl_percent(ale_real value) {
		mpl_value = value;
		mpl_type = 1;
	}

	static void mpu_percent(ale_real value) {
		mpu_value = value;
		mpu_type = 1;
	}

	static void fc(ale_pos fc) {
		front_clip = fc;
	}

	static void rc(ale_pos rc) {
		rear_clip = rc;
	}

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
		 * Rear clip value of 0 is converted to infinity.
		 */

		if (rear_clip == 0) {
			ale_pos one = +1;
			ale_pos zero = +0;

			rear_clip = one / zero;
			assert(isinf(rear_clip) == +1);
		}

		/*
		 * Scale and translate clipping plane depths.
		 */

		ale_pos cp_scalar = d3::align::projective(0).wc(point(0, 0, 0))[2];

		front_clip = front_clip * cp_scalar - cp_scalar;
		rear_clip = rear_clip * cp_scalar - cp_scalar;

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

		bl = cl;

		/*
		 * If there is a file to load scene data from, then use it.
		 */

		if (load_model_name) {
			read_model_file();
			return;
		}

		/*
		 * Otherwise, initialize the model explicitly.
		 */

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

		triangle::fix_vertex(triangle_head[0]->vertices[0]);
		triangle::fix_vertex(triangle_head[0]->vertices[1]);
		triangle::fix_vertex(triangle_head[0]->vertices[2]);

		triangle::fix_vertex(triangle_head[1]->vertices[0]);
		triangle::fix_vertex(triangle_head[1]->vertices[1]);
		triangle::fix_vertex(triangle_head[1]->vertices[2]);
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

	/*
	 * Run a single iteration of the spatial_info update cycle.
	 */
	static void spatial_info_update() {
		/*
		 * Iterate through each frame.
		 */
		for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
			/*
			 * Get transformation data.
			 */
			pt _pt = align::projective(f);

			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));

			/*
			 * Get color data for the frames.
			 */
			const d2::image *im = d2::image_rw::open(f);

			/*
			 * Allocate an image for storing encounter probabilities.
			 */
			d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()), 
					(int) floor(_pt.scaled_width()), 3);

			/*
			 * 
			 */

			/*
			 * Visit spatial_info nodes in order.
			 */

			space_iterate si(_pt);

			do {

				space_traverse st = si.get();

				/*
				 * XXX: This could be more efficient, perhaps.
				 */

				if (spatial_info_map.count(st.get_space()) == 0)
					continue;

				spatial_info *sn = &spatial_info_map[st.get_space()];

				/*
				 * Get information on the subspace.
				 */

				d2::pixel color = sn->get_color();
				ale_real occupancy = sn->get_occupancy();

				/*
				 * Determine the view-local bounding box for the
				 * subspace.
				 */

				point min = d2::point(_pt.scaled_height(), _pt.scaled_width());
				point max = d2::point(0, 0);

				point wbb[2] = { st.get_min(), st.get_max() };

				for (int x = 0; x < 2; x++)
				for (int y = 0; y < 2; y++)
				for (int z = 0; z < 2; z++) {
					point p = _pt.wp_scaled(point(wbb[x][0], wbb[y][1], wbb[z][2]));

					if (p[0] < min[0])
						min[0] = p[0];
					if (p[0] > max[0])
						max[0] = p[0];
					if (p[1] < min[1])
						min[1] = p[1];
					if (p[1] > max[1])
						max[1] = p[1];
				}

				/*
				 * Clip bounding box to image extents.
				 */

				if (min[0] < 0)
					min[0] = 0;
				if (min[1] < 0)
					min[1] = 0;
				if (max[0] > _pt.scaled_height())
					max[0] = _pt.scaled_height();
				if (max[1] > _pt.scaled_width())
					max[1] = _pt.scaled_width();

				/*
				 * Iterate over pixels in the bounding box,
				 * adding new data to the subspace.  XXX:
				 * assume for now that all pixels in the
				 * bounding box intersect the subspace.
				 */

				for (unsigned int i = (unsigned int) floor(min[0]); i < (unsigned int) ceil(max[0]); i++)
				for (unsigned int j = (unsigned int) floor(min[1]); j < (unsigned int) ceil(max[1]); j++) {
					/*
					 * Determine the probability of encounter.
					 */

					d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j)) * occupancy;
					d2::pixel colordiff = color - im->get_pixel(i, j);
					ale_real exp_scale = 256 * 256;

					/*
					 * Update subspace.
					 */

					sn->accumulate_color_1(color, encounter);
					d2::pixel channel_occ = pexp(-colordiff * colordiff * exp_scale);
					for (int k = 0; k < 3; k++)
						sn->accumulate_occupancy_1(channel_occ[k], encounter[k]);

					weights->pix(i, j) += encounter;
				}

			} while (si.next());

			d2::image_rw::close(f);

			delete weights;
		}

		/*
		 * Update all spatial_info structures.
		 */
		for (std::map<space *,spatial_info>::iterator i = spatial_info_map.begin(); i != spatial_info_map.end(); i++) {
			i->second.update_color();
			i->second.update_occupancy();
		}
	}

	/*
	 * Generate an image from a specified view.
	 */
	static const d2::image *view(pt _pt, int n = -1) {
		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		d2::image *im = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

#if 1
		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space_iterate si(_pt);

		do {
			space_traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_space()) == 0)
				continue;

			spatial_info sn = spatial_info_map[st.get_space()];

			/*
			 * Get information on the subspace.
			 */

			d2::pixel color = sn.get_color();
			ale_real occupancy = sn.get_occupancy();

			/*
			 * Determine the view-local bounding box for the
			 * subspace.
			 */

			point min = d2::point(_pt.scaled_height(), _pt.scaled_width());
			point max = d2::point(0, 0);

			point wbb[2] = { st.get_min(), st.get_max() };

			for (int x = 0; x < 2; x++)
			for (int y = 0; y < 2; y++)
			for (int z = 0; z < 2; z++) {
				point p = _pt.wp_scaled(point(wbb[x][0], wbb[y][1], wbb[z][2]));

				if (p[0] < min[0])
					min[0] = p[0];
				if (p[0] > max[0])
					max[0] = p[0];
				if (p[1] < min[1])
					min[1] = p[1];
				if (p[1] > max[1])
					max[1] = p[1];
			}

			/*
			 * Clip bounding box to image extents.
			 */

			if (min[0] < 0)
				min[0] = 0;
			if (min[1] < 0)
				min[1] = 0;
			if (max[0] > _pt.scaled_height())
				max[0] = _pt.scaled_height();
			if (max[1] > _pt.scaled_width())
				max[1] = _pt.scaled_width();

			/*
			 * Iterate over pixels in the bounding box, finding
			 * pixels that intersect the subspace.  XXX: assume
			 * for now that all pixels in the bounding box
			 * intersect the subspace.
			 */

			for (unsigned int i = (unsigned int) floor(min[0]); i < (unsigned int) ceil(max[0]); i++)
			for (unsigned int j = (unsigned int) floor(min[1]); j < (unsigned int) ceil(max[1]); j++) {
				/*
				 * Determine the probability of encounter.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j)) * occupancy;

				/*
				 * Update images.
				 */

				weights->pix(i, j) += encounter;
				im->pix(i, j)      += encounter * color;
			}

		} while (si.next());

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			im->pix(i, j) /= weights->pix(i, j);
		}


#else
		/*
		 * Use hierarchical triangle data.
		 */

		zbuf_elem *zbuf = init_zbuf(_pt);

		zbuffer(_pt, zbuf, triangle_head[0]);
		zbuffer(_pt, zbuf, triangle_head[1]);

		zbuf_elem **zbsu = construct_zbuffers(bl);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			d2::pixel val;
			ale_real div = 0;

			/*
			 * Approximates filter fine:box:1
			 */

			for (unsigned int f = 0; f < d2::image_rw::count(); f++) {

				pt _ptf = align::projective(f);
				_ptf.scale(1 / _ptf.scale_2d());

				point bounds[2];

				bounds[0] = bounds[1] = point::undefined();

				bounds[0][2] = -1;
				bounds[1][2] = -1;

				for (ale_pos ii = -0.5; ii <= 0.5; ii += 1)
				for (ale_pos jj = -0.5; jj <= 0.5; jj += 1) {

					point p = frame_to_frame(d2::point(i + ii, j + jj), _pt, _ptf, zbuf, zbsu[f]);

					if (p[0] < bounds[0][0] || !finite(bounds[0][0]))
						bounds[0][0] = p[0];
					if (p[1] < bounds[0][1] || !finite(bounds[0][1]))
						bounds[0][1] = p[1];
					if (p[0] > bounds[1][0] || !finite(bounds[1][0]))
						bounds[1][0] = p[0];
					if (p[1] > bounds[1][1] || !finite(bounds[1][1]))
						bounds[1][1] = p[1];
				}

				if (!bounds[0].defined() || !bounds[1].defined())
					continue;

				for (int ii = (int) ceil(bounds[0][0]); ii <= (int) floor(bounds[1][0]); ii++)
				for (int jj = (int) ceil(bounds[0][1]); jj <= (int) floor(bounds[1][1]); jj++) {
					if ((n > 0 && f == (unsigned int) n) 
					 || frame_to_frame(d2::point(ii, jj), _ptf, _pt, zbsu[f], zbuf).defined()) {
						val += bl->reference[f]->get_pixel(ii, jj);
						div += 1;
					}
				}
			}

			if (div != 0) {
				im->pix(i, j) = val / div;
				continue;
			}

			/*
			 * Approximates filter triangle:2
			 */

			for (unsigned int f = 0; f < d2::image_rw::count(); f++) {

				pt _ptf = align::projective(f);
				_ptf.scale(1 / _ptf.scale_2d());

				if (n > 0 && f == (unsigned int) n) {
					point p = _ptf.wp_scaled(_pt.pw_scaled(point(i, j, -1)));

					if (!bl->reference[f]->in_bounds(p.xy()))
						continue;

					val += bl->reference[f]->get_bl(p.xy());
					div += 1;

					continue;
				}

				point p = frame_to_frame(d2::point(i, j), _pt, _ptf, zbuf, zbsu[f]);

				if (!p.defined())
					continue;

				d2::pixel v = bl->reference[f]->get_bl(p.xy());

				val += v;
				div += 1;
			}

			im->pix(i, j) = val / div;
		}

		free_zbuffers(zbsu);

		delete[] zbuf;
#endif

		return im;
	}

	static const d2::image *view(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return view(_pt, n);
	}

	static const d2::image *depth(pt _pt, int n = -1) {
		assert (n < 0 || (unsigned int) n < d2::image_rw::count());

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		d2::image *im = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

#if 1
		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space_iterate si(_pt);

		do {
			space_traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_space()) == 0)
				continue;

			spatial_info sn = spatial_info_map[st.get_space()];

			/*
			 * Get information on the subspace.
			 */

			d2::pixel color = sn.get_color();
			ale_real occupancy = sn.get_occupancy();

			/*
			 * Determine the view-local bounding box for the
			 * subspace.
			 */

			point min = d2::point(_pt.scaled_height(), _pt.scaled_width());
			point max = d2::point(0, 0);

			point wbb[2] = { st.get_min(), st.get_max() };

			for (int x = 0; x < 2; x++)
			for (int y = 0; y < 2; y++)
			for (int z = 0; z < 2; z++) {
				point p = _pt.wp_scaled(point(wbb[x][0], wbb[y][1], wbb[z][2]));

				if (p[0] < min[0])
					min[0] = p[0];
				if (p[0] > max[0])
					max[0] = p[0];
				if (p[1] < min[1])
					min[1] = p[1];
				if (p[1] > max[1])
					max[1] = p[1];
			}

			/*
			 * Clip bounding box to image extents.
			 */

			if (min[0] < 0)
				min[0] = 0;
			if (min[1] < 0)
				min[1] = 0;
			if (max[0] > _pt.scaled_height())
				max[0] = _pt.scaled_height();
			if (max[1] > _pt.scaled_width())
				max[1] = _pt.scaled_width();

			/*
			 * Iterate over pixels in the bounding box, finding
			 * pixels that intersect the subspace.  XXX: assume
			 * for now that all pixels in the bounding box
			 * intersect the subspace.
			 */

			for (unsigned int i = (unsigned int) floor(min[0]); i < (unsigned int) ceil(max[0]); i++)
			for (unsigned int j = (unsigned int) floor(min[1]); j < (unsigned int) ceil(max[1]); j++) {
				/*
				 * Determine the probability of encounter.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j)) * occupancy;

				/*
				 * Update images.
				 */

				weights->pix(i, j) += encounter;
				im->pix(i, j)      += encounter * _pt.wp_scaled(st.get_min())[2];
			}

		} while (si.next());

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			im->pix(i, j) /= weights->pix(i, j);
		}

#else

		/*
		 * Use triangle model data.
		 */

		zbuf_elem *zbuf = init_zbuf(_pt);

		zbuffer(_pt, zbuf, triangle_head[0]);
		zbuffer(_pt, zbuf, triangle_head[1]);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			triangle *t = zbuf[i * im->width() + j].nearest(_pt, i, j);

			if (!t)
				continue;
#if 0
			point r = _pt.pc_scaled(point(i, j, -1));
			point vertices[3];

			for (int v = 0; v < 3; v++)
				vertices[v] = _pt.wc(*t->vertices[v]);

			point intersect = rt_intersect(r, vertices);

			im->pix(i, j) = d2::pixel(1, 1, 1) * t->area() * (intersect[0] + intersect[1]);
#elif 1
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

#endif

		return im;
	}

	static const d2::image *depth(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return depth(_pt, n);
	}

	/*
	 * Add specified control points to the model
	 */
	static void add_control_points() {

		/*
		 * Don't add control points to loaded scenes.
		 */
		if (load_model_name)
			return;

		for (unsigned int i = 0; i < cpf::count(); i++) {

			point control_point = cpf::get(i);

			if (!control_point.defined())
				continue;

			if (control_point[2] > front_clip
			 || control_point[2] < rear_clip)
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

		/*
		 * Don't relax loaded scenes
		 */
		if (load_model_name)
			return;

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
	 * Initialize space and identify regions of interest for the adaptive
	 * subspace model.
	 */
	static void make_space() {

		/*
		 * Initialize root space.
		 */

		root_space = new space;

		/*
		 * Subdivide space to resolve intensity matches between pairs
		 * of frames.
		 */

		for (unsigned int f1 = 0; f1 < d2::image_rw::count(); f1++)
		for (unsigned int f2 = 0; f2 < d2::image_rw::count(); f2++) {
			if (f1 == f2)
				continue;

			const d2::image *if1 = d2::image_rw::open(f1);
			const d2::image *if2 = d2::image_rw::open(f2);

			pt _pt1 = align::projective(f1);
			pt _pt2 = align::projective(f2);

			/*
			 * Iterate over all points in the primary frame.
			 */

			for (unsigned int i = 0; i < if1->height(); i++)
			for (unsigned int j = 0; j < if1->width();  j++) {

				/*
				 * Get the pixel color in the primary frame
				 */

				d2::pixel color_primary = if1->get_pixel(i, j);

				/*
				 * Map two depths to the secondary frame.
				 */

				point p1 = _pt2.wp_unscaled(_pt1.pw_unscaled(point(i, j,  1000)));
				point p2 = _pt2.wp_unscaled(_pt1.pw_unscaled(point(i, j, -1000)));

				/*
				 * For cases where the mapped points define a
				 * line and where points on the line fall
				 * within the defined area of the frame,
				 * determine the starting point for inspection.
				 * In other cases, continue to the next pixel.
				 */

				ale_pos diff_i = p2[0] - p1[0];
				ale_pos diff_j = p2[1] - p1[1];
				ale_pos slope = diff_j / diff_i;

				if (isnan(slope))
					continue;

				ale_pos top_intersect = p1[1] - p1[0] * slope;
				ale_pos lef_intersect = p1[0] - p1[1] / slope;
				ale_pos rig_intersect = p1[0] - (p1[1] - if2->width() + 2) / slope;
				ale_pos sp_i, sp_j;

				if (finite(slope) && top_intersect >= 0 && top_intersect < if2->width() - 1) {
					sp_i = 0;
					sp_j = top_intersect;
				} else if (slope > 0 && lef_intersect >= 0 && lef_intersect < if2->height() - 1) {
					sp_i = lef_intersect;
					sp_j = 0;
				} else if (slope < 0 && rig_intersect >= 0 && rig_intersect < if2->height() - 1) {
					sp_i = rig_intersect;
					sp_j = if2->width() - 2;
				} else 
					continue;

				/*
				 * Determine increment values for examining
				 * point, ensuring that incr_i is always
				 * positive.
				 */

				ale_pos incr_i, incr_j;

				if (fabs(diff_i) > fabs(diff_j)) {
					incr_i = 1;
					incr_j = slope;
				} else if (slope > 0) {
					incr_i = 1 / slope;
					incr_j = 1;
				} else {
					incr_i = -1 / slope;
					incr_j = -1;
				}

				/*
				 * Examine regions near the projected line.
				 */

				for (ale_pos ii = sp_i, jj = sp_j; 
					ii < if2->height() - 1 && jj < if2->width() - 1 && ii >= 0 && jj >= 0; 
					ii += incr_i, jj += incr_j) {

					/*
					 * Check for higher, lower, and nearby points.
					 *
					 * 	Red   = 2^0
					 * 	Green = 2^1
					 * 	Blue  = 2^2
					 */

					int higher = 0, lower = 0, nearby = 0;

					for (int iii = 0; iii < 2; iii++)
					for (int jjj = 0; jjj < 2; jjj++) {
						d2::pixel p = if2->get_pixel((int) floor(ii) + iii, (int) floor(jj) + jjj);

						for (int k = 0; k < 3; k++) {
							int bitmask = (int) pow(2, k);

							if (p[k] > color_primary[k])
								higher |= bitmask;
							if (p[k] < color_primary[k])
								lower  |= bitmask;
							if (fabs(p[k] - color_primary[k]) < nearness)
								nearby |= bitmask;
						}
					}

					/*
					 * If this is not a region of interest,
					 * then continue.
					 */

					if (((higher & lower) | nearby) != 0x7)
						continue;

					/*
					 * Create an orthonormal basis to
					 * determine line intersection.
					 */

					point bp0 = _pt1.pw_unscaled(point(i, j, 0));
					point bp1 = _pt1.pw_unscaled(point(i, j, 1));
					point bp2 = _pt2.pw_unscaled(point(ii, jj, 0));

					point b0  = (bp1 - bp0).normalize();
					point b1n = bp2 - bp0;
					point b1  = (b1n - b1n.dproduct(b0) * b0).normalize();

					/*
					 * Select a fourth point to define a second line.
					 */

					point p3  = _pt2.pw_unscaled(point(ii, jj, 1));

					/*
					 * Representation in the new basis.
					 */

					d2::point nbp0 = d2::point(bp0.dproduct(b0), bp0.dproduct(b1));
					d2::point nbp1 = d2::point(bp1.dproduct(b0), bp1.dproduct(b1));
					d2::point nbp2 = d2::point(bp2.dproduct(b0), bp2.dproduct(b1));
					d2::point np3  = d2::point( p3.dproduct(b0),  p3.dproduct(b1));

					/*
					 * Determine intersection of line
					 * (nbp0, nbp1), which is parallel to
					 * b0, with line (nbp2, np3).
					 */

					/*
					 * XXX: a stronger check would be
					 * better here, e.g., involving the
					 * ratio (np3[0] - nbp2[0]) / (np3[1] -
					 * nbp2[1]).  Also, acceptance of these
					 * cases is probably better than
					 * rejection.
					 */

					if (np3[1] - nbp2[1] == 0)
						continue;

					d2::point intersection = d2::point(nbp2[0] 
							+ (nbp0[1] - nbp2[1]) * (np3[0] - nbp2[0]) / (np3[1] - nbp2[1]),
							nbp0[1]);

					/*
					 * Map the intersection back to the world
					 * basis.
					 */

					point iw = intersection[0] * b0 + intersection[1] * b1;

					/*
					 * Reject interesection points behind a
					 * camera.
					 */

					if (_pt1.wc(iw)[2] >= 0 || _pt2.wc(iw)[2] >= 0)
						continue;

					/*
					 * Refine space around the intersection point.
					 */

					space_traverse st = space_traverse::root();

					if (!st.includes(iw))
						continue;

					for(;;) {

						point frame_min[2] = { point::posinf(), point::posinf() },
						      frame_max[2] = { point::neginf(), point::neginf() };

						point p[2] = { st.get_min(), st.get_max() };

						for (int ibit = 0; ibit < 2; ibit++)
						for (int jbit = 0; jbit < 2; jbit++)
						for (int kbit = 0; kbit < 2; kbit++) {
							point pp = point(p[ibit][0], p[jbit][1], p[kbit][2]);

							point ppp[2] = {_pt1.wp_unscaled(pp), _pt2.wp_unscaled(pp)};

							for (int f = 0; f < 2; f++)
							for (int d = 0; d < 3; d++) {
								if (ppp[f][d] < frame_min[f][d])
									frame_min[f][d] = ppp[f][d];
								if (ppp[f][d] > frame_max[f][d])
									frame_max[f][d] = ppp[f][d];
							}
						}

						if (frame_min[0].lengthtosq(frame_max[0]) < 2
						 && frame_min[1].lengthtosq(frame_max[1]) < 2)
							break;

						if (st.positive().includes(iw))
							st = st.positive();
						else if (st.negative().includes(iw))
							st = st.negative();
						else
							assert(0);
					}

					/*
					 * Associate refined space with a
					 * spatial info structure.
					 */

					spatial_info_map[st.get_space()];
				}
			}

			/*
			 * This ordering should ensure that image f1 is cached.
			 */

			d2::image_rw::close(f2);
			d2::image_rw::close(f1);
		}
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
	static void reduce_cost_to_search_depth(const char *d_out[], const char *v_out[], 
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt,
			d2::exposure *exp_out, int inc_bit) {

#if 1
		/*
		 * Subspace model
		 */

		for (int i = 0; i < 10; i++)
			spatial_info_update();

#else

		/*
		 * Triangle model
		 */

		int improved = 1;
		int count = 0;

		/*
		 * Make all perturbation bounds absolute.
		 */

		ale_pos relative_unit = 0;

		for (int v = 0; v < 3; v++) {
			ale_pos distance = triangle_head[0]->vertices[v]->lengthto(*triangle_head[1]->vertices[v]);

			if (distance > relative_unit)
				relative_unit = distance;
		}

		if (mpl_type == 1) {
			mpl_type = 0;
			if (mpl_value > 25)
				mpl_value = 25;
			mpl_value = (mpl_value / 100) * relative_unit;
		}

		if (mpu_type == 1) {
			mpu_type = 0;
			if (mpu_value > 25)
				mpu_value = 25;
			mpu_value = (mpu_value / 100) * relative_unit;
		}

		/*
		 * To start, use the largest perturbation value.
		 */

		perturb = mpu_value;

#if 0
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

		return;
#endif

		/*
		 * Change sets.
		 */
		change_set_t change_set[2];
		int csi = 0;

		/*
		 * Reduce level-of-detail.
		 */

		while (perturb > 2 / cl->sf && reduce_lod());

		/*
		 * Perform cost-reduction.
		 */

		while (perturb >= mpl_value) {

			fprintf(stderr, "begin density test split[%u] \n", time(NULL));

			if (density_test_split())
				continue;

			fprintf(stderr, "end density test split[%u] \n", time(NULL));

			/*
			 * Gather vertex data
			 */

			std::set<point *> vertex_set;
			std::set<point *> fixed_set;
			std::multimap<point *, triangle *> vertex_map;

			fprintf(stderr, "begin collecting vertex data[%u]\n", time(NULL));
			triangle_head[0]->collect_vertex_data(&vertex_set, &vertex_map, &fixed_set);
			triangle_head[1]->collect_vertex_data(&vertex_set, &vertex_map, &fixed_set);
			fprintf(stderr, "end collecting vertex data[%u]\n", time(NULL));

			fprintf(stderr, "Number of vertices: %u\n", vertex_set.size());

			std::set<point *>::iterator aa = vertex_set.begin();
			std::set<point *>::iterator bb = fixed_set.begin();

			while (aa != vertex_set.end() && bb != fixed_set.end()) {
				if (*aa == *bb) {
					fprintf(stderr, "Error: the fixed and free vertex sets have elements in common.\n");
					assert(0);
				} else if (*aa < *bb) {
					aa++;
				} else if (*bb < *aa) {
					bb++;
				} else
					assert(0);
			}

			fprintf(stderr, "Finished checking for common elements.\n");

			/*
			 * Limit the number of iterations at a given LOD. 
			 */

			if (!improved || count > 40) {

				/*
				 * Write output incrementally, if desired.
				 */

				fprintf(stderr, "begin output[%u] \n", time(NULL));

				if (inc_bit) {

					write_model_file();

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

					for (std::map<const char *, pt>::iterator i = d3_depth_pt->begin();
							i != d3_depth_pt->end(); i++) {

						const d2::image *im = depth(i->second);
						d2::image_rw::write_image(i->first, im, exp_out, 1, 1);
						delete im;
					}

					for (std::map<const char *, pt>::iterator i = d3_output_pt->begin();
							i != d3_output_pt->end(); i++) {

						const d2::image *im = view(i->second);
						d2::image_rw::write_image(i->first, im, exp_out);
						delete im;
					}
				}
				
				fprintf(stderr, "end output[%u] \n", time(NULL));

				fprintf(stderr, ".");

				perturb /= 2;

				if (cl->next)
					increase_lod();

				/*
				 * Clear change sets.
				 */

				change_set[0].clear();
				change_set[1].clear();

				count = 0;
				improved = 1;
				continue;
			}

			count++;
			improved = 0;

			/*
			 * Try improving the result by moving existing vertices.
			 */

			fprintf(stderr, "begin adjustment[%u] \n", time(NULL));

			improved |= adjust_vertices(&vertex_set, &vertex_map, &change_set[1 - csi], &change_set[csi]);

			/*
			 * Update change sets.
			 */

			csi = 1 - csi;
			change_set[csi].clear();

			fprintf(stderr, "end adjustment[%u] \n", time(NULL));
		}

#endif
		write_model_file();
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
