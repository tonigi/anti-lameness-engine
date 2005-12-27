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
	 * Occupancy attenuation
	 */

	static double occ_att;

	/*
	 * Normalization of output by weight
	 */

	static int normalize_weights;

	/*
	 * Falloff exponent
	 */

	static double falloff_exponent;

	/*
	 * Third-camera error multiplier
	 */
	static double tc_multiplier;

	/*
	 * Occupancy update iterations
	 */
	static unsigned int ou_iterations;

	/*
	 * Pairwise ambiguity
	 */
	static unsigned int pairwise_ambiguity;

	/*
	 * Pairwise comparisons
	 */
	static const char *pairwise_comparisons;

	/*
	 * 3D Post-exclusion
	 */
	static int d3px_count;
	static double *d3px_parameters;

	/*
	 * Nearness threshold
	 */
	static const ale_real nearness;

	/*
	 * Encounter threshold for defined pixels.
	 */
	static double encounter_threshold;

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
	 * Structure to hold a subdivisible region of space.
	 */

	struct space {
		struct space *positive;
		struct space *negative;

		space() {
			positive = NULL;
			negative = NULL;
		}
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

	public:

		int get_next_split() {

			/*
			 * Double-infinite case.
			 */

			for (int d = 0; d < 3; d++)
			if (isinf(max[d]) && isinf(min[d]))
				return d;

			/*
			 * Finite or single-infinite case
			 */

			if (max[0] - min[0] >= max[1] - min[1]
			 && (max[0] >= max[1] || !isinf(min[1]))
			 && (min[0] <= min[1] || !isinf(max[1]))
			 && max[0] - min[0] >= max[2] - max[2]
			 && (max[0] >= max[2] || !isinf(min[2]))
			 && (min[0] <= min[2] || !isinf(max[2])))
				return 0;

			if (max[1] - min[1] > max[2] - min[2]
			 && (max[1] >= max[2] || !isinf(min[2]))
			 && (min[1] <= min[2] || !isinf(max[2])))
				return 1;

			return 2;
		}

		static space_traverse root() {

			space_traverse result;

			result.current = root_space;
			result.min = point::neginf();
			result.max = point::posinf();

			assert(result.current);

			return result;
		}

		ale_pos split_coordinate(int d) {
			if (isinf(max[d]) && isinf(min[d]))
				return 0;

			if (isinf(max[d]))
				return tan((atan(min[d]) + M_PI/2) / 2);

			if (isinf(min[d]))
				return tan((atan(max[d]) - M_PI/2) / 2);

			return (min[d] + max[d]) / 2;
		}

		ale_pos split_coordinate() {
			int next_split = get_next_split();
			return split_coordinate(next_split);
		}

		int precision_wall() {
			int next_split = get_next_split();
			ale_pos split_point = split_coordinate(next_split);

			assert(split_point <= max[next_split]);
			assert(split_point >= min[next_split]);
			
			if (split_point == min[next_split] || split_point == max[next_split]) 
				return 1;

			return 0;
		}

		space_traverse positive() {

			assert(current);

			int next_split = get_next_split();

			if (current->positive == NULL) {
				current->positive = new space;
				total_divisions++;
			}
			
			space_traverse result;

			result.current = current->positive;
			result.min = min;
			result.max = max;

			result.min[next_split] = split_coordinate(next_split);

			assert(result.current);

			return result;
		}

		space_traverse negative() {

			assert(current);

			int next_split = get_next_split();

			if (current->negative == NULL) {
				current->negative = new space;
				total_divisions++;
			}
			
			space_traverse result;

			result.current = current->negative;
			result.min = min;
			result.max = max;

			result.max[next_split] = split_coordinate(next_split);

			assert(result.current);

			return result;
		}

		point get_min() {
			return min;
		}

		point get_max() {
			return max;
		}

		/*
		 * Get bounding box for projection onto a plane.
		 */

		void get_view_local_bb(pt _pt, point bb[2]) {

			point min = point::posinf();
			point max = point::neginf();

			point wbb[2] = { get_min(), get_max() };


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
				if (p[2] < min[2])
					min[2] = p[2];
				if (p[2] > max[2])
					max[2] = p[2];
			}

			/*
			 * Clip bounding box to image extents.
			 */

			if (min[0] < 0)
				min[0] = 0;
			if (min[1] < 0)
				min[1] = 0;
			if (max[0] > _pt.scaled_height() - 1)
				max[0] = _pt.scaled_height() - 1;
			if (max[1] > _pt.scaled_width() - 1)
				max[1] = _pt.scaled_width() - 1;

			bb[0] = min;
			bb[1] = max;
		}

		int includes(point p) {

			for (int d = 0; d < 3; d++) {
				if (p[d] > max[d])
					return 0;
				if (p[d] < min[d])
					return 0;
				if (isnan(p[d]))
					return 0;
			}

			return 1;
		}

		space *get_space() {
			assert(current);
			return current;
		}

	};

	/*
	 * Class to iterate through subspaces based on proximity to a camera.
	 */

	class space_iterate {
		std::stack<space_traverse> space_stack;
		point camera_origin;

		space_iterate(point co, space_traverse top) {
			camera_origin = co;
			space_stack.push(top);
		}

	public:
		space_iterate(pt _pt, space_traverse top = space_traverse::root()) {
			camera_origin = _pt.cw(point(0, 0, 0));
			space_stack.push(top);
		}

		int next() {
			if (space_stack.empty())
				return 0;

			space_traverse st = space_stack.top();

			int d = st.get_next_split();

			ale_pos split_coordinate = st.split_coordinate();

			space *n = st.get_space()->negative;
			space *p = st.get_space()->positive;

			if (camera_origin[d] > split_coordinate) {
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

			return (!space_stack.empty());
		}

		space_iterate cleave() {
			assert (!space_stack.empty());

			space_iterate result(camera_origin, space_stack.top());
			
			space_stack.pop();

			return result;
		}

		int done() {
			return space_stack.empty();
		}

		space_traverse get() {
			assert (!space_stack.empty());
			return space_stack.top();
		}
	};

	/*
	 * List for calculating weighted median.
	 */
	class wml {
		ale_real *data;
		unsigned int size;
		unsigned int used;

		ale_real &_w(unsigned int i) {
			assert(i <= used);
			return data[i * 2];
		}

		ale_real &_d(unsigned int i) {
			assert(i <= used);
			return data[i * 2 + 1];
		}

		void increase_capacity() {

			if (size > 0)
				size *= 2;
			else
				size = 1;

			data = (ale_real *) realloc(data, sizeof(ale_real) * 2 * (size * 1));

			assert(data);
			assert (size > used);

			if (!data) {
				fprintf(stderr, "Unable to allocate %d bytes of memory\n",
						sizeof(ale_real) * 2 * (size * 1));
				exit(1);
			} 
		}

		void insert_weight(unsigned int i, ale_real v, ale_real w) {
			assert(used < size);
			assert(used >= i);
			for (unsigned int j = used; j > i; j--) {
				_w(j) = _w(j - 1);
				_d(j) = _d(j - 1);
			}

			_w(i) = w;
			_d(i) = v;

			used++;
		}

	public:

		unsigned int get_size() {
			return size;
		}

		unsigned int get_used() {
			return used;
		}

		void print_info() {
			fprintf(stderr, "[st %p sz %u el", this, size);
			for (unsigned int i = 0; i < used; i++)
				fprintf(stderr, " (%f, %f)", _d(i), _w(i));
			fprintf(stderr, "]\n");
		}

		void clear() {
			used = 0;
		}

		void insert_weight(ale_real v, ale_real w) {
			for (unsigned int i = 0; i < used; i++) {
				if (_d(i) == v) {
					_w(i) += w;
					return;
				}
				if (_d(i) > v) {
					if (used == size)
						increase_capacity();
					insert_weight(i, v, w);
					return;
				}
			}
			if (used == size)
				increase_capacity();
			insert_weight(used, v, w);
		}

		/*
		 * Finds the median at half-weight, or between half-weight
		 * and zero-weight, depending on the attenuation value.
		 */

		ale_real find_median(double attenuation = 0) {

			assert(attenuation >= 0);
			assert(attenuation <= 1);

			ale_real zero1 = 0;
			ale_real zero2 = 0;
			ale_real undefined = zero1 / zero2;

			ale_accum weight_sum = 0;

			for (unsigned int i = 0; i < used; i++)
				weight_sum += _w(i);

//			if (weight_sum == 0)
//				return undefined;

			if (used == 0 || used == 1)
				return undefined;

			if (weight_sum == 0) {
				ale_accum data_sum = 0;
				for (unsigned int i = 0; i < used; i++)
					data_sum += _d(i);
				return data_sum / used;
			}
					

			ale_accum midpoint = weight_sum * (0.5 - 0.5 * attenuation);

			ale_accum weight_sum_2 = 0;

			for (unsigned int i = 0; i < used && weight_sum_2 < midpoint; i++) {
				weight_sum_2 += _w(i);

				if (weight_sum_2 > midpoint) {
					return _d(i);
				} else if (weight_sum_2 == midpoint) {
					assert (i + 1 < used);
					return (_d(i) + _d(i + 1)) / 2;
				} 
			}

			return undefined;
		}

		wml(int initial_size = 0) {

//			if (initial_size == 0) {
//				initial_size = (int) (d2::image_rw::count() * 1.5);
//			}

			size = initial_size;
			used = 0;

			if (size > 0) {
				data = (ale_real *) malloc(size * sizeof(ale_real) * 2);
				assert(data);
			} else {
				data = NULL;
			}
		}

		/*
		 * copy constructor.  This is required to avoid undesired frees.
		 */

		wml(const wml &w) {
			size = w.size;
			used = w.used;
			data = (ale_real *) malloc(size * sizeof(ale_real) * 2);
			assert(data);

			memcpy(data, w.data, size * sizeof(ale_real) * 2);
		}

		~wml() {
			free(data);
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
		wml color_weights_1[3];
		wml color_weights_2[3];

		/*
		 * Current color.
		 */
		d2::pixel color;

		/*
		 * Map occupancy value --> weight.
		 */
		wml occupancy_weights_1;
		wml occupancy_weights_2;

		/*
		 * Current occupancy value.
		 */
		ale_real occupancy;

		/*
		 * pocc/socc density
		 */

		unsigned int pocc_density;
		unsigned int socc_density;

		/*
		 * Insert a weight into a list.
		 */
		void insert_weight(wml *m, ale_real v, ale_real w) {
			m->insert_weight(v, w);
		}

		/*
		 * Find the median of a weighted list.  Uses NaN for undefined.
		 */
		ale_real find_median(wml *m, double attenuation = 0) {
			return m->find_median(attenuation);
		}

	public:
		/*
		 * Constructor.
		 */
		spatial_info() {
			color = d2::pixel::zero();
			occupancy = 0;
			pocc_density = 0;
			socc_density = 0;
		}

		/*
		 * Accumulate color; primary data set.
		 */
		void accumulate_color_1(int f, int i, int j, d2::pixel color, d2::pixel weight) {
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
		void accumulate_occupancy_1(int f, int i, int j, ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_1, occupancy, weight);
		}

		/*
		 * Accumulate occupancy; secondary data set.
		 */
		void accumulate_occupancy_2(ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_2, occupancy, weight);
			
			if (occupancy == 0 || occupancy_weights_2.get_size() < 96)
				return;

			// fprintf(stderr, "%p updated socc with: %f %f\n", this, occupancy, weight);
			// occupancy_weights_2.print_info();
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
			ale_real o = find_median(&occupancy_weights_1, occ_att);
			if (isnan(o))
				o = find_median(&occupancy_weights_2, occ_att);
			if (isnan(o))
				o = 0;

			occupancy = o;

			pocc_density = occupancy_weights_1.get_used();
			socc_density = occupancy_weights_2.get_used();

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

		/*
		 * Get primary color density.
		 */

		unsigned int get_pocc_density() {
			return pocc_density;
		}

		unsigned int get_socc_density() {
			return socc_density;
		}
	};

	/*
	 * DEBUG: This is a debugging variable, added for the sake of
	 * tracking particular occupancy registers, when such tracking is
	 * desired.
	 */

//	static spatial_info *tracked_space;

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

public:
	/*
	 * Debugging variables.
	 */

	static unsigned long total_ambiguity;
	static unsigned long total_pixels;
	static unsigned long total_divisions;
	static unsigned long total_tsteps;

	/*
	 * Member functions
	 */

	static void et(double et_parameter) {
		encounter_threshold = et_parameter;
	}

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

#if 0
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
#endif
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
	 * Perform spatial_info updating on a given subspace, for given
	 * parameters.
	 */
	static void subspace_info_update(space_iterate si, int f, d2::image *weights, const d2::image *im, pt _pt) {
		while(!si.done()) {

			space_traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_space()) == 0) {
				si.next();
				continue;
			}

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

			point bb[2];

			st.get_view_local_bb(_pt, bb);

			point min = bb[0];
			point max = bb[1];

//				fprintf(stderr, "frame %d color update space pointer %p, bb (%f, %f) -> (%f, %f)\n", 
//						f, st.get_space(), min[0], min[1], max[0], max[1]);
//
//				fprintf(stderr, "space %p c=[%f %f %f]\n", st.get_space(), color[0], color[1], color[2]);
//				fprintf(stderr, "space %p occ=[%g]\n", st.get_space(), occupancy);

			/*
			 * Use the center of the bounding box to grab interpolation data.
			 */

			d2::point interp((min[0] + max[0]) / 2, (min[1] + max[1]) / 2);

//				fprintf(stderr, "interp=(%f, %f)\n", interp[0], interp[1]);

			/*
			 * For interpolation points, ensure that the
			 * bounding box area is at least 0.25. XXX: Why?
			 * Remove this constraint.
			 *
			 * XXX: Is interpolation useful for anything, given
			 * that we're using spatial info registers at multiple
			 * resolutions?
			 */

			if (/* (max[0] - min[0]) * (max[1] - min[1]) > 0.25
			 && */ max[0] > min[0]
			 && max[1] > min[1]) {
				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_bl(interp));
				d2::pixel pcolor = im->get_bl(interp);
				d2::pixel colordiff = (color - pcolor) * (ale_real) 256;

				if (falloff_exponent != 0) {
					d2::pixel max_diff = im->get_max_diff(interp) * (ale_real) 256;

					for (int k = 0; k < 3; k++)
					if (max_diff[k] > 1)
						colordiff[k] /= pow(max_diff[k], falloff_exponent);
				}

//					fprintf(stderr, "color_interp=(%f, %f, %f)\n", pcolor[0], pcolor[1], pcolor[2]);

//					sn->accumulate_color_2(pcolor, encounter);
				d2::pixel channel_occ = pexp(-colordiff * colordiff);
//					fprintf(stderr, "color_diff=(%f, %f, %f)\n", colordiff[0], colordiff[1], colordiff[2]);
//					fprintf(stderr, "channel_occ=(%g, %g, %g)\n", channel_occ[0], channel_occ[1], channel_occ[2]);

				/*
				 * XXX: the best approach is probably to use 3 separate occupancy
				 * data sets, just as there are 3 separate color data sets.
				 */

				ale_accum occ = channel_occ[0];

				for (int k = 1; k < 3; k++)
					if (channel_occ[k] < occ)
						occ = channel_occ[k];

				sn->accumulate_occupancy_2(occ, encounter[0]);
#if 0
				for (int k = 0; k < 3; k++)
					sn->accumulate_occupancy_2(channel_occ[k], encounter[k]);
#endif
			}

			/*
			 * Data structure to check modification of weights by
			 * higher-resolution subspaces.
			 */

			std::queue<d2::pixel> weight_queue;

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_space()->positive
			 || st.get_space()->negative) {

				/*
				 * Store information about current weights,
				 * so we will know which areas have been
				 * covered by higher-resolution subspaces.
				 */

				for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
				for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {
					if (i < 0 || j < 0)
						continue;
					weight_queue.push(weights->get_pixel(i, j));
				}
				
				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that later.
				 */

				space_iterate cleaved_space = si.cleave();

				cleaved_space.next();

				subspace_info_update(cleaved_space, f, weights, im, _pt);
			} else {
				si.next();
			}
				

			/*
			 * Iterate over pixels in the bounding box,
			 * adding new data to the subspace.  XXX:
			 * assume for now that all pixels in the
			 * bounding box intersect the subspace.
			 */

			for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
			for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {

				if (i < 0 || j < 0)
					continue;

				d2::pixel pcolor = im->get_pixel(i, j);
				d2::pixel colordiff = (color - pcolor) * (ale_real) 256;

				if (falloff_exponent != 0) {
					d2::pixel max_diff = im->get_max_diff(interp) * (ale_real) 256;

					for (int k = 0; k < 3; k++)
					if (max_diff[k] > 1)
						colordiff[k] /= pow(max_diff[k], falloff_exponent);
				}

//					fprintf(stderr, "(i, j) == (%d, %d); c=[%f %f %f]\n",
//							i, j, pcolor[0], pcolor[1], pcolor[2]);

				/*
				 * Determine the probability of
				 * encounter, divided by the occupancy.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j));

				/*
				 * Check for higher-resolution modifications.
				 */

				int high_res_mod = 0;

				if (weight_queue.size()) {
					if (weight_queue.front() != weights->get_pixel(i, j)) {
						high_res_mod = 1;
						encounter = d2::pixel(1, 1, 1) - weight_queue.front();
					}
					weight_queue.pop();
				}

				/*
				 * Update subspace.
				 */

				sn->accumulate_color_1(f, i, j, pcolor, encounter);
				d2::pixel channel_occ = pexp(-colordiff * colordiff);
//					fprintf(stderr, "encounter=(%f, %f, %f)\n", encounter[0], encounter[1], encounter[2]);
//					fprintf(stderr, "color_diff=(%f, %f, %f)\n", colordiff[0], colordiff[1], colordiff[2]);
//					fprintf(stderr, "channel_occ=(%g, %g, %g)\n", channel_occ[0], channel_occ[1], channel_occ[2]);

				ale_accum occ = channel_occ[0];

				for (int k = 1; k < 3; k++)
					if (channel_occ[k] < occ)
						occ = channel_occ[k];

				sn->accumulate_occupancy_1(f, i, j, occ, encounter[0]);

				/*
				 * If weights have not been updated by
				 * higher-resolution cells, then update
				 * weights at the current resolution.
				 */

				if (!high_res_mod)
					weights->pix(i, j) += encounter * occupancy;
			}

		}
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
			_pt.scale(1 / _pt.scale_2d());

			assert((int) floor(d2::align::of(f).unscaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(f).unscaled_width())
			     == (int) floor(_pt.scaled_width()));

			/*
			 * Get color data for the frames.
			 */
			const d2::image *im = d2::image_rw::open(f);

			assert(im);

			/*
			 * Allocate an image for storing encounter probabilities.
			 */
			d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()), 
					(int) floor(_pt.scaled_width()), 3);

			assert(weights);

			/*
			 * Call subspace_info_update for the root space.
			 */

			subspace_info_update(space_iterate(_pt), f, weights, im, _pt);

			d2::image_rw::close(f);

			delete weights;
		}

		/*
		 * Update all spatial_info structures.
		 */
		for (std::map<space *,spatial_info>::iterator i = spatial_info_map.begin(); i != spatial_info_map.end(); i++) {
			i->second.update_color();
			i->second.update_occupancy();

//			d2::pixel color = i->second.get_color();

//			fprintf(stderr, "space p=%p updated to c=[%f %f %f] o=%f\n",
//					i->first, color[0], color[1], color[2], 
//					i->second.get_occupancy());
		}
	}

	/*
	 * Support function for view() and depth().
	 */

	static const void view_recurse(int type, d2::image *im, d2::image *weights, space_iterate si, pt _pt) {
		while (!si.done()) {
			space_traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_space()) == 0) {
				si.next();
				continue;
			}

			spatial_info sn = spatial_info_map[st.get_space()];

			/*
			 * Get information on the subspace.
			 */

			d2::pixel color = sn.get_color();
			// d2::pixel color = d2::pixel(1, 1, 1) * (double) (((unsigned int) (st.get_space()) / sizeof(space)) % 65535);
			ale_real occupancy = sn.get_occupancy();

			/*
			 * Determine the view-local bounding box for the
			 * subspace.
			 */

			point bb[2];

			st.get_view_local_bb(_pt, bb);

			point min = bb[0];
			point max = bb[1];

			/*
			 * Data structure to check modification of weights by
			 * higher-resolution subspaces.
			 */

			std::queue<d2::pixel> weight_queue;

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_space()->positive
			 || st.get_space()->negative) {

				/*
				 * Store information about current weights,
				 * so we will know which areas have been
				 * covered by higher-resolution subspaces.
				 */

				for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
				for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++)
					weight_queue.push(weights->get_pixel(i, j));
				
				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that afterward.
				 */

				space_iterate cleaved_space = si.cleave();

				cleaved_space.next();

				view_recurse(type, im, weights, cleaved_space, _pt);

			} else {
				si.next();
			}
				

			/*
			 * Iterate over pixels in the bounding box, finding
			 * pixels that intersect the subspace.  XXX: assume
			 * for now that all pixels in the bounding box
			 * intersect the subspace.
			 */

			for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
			for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {

				/*
				 * Check for higher-resolution updates.
				 */

				if (weight_queue.size()) {
					if (weight_queue.front() != weights->get_pixel(i, j)) {
						weight_queue.pop();
						continue;
					}
					weight_queue.pop();
				}

				/*
				 * Determine the probability of encounter.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j)) * occupancy;

				/*
				 * Update images.
				 */

				if (type == 0) {

					/*
					 * Color view
					 */

					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * color;

				} else if (type == 1) {

					/*
					 * Weighted (transparent) depth display
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * depth_value;

				} else if (type == 2) {

					/*
					 * Ambiguity (ambivalence) measure.
					 */

					weights->pix(i, j) = d2::pixel(1, 1, 1);
					im->pix(i, j) += 0.1 * d2::pixel(1, 1, 1);

				} else if (type == 3) {

					/*
					 * Closeness measure.
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					if (weights->pix(i, j)[0] == 0) {
						weights->pix(i, j) = d2::pixel(1, 1, 1);
						im->pix(i, j) = d2::pixel(1, 1, 1) * depth_value;
					} else if (im->pix(i, j)[2] < depth_value) {
						im->pix(i, j) = d2::pixel(1, 1, 1) * depth_value;
					} else {
						continue;
					}

				} else if (type == 4) {

					/*
					 * Weighted (transparent) contribution display
					 */

					ale_pos contribution_value = sn.get_pocc_density() + sn.get_socc_density();
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

				} else if (type == 5) {

					/*
					 * Weighted (transparent) occupancy display
					 */

					ale_pos contribution_value = occupancy;
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

				} else 
					assert(0);
			}
		}
	}

	/*
	 * Generate an depth image from a specified view.
	 */
	static const d2::image *depth(pt _pt, int n = -1) {
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

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space_iterate si(_pt);

		view_recurse(1, im, weights, si, _pt);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {

			point w = _pt.pw_scaled(point(i, j, im->pix(i, j)[0] / weights->pix(i, j)[0]));

			if (weights->pix(i, j).min_norm() < encounter_threshold || excluded(w)) {
				im->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
				weights->pix(i, j) = d2::pixel::zero();
			} else if (normalize_weights)
				im->pix(i, j) /= weights->pix(i, j);
		}

		delete weights;

		return im;
	}

	static const d2::image *depth(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return depth(_pt, n);
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

		const d2::image *depths = NULL;

		if (d3px_count > 0)
			depths = depth(_pt, n);

		d2::image *im = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space_iterate si(_pt);

		view_recurse(0, im, weights, si, _pt);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			if (weights->pix(i, j).min_norm() < encounter_threshold
			 || (d3px_count > 0 && isnan(depths->pix(i, j)[0]))) {
				im->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
				weights->pix(i, j) = d2::pixel::zero();
			} else if (normalize_weights)
				im->pix(i, j) /= weights->pix(i, j);
		}

		if (d3px_count > 0)
			delete depths;

		delete weights;

		return im;
	}

	static void tcem(double _tcem) {
		tc_multiplier = _tcem;
	}

	static void oui(unsigned int _oui) {
		ou_iterations = _oui;
	}

	static void pa(unsigned int _pa) {
		pairwise_ambiguity = _pa;
	}

	static void pc(const char *_pc) {
		pairwise_comparisons = _pc;
	}

	static void d3px(int _d3px_count, double *_d3px_parameters) {
		d3px_count = _d3px_count;
		d3px_parameters = _d3px_parameters;
	}

	static void fx(double _fx) {
		falloff_exponent = _fx;
	}

	static void nw() {
		normalize_weights = 1;
	}

	static void no_nw() {
		normalize_weights = 0;
	}

	static int excluded(point p) {
		for (int n = 0; n < d3px_count; n++) {
			double *region = d3px_parameters + (6 * n);
			if (p[0] >= region[0]
			 && p[0] <= region[1]
			 && p[1] >= region[2]
			 && p[1] <= region[3]
			 && p[2] >= region[4]
			 && p[2] <= region[5])
				return 1;
		}

		return 0;
	}

	static const d2::image *view(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return view(_pt, n);
	}

	/*
	 * Add specified control points to the model
	 */
	static void add_control_points() {
	}

	typedef struct {point iw; point ip, is;} analytic;
	typedef std::multimap<ale_real,analytic> score_map;
	typedef std::pair<ale_real,analytic> score_map_element;

	/*
	 * Generate a map from scores to 3D points for various depths at point (i, j) in f1.
	 */
	static score_map p2f_score_map(unsigned int f1, unsigned int f2, pt _pt1, pt _pt2, 
			const d2::image *if1, const d2::image *if2, unsigned int i, unsigned int j) {

		score_map result;

		// fprintf(stderr, "generating score map (i, j) == (%u, %u)\n", i, j);

		// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

		/*
		 * Get the pixel color in the primary frame
		 */

		// d2::pixel color_primary = if1->get_pixel(i, j);

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

		if (isnan(slope)) {
			fprintf(stderr, "%d->%d (%d, %d) has undefined slope\n",
					f1, f2, i, j);
			return result;
		}

		/*
		 * Make absurdly large/small slopes either infinity, negative infinity, or zero.
		 */

		if (fabs(slope) > if2->width() * 100) {
			double zero = 0;
			double one  = 1;
			double inf  = one / zero;
			slope = inf;
		} else if (slope < 1 / (double) if2->height() / 100
		        && slope > -1/ (double) if2->height() / 100) {
			slope = 0;
		}

		ale_pos top_intersect = p1[1] - p1[0] * slope;
		ale_pos lef_intersect = p1[0] - p1[1] / slope;
		ale_pos rig_intersect = p1[0] - (p1[1] - if2->width() + 2) / slope;
		ale_pos sp_i, sp_j;

		if (slope == 0) {
			sp_i = lef_intersect;
			sp_j = 0;
		} else if (finite(slope) && top_intersect >= 0 && top_intersect < if2->width() - 1) {
			sp_i = 0;
			sp_j = top_intersect;
		} else if (slope > 0 && lef_intersect >= 0 && lef_intersect <= if2->height() - 1) {
			sp_i = lef_intersect;
			sp_j = 0;
		} else if (slope < 0 && rig_intersect >= 0 && rig_intersect <= if2->height() - 1) {
			sp_i = rig_intersect;
			sp_j = if2->width() - 2;
		} else {
			return result;
		}

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
			ii <= if2->height() - 1 && jj <= if2->width() - 1 && ii >= 0 && jj >= 0; 
			ii += incr_i, jj += incr_j) {

			/*
			 * Create an orthonormal basis to
			 * determine line intersection.
			 */

			point bp0 = _pt1.pw_unscaled(point(i, j, 0));
			point bp1 = _pt1.pw_unscaled(point(i, j, 10));
			point bp2 = _pt2.pw_unscaled(point(ii, jj, 0));

			point b0  = (bp1 - bp0).normalize();
			point b1n = bp2 - bp0;
			point b1  = (b1n - b1n.dproduct(b0) * b0).normalize();
			point b2  = point(0, 0, 0).xproduct(b0, b1).normalize(); // Should already have norm=1
			
			/*
			 * Select a fourth point to define a second line.
			 */

			point p3  = _pt2.pw_unscaled(point(ii, jj, 10));

			/*
			 * Representation in the new basis.
			 */

			d2::point nbp0 = d2::point(bp0.dproduct(b0), bp0.dproduct(b1));
			// d2::point nbp1 = d2::point(bp1.dproduct(b0), bp1.dproduct(b1));
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

			ale_pos b2_offset = b2.dproduct(bp0);

			/*
			 * Map the intersection back to the world
			 * basis.
			 */

			point iw = intersection[0] * b0 + intersection[1] * b1 + b2_offset * b2;

			/*
			 * Reject intersection points behind a
			 * camera.
			 */

			point icp = _pt1.wc(iw);
			point ics = _pt2.wc(iw);


			if (icp[2] >= 0 || ics[2] >= 0)
				continue;

			/*
			 * Reject clipping plane violations.
			 */

			if (iw[2] > front_clip
			 || iw[2] < rear_clip)
				continue;

			/*
			 * Score the point.
			 */

			point ip = _pt1.wp_unscaled(iw);

			point is = _pt2.wp_unscaled(iw);

			analytic _a = { iw, ip, is };

			/*
			 * Calculate score from color match.  Assume for now
			 * that the transformation can be approximated locally
			 * with a translation.
			 */

			ale_pos score = 0;
			ale_pos divisor = 0;
			ale_pos l1_multiplier = 0.125;

			if (if1->in_bounds(ip.xy())
			 && if2->in_bounds(is.xy())) {
				divisor += 1 - l1_multiplier;
				score += (1 - l1_multiplier)
				       * (if1->get_bl(ip.xy()) - if2->get_bl(is.xy())).normsq();
			}

			for (int iii = -1; iii <= 1; iii++)
			for (int jjj = -1; jjj <= 1; jjj++) {
				d2::point t(iii, jjj);

				if (!if1->in_bounds(ip.xy() + t)
				 || !if2->in_bounds(is.xy() + t))
					continue;

				divisor += l1_multiplier;
				score   += l1_multiplier
					 * (if1->get_bl(ip.xy() + t) - if2->get_bl(is.xy() + t)).normsq();
					 
			}

			/*
			 * Include third-camera contributions in the score.
			 */

			if (tc_multiplier != 0)
			for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
				if (f == f1 || f == f2)
					continue;

				assert(cl);
				assert(cl->reference);
				assert(cl->reference[f]);

				pt _ptf = align::projective(f);

				point p = _ptf.wp_unscaled(iw);

				if (!cl->reference[f]->in_bounds(p.xy())
				 || !if2->in_bounds(ip.xy()))
					continue;

				divisor += tc_multiplier;
				score   += tc_multiplier
					 * (if1->get_bl(ip.xy()) - cl->reference[f]->get_bl(p.xy())).normsq();
			}

			/*
			 * Reject points with undefined score.
			 */

			if (!finite(score / divisor))
				continue;

			/*
			 * Add the point to the score map.
			 */

			result.insert(score_map_element(score / divisor, _a));
		}

		return result;
	}

	/*
	 * Attempt to refine space around a point, to high and low resolutions
	 * for two cameras, resulting in four resolutions in total.
	 */

	static void refine_space(point iw, pt _pt1, pt _pt2) {

		space_traverse st = space_traverse::root();

		if (!st.includes(iw)) {
			assert(0);
			return;
		}

		int camera_highres[2] = {0, 0};
		int camera_lowres[2] = {0, 0};

		/*
		 * Loop until all resolutions of interest have been generated.
		 */
		
		for(;;) {

			point frame_min[2] = { point::posinf(), point::posinf() },
			      frame_max[2] = { point::neginf(), point::neginf() };

			point p[2] = { st.get_min(), st.get_max() };

			/*
			 * Cycle through the corner points bounding the
			 * subspace to determine a bounding box.  
			 *
			 * NB: This code is not identical to
			 * get_view_local_bb(), as it does not clip the
			 * results.
			 */

			for (int ibit = 0; ibit < 2; ibit++)
			for (int jbit = 0; jbit < 2; jbit++)
			for (int kbit = 0; kbit < 2; kbit++) {
				point pp = point(p[ibit][0], p[jbit][1], p[kbit][2]);

				point ppp[2] = {_pt1.wp_unscaled(pp), _pt2.wp_unscaled(pp)};

				for (int f = 0; f < 2; f++)
				for (int d = 0; d < 3; d++) {
					if (ppp[f][d] < frame_min[f][d] || isnan(ppp[f][d]))
						frame_min[f][d] = ppp[f][d];
					if (ppp[f][d] > frame_max[f][d] || isnan(ppp[f][d]))
						frame_max[f][d] = ppp[f][d];
				}
			}

			/*
			 * Generate any new desired spatial registers.
			 */

			for (int f = 0; f < 2; f++) {

				/*
				 * Low resolution
				 */

				if (frame_max[f][0] - frame_min[f][0] < 2
				 && frame_max[f][1] - frame_min[f][1] < 2
				 && camera_lowres[f] == 0) {
					spatial_info_map[st.get_space()];
					camera_lowres[f] = 1;
				}

				/*
				 * High resolution.
				 */

				if (frame_max[f][0] - frame_min[f][0] < 1
				 && frame_max[f][1] - frame_min[f][1] < 1
				 && camera_highres[f] == 0) {
					spatial_info_map[st.get_space()];
					camera_highres[f] = 1;
				}
			}

			/*
			 * Check for completion
			 */

			if (camera_highres[0]
			 && camera_highres[1]
			 && camera_lowres[0]
			 && camera_lowres[1])
				return;

			/*
			 * Check precision before analyzing space further.
			 */

			if (st.precision_wall()) {
				fprintf(stderr, "\n\n*** Error: reached subspace precision wall ***\n\n");
				assert(0);
				return;
			}

			if (st.positive().includes(iw)) {
				st = st.positive();
				total_tsteps++;
			} else if (st.negative().includes(iw)) {
				st = st.negative();
				total_tsteps++;
			} else {
				fprintf(stderr, "failed iw = (%f, %f, %f)\n", 
						iw[0], iw[1], iw[2]);
				assert(0);
			}
		}
	}

	/*
	 * Analyze space in a manner dependent on the score map.
	 */

	static void analyze_space_from_map(unsigned int f1, unsigned int f2, unsigned int i, 
			unsigned int j, pt _pt1, pt _pt2, score_map _sm) {

		int accumulated_ambiguity = 0;
		int max_acc_amb = pairwise_ambiguity;

		for(score_map::iterator smi = _sm.begin(); smi != _sm.end(); smi++) {

			point iw = smi->second.iw;
			point ip = smi->second.ip;
			// point is = smi->second.is;

			if (accumulated_ambiguity++ >= max_acc_amb)
				break;

			total_ambiguity++;

			/*
			 * Attempt to refine space around the intersection point.
			 */

			refine_space(iw, _pt1, _pt2);
		}
	}

	static void refine_space_for_output(pt _pt, space_traverse st = space_traverse::root()) {
		if (!spatial_info_map.count(st.get_space())) {
			if (st.get_space()->positive)
				refine_space_for_output(_pt, st.positive());
			if (st.get_space()->negative)
				refine_space_for_output(_pt, st.negative());
			return;
		}

		point bb[2];
		st.get_view_local_bb(_pt, bb);

		if (bb[0].xy().lengthtosq(bb[1].xy()) < 1)
			return;

		if (!_pt.scaled_in_bounds(bb[0]) || !_pt.scaled_in_bounds(bb[1]))
			return;

		spatial_info_map[st.positive().get_space()];
		spatial_info_map[st.negative().get_space()];

		refine_space_for_output(_pt, st.positive());
		refine_space_for_output(_pt, st.negative());
	}

	static void refine_space_for_output(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		for (unsigned int f = 0; f < d2::image_rw::count(); f++) 
			if (d_out[f] || v_out[f])
				refine_space_for_output(d3::align::projective(f));

		for (std::map<const char *, pt>::iterator i = d3_depth_pt->begin(); i != d3_depth_pt->end(); i++)
			refine_space_for_output(i->second);

		for (std::map<const char *, pt>::iterator i = d3_output_pt->begin(); i != d3_output_pt->end(); i++)
			refine_space_for_output(i->second);
	}


	/*
	 * Initialize space and identify regions of interest for the adaptive
	 * subspace model.
	 */
	static void make_space(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		fprintf(stderr, "Subdividing 3D space");

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

			if (!d_out[f1] && !v_out[f1] && !d3_depth_pt->size()
			 && !d3_output_pt->size() && strcmp(pairwise_comparisons, "all"))
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

				total_pixels++;

				/*
				 * Generate a map from scores to 3D points for
				 * various depths in f1.
				 */

				score_map _sm = p2f_score_map(f1, f2, _pt1, _pt2, if1, if2, i, j);

				/*
				 * Analyze space in a manner dependent on the score map.
				 */

				analyze_space_from_map(f1, f2, i, j, _pt1, _pt2, _sm);

			}

			/*
			 * This ordering should ensure that image f1 is cached.
			 */

			d2::image_rw::close(f2);
			d2::image_rw::close(f1);
		}

		/*
		 * This is expensive.
		 */

		// refine_space_for_output(d_out, v_out, d3_depth_pt, d3_output_pt);

		fprintf(stderr, ".\n");
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
	static void reduce_cost_to_search_depth(d2::exposure *exp_out, int inc_bit) {

		/*
		 * Subspace model
		 */

		for (unsigned int i = 0; i < ou_iterations; i++)
			spatial_info_update();
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
