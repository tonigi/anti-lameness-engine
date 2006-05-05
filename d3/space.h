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
 * d3/space.h: Representation of 3D space.
 */

#ifndef __space_h__
#define __space_h__

#include "point.h"

class space {
public:

	/*
	 * Structure to hold a subdivisible region of space.
	 */
	struct node {
		struct node *positive;
		struct node *negative;

		node() {
			positive = NULL;
			negative = NULL;
		}
	};

private:
	/*
	 * Space root pointer
	 */
	static node *root_node;

public:

	static void init_root() {
		root_node = new node;
	}

	/*
	 * Space traversal and navigation class.
	 */

	class traverse {
		node *current;
		point bounds[2];
	
	public:

		static int get_next_split(point min, point max) {
			assert(min[0] < max[0]);
			assert(min[1] < max[1]);
			assert(min[2] < max[2]);

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
			 && max[0] - min[0] >= max[2] - min[2]
			 && (max[0] >= max[2] || !isinf(min[2]))
			 && (min[0] <= min[2] || !isinf(max[2])))
				return 0;

			if (max[1] - min[1] > max[2] - min[2]
			 && (max[1] >= max[2] || !isinf(min[2]))
			 && (min[1] <= min[2] || !isinf(max[2])))
				return 1;

			return 2;
		}

		static ale_pos split_coordinate(int d, point min, point max) {
			if (isinf(max[d]) && isinf(min[d]))
				return 0;

			if (isinf(max[d]))
				return tan((atan(min[d]) + M_PI/2) / 2);

			if (isinf(min[d]))
				return tan((atan(max[d]) - M_PI/2) / 2);

			return (min[d] + max[d]) / 2;
		}

		static int get_next_cells(int d, point min, point max, point cells[2][2]) {
			cells[0][0] = min;
			cells[0][1] = max;
			cells[1][0] = min;
			cells[1][1] = max;

			ale_pos split_point = split_coordinate(d, min, max);

			if (split_point == min[d]
			 || split_point == max[d]
			 || !finite(split_point))
				return 0;

			cells[0][1][d] = split_point;
			cells[1][0][d] = split_point;

			return 1;
		}

		int get_next_split() {
			return get_next_split(bounds[0], bounds[1]);
		}

		ale_pos split_coordinate(int d) {
			return split_coordinate(d, bounds[0], bounds[1]);
		}

		ale_pos split_coordinate() {
			int next_split = get_next_split();
			return split_coordinate(next_split);
		}

		static traverse root() {

			traverse result;

			result.current = root_node;
			result.bounds[0] = point::neginf();
			result.bounds[1] = point::posinf();

			assert(result.current);

			return result;
		}

		int precision_wall() {
			int next_split = get_next_split();
			ale_pos split_point = split_coordinate(next_split);

			point &min = bounds[0];
			point &max = bounds[1];

			assert(split_point <= max[next_split]);
			assert(split_point >= min[next_split]);
			
			if (split_point == min[next_split] || split_point == max[next_split]) 
				return 1;

			return 0;
		}

		traverse positive() {

			assert(current);

			int next_split = get_next_split();

			if (current->positive == NULL) {
				current->positive = new node;
			}
			
			traverse result;

			result.current = current->positive;
			result.bounds[0] = bounds[0];
			result.bounds[1] = bounds[1];

			result.bounds[0][next_split] = split_coordinate(next_split);

			assert(result.current);

			return result;
		}

		traverse negative() {

			assert(current);

			int next_split = get_next_split();

			if (current->negative == NULL) {
				current->negative = new node;
			}
			
			traverse result;

			result.current = current->negative;
			result.bounds[0] = bounds[0];
			result.bounds[1] = bounds[1];

			result.bounds[1][next_split] = split_coordinate(next_split);

			assert(result.current);

			return result;
		}

		point get_min() const {
			return bounds[0];
		}

		point get_max() const {
			return bounds[1];
		}

		const point *get_bounds() const {
			return bounds;
		}

		point get_centroid() const {
			return (bounds[0] + bounds[1]) / 2;
		}

		int includes(point p) {

			for (int d = 0; d < 3; d++) {
				if (p[d] > bounds[1][d])
					return 0;
				if (p[d] < bounds[0][d])
					return 0;
				if (isnan(p[d]))
					return 0;
			}

			return 1;
		}

		node *get_node() {
			assert(current);
			return current;
		}

	};

	/*
	 * Class to iterate through subspaces based on proximity to a camera.
	 */

	class iterate {
		std::stack<traverse> node_stack;
		point camera_origin;

	public:
		iterate(point co, traverse top = traverse::root()) {
			camera_origin = co;
			node_stack.push(top);
		}

		int next() {
			if (node_stack.empty())
				return 0;

			traverse st = node_stack.top();

			int d = st.get_next_split();

			ale_pos split_coordinate = st.split_coordinate();

			node *n = st.get_node()->negative;
			node *p = st.get_node()->positive;

			if (camera_origin[d] > split_coordinate) {
				if (n) {
					node_stack.top() = st.negative();
					if (p)
						node_stack.push(st.positive());
				} else {
					if (p)
						node_stack.top() = st.positive();
					else
						node_stack.pop();
				}
			} else {
				if (p) {
					node_stack.top() = st.positive();
					if (n)
						node_stack.push(st.negative());
				} else {
					if (n)
						node_stack.top() = st.negative();
					else
						node_stack.pop();
				}
			}

			return (!node_stack.empty());
		}

		iterate cleave() {
			assert (!node_stack.empty());

			iterate result(camera_origin, node_stack.top());
			
			node_stack.pop();

			return result;
		}

		int done() {
			return node_stack.empty();
		}

		traverse get() {
			assert (!node_stack.empty());
			return node_stack.top();
		}
	};
};

#endif
