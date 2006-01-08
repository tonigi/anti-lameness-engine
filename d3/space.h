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

		static traverse root() {

			traverse result;

			result.current = root_node;
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

		traverse positive() {

			assert(current);

			int next_split = get_next_split();

			if (current->positive == NULL) {
				current->positive = new node;
			}
			
			traverse result;

			result.current = current->positive;
			result.min = min;
			result.max = max;

			result.min[next_split] = split_coordinate(next_split);

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

		iterate(point co, traverse top) {
			camera_origin = co;
			node_stack.push(top);
		}

	public:
		iterate(pt _pt, traverse top = traverse::root()) {
			camera_origin = _pt.cw(point(0, 0, 0));
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
