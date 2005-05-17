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

#include "scene.h"

/*
 * See scene.h for details on these variables.
 */

scene::lod *scene::cl;
scene::triangle *scene::triangle_head[2];
std::map<point *, scene::vertex_aux> scene::vam;
ale_pos scene::edge_cost_multiplier = 0.001;
ale_pos scene::angle_cost_multiplier = 0.001;

/*
 * Functions.
 */

point scene::frame_to_frame(d2::point p, int f1, int f2, zbuf_elem *z1, zbuf_elem *z2) {
	pt _pt1 = align::of(f1);
	pt _pt2 = align::of(f2);

	int i = round(p[0]);
	int j = round(p[1]);

	if (i < 0 || i >= _pt1.scaled_height())
		return point::undefined();
	if (j < 0 || j >= _pt1.scaled_width())
		return point::undefined();

	z1[i * _pt1.scaled_width() + j]->find_nearest();

	triangle *t = z1[i * _pt1.scaled_width() + j]->nearest();

	if (!t)
		return point::undefined();

	point c_ray = _pt1.pc_scaled(point(i, j, -1));
	point c[3];

	for (int v = 0; v < 3; v++)
		c[v] = _pt1.wc(*t->vertices[v]);

	point multipliers = rt_intersect(c, c_ray);

	c_ray *= multipliers[2];

	point p2 = _pt2.wp_scaled(c_ray);

	int i2 = round(p2[0]);
	int j2 = round(p2[1]);

	if (i2 < 0 || i2 >= _pt2.scaled_height())
		return point::undefined();
	if (j2 < 0 || j2 >= _pt2.scaled_width())
		return point::undefined();

	z2[i2 * _pt2.scaled_width() + j2].find_nearest();

	triangle *t2 = z2[i2 * _pt2.scaled_width() + j2].nearest();

	if (t != t2)
		return point::undefined();

	return p2;
}

void scene::zbuf_elem::find_nearest(pt _pt, int i, int j) {
	clear_nearest();

	point c_ray = _pt.pc_scaled(point(i, j, -1));

	ale_pos least_distance = +1;
	ale_pos zero = +0;
	least_distance = least_distance / zero;

	assert (isinf(least_distance) == 1);

	if (array == NULL)
		return;

	for (int t_index = 0; t_index < size; t_index++) {
		triangle *t = array[t_index];
		point c[3];

		for (int v = 0; v < 3; v++)
			c[v] = _pt.wc(*t->vertices[v]);

		if (!is_interior_c(c, c_ray))
			continue;

		point multipliers = rt_intersect(c, c_ray);
		ale_pos ray_multiplier = fabs(multipliers[2]);

		if (ray_multiplier >= least_distance)
			continue;

		least_distance = ray_multiplier;
		_nearest = t;
	}
}

ale_accum scene::vertex_movement_cost(scene::triangle *t, point *vertex, point new_position, zbuf_elem **z) {
	ale_accum orig_color_cost = 0, new_color_cost = 0;
	ale_accum orig_color_div  = 0, new_color_div  = 0;
	ale_accum orig_geom_cost  = 0, new_geom_cost  = 0;

	point original_position = *vertex;

	/*
	 * Determine the triangles surrounding the vertex.
	 */

	std::set<triangle *> *t_set = new std::set<triangle *>;

	t->triangles_around_vertex(vertex, t_set);

	/*
	 * Determine the set of vertices associated with the triangles
	 * surrounding the specified vertex.
	 */

	std::set<point *> *v_set = new std::set<point *>;

	for (std::forward_iterator i = t_set->begin(); i != t_set->end(); i++)
	for (int v = 0; v < 3; v++)
		v_set.insert((*i)->vertices[v]);

	/*
	 * Determine geometric costs.
	 */

	for (std::forward_iterator i = t_set->begin(); i != t_set->end(); i++)
	for (int v = 0; v < 3; v++)
		orig_geom_cost += (*i)->edge_cost() + (*i)->angle_cost();

	/*
	 * Determine color costs
	 */

	assert(0);
	
	/*
	 * Move the target vertex
	 */

	*vertex = new_position;

	/*
	 * Determine geometric costs.
	 */

	for (std::forward_iterator i = t_set->begin(); i != t_set->end(); i++)
	for (int v = 0; v < 3; v++)
		new_geom_cost += (*i)->edge_cost() + (*i)->angle_cost();

	/*
	 * Determine color costs.
	 */
	 
	assert(0);

	/*
	 * Free sets.
	 */
	
	delete t_set;
	delete v_set;

	/*
	 * Restore the target vertex
	 */

	*vertex = original_position;

	/*
	 * Apply divisors.
	 */

	orig_color_cost /= orig_color_div;
	new_color_cost  /= new_color_div;

	if (!finite(orig_color_cost) || !finite(new_color_cost))
		orig_color_cost = new_color_cost = 0;

	/*
	 * Return the error difference
	 */

	return (new_color_cost + new_geom_cost)
	     - (orig_color_cost + orig_geom_cost);
}
