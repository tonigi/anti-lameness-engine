// Copyright 2003, 2004, 2005 David Hilvert <dhilvert@gmail.com>,
//                                          <dhilvert@auricle.dyndns.org>,
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

point scene::frame_to_frame(d2::point p, const pt &_pt1, const pt &_pt2, zbuf_elem *z1, zbuf_elem *z2) {

	int i = (int) round(p[0]);
	int j = (int) round(p[1]);

	if (i < 0 || i > _pt1.scaled_height() - 1)
		return point::undefined();
	if (j < 0 || j > _pt1.scaled_width() - 1)
		return point::undefined();

	triangle *t = z1[i * (int) floor(_pt1.scaled_width()) + j].nearest(_pt1, i, j);

	if (!t)
		return point::undefined();

	point c_ray = _pt1.pc_scaled(point(i, j, -1));
	point c[3];

	for (int v = 0; v < 3; v++)
		c[v] = _pt1.wc(*t->vertices[v]);

	point multipliers = rt_intersect(c_ray, c);

	c_ray *= multipliers[2];

	point p2 = _pt2.wp_scaled(_pt1.cw(c_ray));

	if (p2[0] < 0 || p2[0] > floor(_pt2.scaled_height() - 1))
		return point::undefined();
	if (p2[1] < 0 || p2[1] > floor(_pt2.scaled_width() - 1))
		return point::undefined();

	int i2 = (int) round(p2[0]);
	int j2 = (int) round(p2[1]);

	triangle *t2 = z2[i2 * (int) floor(_pt2.scaled_width()) + j2].nearest(_pt2, i2, j2);

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

	assert(tset);

	for (std::set<triangle *>::iterator ii = tset->begin(); ii != tset->end(); ii++) {
		triangle *t = (*ii);
		point c[3];

		for (int v = 0; v < 3; v++)
			c[v] = _pt.wc(*t->vertices[v]);

		if (!is_interior_c(c, c_ray))
			continue;

		point multipliers = rt_intersect(c_ray, c);
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

	for (std::set<triangle *>::iterator i = t_set->begin(); i != t_set->end(); i++)
	for (int v = 0; v < 3; v++)
		v_set->insert((*i)->vertices[v]);

	/*
	 * Determine bounding boxes for calculating color costs
	 */

	point *bb = new point[2 * d2::image_rw::count()];
	for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
		point *bbp = bb + 2 * f;
		pt _pt = align::projective(f);
		_pt.scale(cl->sf / _pt.scale_2d());

		ale_accum inf = +1;
		ale_accum zero = 0;
		inf /= zero;
		assert(isinf(inf)  ==  1);
		assert(isinf(-inf) == -1);

		bbp[0][0] = bbp[0][1] =  inf;
		bbp[1][0] = bbp[1][1] = -inf;

		for (std::set<point *>::iterator i = v_set->begin(); i != v_set->end(); i++) {
			point p = _pt.wp_scaled(**i);

			if (p[0] < bbp[0][0])
				bbp[0][0] = p[0];
			if (p[1] < bbp[0][1])
				bbp[0][1] = p[1];
			if (p[0] > bbp[1][0])
				bbp[1][0] = p[0];
			if (p[1] > bbp[1][1])
				bbp[1][1] = p[1];
		}

		point np = _pt.wp_scaled(new_position);

		if (np[0] < bbp[0][0])
			bbp[0][0] = np[0];
		if (np[1] < bbp[0][1])
			bbp[0][1] = np[1];
		if (np[0] > bbp[1][0])
			bbp[1][0] = np[0];
		if (np[1] > bbp[1][1])
			bbp[1][1] = np[1];

		for (int d = 0; d < 2; d++) {
			if (bbp[d][0] < 0)
				bbp[d][0] = 0;
			if (bbp[d][0] > floor(_pt.scaled_height()) - 1)
				bbp[d][0] = floor(_pt.scaled_height()) - 1;
			if (bbp[d][1] < 0)
				bbp[d][1] = 0;
			if (bbp[d][1] > floor(_pt.scaled_width()) - 1)
				bbp[d][1] = floor(_pt.scaled_width()) - 1;
		}
	}

	/*
	 * Add triangles to bounding box areas.
	 */

	for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
		pt _pt = align::projective(f);

		_pt.scale(cl->sf / _pt.scale_2d());
			
		for (int i = (int) floor(bb[f * 2 + 0][0]); i <= (int) ceil(bb[f * 2 + 1][0]); i++)
		for (int j = (int) floor(bb[f * 2 + 0][1]); j <= (int) ceil(bb[f * 2 + 1][1]); j++) {
			int n = i * (int) floor(_pt.scaled_width()) + j;
			z[f][n].insert(t_set->begin(), t_set->end());
		}
	}

	/*
	 * Determine geometric costs.
	 */

	for (std::set<triangle *>::iterator i = t_set->begin(); i != t_set->end(); i++)
	for (int v = 0; v < 3; v++)
		orig_geom_cost += (*i)->edge_cost() + (*i)->angle_cost();

	/*
	 * Determine color costs
	 */

	for (unsigned int f1 = 0; f1 < d2::image_rw::count(); f1++)
	for (unsigned int f2 = 0; f2 < d2::image_rw::count(); f2++) {
		pt _pt1 = align::projective(f1);
		pt _pt2 = align::projective(f2);

		_pt1.scale(cl->sf / _pt1.scale_2d());
		_pt2.scale(cl->sf / _pt2.scale_2d());

		if (f1 == f2)
			continue;
			
		for (int i = (int) floor(bb[f1 * 2 + 0][0]); i <= (int) ceil(bb[f1 * 2 + 1][0]); i++)
		for (int j = (int) floor(bb[f1 * 2 + 0][1]); j <= (int) ceil(bb[f1 * 2 + 1][1]); j++) {
			d2::point p1(i, j);
			d2::point p2 = frame_to_frame(p1, _pt1, _pt2, z[f1], z[f2]).xy();

			if (!p2.defined())
				continue;

			orig_color_cost += (cl->reference[f1]->get_bl(p1)
			                  - cl->reference[f2]->get_bl(p2)).normsq();

			orig_color_div  += 1;
		}
	}
	
	/*
	 * Move the target vertex
	 */

	*vertex = new_position;

	/*
	 * Determine geometric costs.
	 */

	for (std::set<triangle *>::iterator i = t_set->begin(); i != t_set->end(); i++)
	for (int v = 0; v < 3; v++)
		new_geom_cost += (*i)->edge_cost() + (*i)->angle_cost();

	/*
	 * Determine color costs.
	 */
	 
	for (unsigned int f1 = 0; f1 < d2::image_rw::count(); f1++)
	for (unsigned int f2 = 0; f2 < d2::image_rw::count(); f2++) {
		pt _pt1 = align::projective(f1);
		pt _pt2 = align::projective(f2);

		_pt1.scale(cl->sf / _pt1.scale_2d());
		_pt2.scale(cl->sf / _pt2.scale_2d());

		if (f1 == f2)
			continue;
		for (int i = (int) floor(bb[f1 * 2 + 0][0]); i <= (int) ceil(bb[f1 * 2 + 1][0]); i++)
		for (int j = (int) floor(bb[f1 * 2 + 0][1]); j <= (int) ceil(bb[f1 * 2 + 1][1]); j++) {
			d2::point p1(i, j);
			d2::point p2 = frame_to_frame(p1, _pt1, _pt2, z[f1], z[f2]).xy();

			if (!p2.defined())
				continue;

			new_color_cost += (cl->reference[f1]->get_bl(p1)
			                 - cl->reference[f2]->get_bl(p2)).normsq();

			new_color_div  += 1;
		}
	}

	/*
	 * Free allocated memory.
	 */
	
	delete[] bb;
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

	if (sqrt(new_color_cost) + new_geom_cost
	  - sqrt(orig_color_cost) - orig_geom_cost < 0)
		fprintf(stderr, "[ncc=%f ngc=%f occ=%f ogc=%f] ", 
			sqrt(new_color_cost), new_geom_cost,
			sqrt(orig_color_cost), orig_geom_cost);

	return (sqrt(new_color_cost) + new_geom_cost)
	     - (sqrt(orig_color_cost) + orig_geom_cost);
}
