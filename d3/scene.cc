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

/*
 * Function bodies for recursive cost adjustment.
 */

void scene::local_puts(char *s, unsigned int search_depth, unsigned int f, 
                       unsigned int x, unsigned int y, double value) {

	static int success_bit = 0;
	static int track_bit = 0;
	static unsigned int last_sum = 0;


	// if (success_bit == 0 && x != 201)
	//	return;

	if (search_depth == 0 && (x == 59 || x == 60)) 
		track_bit = 0;
	else if (search_depth == 0) {
		track_bit = 0;
	}

	// if (s != "success" && s != "early_success" && success_bit != 1)
	// 	return;

	if (success_bit != 1 && (s == "success" || s == "early_success"))
		success_bit = 1;

	if (search_depth == 0 && (s == "success" || s == "early_success"))
		success_bit = 0;

	if (!success_bit && !track_bit)
		return;

	if (search_depth + f + x + y != last_sum) {
		last_sum = search_depth + f + x + y;
	} else if (s != "max_struct.cost") {
		return;
	}

	return;

	if (search_depth > 2)
	     return;

	for (unsigned int i = 0; i < search_depth; i++)
		fprintf(stderr, " ");
	
	// fprintf(stderr, "[d=%d f=%d x=%d y=%d color=(%.10f, %.10f, %.10f) depth=%f %s=%f]\n",
	fprintf(stderr, "[d=%d f=%d x=%d y=%d %s=%f]\n",
			search_depth,
			f, x, y, 
//			cl->pz_color[f]->get_pixel(x, y)[0],
//			cl->pz_color[f]->get_pixel(x, y)[1],
//			cl->pz_color[f]->get_pixel(x, y)[2],
//			cl->pz_depth[f]->get_pixel(x, y)[0],
			s, value);

}

