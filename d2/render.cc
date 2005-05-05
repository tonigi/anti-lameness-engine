// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                              <dhilvert@ugcs.caltech.edu>

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

#include "render.h"

/*
 * See align.h for details on these variables.
 */

unsigned int render::rx_count;
ale_pos *render::rx_parameters;
int render::rx_show;
render *render::directory[ACTIVE_RENDERER_COUNT];
int render::directory_length;
int render::extend;
double render::scale_factor;
double render::wt = 0.1;
