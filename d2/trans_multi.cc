// Copyright 2008 David Hilvert <dhilvert@auricle.dyndns.org>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "trans_multi.h"

/*
 * See trans_multi.h for details on these variables.
 */

unsigned int trans_multi::_multi = 2;
ale_pos trans_multi::_multi_decomp = 100;
ale_pos trans_multi::_multi_improvement = 0;
unsigned int trans_multi::_track = 0;  /* None */
ale_pos trans_multi::track_x = 0;
ale_pos trans_multi::track_y = 0;
