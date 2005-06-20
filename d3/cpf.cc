// Copyright 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
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

/*
 * d3/cpf.cc: Control point file handler static variables.
 */

#include "cpf.h"

FILE *cpf::load_f = NULL;
FILE *cpf::save_f = NULL;
int cpf::load_version = -1;

const char *cpf::load_n = NULL;
const char *cpf::save_n = NULL;
int cpf::save_version = 0;

struct cpf::control_point *cpf::cp_array = NULL;
unsigned int cpf::cp_array_max = 0;
unsigned int cpf::cp_index = 0;

ale_pos cpf::cpp_lower = 0.125;
ale_pos cpf::cpp_upper = 32;
ale_pos cpf::va_upper  = 32;

ale_pos cpf::stereo_threshold = 4;

unsigned int cpf::systems_solved = 0;

int cpf::total_error_mean = 1;
