// Copyright 2006 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                              <dhilvert@ugcs.caltech.edu>

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

#include "input.h"

/*
 * See input.h for details on these variables.
 */
std::stack<input::environment *> input::environment::environment_stack;
std::set<input::environment *> input::environment::environment_set;

int input::global_options = 1;

input::environment *input::genv = NULL;

/*
 * List of options that can be used as nonglobals.
 */

const char *input::supported_nonglobal_option_table[] = {
	"threads",
	"per-cpu",
	"perturb-upper",
	"ev",
	"gs",
	"gs-mo",
	"black",
	"ma-card",
	"ma-cont",
	NULL
};

const char *input::focus_prefixes[] = {
	"ci=",
	"fr=",
	"ht=",
	"vt=",
	"sy=",
	"ey=",
	"sx=",
	"ex=",
	"sd=",
	"ed=",
	"ap=",
	"sc=",
	"sr=",
	"fs=",
	NULL
};

/*
 * Entries in this table are:
 *
 *   name, map_name, map_value, arg_count, multi
 *
 * The table must be terminated with name of NULL.
 */

input::simple_option input::simple_option_table[] = {

	{ "8bpc", "bpc" }, 
	{ "16bpc", "bpc" }, 

	{ "plain", "format" }, 
	{ "raw", "format" },
	{ "auto", "format" },

	{ "align-all", "align" },
	{ "align-green", "align" },
	{ "align-sum", "align" },

	{ "translation", "transformation" },
	{ "euclidean", "transformation" },
	{ "projective", "transformation" },

	{ "identity", "transformation-default" },
	{ "follow", "transformation-default" },

	{ "perturb-output", "perturb" },
	{ "perturb-source", "perturb" },

	{ "fail-optimal", "fail" },
	{ "fail-default", "fail" },

	{ "profile" },

	{ "extend" },
	{ "no-extend", "extend", "0" },

	{ "cache", NULL, NULL, 1 },

	{ "ev", NULL, NULL, 1 },
	{ "black", NULL, NULL, 1 },

	{ "threads", NULL, NULL, 1 },
	{ "per-cpu", NULL, NULL, 1 },

	{ "oc" }, 
	{ "no-oc", "oc", "0" },

	{ "gs", NULL, NULL, 1 },
	{ "gs-mo", NULL, NULL, 1 },

	{ "ma-card", NULL, NULL, 1 },
	{ "ma-cont", NULL, NULL, 1 },

	{ "focus", "error" },

	{ "3ddp", NULL, NULL, 10, 1 },
	{ "3dvp", NULL, NULL, 10, 1 },

	{ "3dv", NULL, NULL, 2, 1 },
	{ "3dd", NULL, NULL, 2, 1 },

	{ "view-angle", NULL, NULL, 1 },

	{ "cpf-load", NULL, NULL, 1 },

	{ "ui", NULL, NULL, 1 },

	{ "3d-fmr", NULL, NULL, 1 },

	{ "3d-dmr", NULL, NULL, 1 },

	{ "et", NULL, NULL, 1 },
	
	{ "st", NULL, NULL, 1 },

	{ "di-lower", NULL, NULL, 1 },

	{ "rc", NULL, NULL, 1 },

	{ "do-try", NULL, NULL, 1 },
	{ "di-upper", NULL, NULL, 1 },

	{ "fc", NULL, NULL, 1 },

	{ "ecm", NULL, NULL, 1 },
	{ "acm", NULL, NULL, 1 },

	{ "def-nn", NULL, NULL, 1 },
	
	{ "fx", NULL, NULL, 1 },

	{ "tcem", NULL, NULL, 1 },

	{ "oui", NULL, NULL, 1 },

	{ "pa", NULL, NULL, 1 },

	{ "pc", NULL, NULL, 1 },

	{ "cw" },
	{ "no-cw", "cw", "0" },

	{ "wm", NULL, NULL, 3 },

	{ "fl", NULL, NULL, 3 },

	{ "wmx", NULL, NULL, 3 },

	{ "flshow", NULL, NULL, 1 },

	{ "3dpx", NULL, NULL, 6 },

	{ "ex", NULL, NULL, 6, 1 },
	{ "crop", NULL, NULL, 6, 1 },
	{ "fex", NULL, NULL, 6, 1 },
	{ "fcrop", NULL, NULL, 6, 1 },
	{ "exshow" },

	{ "wt", NULL, NULL, 1 },
	{ "3d-chain", NULL, NULL, 1 },
	{ "dchain", NULL, NULL, 1 },
	{ "achain", NULL, NULL, 1 },
	{ "afilter", NULL, NULL, 1 },
	{ "ochain", NULL, NULL, 2 },

	{ "visp", NULL, NULL, 5 },
	{ "cx", NULL, NULL, 1 },
	{ "no-cx", "cx", "0" },

	{ "ip", NULL, NULL, 0 },

	{ "bayer", NULL, NULL, 1 },

	{ "lpsf", NULL, NULL, 1 },
	{ "nlpsf", NULL, NULL, 1 }, 
	{ "psf-match", NULL, NULL, 6 },

	{ "device", NULL, NULL, 1 },

	{ "usm", NULL, NULL, 1 },

	{ "ipr", NULL, NULL, 1 },

	{ "cpp-err-median", "cpp-err", "median" },
	{ "cpp-err-mean", "cpp-err", "mean" },

	{ "vp-adjust" },
	{ "vp-noadjust", "vp-adjust", "0" },

	{ "vo-adjust" },
	{ "vo-noadjust", "vo-adjust", "0" },

	{ "ip-mean", "ip-statistic", "mean" },
	{ "ip-median", "ip-statistic", "median" },

	{ "ip-wl", "ip-wl", "1", 1 },
	{ "ip-nowl", "ip-wl", "0" },

	{ "ips", NULL, NULL, 1 },

	{ "ipc", NULL, NULL, 2 },

	{ "exp-extend" },
	{ "exp-noextend", "exp-extend", "0" },

	{ "exp-register" },
	{ "exp-noregister", "exp-register", "0" },
	{ "exp-meta-only", "exp-register", "2" },

	{ "drizzle-only" },

	{ "subspace-traverse" },

	{ "3d-filter" },
	{ "3d-nofilter", "3d-filter", "0" },

	{ "occ-norm" },
	{ "occ-nonorm", "occ-norm", "0" },

	{ "inc" },
	{ "no-inc", "inc", "0" },

	{ "exp-mult", NULL, NULL, 3 },
	
	{ "visp-scale", NULL, NULL, 1 },

	{ "scale", NULL, NULL, 1 },

	{ "metric", NULL, NULL, 1 },

	{ "threshold", NULL, NULL, 1 },

	{ "drizzle-diam", NULL, NULL, 1 },

	{ "perturb-upper", NULL, NULL, 1 },
	{ "perturb-lower", NULL, NULL, 1 },

	{ "stepsize", NULL, NULL, 1 },

	{ "va-upper", NULL, NULL, 1 },
	{ "cpp-upper", NULL, NULL, 1 },
	{ "cpp-lower", NULL, NULL, 1 }, 
	{ "hf-enhance", NULL, NULL, 1 },
	{ "rot-upper", NULL, NULL, 1 },

	{ "bda-mult", NULL, NULL, 1 },
	{ "bda-rate", NULL, NULL, 1 },

	{ "lod-max", NULL, NULL, 1 },

	{ "cpf-load", NULL, NULL, 1 },
	{ "model-load", NULL, NULL, 1 },
	{ "model-save", NULL, NULL, 1 },
	{ "trans-load", NULL, NULL, 1 },
	{ "trans-save", NULL, NULL, 1 },
	{ "3d-trans-load", NULL, NULL, 1 },
	{ "3d-trans-save", NULL, NULL, 1 },

	/*
	 * End of table.
	 */

	{ NULL }
};
