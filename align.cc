#include "align.h"

/*
 * See align.h for details on these variables.
 */

double align::scale_factor;

int align::_keep = 0;
transformation *align::kept_t = NULL;
int *align::kept_ok = NULL;

tload_t *align::tload = NULL;
tsave_t *align::tsave = NULL;
render *align::reference = NULL;
const image *align::reference_image = NULL;
const image *align::input_frame = NULL;
const image_weights *align::reference_defined = NULL;
transformation align::latest_t;
int align::latest_ok;
int align::latest = 0;
int align::extend = 0;
int align::extend_orig_height = 0;
int align::extend_orig_width = 0;

int align::alignment_class = 1;
int align::default_initial_alignment_type = 0;
transformation align::default_initial_alignment;
int align::is_default = 1;
int align::old_is_default;
transformation align::old_initial_alignment;
transformation align::old_final_alignment;
int align::channel_alignment_type = 2;
double align::metric_exponent = 2;
double align::match_threshold = 0;
double align::perturb_lower = 0.125;
double align::perturb_upper = 32;

/*
 * Empirically, it's okay to use a level-of-detail equal to twice the
 * resolution of the perturbation, so we set the default lod_max to 1, as
 * 2^1==2.  lod_max of zero seems okay also, but lower values seem to cause
 * problems.
 */

int align::lod_max = 1;

double align::rot_max = 32.0;
double align::match_sum = 0;
int align::match_count = 0;

double align::_mc = 0;
