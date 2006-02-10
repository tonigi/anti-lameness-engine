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
 * d3/scene.h: Representation of a 3D scene.
 */

#ifndef __scene_h__
#define __scene_h__

#include "point.h"

/*
 * View angle multiplier.  
 *
 * Setting this to a value larger than one can be useful for debugging.
 */

#define VIEW_ANGLE_MULTIPLIER 1

class scene {

	/*
	 * Clipping planes
	 */
	static ale_pos front_clip;
	static ale_pos rear_clip;

	/*
	 * Decimation exponents for geometry
	 */
	static int primary_decimation_upper;
	static int input_decimation_lower;
	static int output_decimation_preferred;

	/*
	 * Output clipping
	 */
	static int output_clip;

	/*
	 * Model files
	 */
	static const char *load_model_name;
	static const char *save_model_name;

	/*
	 * Occupancy attenuation
	 */

	static double occ_att;

	/*
	 * Normalization of output by weight
	 */

	static int normalize_weights;

	/*
	 * Falloff exponent
	 */

	static double falloff_exponent;

	/*
	 * Third-camera error multiplier
	 */
	static double tc_multiplier;

	/*
	 * Occupancy update iterations
	 */
	static unsigned int ou_iterations;

	/*
	 * Pairwise ambiguity
	 */
	static unsigned int pairwise_ambiguity;

	/*
	 * Pairwise comparisons
	 */
	static const char *pairwise_comparisons;

	/*
	 * 3D Post-exclusion
	 */
	static int d3px_count;
	static double *d3px_parameters;

	/*
	 * Nearness threshold
	 */
	static const ale_real nearness;

	/*
	 * Encounter threshold for defined pixels.
	 */
	static double encounter_threshold;

	/*
	 * Structure to hold input frame information at levels of 
	 * detail between full detail and full decimation.
	 */
	class lod_image {
		unsigned int f;
		unsigned int entries;
		std::vector<const d2::image *> im;
		std::vector<pt> transformation;

	public:
		/*
		 * Constructor
		 */
		lod_image(unsigned int _f) {

			pt _pt;
			
			f = _f;
			im.push_back(d2::image_rw::copy(f, "3D reference image"));
			assert(im.back());
			entries = 1;
			_pt = d3::align::projective(f);
			_pt.scale(1 / _pt.scale_2d());
			transformation.push_back(_pt);

			while (im.back()->height() > 4
			    && im.back()->width() > 4) {

				im.push_back(im.back()->scale_by_half("3D, reduced LOD"));
				assert(im.back());

				_pt.scale(1 / _pt.scale_2d() / pow(2, entries));
				transformation.push_back(_pt);

				entries += 1;
			}
		}

		/*
		 * Get the number of scales
		 */
		unsigned int count() {
			return entries;
		}

		/*
		 * Get an image.
		 */
		const d2::image *get_image(unsigned int i) {
			assert(i >= 0);
			assert(i < entries);
			return im[i];
		}

		int in_bounds(d2::point p) {
			return im[0]->in_bounds(p);
		}

		/*
		 * Get a 'trilinear' color.  We currently don't do interpolation
		 * between levels of detail; hence, it's discontinuous in tl_coord.
		 */
		d2::pixel get_tl(d2::point p, ale_pos tl_coord) {

			assert(in_bounds(p));

			tl_coord = round(tl_coord);

			if (tl_coord >= entries)
				tl_coord = entries;
			if (tl_coord < 0)
				tl_coord = 0;

			p = p / (ale_pos) pow(2, tl_coord);

			unsigned int itlc = (unsigned int) tl_coord;

			if (p[0] > im[itlc]->height() - 1)
				p[0] = im[itlc]->height() - 1;
			if (p[1] > im[itlc]->width() - 1)
				p[1] = im[itlc]->width() - 1;

			assert(p[0] >= 0);
			assert(p[1] >= 0);

			return im[itlc]->get_bl(p);
		}

		d2::pixel get_max_diff(d2::point p, ale_pos tl_coord) {
			assert(in_bounds(p));

			tl_coord = round(tl_coord);

			if (tl_coord >= entries)
				tl_coord = entries;
			if (tl_coord < 0)
				tl_coord = 0;

			p = p / (ale_pos) pow(2, tl_coord);

			unsigned int itlc = (unsigned int) tl_coord;

			if (p[0] > im[itlc]->height() - 1)
				p[0] = im[itlc]->height() - 1;
			if (p[1] > im[itlc]->width() - 1)
				p[1] = im[itlc]->width() - 1;

			assert(p[0] >= 0);
			assert(p[1] >= 0);

			return im[itlc]->get_max_diff(p);
		}

		/*
		 * Get the transformation
		 */
		pt get_t(unsigned int i) {
			assert(i >= 0);
			assert(i < entries);
			return transformation[i];
		}

		/*
		 * Get the camera origin in world coordinates
		 */
		point origin() {
			return transformation[0].origin();
		}

		/*
		 * Destructor
		 */
		~lod_image() {
			for (unsigned int i = 0; i < entries; i++) {
				delete im[i];
			}
		}
	};

	/*
	 * Structure to hold weight information for reference images.
	 */
	class ref_weights {
		unsigned int f;
		std::vector<d2::image *> weights;
		pt transformation;
		int tc_low;
		int tc_high;
		int initialized;

		void set_image(d2::image *im, ale_real value) {
			assert(im);
			for (unsigned int i = 0; i < im->height(); i++)
			for (unsigned int j = 0; j < im->width(); j++)
				im->pix(i, j) = d2::pixel(value, value, value);
		}

		d2::image *make_image(ale_pos sf, ale_real init_value = 0) {
			d2::image *result = new d2::image_ale_real(
					(unsigned int) ceil(transformation.unscaled_height() * sf),
					(unsigned int) ceil(transformation.unscaled_width() * sf), 3);
			assert(result);

			if (init_value != 0)
				set_image(result, init_value);

			return result;
		}

	public:

		/*
		 * Explicit weight subtree
		 */
		struct subtree {
			ale_real node_value;
			subtree *children[2][2];

			subtree(ale_real nv, subtree *a, subtree *b, subtree *c, subtree *d) {
				node_value = nv;
				children[0][0] = a;
				children[0][1] = b;
				children[1][0] = c;
				children[1][1] = d;
			}

			~subtree() {
				for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					delete children[i][j];
			}
		};

		/*
		 * Constructor
		 */
		ref_weights(unsigned int _f) {
			f = _f;
			transformation = d3::align::projective(f);
			initialized = 0;
		}

		/*
		 * Check spatial bounds.
		 */
		int in_spatial_bounds(point p) {

			if (!p.defined())
				return 0;

			if (p[0] < 0)
				return 0;
			if (p[1] < 0)
				return 0;
			if (p[0] > transformation.unscaled_height() - 1)
				return 0;
			if (p[1] > transformation.unscaled_width() - 1)
				return 0;
			if (p[2] >= 0)
				return 0;

			return 1;
		}

		int in_spatial_bounds(const space::traverse &t) {
			point p = transformation.centroid(t);
			return in_spatial_bounds(p);
		}

		/*
		 * Increase resolution to the given level.
		 */
		void increase_resolution(int tc, unsigned int i, unsigned int j) {
			d2::image *im = weights[tc - tc_low];
			assert(im);
			assert(i <= im->height() - 1);
			assert(j <= im->width() - 1);

			/*
			 * Check for the cases known to have no lower level of detail.
			 */

			if (im->pix(i, j)[0] == -1)
				return;

			if (tc == tc_high)
				return;

			increase_resolution(tc + 1, i / 2, j / 2);

			/*
			 * Load the lower-level image structure.
			 */

			d2::image *iim = weights[tc + 1 - tc_low];

			assert(iim);
			assert(i / 2 <= iim->height() - 1);
			assert(j / 2 <= iim->width() - 1);

			/*
			 * Check for the case where no lower level of detail is
			 * available.
			 */

			if (iim->pix(i / 2, j / 2)[0] == -1)
				return;

			/*
			 * Spread out the lower level of detail among (uninitialized)
			 * peer values.
			 */

			for (unsigned int ii = (i / 2) * 2; ii < (i / 2) * 2 + 1; ii++)
			for (unsigned int jj = (j / 2) * 2; jj < (j / 2) * 2 + 1; jj++) {
				assert(ii <= im->height() - 1);
				assert(jj <= im->width() - 1);
				assert(im->pix(ii, jj)[0] == 0);

				im->pix(ii, jj) = iim->pix(i / 2, j / 2);
			}

			/*
			 * Set the lower level of detail to point here.
			 */

			iim->pix(i / 2, j / 2)[0] = -1;
		}

		/*
		 * Add weights to positive higher-resolution pixels where
		 * found when their current values match the given subtree
		 * values; set negative pixels to zero and return 0 if no
		 * positive higher-resolution pixels are found.  
		 */
		int add_partial(int tc, unsigned int i, unsigned int j, ale_real weight, subtree *st) {
			d2::image *im = weights[tc - tc_low];
			assert(im);
			assert(i <= im->height() - 1);
			assert(j <= im->width() - 1);

			/*
			 * Check for positive values.
			 */

			if (im->pix(i, j)[0] > 0) {
				if (st && st->node_value == im->pix(i, j)[0])
					im->pix(i, j)[0] += weight;
				return 1;
			}

			/*
			 * Handle the case where there are no higher levels of detail.
			 */

			if (tc == tc_low) {
				assert(im->pix(i, j)[0] == 0);
				return 0;
			}

			/*
			 * Handle the case where higher levels of detail are available.
			 */

			int success[2][2];

			for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 2; jj++)
				success[ii][jj] = add_partial(tc - 1, i * 2 + ii, j * 2 + jj, weight, 
						st ? st->children[ii][jj] : NULL);

			if (!success[0][0]
			 && !success[0][1]
			 && !success[1][0]
			 && !success[1][1]) {
				im->pix(i, j)[0] = 0;
				return 0;
			}

			d2::image *iim = weights[tc - 1 - tc_low];
			assert(iim);

			for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 2; jj++)
				if (success[ii][jj] == 0) {
					assert(i * 2 + ii < iim->height());
					assert(j * 2 + jj < iim->width());

					im->pix(i * 2 + ii, j * 2 + jj)[0] = weight;
				}

			im->pix(i, j)[0] = -1;

			return 1;
		}

		/*
		 * Add weight.
		 */
		void add_weight(int tc, unsigned int i, unsigned int j, ale_real weight, subtree *st) {
			d2::image *im = weights[tc - tc_low];
			assert(im);
			assert(i <= im->height() - 1);
			assert(j <= im->width() - 1);

			/*
			 * Increase resolution, if necessary
			 */

			increase_resolution(tc, i, j);

			/*
			 * Attempt to add the weight at levels of detail
			 * where weight is defined.
			 */

			if (add_partial(tc, i, j, weight, st))
				return;

			/*
			 * If no weights are defined at any level of detail,
			 * then set the weight here.
			 */

			im->pix(i, j)[0] = weight;
		}

		void add_weight(int tc, d2::point p, ale_real weight, subtree *st) {

			p *= pow(2, -tc);

			unsigned int i = (unsigned int) floor(p[0]);
			unsigned int j = (unsigned int) floor(p[1]);

			add_weight(tc, i, j, weight, st);
		}

		void add_weight(const space::traverse &t, ale_real weight, subtree *st) {

			if (weight == 0)
				return;
			
			ale_pos tc = transformation.trilinear_coordinate(t);
			point p = transformation.centroid(t);
			assert(in_spatial_bounds(p));

			tc = round(tc);

			/*
			 * Establish a reasonable (?) upper bound on resolution.
			 */

			if (tc < input_decimation_lower) {
				weight /= pow(4, (input_decimation_lower - tc));
				tc = input_decimation_lower;
			}

			/*
			 * Initialize, if necessary.
			 */

			if (!initialized) {
				tc_low = tc_high = (int) tc;

				ale_real sf = pow(2, -tc);

				weights.push_back(make_image(sf));
			}

			/*
			 * Check resolution bounds
			 */

			assert (tc_low <= tc_high);

			/*
			 * Generate additional levels of detail, if necessary.
			 */

			while (tc < tc_low) {
				tc_low--;

				ale_real sf = pow(2, -tc_low);

				weights.insert(weights.begin(), make_image(sf));
			}

			while (tc > tc_high) {
				tc_high++;

				ale_real sf = pow(2, -tc_high);

				weights.push_back(make_image(sf, -1));
			}

			add_weight((int) tc, p.xy(), weight, st);
		}

		/*
		 * Get weight
		 */
		ale_real get_weight(int tc, unsigned int i, unsigned int j) {
			d2::image *im = weights[tc - tc_low];
			assert(im);
			assert(i < im->height());
			assert(j < im->width());

			if (im->pix(i, j)[0] == -1) {
				return (get_weight(tc - 1, i * 2 + 0, j * 2 + 0)
				      + get_weight(tc - 1, i * 2 + 1, j * 2 + 0)
				      + get_weight(tc - 1, i * 2 + 1, j * 2 + 1)
				      + get_weight(tc - 1, i * 2 + 0, j * 2 + 1)) / 4;
			}

			if (im->pix(i, j)[0] == 0) {
				return get_weight(tc + 1, i / 2, j / 2);
			}

			return im->pix(i, j)[0];
		}

		ale_real get_weight(int tc, d2::point p) {

			p *= pow(2, -tc);

			unsigned int i = (unsigned int) floor(p[0]);
			unsigned int j = (unsigned int) floor(p[1]);

			return get_weight(tc, i, j);
		}

		ale_real get_weight(const space::traverse &t) {
			ale_pos tc = transformation.trilinear_coordinate(t);
			point p = transformation.centroid(t);
			assert(in_spatial_bounds(p));

			if (!initialized)
				return 0;

			tc = round(tc);

			if (tc < tc_low) {
				tc = tc_low;
			}

			if (tc <= tc_high) {
				return get_weight((int) tc, p.xy());
			}

			assert(tc > tc_high);

			int multiplier = (int) pow(2, (tc - tc_high));
			ale_real result = 0;

			for (int i = 0; i < multiplier; i++)
			for (int j = 0; j < multiplier; j++) {
				result += get_weight(tc_high,
						(unsigned int) floor(p[0]) * multiplier + i,
						(unsigned int) floor(p[1]) * multiplier + j);
			}

			return result;
		}

		/*
		 * Check whether a subtree is simple.
		 */
		int is_simple(subtree *s) {
			assert (s);

			if (s->node_value == -1
			 && s->children[0][0] == NULL
			 && s->children[0][1] == NULL
			 && s->children[1][0] == NULL
			 && s->children[1][1] == NULL)
				return 1;

			return 0;
		}

		/*
		 * Get a weight subtree.
		 */
		subtree *get_subtree(int tc, unsigned int i, unsigned int j) {

			/*
			 * tc > tc_high is handled recursively.
			 */

			if (tc > tc_high) {
				subtree *result = new subtree(-1, 
						get_subtree(tc - 1, i * 2 + 0, j * 2 + 0),
						get_subtree(tc - 1, i * 2 + 0, j * 2 + 1),
						get_subtree(tc - 1, i * 2 + 1, j * 2 + 0),
						get_subtree(tc - 1, i * 2 + 1, j * 2 + 1));

				if (is_simple(result)) {
					delete result;
					return NULL;
				}

				return result;
			}

			assert(tc >= tc_low);
			assert(weights[tc - tc_low]);

			d2::image *im = weights[tc - tc_low];

			/*
			 * Rectangular images will, in general, have
			 * out-of-bounds tree sections.  Handle this case.
			 */

			if (i >= im->height())
				return NULL;
			if (j >= im->width())
				return NULL;

			/*
			 * -1 weights are handled recursively
			 */

			if (im->pix(i, j)[0] == -1) {
				subtree *result = new subtree(-1, 
						get_subtree(tc - 1, i * 2 + 0, j * 2 + 0),
						get_subtree(tc - 1, i * 2 + 0, j * 2 + 1),
						get_subtree(tc - 1, i * 2 + 1, j * 2 + 0),
						get_subtree(tc - 1, i * 2 + 1, j * 2 + 1));

				if (is_simple(result)) {
					im->pix(i, j)[0] = 0;
					delete result;
					return NULL;
				}

				return result;
			}

			/*
			 * Zero weights have NULL subtrees.
			 */

			if (im->pix(i, j)[0] == 0)
				return NULL;

			/*
			 * Handle the remaining case.
			 */

			return new subtree(im->pix(i, j)[0], NULL, NULL, NULL, NULL);
		}

		subtree *get_subtree(int tc, d2::point p) {
			p *= pow(2, -tc);

			unsigned int i = (unsigned int) floor(p[0]);
			unsigned int j = (unsigned int) floor(p[1]);

			return get_subtree(tc, i, j);
		}

		subtree *get_subtree(const space::traverse &t) {
			ale_pos tc = transformation.trilinear_coordinate(t);
			point p = transformation.centroid(t);
			assert(in_spatial_bounds(p));

			if (!initialized)
				return NULL;

			if (tc < input_decimation_lower)
				tc = input_decimation_lower;

			tc = round(tc);

			if (tc < tc_low)
				return NULL;

			return get_subtree((int) tc, p.xy());
		}

		/*
		 * Destructor
		 */
		~ref_weights() {
			for (unsigned int i = 0; i < weights.size(); i++) {
				delete weights[i];
			}
		}
	};

	/*
	 * Resolution check.
	 */
	static int resolution_ok(pt transformation, ale_pos tc) {

		if (pow(2, tc) > transformation.unscaled_height()
		 || pow(2, tc) > transformation.unscaled_width())
			return 0;

		if (tc < input_decimation_lower)
			return 0;

		return 1;
	}

	/*
	 * Structure to hold input frame information at all levels of detail.
	 */
	class lod_images {

		/*
		 * All images.
		 */

		std::vector<lod_image *> images;

	public:

		lod_images() {
			images.resize(d2::image_rw::count(), NULL);
		}

		void open(unsigned int f) {
			assert (images[f] == NULL);

			if (images[f] == NULL)
				images[f] = new lod_image(f);
		}

		void open_all() {
			for (unsigned int f = 0; f < d2::image_rw::count(); f++)
				open(f);
		}

		lod_image *get(unsigned int f) {
			assert (images[f] != NULL);
			return images[f];
		}

		void close(unsigned int f) {
			assert (images[f] != NULL);
			delete images[f];
			images[f] = NULL;
		}

		void close_all() {
			for (unsigned int f = 0; f < d2::image_rw::count(); f++)
				close(f);
		}

		~lod_images() {
			close_all();
		}
	};

	/*
	 * All levels-of-detail
	 */

	static struct lod_images *al;

	/*
	 * Data structure for storing best encountered subspace candidates.
	 */
	class candidates {
		std::vector<std::vector<std::pair<ale_pos, ale_real> > > levels;
		int image_index;
		unsigned int height;
		unsigned int width;

		/*
		 * Point p is in world coordinates.
		 */
		void generate_subspace(point iw, ale_pos diagonal) {

			fprintf(stderr, "*");

			space::traverse st = space::traverse::root();

			if (!st.includes(iw)) {
				assert(0);
				return;
			}

			int highres = 0;
			int lowres = 0;

			/*
			 * Loop until resolutions of interest have been generated.
			 */
			
			for(;;) {

				ale_pos current_diagonal = (st.get_max() - st.get_min()).norm();

				assert(!isnan(current_diagonal));

				/*
				 * Generate any new desired spatial registers.
				 */

				/*
				 * Inputs
				 */

				for (int f = 0; f < 2; f++) {

					/*
					 * Low resolution
					 */

					if (current_diagonal < 2 * diagonal
					 && lowres == 0) {
						spatial_info_map[st.get_node()];
						lowres = 1;
					}

					/*
					 * High resolution.
					 */

					if (current_diagonal < diagonal
					 && highres == 0) {
						spatial_info_map[st.get_node()];
						highres = 1;
					}
				}

				/*
				 * Check for completion
				 */

				if (highres && lowres)
					return;

				/*
				 * Check precision before analyzing space further.
				 */

				if (st.precision_wall()) {
					fprintf(stderr, "\n\n*** Error: reached subspace precision wall ***\n\n");
					assert(0);
					return;
				}

				if (st.positive().includes(iw)) {
					st = st.positive();
					total_tsteps++;
				} else if (st.negative().includes(iw)) {
					st = st.negative();
					total_tsteps++;
				} else {
					fprintf(stderr, "failed iw = (%f, %f, %f)\n", 
							iw[0], iw[1], iw[2]);
					assert(0);
				}
			}
		}

	public:
		candidates(int f) {

			image_index = f;
			height = (unsigned int) al->get(f)->get_t(0).unscaled_height();
			width = (unsigned int) al->get(f)->get_t(0).unscaled_width();

			/*
			 * Is this necessary?
			 */

			levels.resize(primary_decimation_upper - input_decimation_lower + 1);
			for (int l = input_decimation_lower; l <= primary_decimation_upper; l++) {
				levels[l].resize((unsigned int) (floor(height / pow(2, l))
					                       * floor(width / pow(2, l))
							       * pairwise_ambiguity),
						 std::pair<ale_pos, ale_real>(0, 0));
			}
		}

		/*
		 * Point p is expected to be in local projective coordinates.
		 */

		void add_candidate(point p, int tc, ale_real score) {
			assert(tc <= primary_decimation_upper);
			assert(tc >= input_decimation_lower);
			assert(p[2] < 0);
			assert(score >= 0);

			int i = (unsigned int) floor(p[0] / pow(2, tc));
			int j = (unsigned int) floor(p[1] / pow(2, tc));

			assert(i < floor(height / pow(2, tc)));
			assert(j < floor(width / pow(2, tc)));

			for (unsigned int k = 0; k < pairwise_ambiguity; k++) {
				std::pair<ale_pos, ale_real> *pk =
					&(levels[tc][i * width * pairwise_ambiguity + j * pairwise_ambiguity + k]);

				if (pk->first != 0 && score >= pk->second)
					continue;

				ale_pos tp = pk->first;
				ale_real tr = pk->second;

				pk->first = p[2];
				pk->second = score;

				p[2] = tp;
				score = tr;

				if (p[2] == 0)
					break;
			}
		}

		/*
		 * Generate subspaces for candidates.
		 */

		void generate_subspaces() {
			fprintf(stderr, "+");
			for (int l = input_decimation_lower; l <= primary_decimation_upper; l++) {
				for (unsigned int i = 0; i < (unsigned int) floor(height / pow(2, l)); i++)
				for (unsigned int j = 0; j < (unsigned int) floor(width / pow(2, l)); j++)
				for (unsigned int k = 0; k < pairwise_ambiguity; k++) {
					std::pair<ale_pos, ale_real> *pk =
						&(levels[l - input_decimation_lower]
						        [i * width * pairwise_ambiguity + j * pairwise_ambiguity + k]);

					if (pk->first == 0) {
						fprintf(stderr, "o");
						continue;
					} else {
						fprintf(stderr, "|");
					}

					ale_pos si = i * pow(2, l) + ((l > 0) ? pow(2, l - 1) : 0);
					ale_pos sj = j * pow(2, l) + ((l > 0) ? pow(2, l - 1) : 0);

					point p = al->get(image_index)->get_t(0).pw_unscaled(point(si, sj, pk->first));

					generate_subspace(p, 
							al->get(image_index)->get_t(0).diagonal_distance_3d(pk->first, l));
				}
			}
		}
	};

	/*
	 * List for calculating weighted median.
	 */
	class wml {
		ale_real *data;
		unsigned int size;
		unsigned int used;

		ale_real &_w(unsigned int i) {
			assert(i <= used);
			return data[i * 2];
		}

		ale_real &_d(unsigned int i) {
			assert(i <= used);
			return data[i * 2 + 1];
		}

		void increase_capacity() {

			if (size > 0)
				size *= 2;
			else
				size = 1;

			data = (ale_real *) realloc(data, sizeof(ale_real) * 2 * (size * 1));

			assert(data);
			assert (size > used);

			if (!data) {
				fprintf(stderr, "Unable to allocate %d bytes of memory\n",
						sizeof(ale_real) * 2 * (size * 1));
				exit(1);
			} 
		}

		void insert_weight(unsigned int i, ale_real v, ale_real w) {
			assert(used < size);
			assert(used >= i);
			for (unsigned int j = used; j > i; j--) {
				_w(j) = _w(j - 1);
				_d(j) = _d(j - 1);
			}

			_w(i) = w;
			_d(i) = v;

			used++;
		}

	public:

		unsigned int get_size() {
			return size;
		}

		unsigned int get_used() {
			return used;
		}

		void print_info() {
			fprintf(stderr, "[st %p sz %u el", this, size);
			for (unsigned int i = 0; i < used; i++)
				fprintf(stderr, " (%f, %f)", _d(i), _w(i));
			fprintf(stderr, "]\n");
		}

		void clear() {
			used = 0;
		}

		void insert_weight(ale_real v, ale_real w) {
			for (unsigned int i = 0; i < used; i++) {
				if (_d(i) == v) {
					_w(i) += w;
					return;
				}
				if (_d(i) > v) {
					if (used == size)
						increase_capacity();
					insert_weight(i, v, w);
					return;
				}
			}
			if (used == size)
				increase_capacity();
			insert_weight(used, v, w);
		}

		/*
		 * Finds the median at half-weight, or between half-weight
		 * and zero-weight, depending on the attenuation value.
		 */

		ale_real find_median(double attenuation = 0) {

			assert(attenuation >= 0);
			assert(attenuation <= 1);

			ale_real zero1 = 0;
			ale_real zero2 = 0;
			ale_real undefined = zero1 / zero2;

			ale_accum weight_sum = 0;

			for (unsigned int i = 0; i < used; i++)
				weight_sum += _w(i);

//			if (weight_sum == 0)
//				return undefined;

			if (used == 0 || used == 1)
				return undefined;

			if (weight_sum == 0) {
				ale_accum data_sum = 0;
				for (unsigned int i = 0; i < used; i++)
					data_sum += _d(i);
				return data_sum / used;
			}
					

			ale_accum midpoint = weight_sum * (0.5 - 0.5 * attenuation);

			ale_accum weight_sum_2 = 0;

			for (unsigned int i = 0; i < used && weight_sum_2 < midpoint; i++) {
				weight_sum_2 += _w(i);

				if (weight_sum_2 > midpoint) {
					return _d(i);
				} else if (weight_sum_2 == midpoint) {
					assert (i + 1 < used);
					return (_d(i) + _d(i + 1)) / 2;
				} 
			}

			return undefined;
		}

		wml(int initial_size = 0) {

//			if (initial_size == 0) {
//				initial_size = (int) (d2::image_rw::count() * 1.5);
//			}

			size = initial_size;
			used = 0;

			if (size > 0) {
				data = (ale_real *) malloc(size * sizeof(ale_real) * 2);
				assert(data);
			} else {
				data = NULL;
			}
		}

		/*
		 * copy constructor.  This is required to avoid undesired frees.
		 */

		wml(const wml &w) {
			size = w.size;
			used = w.used;
			data = (ale_real *) malloc(size * sizeof(ale_real) * 2);
			assert(data);

			memcpy(data, w.data, size * sizeof(ale_real) * 2);
		}

		~wml() {
			free(data);
		}
	};

	/*
	 * Class for information regarding spatial regions of interest.
	 *
	 * This class is configured for convenience in cases where sampling is
	 * performed using an approximation of the fine:box:1,triangle:2 chain.
	 * In this case, the *_1 variables would store the fine data and the
	 * *_2 variables would store the coarse data.  Other uses are also
	 * possible.
	 */

	class spatial_info {
		/*
		 * Map channel value --> weight.
		 */
		wml color_weights_1[3];
		wml color_weights_2[3];

		/*
		 * Current color.
		 */
		d2::pixel color;

		/*
		 * Map occupancy value --> weight.
		 */
		wml occupancy_weights_1;
		wml occupancy_weights_2;

		/*
		 * Current occupancy value.
		 */
		ale_real occupancy;

		/*
		 * pocc/socc density
		 */

		unsigned int pocc_density;
		unsigned int socc_density;

		/*
		 * Insert a weight into a list.
		 */
		void insert_weight(wml *m, ale_real v, ale_real w) {
			m->insert_weight(v, w);
		}

		/*
		 * Find the median of a weighted list.  Uses NaN for undefined.
		 */
		ale_real find_median(wml *m, double attenuation = 0) {
			return m->find_median(attenuation);
		}

	public:
		/*
		 * Constructor.
		 */
		spatial_info() {
			color = d2::pixel::zero();
			occupancy = 0;
			pocc_density = 0;
			socc_density = 0;
		}

		/*
		 * Accumulate color; primary data set.
		 */
		void accumulate_color_1(int f, d2::pixel color, d2::pixel weight) {
			for (int k = 0; k < 3; k++)
				insert_weight(&color_weights_1[k], color[k], weight[k]);
		}

		/*
		 * Accumulate color; secondary data set.
		 */
		void accumulate_color_2(d2::pixel color, d2::pixel weight) {
			for (int k = 0; k < 3; k++)
				insert_weight(&color_weights_2[k], color[k], weight[k]);
		}

		/*
		 * Accumulate occupancy; primary data set.
		 */
		void accumulate_occupancy_1(int f, ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_1, occupancy, weight);
		}

		/*
		 * Accumulate occupancy; secondary data set.
		 */
		void accumulate_occupancy_2(ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_2, occupancy, weight);
			
			if (occupancy == 0 || occupancy_weights_2.get_size() < 96)
				return;

			// fprintf(stderr, "%p updated socc with: %f %f\n", this, occupancy, weight);
			// occupancy_weights_2.print_info();
		}

		/*
		 * Update color (and clear accumulation structures).
		 */
		void update_color() {
			for (int d = 0; d < 3; d++) {
				ale_real c = find_median(&color_weights_1[d]);
				if (isnan(c))
					c = find_median(&color_weights_2[d]);
				if (isnan(c))
					c = 0;

				color[d] = c;

				color_weights_1[d].clear();
				color_weights_2[d].clear();
			}
		}

		/*
		 * Update occupancy (and clear accumulation structures).
		 */
		void update_occupancy() {
			ale_real o = find_median(&occupancy_weights_1, occ_att);
			if (isnan(o))
				o = find_median(&occupancy_weights_2, occ_att);
			if (isnan(o))
				o = 0;

			occupancy = o;

			pocc_density = occupancy_weights_1.get_used();
			socc_density = occupancy_weights_2.get_used();

			occupancy_weights_1.clear();
			occupancy_weights_2.clear();

		}

		/*
		 * Get current color.
		 */
		d2::pixel get_color() {
			return color;
		}

		/*
		 * Get current occupancy.
		 */
		ale_real get_occupancy() {
			return occupancy;
		}

		/*
		 * Get primary color density.
		 */

		unsigned int get_pocc_density() {
			return pocc_density;
		}

		unsigned int get_socc_density() {
			return socc_density;
		}
	};

	/*
	 * Map spatial regions of interest to spatial info structures.  XXX:
	 * This may get very poor cache behavior in comparison with, say, an
	 * array.  Unfortunately, there is no immediately obvious array
	 * representation.  If some kind of array representation were adopted,
	 * it would probably cluster regions of similar depth from the
	 * perspective of the typical camera.  In particular, for a
	 * stereoscopic view, depth ordering for two random points tends to be
	 * similar between cameras, I think.  Unfortunately, it is never
	 * identical for all points (unless cameras are co-located).  One
	 * possible approach would be to order based on, say, camera 0's idea
	 * of depth.
	 */

	typedef std::map<struct space::node *, spatial_info> spatial_info_map_t;

	static spatial_info_map_t spatial_info_map;

public:

	/*
	 * Debugging variables.
	 */

	static unsigned long total_ambiguity;
	static unsigned long total_pixels;
	static unsigned long total_divisions;
	static unsigned long total_tsteps;

	/*
	 * Member functions
	 */

	static void et(double et_parameter) {
		encounter_threshold = et_parameter;
	}

	static void load_model(const char *name) {
		load_model_name = name;
	}

	static void save_model(const char *name) {
		save_model_name = name;
	}

	static void fc(ale_pos fc) {
		front_clip = fc;
	}

	static void di_upper(ale_pos _dgi) {
		primary_decimation_upper = (int) round(_dgi);
	}

	static void do_try(ale_pos _dgo) {
		output_decimation_preferred = (int) round(_dgo);
	}

	static void di_lower(ale_pos _idiv) {
		input_decimation_lower = (int) round(_idiv);
	}

	static void oc() {
		output_clip = 1;
	}

	static void no_oc() {
		output_clip = 0;
	}

	static void rc(ale_pos rc) {
		rear_clip = rc;
	}

	/*
	 * Initialize 3D scene from 2D scene, using 2D and 3D alignment
	 * information.
	 */
	static void init_from_d2() {

		/*
		 * Rear clip value of 0 is converted to infinity.
		 */

		if (rear_clip == 0) {
			ale_pos one = +1;
			ale_pos zero = +0;

			rear_clip = one / zero;
			assert(isinf(rear_clip) == +1);
		}

		/*
		 * Scale and translate clipping plane depths.
		 */

		ale_pos cp_scalar = d3::align::projective(0).wc(point(0, 0, 0))[2];

		front_clip = front_clip * cp_scalar - cp_scalar;
		rear_clip = rear_clip * cp_scalar - cp_scalar;


		/*
		 * Allocate image structures.
		 */

		al = new lod_images;

		if (tc_multiplier != 0) {
			al->open_all();
		}
	}

	/*
	 * Perform spatial_info updating on a given subspace, for given
	 * parameters.
	 */
	static void subspace_info_update(space::iterate si, int f, ref_weights *weights) {
		while(!si.done()) {

			space::traverse st = si.get();

			/*
			 * Skip spaces with no color information.
			 *
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_node()) == 0) {
				si.next();
				continue;
			}

			/*
			 * Get in-bounds centroid, if one exists.
			 */

			point p = al->get(f)->get_t(0).centroid(st);

			if (!p.defined()) {
				si.next();
				continue;
			}
				
			/*
			 * Get information on the subspace.
			 */

			spatial_info *sn = &spatial_info_map[st.get_node()];
			d2::pixel color = sn->get_color();
			ale_real occupancy = sn->get_occupancy();

			/*
			 * Store current weight so we can later check for
			 * modification by higher-resolution subspaces.
			 */

			ref_weights::subtree *tree = weights->get_subtree(st);

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_node()->positive
			 || st.get_node()->negative) {

				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that later.
				 */

				space::iterate cleaved_space = si.cleave();

				cleaved_space.next();

				subspace_info_update(cleaved_space, f, weights);

			} else {
				si.next();
			}

			/*
			 * Add new data on the subspace and update weights.
			 */

			ale_pos tc = al->get(f)->get_t(0).trilinear_coordinate(st);
			d2::pixel pcolor = al->get(f)->get_tl(p.xy(), tc);
			d2::pixel colordiff = (color - pcolor) * (ale_real) 256;

			if (falloff_exponent != 0) {
				d2::pixel max_diff = al->get(f)->get_max_diff(p.xy(), tc) * (ale_real) 256;

				for (int k = 0; k < 3; k++)
				if (max_diff[k] > 1)
					colordiff[k] /= pow(max_diff[k], falloff_exponent);
			}

			/*
			 * Determine the probability of
			 * encounter, divided by the occupancy.
			 */

			d2::pixel encounter = d2::pixel(1, 1, 1) * (1 - weights->get_weight(st));

			/*
			 * Update weights
			 */

			weights->add_weight(st, (encounter * occupancy)[0], tree);

			/*
			 * Delete the subtree, if necessary.
			 */

			delete tree;

			/*
			 * Check for cases in which the subspace should not be
			 * updated.
			 */

			if (!resolution_ok(al->get(f)->get_t(0), tc))
				return;

			/*
			 * Update subspace.
			 */

			sn->accumulate_color_1(f, pcolor, encounter);
			d2::pixel channel_occ = pexp(-colordiff * colordiff);

			ale_accum occ = channel_occ[0];

			for (int k = 1; k < 3; k++)
				if (channel_occ[k] < occ)
					occ = channel_occ[k];

			sn->accumulate_occupancy_1(f, occ, encounter[0]);

		}
	}

	/*
	 * Run a single iteration of the spatial_info update cycle.
	 */
	static void spatial_info_update() {
		/*
		 * Iterate through each frame.
		 */
		for (unsigned int f = 0; f < d2::image_rw::count(); f++) {

			/*
			 * Open the frame and transformation.
			 */

			if (tc_multiplier == 0)
				al->open(f);

			/*
			 * Allocate weights data structure for storing encounter
			 * probabilities.
			 */

			ref_weights *weights = new ref_weights(f);

			/*
			 * Call subspace_info_update for the root space.
			 */

			subspace_info_update(space::iterate(al->get(f)->origin()), f, weights);

			/*
			 * Free weights.
			 */

			delete weights;

			/*
			 * Close the frame and transformation.
			 */

			if (tc_multiplier == 0)
				al->close(f);
		}

		/*
		 * Update all spatial_info structures.
		 */
		for (spatial_info_map_t::iterator i = spatial_info_map.begin(); i != spatial_info_map.end(); i++) {
			i->second.update_color();
			i->second.update_occupancy();

//			d2::pixel color = i->second.get_color();

//			fprintf(stderr, "space p=%p updated to c=[%f %f %f] o=%f\n",
//					i->first, color[0], color[1], color[2], 
//					i->second.get_occupancy());
		}
	}

	/*
	 * Support function for view() and depth().
	 */

	static const void view_recurse(int type, d2::image *im, d2::image *weights, space::iterate si, pt _pt) {
		while (!si.done()) {
			space::traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_node()) == 0) {
				si.next();
				continue;
			}

			spatial_info sn = spatial_info_map[st.get_node()];

			/*
			 * Get information on the subspace.
			 */

			d2::pixel color = sn.get_color();
			// d2::pixel color = d2::pixel(1, 1, 1) * (double) (((unsigned int) (st.get_node()) / sizeof(space)) % 65535);
			ale_real occupancy = sn.get_occupancy();

			/*
			 * Determine the view-local bounding box for the
			 * subspace.
			 */

			point bb[2];

			_pt.get_view_local_bb_scaled(st, bb);

			point min = bb[0];
			point max = bb[1];

			/*
			 * Data structure to check modification of weights by
			 * higher-resolution subspaces.
			 */

			std::queue<d2::pixel> weight_queue;

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_node()->positive
			 || st.get_node()->negative) {

				/*
				 * Store information about current weights,
				 * so we will know which areas have been
				 * covered by higher-resolution subspaces.
				 */

				for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
				for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++)
					weight_queue.push(weights->get_pixel(i, j));
				
				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that afterward.
				 */

				space::iterate cleaved_space = si.cleave();

				cleaved_space.next();

				view_recurse(type, im, weights, cleaved_space, _pt);

			} else {
				si.next();
			}
				

			/*
			 * Iterate over pixels in the bounding box, finding
			 * pixels that intersect the subspace.  XXX: assume
			 * for now that all pixels in the bounding box
			 * intersect the subspace.
			 */

			for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
			for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {

				/*
				 * Check for higher-resolution updates.
				 */

				if (weight_queue.size()) {
					if (weight_queue.front() != weights->get_pixel(i, j)) {
						weight_queue.pop();
						continue;
					}
					weight_queue.pop();
				}

				/*
				 * Determine the probability of encounter.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j)) * occupancy;

				/*
				 * Update images.
				 */

				if (type == 0) {

					/*
					 * Color view
					 */

					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * color;

				} else if (type == 1) {

					/*
					 * Weighted (transparent) depth display
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * depth_value;

				} else if (type == 2) {

					/*
					 * Ambiguity (ambivalence) measure.
					 */

					weights->pix(i, j) = d2::pixel(1, 1, 1);
					im->pix(i, j) += 0.1 * d2::pixel(1, 1, 1);

				} else if (type == 3) {

					/*
					 * Closeness measure.
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					if (weights->pix(i, j)[0] == 0) {
						weights->pix(i, j) = d2::pixel(1, 1, 1);
						im->pix(i, j) = d2::pixel(1, 1, 1) * depth_value;
					} else if (im->pix(i, j)[2] < depth_value) {
						im->pix(i, j) = d2::pixel(1, 1, 1) * depth_value;
					} else {
						continue;
					}

				} else if (type == 4) {

					/*
					 * Weighted (transparent) contribution display
					 */

					ale_pos contribution_value = sn.get_pocc_density() + sn.get_socc_density();
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

				} else if (type == 5) {

					/*
					 * Weighted (transparent) occupancy display
					 */

					ale_pos contribution_value = occupancy;
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

				} else if (type == 6) {
					
					/*
					 * (Depth, xres, yres) triple
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					weights->pix(i, j)[0] += encounter[0];
					if (weights->pix(i, j)[1] < encounter[0]) {
						weights->pix(i, j)[1] = encounter[0];
						im->pix(i, j)[0] = weights->pix(i, j)[1] * depth_value;
						im->pix(i, j)[1] = max[0] - min[0];
						im->pix(i, j)[2] = max[1] - min[1];
					}

				} else if (type == 7) {
					
					/*
					 * (xoff, yoff, 0) triple
					 */

					weights->pix(i, j)[0] += encounter[0];
					if (weights->pix(i, j)[1] < encounter[0]) {
						weights->pix(i, j)[1] = encounter[0];
						im->pix(i, j)[0] = i - min[0];
						im->pix(i, j)[1] = j - min[1];
						im->pix(i, j)[2] = 0;
					}

				} else 
					assert(0);
			}
		}
	}

	/*
	 * Generate an depth image from a specified view.
	 */
	static const d2::image *depth(pt _pt, int n = -1) {
		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		d2::image *im1 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		d2::image *im2 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		d2::image *im3 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space::iterate si(_pt.origin());

		view_recurse(6, im1, weights, si, _pt);

		delete weights;
		weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		view_recurse(7, im2, weights, si, _pt);

		/*
		 * Normalize depths by weights
		 */

		if (normalize_weights)
		for (unsigned int i = 0; i < im1->height(); i++)
		for (unsigned int j = 0; j < im1->width();  j++)
			im1->pix(i, j)[0] /= weights->pix(i, j)[1];

	
		for (unsigned int i = 0; i < im1->height(); i++)
		for (unsigned int j = 0; j < im1->width();  j++) {

			/*
			 * Handle interpolation.
			 */

			d2::point x;
			d2::point blx;
			d2::point res(im1->pix(i, j)[1], im1->pix(i, j)[2]);

			for (int d = 0; d < 2; d++) {

				if (im2->pix(i, j)[d] < res[d] / 2)
					x[d] = (ale_pos) (d?j:i) - res[d] / 2 - im2->pix(i, j)[d];
				else
					x[d] = (ale_pos) (d?j:i) + res[d] / 2 - im2->pix(i, j)[d];

				blx[d] = 1 - ((d?j:i) - x[d]) / res[d];
			}

			ale_real depth_val = 0;
			ale_real depth_weight = 0;

			for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 2; jj++) {
				d2::point p = x + d2::point(ii, jj) * res;
				if (im1->in_bounds(p)) {

					ale_real d = im1->get_bl(p)[0];

					if (isnan(d))
						continue;

					ale_real w = ((ii ? (1 - blx[0]) : blx[0]) * (jj ? (1 - blx[1]) : blx[1]));
					depth_weight += w;
					depth_val += w * d;
				}
			}

			ale_real depth = depth_val / depth_weight;

			/*
			 * Handle exclusions and encounter thresholds
			 */

			point w = _pt.pw_scaled(point(i, j, depth));

			if (weights->pix(i, j)[0] < encounter_threshold || excluded(w)) {
				im3->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
			} else {
				im3->pix(i, j) = d2::pixel(1, 1, 1) * depth;
			}
		}

		delete weights;
		delete im1;
		delete im2;

		return im3;
	}

	static const d2::image *depth(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return depth(_pt, n);
	}

	/*
	 * Generate an image from a specified view.
	 */
	static const d2::image *view(pt _pt, int n = -1) {
		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		const d2::image *depths = NULL;

		if (d3px_count > 0)
			depths = depth(_pt, n);

		d2::image *im = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space::iterate si(_pt.origin());

		view_recurse(0, im, weights, si, _pt);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			if (weights->pix(i, j).min_norm() < encounter_threshold
			 || (d3px_count > 0 && isnan(depths->pix(i, j)[0]))) {
				im->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
				weights->pix(i, j) = d2::pixel::zero();
			} else if (normalize_weights)
				im->pix(i, j) /= weights->pix(i, j);
		}

		if (d3px_count > 0)
			delete depths;

		delete weights;

		return im;
	}

	static void tcem(double _tcem) {
		tc_multiplier = _tcem;
	}

	static void oui(unsigned int _oui) {
		ou_iterations = _oui;
	}

	static void pa(unsigned int _pa) {
		pairwise_ambiguity = _pa;
	}

	static void pc(const char *_pc) {
		pairwise_comparisons = _pc;
	}

	static void d3px(int _d3px_count, double *_d3px_parameters) {
		d3px_count = _d3px_count;
		d3px_parameters = _d3px_parameters;
	}

	static void fx(double _fx) {
		falloff_exponent = _fx;
	}

	static void nw() {
		normalize_weights = 1;
	}

	static void no_nw() {
		normalize_weights = 0;
	}

	static int excluded(point p) {
		for (int n = 0; n < d3px_count; n++) {
			double *region = d3px_parameters + (6 * n);
			if (p[0] >= region[0]
			 && p[0] <= region[1]
			 && p[1] >= region[2]
			 && p[1] <= region[3]
			 && p[2] >= region[4]
			 && p[2] <= region[5])
				return 1;
		}

		return 0;
	}

	static const d2::image *view(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return view(_pt, n);
	}

	/*
	 * Add specified control points to the model
	 */
	static void add_control_points() {
	}

	typedef struct {point iw; point ip, is;} analytic;
	typedef std::multimap<ale_real,analytic> score_map;
	typedef std::pair<ale_real,analytic> score_map_element;

	/*
	 * Make pt list.
	 */
	static std::vector<pt> make_pt_list(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		std::vector<pt> result;

		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
			if (d_out[n] || v_out[n]) {
				result.push_back(align::projective(n));
			}
		}

		for (std::map<const char *, pt>::iterator i = d3_depth_pt->begin(); i != d3_depth_pt->end(); i++) {
			result.push_back(i->second);
		}

		for (std::map<const char *, pt>::iterator i = d3_output_pt->begin(); i != d3_output_pt->end(); i++) {
			result.push_back(i->second);
		}

		return result;
	}

	/*
	 * Get a trilinear coordinate for an anisotropic candidate cell.
	 */
	static ale_pos get_trilinear_coordinate(point min, point max, pt _pt) {

		d2::point local_min, local_max;

		local_min = _pt.wp_unscaled(min).xy();
		local_max = _pt.wp_unscaled(min).xy();

		point cell[2] = {min, max};

		/*
		 * Determine the view-local extrema in 2 dimensions.
		 */

		for (int r = 1; r < 8; r++) {
			point local = _pt.wp_unscaled(point(cell[r>>2][0], cell[(r>>1)%2][1], cell[r%2][2]));
			
			for (int d = 0; d < 2; d++) {
				if (local[d] < local_min[d])
					local_min[d] = local[d];
				if (local[d] > local_max[d])
					local_max[d] = local[d];
			}
		}

		ale_pos diameter = (local_max - local_min).norm();

		return log(diameter / sqrt(2)) / log(2);
	}

	/*
	 * Check whether a cell is visible from a given viewpoint.  This function
	 * is guaranteed to return 1 when a cell is visible, but it is not guaranteed
	 * to return 0 when a cell is invisible.
	 */
	static int pt_might_be_visible(const pt &viewpoint, point min, point max) {

		int doc = (rand() % 100000) ? 0 : 1;

		if (doc)
			fprintf(stderr, "checking visibility:\n");

		point cell[2] = {min, max};

		/*
		 * Cycle through all vertices of the cell to check certain
		 * properties.
		 */
		int pos[3] = {0, 0, 0};
		int neg[3] = {0, 0, 0};
		for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++) {
			point p = viewpoint.wp_unscaled(point(cell[i][0], cell[j][1], cell[k][2]));

			if (p[2] < 0 && viewpoint.unscaled_in_bounds(p))
				return 1;

			if (isnan(p[0])
			 || isnan(p[1])
			 || isnan(p[2]))
				return 1;

			if (p[2] > 0)
				for (int d = 0; d < 2; d++)
					p[d] *= -1;

			if (doc)
				fprintf(stderr, "\t[%f %f %f] --> [%f %f %f]\n", 
						cell[i][0], cell[j][1], cell[k][2],
						p[0], p[1], p[2]);

			for (int d = 0; d < 3; d++)
				if (p[d] >= 0)
					pos[d] = 1;

			if (p[0] <= viewpoint.unscaled_height() - 1)
				neg[0] = 1;

			if (p[1] <= viewpoint.unscaled_width() - 1)
				neg[1] = 1;

			if (p[2] <= 0)
				neg[2] = 1;
		}
		
		if (!neg[2])
			return 0;

		if (!pos[0]
		 || !neg[0]
		 || !pos[1]
		 || !neg[1])
			return 0;

		return 1;
	}

	/*
	 * Check whether a cell is output-visible.
	 */
	static int output_might_be_visible(const std::vector<pt> &pt_outputs, point min, point max) {
		for (unsigned int n = 0; n < pt_outputs.size(); n++)
			if (pt_might_be_visible(pt_outputs[n], min, max))
				return 1;
		return 0;
	}

	/*
	 * Check whether a cell is input-visible.
	 */
	static int input_might_be_visible(unsigned int f, point min, point max) {
		return pt_might_be_visible(align::projective(f), min, max);
	}

	/*
	 * Return true if a cell fails an output resolution bound.
	 */
	static int fails_output_resolution_bound(point min, point max, const std::vector<pt> &pt_outputs) {
		for (unsigned int n = 0; n < pt_outputs.size(); n++) {

			point p = pt_outputs[n].centroid(min, max);

			if (!p.defined())
				continue;

			if (get_trilinear_coordinate(min, max, pt_outputs[n]) < output_decimation_preferred)
				return 1;
		}
		
		return 0;
	}

	/*
	 * Check lower-bound resolution constraints 
	 */
	static int exceeds_resolution_lower_bounds(unsigned int f1, unsigned int f2,
			point min, point max, const std::vector<pt> &pt_outputs) {

		pt _pt = al->get(f1)->get_t(0);
		point p = _pt.centroid(min, max);

		if (get_trilinear_coordinate(min, max, _pt) < input_decimation_lower)
			return 1;

		if (fails_output_resolution_bound(min, max, pt_outputs))
			return 0;

		if (get_trilinear_coordinate(min, max, _pt) < primary_decimation_upper)
			return 1;

		return 0;
	}

	/*
	 * Try the candidate nearest to the specified cell.
	 */
	static void try_nearest_candidate(unsigned int f1, unsigned int f2, candidates *c, point min, point max) {
		point centroid = (max + min) / 2;
		pt _pt[2] = { al->get(f1)->get_t(0), al->get(f2)->get_t(0) };
		point p[2];

		/*
		 * Reject clipping plane violations.
		 */

		if (centroid[2] > front_clip
		 || centroid[2] < rear_clip)
			return;

		/*
		 * Calculate projections.
		 */

		for (int n = 0; n < 2; n++) {

			p[n] = _pt[n].wp_unscaled(centroid);

			if (!_pt[n].unscaled_in_bounds(p[n]))
				return;

			if (p[n][2] >= 0)
				return;
		}

		int tc = (int) round(get_trilinear_coordinate(min, max, _pt[0]));
		int stc = (int) round(get_trilinear_coordinate(min, max, _pt[1]));

		while (tc < input_decimation_lower || stc < input_decimation_lower) {
			tc++;
			stc++;
		}

		if (tc > primary_decimation_upper)
			return;

		/*
		 * Calculate score from color match.  Assume for now
		 * that the transformation can be approximated locally
		 * with a translation.
		 */

		ale_pos score = 0;
		ale_pos divisor = 0;
		ale_pos l1_multiplier = 0.125;
		lod_image *if1 = al->get(f1);
		lod_image *if2 = al->get(f2);

		if (if1->in_bounds(p[0].xy())
		 && if2->in_bounds(p[1].xy())) {
			divisor += 1 - l1_multiplier;
			score += (1 - l1_multiplier)
			       * (if1->get_tl(p[0].xy(), tc) - if2->get_tl(p[1].xy(), stc)).normsq();
		}

		for (int iii = -1; iii <= 1; iii++)
		for (int jjj = -1; jjj <= 1; jjj++) {
			d2::point t(iii, jjj);

			if (!if1->in_bounds(p[0].xy() + t)
			 || !if2->in_bounds(p[1].xy() + t))
				continue;

			divisor += l1_multiplier;
			score   += l1_multiplier
				 * (if1->get_tl(p[0].xy() + t, tc) - if2->get_tl(p[1].xy() + t, tc)).normsq();
				 
		}

		/*
		 * Include third-camera contributions in the score.
		 */

		if (tc_multiplier != 0)
		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
			if (n == f1 || n == f2)
				continue;

			lod_image *ifn = al->get(n);
			pt _ptn = ifn->get_t(0);
			point pn = _ptn.wp_unscaled(centroid);

			if (!_ptn.unscaled_in_bounds(pn))
				continue;

			if (pn[2] >= 0)
				continue;

			ale_pos ttc = get_trilinear_coordinate(min, max, _ptn);

			divisor += tc_multiplier;
			score   += tc_multiplier
				 * (if1->get_tl(p[0].xy(), tc) - ifn->get_tl(pn.xy(), ttc)).normsq();
		}

		c->add_candidate(p[0], tc, score / divisor);
	}

	/*
	 * Check for cells that are completely clipped.
	 */
	static int completely_clipped(point min, point max) {
		return (min[2] > front_clip
		     || max[2] < rear_clip);
	}

	/*
	 * Update extremum variables for cell points mapped to a particular view.
	 */
	static void update_extrema(point min, point max, pt _pt, int *extreme_dim, ale_pos *extreme_ratio) {

		point local_min, local_max;

		local_min = _pt.wp_unscaled(min);
		local_max = _pt.wp_unscaled(min);

		point cell[2] = {min, max};

		int near_vertex = 0;

		/*
		 * Determine the view-local extrema in all dimensions, and
		 * determine the vertex of closest z coordinate.
		 */

		for (int r = 1; r < 8; r++) {
			point local = _pt.wp_unscaled(point(cell[r>>2][0], cell[(r>>1)%2][1], cell[r%2][2]));
			
			for (int d = 0; d < 3; d++) {
				if (local[d] < local_min[d])
					local_min[d] = local[d];
				if (local[d] > local_max[d])
					local_max[d] = local[d];
			}

			if (local[2] == local_max[2])
				near_vertex = r;
		}

		ale_pos diameter = (local_max.xy() - local_min.xy()).norm();

		/*
		 * Update extrema as necessary for each dimension.
		 */

		for (int d = 0; d < 3; d++) {

			int r = near_vertex;

			int p1[3] = {r>>2, (r>>1)%2, r%2};
			int p2[3] = {r>>2, (r>>1)%2, r%2};

			p2[d] = 1 - p2[d];

			ale_pos local_distance = (_pt.wp_unscaled(point(cell[p1[0]][0], cell[p1[1]][1], cell[p1[2]][2])).xy()
					        - _pt.wp_unscaled(point(cell[p2[0]][0], cell[p2[1]][1], cell[p2[2]][2])).xy()).norm();

			if (local_distance / diameter > *extreme_ratio) {
				*extreme_ratio = local_distance / diameter;
				*extreme_dim = d;
			}
		}
	}

	/*
	 * Get the next split dimension.
	 */
	static int get_next_split(int f1, int f2, point min, point max, const std::vector<pt> &pt_outputs) {
		for (int d = 0; d < 3; d++)
			if (isinf(min[d]) || isinf(max[d]))
				return space::traverse::get_next_split(min, max);

		int extreme_dim = 0;
		ale_pos extreme_ratio = 0;

		update_extrema(min, max, al->get(f1)->get_t(0), &extreme_dim, &extreme_ratio);
		update_extrema(min, max, al->get(f2)->get_t(0), &extreme_dim, &extreme_ratio);

		for (unsigned int n = 0; n < pt_outputs.size(); n++) {
			update_extrema(min, max, pt_outputs[n], &extreme_dim, &extreme_ratio);
		}

		return extreme_dim;
	}

	/*
	 * Find candidates for subspace creation.
	 */
	static void find_candidates(unsigned int f1, unsigned int f2, candidates *c, point min, point max,
			const std::vector<pt> &pt_outputs) {

		if (completely_clipped(min, max))
			return;

		if (!input_might_be_visible(f1, min, max)
		 || !input_might_be_visible(f2, min, max))
			return;

		if (output_clip && !output_might_be_visible(pt_outputs, min, max))
			return;

		if (exceeds_resolution_lower_bounds(f1, f2, min, max, pt_outputs)) {
			if (!(rand() % 100000))
				fprintf(stderr, "([%f %f %f], [%f %f %f]) at %d\n", 
						min[0], min[1], min[2],
						max[0], max[1], max[2],
						__LINE__);

			try_nearest_candidate(f1, f2, c, min, max);
			return;
		}

		point new_cells[2][2];

		if (!space::traverse::get_next_cells(get_next_split(f1, f2, min, max, pt_outputs), min, max, new_cells))
			return;

		find_candidates(f1, f2, c, new_cells[0][0], new_cells[0][1], pt_outputs);
		find_candidates(f1, f2, c, new_cells[1][0], new_cells[1][1], pt_outputs);
	}

	/*
	 * Initialize space and identify regions of interest for the adaptive
	 * subspace model.
	 */
	static void make_space(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		fprintf(stderr, "Subdividing 3D space");

		std::vector<pt> pt_outputs = make_pt_list(d_out, v_out, d3_depth_pt, d3_output_pt);

		/*
		 * Initialize root space.
		 */

		space::init_root();

		/*
		 * Subdivide space to resolve intensity matches between pairs
		 * of frames.
		 */

		for (unsigned int f1 = 0; f1 < d2::image_rw::count(); f1++) {

			if (d3_depth_pt->size() == 0
			 && d3_output_pt->size() == 0
			 && d_out[f1] == NULL
			 && v_out[f1] == NULL)
				continue;

			if (tc_multiplier == 0)
				al->open(f1);

			for (unsigned int f2 = 0; f2 < d2::image_rw::count(); f2++) {

				if (f1 == f2)
					continue;

				if (tc_multiplier == 0)
					al->open(f2);

				candidates *c = new candidates(f1);

				find_candidates(f1, f2, c, point::neginf(), point::posinf(), pt_outputs);

				c->generate_subspaces();

				if (tc_multiplier == 0)
					al->close(f2);
			}

			if (tc_multiplier == 0)
				al->close(f1);
		}

		fprintf(stderr, ".\n");
	}


	/*
	 * Update spatial information structures.
	 *
	 * XXX: the name of this function is horribly misleading.  There isn't
	 * even a 'search depth' any longer, since there is no longer any
	 * bounded DFS occurring.
	 */
	static void reduce_cost_to_search_depth(d2::exposure *exp_out, int inc_bit) {

		/*
		 * Subspace model
		 */

		for (unsigned int i = 0; i < ou_iterations; i++)
			spatial_info_update();

		fprintf(stderr, "Final spatial map size: %u\n", spatial_info_map.size());
	}

#if 0
	/*
	 * Describe a scene to a renderer
	 */
	static void describe(render *r) {
	}
#endif
};

#endif
