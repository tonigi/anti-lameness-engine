// Copyright 2003, 2004, 2005, 2006 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                                <dhilvert@ugcs.caltech.edu>

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
	 * Filtering data
	 */

	static int use_filter;
	static const char *d3chain_type;

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
	 * Median calculation radii.
	 */
	static double depth_median_radius;
	static double diff_median_radius;

	/*
	 * Flag for subspace traversal.
	 */
	static int subspace_traverse;

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

			if (i == im->height() - 1
			 || j == im->width() - 1) {
				return 1;
			}

			assert(i <= im->height() - 1);
			assert(j <= im->width() - 1);

			/*
			 * Check for positive values.
			 */

			if (im->pix(i, j)[0] > 0) {
				if (st && st->node_value == im->pix(i, j)[0])
					im->pix(i, j)[0] += weight * (1 - im->pix(i, j)[0]);
				return 1;
			}

			/*
			 * Handle the case where there are no higher levels of detail.
			 */

			if (tc == tc_low) {
				if (im->pix(i, j)[0] != 0) {
					fprintf(stderr, "failing assertion: im[%d]->pix(%d, %d)[0] == %g\n", tc, i, j, 
							im->pix(i, j)[0]);
				}
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

					iim->pix(i * 2 + ii, j * 2 + jj)[0] = weight;
				}

			im->pix(i, j)[0] = -1;

			return 1;
		}

		/*
		 * Add weight.
		 */
		void add_weight(int tc, unsigned int i, unsigned int j, ale_real weight, subtree *st) {

			assert (weight >= 0);

			d2::image *im = weights[tc - tc_low];
			assert(im);

//			fprintf(stderr, "[aw tc=%d i=%d j=%d imax=%d jmax=%d]\n",
//					tc, i, j, im->height(), im->width());
			
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

			assert (weight >= 0);

			p *= pow(2, -tc);

			unsigned int i = (unsigned int) floor(p[0]);
			unsigned int j = (unsigned int) floor(p[1]);

			add_weight(tc, i, j, weight, st);
		}

		void add_weight(const space::traverse &t, ale_real weight, subtree *st) {

			assert (weight >= 0);

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

				initialized = 1;
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

//			fprintf(stderr, "[gw tc=%d i=%u j=%u tclow=%d tchigh=%d]\n", 
//					tc, i, j, tc_low, tc_high);

			if (tc < tc_low || !initialized)
				return 0;

			if (tc > tc_high) {
				return (get_weight(tc - 1, i * 2 + 0, j * 2 + 0)
				      + get_weight(tc - 1, i * 2 + 1, j * 2 + 0)
				      + get_weight(tc - 1, i * 2 + 1, j * 2 + 1)
				      + get_weight(tc - 1, i * 2 + 0, j * 2 + 1)) / 4;
			}

			assert(weights.size() > (unsigned int) (tc - tc_low));

			d2::image *im = weights[tc - tc_low];
			assert(im);

			if (i == im->height())
				return 1;
			if (j == im->width())
				return 1;

			assert(i < im->height());
			assert(j < im->width());

			if (im->pix(i, j)[0] == -1) {
				return (get_weight(tc - 1, i * 2 + 0, j * 2 + 0)
				      + get_weight(tc - 1, i * 2 + 1, j * 2 + 0)
				      + get_weight(tc - 1, i * 2 + 1, j * 2 + 1)
				      + get_weight(tc - 1, i * 2 + 0, j * 2 + 1)) / 4;
			}

			if (im->pix(i, j)[0] == 0) {
				if (tc == tc_high)
					return 0;
				if (weights[tc - tc_low + 1]->pix(i / 2, j / 2)[0] == -1)
					return 0;
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

			return get_weight((int) tc, p.xy());
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

		if (tc < input_decimation_lower - 1.5)
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

		unsigned int count() {
			return d2::image_rw::count();
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

//			fprintf(stderr, "[gs iw=%f %f %f d=%f]\n", 
//					iw[0], iw[1], iw[2], diagonal);

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
				levels[l - input_decimation_lower].resize((unsigned int) (floor(height / pow(2, l))
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

			int sheight = (int) floor(height / pow(2, tc));
			int swidth  = (int) floor(width / pow(2, tc));

			assert(i < sheight);
			assert(j < swidth);

			for (unsigned int k = 0; k < pairwise_ambiguity; k++) {
				std::pair<ale_pos, ale_real> *pk =
					&(levels[tc - input_decimation_lower][i * swidth * pairwise_ambiguity + j * pairwise_ambiguity + k]);

				if (pk->first != 0 && score >= pk->second)
					continue;

				if (i == 1 && j == 1 && tc == 4)
					fprintf(stderr, "[ac p2=%f score=%f]\n", p[2], score);

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
				unsigned int sheight = (unsigned int) floor(height / pow(2, l));
				unsigned int swidth  = (unsigned int) floor(width  / pow(2, l));

				for (unsigned int i = 0; i < sheight; i++)
				for (unsigned int j = 0; j < swidth; j++)
				for (unsigned int k = 0; k < pairwise_ambiguity; k++) {
					std::pair<ale_pos, ale_real> *pk =
						&(levels[l - input_decimation_lower]
							[i * swidth * pairwise_ambiguity + j * pairwise_ambiguity + k]);

					if (pk->first == 0) {
						fprintf(stderr, "o");
						continue;
					} else {
						fprintf(stderr, "|");
					}

					ale_pos si = i * pow(2, l) + ((l > 0) ? pow(2, l - 1) : 0);
					ale_pos sj = j * pow(2, l) + ((l > 0) ? pow(2, l - 1) : 0);

//					fprintf(stderr, "[gss l=%d i=%d j=%d d=%g]\n", l, i, j, pk->first);

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
			assert (finite(occupancy));
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

#if !defined(HASH_MAP_GNU) && !defined(HASH_MAP_STD)
	typedef std::map<struct space::node *, spatial_info> spatial_info_map_t;
#elif defined(HASH_MAP_GNU)
        struct node_hash
        {
                size_t operator()(struct space::node *n) const
                {
                        return __gnu_cxx::hash<long>()((long) n);
                }
        };
	typedef __gnu_cxx::hash_map<struct space::node *, spatial_info, node_hash > spatial_info_map_t;
#elif defined(HASH_MAP_STD)
	typedef std::hash_map<struct space::node *, spatial_info> spatial_info_map_t;
#endif

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

	static void dmr(double dmr_parameter) {
		depth_median_radius = dmr_parameter;
	}

	static void fmr(double fmr_parameter) {
		diff_median_radius = fmr_parameter;
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

		fprintf(stderr, "front_clip=%f rear_clip=%f\n", front_clip, rear_clip);

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
			 * Determine the probability of encounter.
			 */

			d2::pixel encounter = d2::pixel(1, 1, 1) * (1 - weights->get_weight(st));

			/*
			 * Update weights
			 */

			weights->add_weight(st, occupancy, tree);

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
	 * Support function for view() and depth().  This function
	 * always performs exclusion.
	 */

	static const void view_recurse(int type, d2::image *im, d2::image *weights, space::iterate si, pt _pt, 
			int prune = 0, d2::point pl = d2::point(0, 0), d2::point ph = d2::point(0, 0)) {
		while (!si.done()) {
			space::traverse st = si.get();

			/*
			 * Remove excluded regions.
			 */

			if (excluded(st)) {
				si.cleave();
				continue;
			}

			/*
			 * Prune.
			 */

			if (prune && !_pt.check_inclusion_scaled(st, pl, ph)) {
				si.cleave();
				continue;
			}

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

			if (prune) {
				if (min[0] > ph[0]
				 || min[1] > ph[1]
				 || max[0] < pl[0]
				 || max[1] < pl[1]) {
					si.next();
					continue;
				}

				if (min[0] < pl[0])
					min[0] = pl[0];
				if (min[1] < pl[1])
					min[1] = pl[1];
				if (max[0] > ph[0])
					max[0] = ph[0];
				if (max[1] > ph[1])
					max[1] = ph[1];

				min[0] -= pl[0];
				min[1] -= pl[1];
				max[0] -= pl[0];
				max[1] -= pl[1];
			}

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

				view_recurse(type, im, weights, cleaved_space, _pt, prune, pl, ph);

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

				d2::pixel encounter = (d2::pixel(1, 1, 1) 
						     - weights->get_pixel(i, j)) 
					            * occupancy;

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

					ale_pos contribution_value = sn.get_pocc_density() /* + sn.get_socc_density() */;
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

					assert (finite(encounter[0]));
					assert (finite(contribution_value));

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
	static const d2::image *depth(pt _pt, int n = -1, int prune = 0, 
			d2::point pl = d2::point(0, 0), d2::point ph = d2::point(0, 0)) {
		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		d2::image *im1, *im2, *im3, *weights;;

		if (prune) {

			im1 = new d2::image_ale_real((int) floor(ph[0] - pl[0]) + 1,
					(int) floor(ph[1] - pl[1]) + 1, 3);

			im2 = new d2::image_ale_real((int) floor(ph[0] - pl[0]) + 1,
					(int) floor(ph[1] - pl[1]) + 1, 3);

			im3 = new d2::image_ale_real((int) floor(ph[0] - pl[0]) + 1,
					(int) floor(ph[1] - pl[1]) + 1, 3);

			weights = new d2::image_ale_real((int) floor(ph[0] - pl[0]) + 1,
					(int) floor(ph[1] - pl[1]) + 1, 3);

		} else {

			im1 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
			     	       (int) floor(_pt.scaled_width()), 3);

			im2 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
			     	       (int) floor(_pt.scaled_width()), 3);

			im3 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
			     	       (int) floor(_pt.scaled_width()), 3);

			weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
							(int) floor(_pt.scaled_width()), 3);
		}

		/*
		 * Iterate through subspaces.
		 */

		space::iterate si(_pt.origin());

		view_recurse(6, im1, weights, si, _pt, prune, pl, ph);

		delete weights;

		if (prune) {
			weights = new d2::image_ale_real((int) floor(ph[0] - pl[0]) + 1,
					(int) floor(ph[1] - pl[1]) + 1, 3);
		} else {
			weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
							(int) floor(_pt.scaled_width()), 3);
		}

#if 1
		view_recurse(7, im2, weights, si, _pt, prune, pl, ph);
#else
		view_recurse(4, im2, weights, si, _pt, prune, pl, ph);
		return im2;
#endif

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
			 * Handle encounter thresholds
			 */

			point w = _pt.pw_scaled(point(i + pl[0], j + pl[1], depth));

			if (weights->pix(i, j)[0] < encounter_threshold) {
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

	/*
	 * Unfiltered function
	 */
	static const d2::image *view_nofilter_focus(pt _pt, int n) {

		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		const d2::image *depths = depth(_pt, n);

		d2::image *im = new d2::image_ale_real((int) floor(_pt.scaled_height()),
					       (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {

			focus::result _focus = focus::get(depths, i, j);

			if (!finite(_focus.focal_distance))
				continue;

			/*
			 * Iterate over views for this focus region.
			 */

			for (unsigned int v = 0; v < _focus.sample_count; v++) {

				/*
				 * Determine the (x, y) offset for this view.
				 */

				ale_pos ofx = _focus.aperture;
				ale_pos ofy = _focus.aperture;

				while (ofx * ofx + ofy * ofy > _focus.aperture * _focus.aperture / 4) {
					ofx = (rand() * _focus.aperture) / RAND_MAX - _focus.aperture / 2;
					ofy = (rand() * _focus.aperture) / RAND_MAX - _focus.aperture / 2;
				}

				// fprintf(stderr, "[vnff ofx=%f ofy=%f]\n", ofx, ofy);

				/*
				 * Generate a new view from the given offset.
				 */

				point new_view = _pt.cw(point(ofx, ofy, 0));
				pt _pt_new = _pt;
				for (int d = 0; d < 3; d++)
					_pt_new.e().set_translation(d, -new_view[d]);

				// fprintf(stderr, "[vnff nv=[%f %f %f]]\n", new_view[0], new_view[1], new_view[2]);

				/*
				 * Map the focused point to the new view.
				 */

				point p = _pt_new.wp_scaled(_pt.pw_scaled(point(i, j, _focus.focal_distance)));

				/*
				 * Determine weight and color for the given point.
				 */

				d2::image *im_point = new d2::image_ale_real(1, 1, 3);
				d2::image *wt_point = new d2::image_ale_real(1, 1, 3);

				view_recurse(0, im_point, wt_point, space::iterate(_pt_new.origin()),
					_pt_new, 1, p.xy(), p.xy());

				im->pix(i, j) += im_point->pix(0, 0);
				weights->pix(i, j) += wt_point->pix(0, 0);

				delete im_point;
				delete wt_point;
			}

			if (weights->pix(i, j).min_norm() < encounter_threshold
			 || (d3px_count > 0 && isnan(depths->pix(i, j)[0]))) {
				im->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
				weights->pix(i, j) = d2::pixel::zero();
			} else if (normalize_weights)
				im->pix(i, j) /= weights->pix(i, j);
		}

		delete weights;

		delete depths;

		return im;
	}

	/*
	 * Unfiltered function
	 */
	static const d2::image *view_nofilter(pt _pt, int n) {

		if (!focus::is_trivial())
			return view_nofilter_focus(_pt, n);

		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		const d2::image *depths = depth(_pt, n);

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

		delete weights;

		delete depths;

		return im;
	}

	/*
	 * This function always performs exclusion.
	 */

	static space::node *most_visible_pointwise(d2::pixel *weight, space::iterate si, pt _pt, d2::point p) {

		space::node *result = NULL;

		while (!si.done()) {
			space::traverse st = si.get();

			/*
			 * Prune certain regions known to be uninteresting.
			 */

			if (excluded(st) || !_pt.check_inclusion_scaled(st, p)) {
				si.cleave();
				continue;
			}

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

			ale_real occupancy = sn.get_occupancy();

			/*
			 * Preserve current weight in order to check for
			 * modification by higher-resolution subspaces.
			 */

			d2::pixel old_weight = *weight;

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_node()->positive
			 || st.get_node()->negative) {

				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that afterward.
				 */

				space::iterate cleaved_space = si.cleave();

				cleaved_space.next();

				space::node *r = most_visible_pointwise(weight, cleaved_space, _pt, p);

				if (old_weight[1] != (*weight)[1])
					result = r;

			} else {
				si.next();
			}
				

			/*
			 * Check for higher-resolution updates.
			 */

			if (old_weight != *weight)
				continue;

			/*
			 * Determine the probability of encounter.
			 */

			ale_pos encounter = (1 - (*weight)[0]) * occupancy;

			/*
			 * (*weight)[0] stores the cumulative weight; (*weight)[1] stores the maximum.
			 */

			if (encounter > (*weight)[1]) {
				result = st.get_node();
				(*weight)[1] = encounter;
			}

			(*weight)[0] += encounter;
		}

		return result;
	}

	/*
	 * This function performs exclusion iff SCALED is true.
	 */
	static void  most_visible_generic(std::vector<space::node *> &results, d2::image *weights, 
			space::iterate si, pt _pt, int scaled) {

		assert (results.size() == weights->height() * weights->width());

		while (!si.done()) {
			space::traverse st = si.get();

			if (scaled && excluded(st)) {
				si.cleave();
				continue;
			}

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

				most_visible_generic(results, weights, cleaved_space, _pt, scaled);

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

				ale_pos encounter = (1 - weights->get_pixel(i, j)[0]) * occupancy;

				/*
				 * weights[0] stores the cumulative weight; weights[1] stores the maximum.
				 */

				if (encounter > weights->get_pixel(i, j)[1]
				 || results[i * weights->width() + j] == NULL) {
					results[i * weights->width() + j] = st.get_node();
					weights->chan(i, j, 1) = encounter;
				}

				weights->chan(i, j, 0) += encounter;
			}
		}
	}

	static std::vector<space::node *> most_visible_scaled(pt _pt) {
		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()), 
				(int) floor(_pt.scaled_width()), 3);
		std::vector<space::node *> results;

		results.resize(weights->height() * weights->width(), 0);
	
		most_visible_generic(results, weights, space::iterate(_pt.origin()), _pt, 1);
		
		return results;
	}

	static std::vector<space::node *> most_visible_unscaled(pt _pt) {
		d2::image *weights = new d2::image_ale_real((int) floor(_pt.unscaled_height()), 
				(int) floor(_pt.unscaled_width()), 3);
		std::vector<space::node *> results;
		
		results.resize(weights->height() * weights->width(), 0);

		most_visible_generic(results, weights, space::iterate(_pt.origin()), _pt, 0);
		
		return results;
	}

	static const int visibility_search(const std::vector<space::node *> &fmv, space::node *mv) {

		if (mv == NULL)
			return 0;

		if (std::binary_search(fmv.begin(), fmv.end(), mv))
			return 1;

		return (visibility_search(fmv, mv->positive)
		     || visibility_search(fmv, mv->negative));

	}

	/*
	 * Filtered function.
	 */
	static const d2::image *view_filter_focus(pt _pt, int n) {

		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		/*
		 * Generate a new 2D renderer for filtering.
		 */

		d2::render::reset();
		d2::render *renderer = d2::render_parse::get(d3chain_type);

		/*
		 * Get depth image for focus region determination.
		 */

		const d2::image *depths = depth(_pt, n);

		unsigned int height = (unsigned int) floor(_pt.scaled_height());
		unsigned int width = (unsigned int) floor(_pt.scaled_width());

		renderer->init_point_renderer(height, width, 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Prepare view data.
		 */

		if (tc_multiplier == 0)
			al->open_all();

		pt *_ptf = new pt[al->count()];
		std::vector<space::node *> *fmv = new std::vector<space::node *>[al->count()];

		for (unsigned int f = 0; f < al->count(); f++) {
			_ptf[f] = al->get(f)->get_t(0);
			fmv[f] = most_visible_unscaled(_ptf[f]);
			std::sort(fmv[f].begin(), fmv[f].end());
		}

		if (tc_multiplier == 0)
			al->close_all();

		/*
		 * Open all files for rendering.
		 */

		d2::image_rw::open_all();

		/*
		 * Iterate over output points.
		 */
		
		for (unsigned int i = 0; i < height; i++)
		for (unsigned int j = 0; j < width; j++) {

			focus::result _focus = focus::get(depths, i, j);

			if (!finite(_focus.focal_distance))
				continue;

			/*
			 * Iterate over views for this focus region.
			 */

			for (unsigned int v = 0; v < _focus.sample_count; v++) {

				/*
				 * Determine the (x, y) offset for this view.
				 */

				ale_pos ofx = _focus.aperture;
				ale_pos ofy = _focus.aperture;

				while (ofx * ofx + ofy * ofy > _focus.aperture * _focus.aperture / 4) {
					ofx = (rand() * _focus.aperture) / RAND_MAX - _focus.aperture / 2;
					ofy = (rand() * _focus.aperture) / RAND_MAX - _focus.aperture / 2;
				}

				/*
				 * Generate a new view from the given offset.
				 */

				point new_view = _pt.cw(point(ofx, ofy, 0));
				pt _pt_new = _pt;
				for (int d = 0; d < 3; d++)
					_pt_new.e().set_translation(d, -new_view[d]);

				/*
				 * Map the focused point to the new view.
				 */

				point p = _pt_new.wp_scaled(_pt.pw_scaled(point(i, j, _focus.focal_distance)));

//					fprintf(stderr, "[vff i=%d j=%d p=[%f %f %f]]\n",
//							i, j, p[0], p[1], p[2]);

				/*
				 * Determine the most-visible subspace.
				 */

				d2::pixel weight(0, 0, 0);
				space::node *mv = most_visible_pointwise(&weight, space::iterate(_pt.origin()), 
					_pt_new, p.xy());

				if (mv == NULL)
					continue;

				/*
				 * Generate a local depth image of required radius.
				 */

				ale_pos radius = 1;

				if (diff_median_radius + 1 > radius)
					radius = diff_median_radius + 1;
				if (depth_median_radius > radius)
					radius = depth_median_radius;

				d2::point pl = p.xy() - d2::point(radius, radius);
				d2::point ph = p.xy() + d2::point(radius, radius);
				const d2::image *local_depth = depth(_pt_new, -1, 1, pl, ph);

				/*
				 * Find depth and diff at this point, check for
				 * undefined values, and generate projections
				 * of the image corners on the estimated normal
				 * surface.
				 */

				d2::image *median_diffs = local_depth->fcdiff_median((int) floor(diff_median_radius));
				d2::image *median_depths = local_depth->medians((int) floor(depth_median_radius));

				d2::pixel depth = median_depths->pix((int) radius, (int) radius);
				d2::pixel diff = median_diffs->pix((int) radius, (int) radius);

				delete median_diffs;
				delete median_depths;
				delete local_depth;

				if (!depth.finite() || !diff.finite())
					continue;

				point local_points[3] = { 
					point(p[0],     p[1],     depth[0]),
					point(p[0] + 1, p[1],     depth[0] + diff[0]),
					point(p[0],     p[1] + 1, depth[0] + diff[1])
				};

				/*
				 * Iterate over files.
				 */

				for (unsigned int f = 0; f < d2::image_rw::count(); f++) {

					if (!visibility_search(fmv[f], mv))
						continue;

					/*
					 * Determine transformation at (i, j).  First
					 * determine transformation from the output to
					 * the input, then invert this, as we need the
					 * inverse transformation for filtering.
					 */

					d2::point remote_points[3] = {
						_ptf[f].wp_unscaled(_pt_new.pw_scaled(point(local_points[0]))).xy(),
						_ptf[f].wp_unscaled(_pt_new.pw_scaled(point(local_points[1]))).xy(),
						_ptf[f].wp_unscaled(_pt_new.pw_scaled(point(local_points[2]))).xy()
					};

					/*
					 * Forward matrix for the linear component of the 
					 * transformation.
					 */

					d2::point forward_matrix[2] = {
						remote_points[1] - remote_points[0],
						remote_points[2] - remote_points[0]
					};

					/*
					 * Inverse matrix for the linear component of
					 * the transformation.  Calculate using the
					 * determinant D.
					 */

					ale_pos D = forward_matrix[0][0] * forward_matrix[1][1]
						  - forward_matrix[0][1] * forward_matrix[1][0];

					if (D == 0)
						continue;

					d2::point inverse_matrix[2] = {
						d2::point( forward_matrix[1][1] / D, -forward_matrix[1][0] / D),
						d2::point(-forward_matrix[0][1] / D,  forward_matrix[0][0] / D)
					};

					/*
					 * Determine the projective transformation parameters for the
					 * inverse transformation.
					 */
					
					const d2::image *imf = d2::image_rw::get_open(f);

					d2::transformation inv_t = d2::transformation::gpt_identity(imf, 1);

					d2::point local_bounds[4];

					for (int n = 0; n < 4; n++) {
						d2::point remote_bound = d2::point((n == 1 || n == 2) ? imf->height() : 0,
										   (n == 2 || n == 3) ? imf->width()  : 0)
								       - remote_points[0];

						local_bounds[n] = d2::point(i, j)
								+ d2::point(remote_bound[0] * inverse_matrix[0][0]
									  + remote_bound[1] * inverse_matrix[1][0],
									    remote_bound[0] * inverse_matrix[0][1]
									  + remote_bound[1] * inverse_matrix[1][1]);
					}

					inv_t.gpt_set(local_bounds);

					/*
					 * Perform render step for the given frame,
					 * transformation, and point.
					 */

					renderer->point_render(i, j, f, inv_t);
				}
			}
		}

		/*
		 * Close all files and finish rendering.
		 */

		d2::image_rw::close_all();

		renderer->finish_point_rendering();

		return renderer->get_image();
	}

	static const d2::image *view_filter(pt _pt, int n) {

		if (!focus::is_trivial())
			return view_filter_focus(_pt, n);

		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		/*
		 * Generate a new 2D renderer for filtering.
		 */

		d2::render::reset();
		d2::render *renderer = d2::render_parse::get(d3chain_type);

		/*
		 * Get depth image in order to estimate normals (and hence
		 * transformations).
		 */

		const d2::image *depths = depth(_pt, n);

		d2::image *median_diffs = depths->fcdiff_median((int) floor(diff_median_radius));
		d2::image *median_depths = depths->medians((int) floor(depth_median_radius));

		unsigned int height = (unsigned int) floor(_pt.scaled_height());
		unsigned int width = (unsigned int) floor(_pt.scaled_width());

		renderer->init_point_renderer(height, width, 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		std::vector<space::node *> mv = most_visible_scaled(_pt);

		for (unsigned int f = 0; f < d2::image_rw::count(); f++) {

			if (tc_multiplier == 0)
				al->open(f);

			pt _ptf = al->get(f)->get_t(0);

			std::vector<space::node *> fmv = most_visible_unscaled(_ptf);
			std::sort(fmv.begin(), fmv.end());

			for (unsigned int i = 0; i < height; i++)
			for (unsigned int j = 0; j < width; j++) {

				/*
				 * Check visibility.
				 */

				int n = i * width + j;

				if (!visibility_search(fmv, mv[n]))
					continue;

				/*
				 * Find depth and diff at this point, check for
				 * undefined values, and generate projections
				 * of the image corners on the estimated normal
				 * surface.
				 */

				d2::pixel depth = median_depths->pix(i, j);
				d2::pixel diff = median_diffs->pix(i, j);
				// d2::pixel diff = d2::pixel(0, 0, 0);

				if (!depth.finite() || !diff.finite())
					continue;

				point local_points[3] = { 
					point(i,     j,     depth[0]),
				        point(i + 1, j,     depth[0] + diff[0]),
				        point(i    , j + 1, depth[0] + diff[1])
				};

				/*
				 * Determine transformation at (i, j).  First
				 * determine transformation from the output to
				 * the input, then invert this, as we need the
				 * inverse transformation for filtering.
				 */

				d2::point remote_points[3] = {
					_ptf.wp_unscaled(_pt.pw_scaled(point(local_points[0]))).xy(),
					_ptf.wp_unscaled(_pt.pw_scaled(point(local_points[1]))).xy(),
					_ptf.wp_unscaled(_pt.pw_scaled(point(local_points[2]))).xy()
				};

				/*
				 * Forward matrix for the linear component of the 
				 * transformation.
				 */

				d2::point forward_matrix[2] = {
					remote_points[1] - remote_points[0],
					remote_points[2] - remote_points[0]
				};

				/*
				 * Inverse matrix for the linear component of
				 * the transformation.  Calculate using the
				 * determinant D.
				 */

				ale_pos D = forward_matrix[0][0] * forward_matrix[1][1]
					  - forward_matrix[0][1] * forward_matrix[1][0];

				if (D == 0)
					continue;

				d2::point inverse_matrix[2] = {
					d2::point( forward_matrix[1][1] / D, -forward_matrix[1][0] / D),
					d2::point(-forward_matrix[0][1] / D,  forward_matrix[0][0] / D)
				};

				/*
				 * Determine the projective transformation parameters for the
				 * inverse transformation.
				 */
				
				const d2::image *imf = d2::image_rw::open(f);

				d2::transformation inv_t = d2::transformation::gpt_identity(imf, 1);

				d2::point local_bounds[4];

				for (int n = 0; n < 4; n++) {
					d2::point remote_bound = d2::point((n == 1 || n == 2) ? imf->height() : 0,
							                   (n == 2 || n == 3) ? imf->width()  : 0)
							       - remote_points[0];

					local_bounds[n] = local_points[0].xy()
						        + d2::point(remote_bound[0] * inverse_matrix[0][0]
							          + remote_bound[1] * inverse_matrix[1][0],
								    remote_bound[0] * inverse_matrix[0][1]
								  + remote_bound[1] * inverse_matrix[1][1]);
				}

				inv_t.gpt_set(local_bounds);

				d2::image_rw::close(f);

				/*
				 * Perform render step for the given frame,
				 * transformation, and point.
				 */

				d2::image_rw::open(f);
				renderer->point_render(i, j, f, inv_t);
				d2::image_rw::close(f);
			}

			if (tc_multiplier == 0) 
				al->close(f);
		}

		renderer->finish_point_rendering();

		return renderer->get_image();
	}

	/*
	 * Generic function.
	 */
	static const d2::image *view(pt _pt, int n = -1) {

		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (use_filter) {
			return view_filter(_pt, n);
		} else {
			return view_nofilter(_pt, n);
		}
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

	static void nofilter() {
		use_filter = 0;
	}

	static void filter() {
		use_filter = 1;
	}

	static void set_filter_type(const char *type) {
		d3chain_type = type;
	}

	static void set_subspace_traverse() {
		subspace_traverse = 1;
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

	/*
	 * This function returns true if a space is completely excluded.
	 */
	static int excluded(const space::traverse &st) {
		for (int n = 0; n < d3px_count; n++) {
			double *region = d3px_parameters + (6 * n);
			if (st.get_min()[0] >= region[0]
			 && st.get_max()[0] <= region[1]
			 && st.get_min()[1] >= region[2]
			 && st.get_max()[1] <= region[3]
			 && st.get_min()[2] >= region[4]
			 && st.get_max()[2] <= region[5])
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
				if (isnan(local[d])) 
					return local[d];
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

		// fprintf(stderr, "[tnc n=%f %f %f x=%f %f %f]\n", min[0], min[1], min[2], max[0], max[1], max[2]);

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

			// fprintf(stderr, ":");

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
			const std::vector<pt> &pt_outputs, int depth = 0) {

		int print = 0;

		if (min[0] < 20.0001 && max[0] > 20.0001
		 && min[1] < 20.0001 && max[1] > 20.0001
		 && min[2] < 0.0001 && max[2] > 0.0001)
			print = 1;

		if (print) {
			for (int i = depth; i > 0; i--) {
				fprintf(stderr, "+");
			}
			fprintf(stderr, "[fc n=%f %f %f x=%f %f %f]\n",
					min[0], min[1], min[2], max[0], max[1], max[2]);
		}

		if (completely_clipped(min, max)) {
			if (print)
				fprintf(stderr, "c");
			return;
		}

		if (!input_might_be_visible(f1, min, max)
		 || !input_might_be_visible(f2, min, max)) {
			if (print)
				fprintf(stderr, "v");
			return;
		}

		if (output_clip && !output_might_be_visible(pt_outputs, min, max)) {
			if (print)
				fprintf(stderr, "o");
			return;
		}

		if (exceeds_resolution_lower_bounds(f1, f2, min, max, pt_outputs)) {
			if (!(rand() % 100000))
				fprintf(stderr, "([%f %f %f], [%f %f %f]) at %d\n", 
						min[0], min[1], min[2],
						max[0], max[1], max[2],
						__LINE__);

			if (print)
				fprintf(stderr, "t");

			try_nearest_candidate(f1, f2, c, min, max);
			return;
		}

		point new_cells[2][2];

		if (!space::traverse::get_next_cells(get_next_split(f1, f2, min, max, pt_outputs), min, max, new_cells)) {
			if (print)
				fprintf(stderr, "n");
			return;
		}

		if (print) {
			fprintf(stderr, "nc[0][0]=%f %f %f nc[0][1]=%f %f %f nc[1][0]=%f %f %f nc[1][1]=%f %f %f\n",
					new_cells[0][0][0],
					new_cells[0][0][1],
					new_cells[0][0][2],
					new_cells[0][1][0],
					new_cells[0][1][1],
					new_cells[0][1][2],
					new_cells[1][0][0],
					new_cells[1][0][1],
					new_cells[1][0][2],
					new_cells[1][1][0],
					new_cells[1][1][1],
					new_cells[1][1][2]);
		}

		find_candidates(f1, f2, c, new_cells[0][0], new_cells[0][1], pt_outputs, depth + 1);
		find_candidates(f1, f2, c, new_cells[1][0], new_cells[1][1], pt_outputs, depth + 1);
	}

	/*
	 * Generate a map from scores to 3D points for various depths at point (i, j) in f1, at 
	 * lowest resolution.
	 */
	static score_map p2f_score_map(unsigned int f1, unsigned int f2, unsigned int i, unsigned int j) {

		score_map result;

		pt _pt1 = al->get(f1)->get_t(primary_decimation_upper);
		pt _pt2 = al->get(f2)->get_t(primary_decimation_upper);

		const d2::image *if1 = al->get(f1)->get_image(primary_decimation_upper);
		const d2::image *if2 = al->get(f2)->get_image(primary_decimation_upper);

		/*
		 * Get the pixel color in the primary frame
		 */

		// d2::pixel color_primary = if1->get_pixel(i, j);

		/*
		 * Map two depths to the secondary frame.
		 */

		point p1 = _pt2.wp_unscaled(_pt1.pw_unscaled(point(i, j, 1000)));
		point p2 = _pt2.wp_unscaled(_pt1.pw_unscaled(point(i, j, -1000)));

//		fprintf(stderr, "%d->%d (%d, %d) point pair: (%d, %d, %d -> %f, %f), (%d, %d, %d -> %f, %f)\n",
//				f1, f2, i, j, i, j, 1000, p1[0], p1[1], i, j, -1000, p2[0], p2[1]);
//		_pt1.debug_output();
//		_pt2.debug_output();


		/*
		 * For cases where the mapped points define a
		 * line and where points on the line fall
		 * within the defined area of the frame,
		 * determine the starting point for inspection.
		 * In other cases, continue to the next pixel.
		 */

		ale_pos diff_i = p2[0] - p1[0];
		ale_pos diff_j = p2[1] - p1[1];
		ale_pos slope = diff_j / diff_i;

		if (isnan(slope)) {
			assert(0);
			fprintf(stderr, "%d->%d (%d, %d) has undefined slope\n",
					f1, f2, i, j);
			return result;
		}

		/*
		 * Make absurdly large/small slopes either infinity, negative infinity, or zero.
		 */

		if (fabs(slope) > if2->width() * 100) {
			double zero = 0;
			double one  = 1;
			double inf  = one / zero;
			slope = inf;
		} else if (slope < 1 / (double) if2->height() / 100
			&& slope > -1/ (double) if2->height() / 100) {
			slope = 0;
		}

		// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

		ale_pos top_intersect = p1[1] - p1[0] * slope;
		ale_pos lef_intersect = p1[0] - p1[1] / slope;
		ale_pos rig_intersect = p1[0] - (p1[1] - if2->width() + 2) / slope;
		ale_pos sp_i, sp_j;

		// fprintf(stderr, "slope == %f\n", slope);


		if (slope == 0) {
			// fprintf(stderr, "case 0\n");
			sp_i = lef_intersect;
			sp_j = 0;
		} else if (finite(slope) && top_intersect >= 0 && top_intersect < if2->width() - 1) {
			// fprintf(stderr, "case 1\n");
			sp_i = 0;
			sp_j = top_intersect;
		} else if (slope > 0 && lef_intersect >= 0 && lef_intersect <= if2->height() - 1) {
			// fprintf(stderr, "case 2\n");
			sp_i = lef_intersect;
			sp_j = 0;
		} else if (slope < 0 && rig_intersect >= 0 && rig_intersect <= if2->height() - 1) {
			// fprintf(stderr, "case 3\n");
			sp_i = rig_intersect;
			sp_j = if2->width() - 2;
		} else {
			// fprintf(stderr, "case 4\n");
			// fprintf(stderr, "%d->%d (%d, %d) does not intersect the defined area\n",
			//		f1, f2, i, j);
			return result;
		}


		// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

		/*
		 * Determine increment values for examining
		 * point, ensuring that incr_i is always
		 * positive.
		 */

		ale_pos incr_i, incr_j;

		if (fabs(diff_i) > fabs(diff_j)) {
			incr_i = 1;
			incr_j = slope;
		} else if (slope > 0) {
			incr_i = 1 / slope;
			incr_j = 1;
		} else {
			incr_i = -1 / slope;
			incr_j = -1;
		}
		
		// fprintf(stderr, "%d->%d (%d, %d) increments are (%f, %f)\n",
		//		f1, f2, i, j, incr_i, incr_j);

		/*
		 * Examine regions near the projected line.
		 */

		for (ale_pos ii = sp_i, jj = sp_j; 
			ii <= if2->height() - 1 && jj <= if2->width() - 1 && ii >= 0 && jj >= 0; 
			ii += incr_i, jj += incr_j) {

			// fprintf(stderr, "%d->%d (%d, %d) checking (%f, %f)\n", 
			//		f1, f2, i, j, ii, jj);

#if 0
			/*
			 * Check for higher, lower, and nearby points.
			 *
			 *	Red   = 2^0
			 *	Green = 2^1
			 *	Blue  = 2^2
			 */

			int higher = 0, lower = 0, nearby = 0;

			for (int iii = 0; iii < 2; iii++)
			for (int jjj = 0; jjj < 2; jjj++) {
				d2::pixel p = if2->get_pixel((int) floor(ii) + iii, (int) floor(jj) + jjj);

				for (int k = 0; k < 3; k++) {
					int bitmask = (int) pow(2, k);

					if (p[k] > color_primary[k])
						higher |= bitmask;
					if (p[k] < color_primary[k])
						lower  |= bitmask;
					if (fabs(p[k] - color_primary[k]) < nearness)
						nearby |= bitmask;
				}
			}

			/*
			 * If this is not a region of interest,
			 * then continue.
			 */


			fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			// if (((higher & lower) | nearby) != 0x7)
			//	continue;
#endif
			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			// fprintf(stderr, "%d->%d (%d, %d) accepted (%f, %f)\n", 
			//		f1, f2, i, j, ii, jj);

			/*
			 * Create an orthonormal basis to
			 * determine line intersection.
			 */

			point bp0 = _pt1.pw_unscaled(point(i, j, 0));
			point bp1 = _pt1.pw_unscaled(point(i, j, 10));
			point bp2 = _pt2.pw_unscaled(point(ii, jj, 0));

			point foo = _pt1.wp_unscaled(bp0);
			// fprintf(stderr, "(%d, %d, 0) transformed to world and back is: (%f, %f, %f)\n",
			//		i, j, foo[0], foo[1], foo[2]);

			foo = _pt1.wp_unscaled(bp1);
			// fprintf(stderr, "(%d, %d, 10)  transformed to world and back is: (%f, %f, %f)\n",
			//		i, j, foo[0], foo[1], foo[2]);

			point b0  = (bp1 - bp0).normalize();
			point b1n = bp2 - bp0;
			point b1  = (b1n - b1n.dproduct(b0) * b0).normalize();
			point b2  = point(0, 0, 0).xproduct(b0, b1).normalize(); // Should already have norm=1
			

			foo = _pt1.wp_unscaled(bp0 + 30 * b0);

			/*
			 * Select a fourth point to define a second line.
			 */

			point p3  = _pt2.pw_unscaled(point(ii, jj, 10));

			/*
			 * Representation in the new basis.
			 */

			d2::point nbp0 = d2::point(bp0.dproduct(b0), bp0.dproduct(b1));
			// d2::point nbp1 = d2::point(bp1.dproduct(b0), bp1.dproduct(b1));
			d2::point nbp2 = d2::point(bp2.dproduct(b0), bp2.dproduct(b1));
			d2::point np3  = d2::point( p3.dproduct(b0),  p3.dproduct(b1));

			/*
			 * Determine intersection of line
			 * (nbp0, nbp1), which is parallel to
			 * b0, with line (nbp2, np3).
			 */

			/*
			 * XXX: a stronger check would be
			 * better here, e.g., involving the
			 * ratio (np3[0] - nbp2[0]) / (np3[1] -
			 * nbp2[1]).  Also, acceptance of these
			 * cases is probably better than
			 * rejection.
			 */


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			if (np3[1] - nbp2[1] == 0)
				continue;


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			d2::point intersection = d2::point(nbp2[0] 
					+ (nbp0[1] - nbp2[1]) * (np3[0] - nbp2[0]) / (np3[1] - nbp2[1]),
					nbp0[1]);

			ale_pos b2_offset = b2.dproduct(bp0);

			/*
			 * Map the intersection back to the world
			 * basis.
			 */

			point iw = intersection[0] * b0 + intersection[1] * b1 + b2_offset * b2;

			/*
			 * Reject intersection points behind a
			 * camera.
			 */

			point icp = _pt1.wc(iw);
			point ics = _pt2.wc(iw);


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			if (icp[2] >= 0 || ics[2] >= 0)
				continue;


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			/*
			 * Reject clipping plane violations.
			 */


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			if (iw[2] > front_clip
			 || iw[2] < rear_clip)
				continue;


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			/*
			 * Score the point.
			 */

			point ip = _pt1.wp_unscaled(iw);

			point is = _pt2.wp_unscaled(iw);

			analytic _a = { iw, ip, is };

			/*
			 * Calculate score from color match.  Assume for now
			 * that the transformation can be approximated locally
			 * with a translation.
			 */

			ale_pos score = 0;
			ale_pos divisor = 0;
			ale_pos l1_multiplier = 0.125;

			if (if1->in_bounds(ip.xy())
			 && if2->in_bounds(is.xy())) {
				divisor += 1 - l1_multiplier;
				score += (1 - l1_multiplier)
				       * (if1->get_bl(ip.xy()) - if2->get_bl(is.xy())).normsq();
			}


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			for (int iii = -1; iii <= 1; iii++)
			for (int jjj = -1; jjj <= 1; jjj++) {
				d2::point t(iii, jjj);

				if (!if1->in_bounds(ip.xy() + t)
				 || !if2->in_bounds(is.xy() + t))
					continue;

				divisor += l1_multiplier;
				score   += l1_multiplier
					 * (if1->get_bl(ip.xy() + t) - if2->get_bl(is.xy() + t)).normsq();
					 
			}

			/*
			 * Include third-camera contributions in the score.
			 */

			if (tc_multiplier != 0)
			for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
				if (f == f1 || f == f2)
					continue;

				const d2::image *if3 = al->get(f)->get_image(primary_decimation_upper);
				pt _pt3 = al->get(f)->get_t(primary_decimation_upper);

				point p = _pt3.wp_unscaled(iw);

				if (!if3->in_bounds(p.xy())
				 || !if1->in_bounds(ip.xy()))
					continue;

				divisor += tc_multiplier;
				score   += tc_multiplier
					 * (if1->get_bl(ip.xy()) - if3->get_bl(p.xy())).normsq();
			}

			


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			/*
			 * Reject points with undefined score.
			 */


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			if (!finite(score / divisor))
				continue;


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

#if 0
			/*
			 * XXX: reject points not on the z=-27.882252 plane.
			 */


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			if (_a.ip[2] > -27 || _a.ip[2] < -28)
				continue;
#endif


			// fprintf(stderr, "score map (%u, %u) line %u\n", i, j, __LINE__);

			/*
			 * Add the point to the score map.
			 */

//			d2::pixel c_ip = if1->in_bounds(ip.xy()) ? if1->get_bl(ip.xy())
//								 : d2::pixel();
//			d2::pixel c_is = if2->in_bounds(is.xy()) ? if2->get_bl(is.xy())
//								 : d2::pixel();

//			fprintf(stderr, "Candidate subspace: f1=%u f2=%u i=%u j=%u ii=%f jj=%f"
//					"cp=[%f %f %f] cs=[%f %f %f]\n",
//					f1, f2, i, j, ii, jj, c_ip[0], c_ip[1], c_ip[2],
//					c_is[0], c_is[1], c_is[2]);

			result.insert(score_map_element(score / divisor, _a));
		}

//		fprintf(stderr, "Iterating through the score map:\n");
//
//		for (score_map::iterator smi = result.begin(); smi != result.end(); smi++) {
//			fprintf(stderr, "%f ", smi->first);
//		}
//
//		fprintf(stderr, "\n");

		return result;
	}


	/*
	 * Attempt to refine space around a point, to high and low resolutions
	 * resulting in two resolutions in total.
	 */

	static space::traverse refine_space(point iw, ale_pos target_dim, int use_filler) {

		space::traverse st = space::traverse::root();

		if (!st.includes(iw)) {
			assert(0);
			return st;
		}

		int lr_done = !use_filler;

		/*
		 * Loop until all resolutions of interest have been generated.
		 */
		
		for(;;) {

			point diff = st.get_max() - st.get_min();

			point p[2] = { st.get_min(), st.get_max() };

			ale_pos dim_max = 0;

			for (int d = 0; d < 3; d++) {
				ale_pos d_value = fabs(p[0][d] - p[1][d]);
				if (d_value > dim_max)
					dim_max = d_value;
			}

			/*
			 * Generate any new desired spatial registers.
			 */

			for (int f = 0; f < 2; f++) {

				/*
				 * Low resolution
				 */

				if (dim_max < 2 * target_dim
				 && lr_done == 0) {
					spatial_info_map[st.get_node()];
					lr_done = 1;
				}

				/*
				 * High resolution.
				 */

				if (dim_max < target_dim) {
					spatial_info_map[st.get_node()];
					return st;
				}
			}

			/*
			 * Check precision before analyzing space further.
			 */

			if (st.precision_wall()) {
				fprintf(stderr, "\n\n*** Error: reached subspace precision wall ***\n\n");
				assert(0);
				return st;
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

	/*
	 * Calculate target dimension
	 */

	static ale_pos calc_target_dim(point iw, pt _pt, const char *d_out[], const char *v_out[], 
			std::map<const char *, pt> *d3_depth_pt, 
			std::map<const char *, pt> *d3_output_pt) {

		ale_pos result = _pt.distance_1d(iw, primary_decimation_upper);

		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
			if (d_out[n] && align::projective(n).distance_1d(iw, 0) < result)
				result = align::projective(n).distance_1d(iw, 0);
			if (v_out[n] && align::projective(n).distance_1d(iw, 0) < result)
				result = align::projective(n).distance_1d(iw, 0);
		}

		for (std::map<const char *, pt>::iterator i = d3_output_pt->begin(); i != d3_output_pt->end(); i++) {
			if (i->second.distance_1d(iw, 0) < result)
				result = i->second.distance_1d(iw, 0);
		}

		for (std::map<const char *, pt>::iterator i = d3_depth_pt->begin(); i != d3_depth_pt->end(); i++) {
			if (i->second.distance_1d(iw, 0) < result)
				result = i->second.distance_1d(iw, 0);
		}

		assert (result > 0);

		return result;
	}

	/*
	 * Calculate level of detail for a given viewpoint.
	 */

	static int calc_lod(ale_pos depth1, pt _pt, ale_pos target_dim) {
		return (int) round(_pt.trilinear_coordinate(depth1, target_dim * sqrt(2)));
	}

	/*
	 * Calculate depth range for a given pair of viewpoints.
	 */

	static ale_pos calc_depth_range(point iw, pt _pt1, pt _pt2) {

		point ip = _pt1.wp_unscaled(iw);

		ale_pos reference_change = fabs(ip[2] / 1000);

		point iw1 = _pt1.pw_scaled(ip + point(0, 0, reference_change));
		point iw2 = _pt1.pw_scaled(ip - point(0, 0, reference_change));

		point is = _pt2.wc(iw);
		point is1 = _pt2.wc(iw1);
		point is2 = _pt2.wc(iw2);

		assert(is[2] < 0);

		ale_pos d1 = (is1.xy() - is.xy()).norm();
		ale_pos d2 = (is2.xy() - is.xy()).norm();

		if (is1[2] < 0 && is2[2] < 0) {

			if (d1 > d2)
				return reference_change / d1;
			else
				return reference_change / d2;
		}

		if (is1[2] < 0)
			return reference_change / d1;

		if (is2[2] < 0)
			return reference_change / d2;

		return 0;
	}

	/*
	 * Calculate a refined point for a given set of parameters.
	 */

	static point get_refined_point(pt _pt1, pt _pt2, int i, int j, 
			int f1, int f2, int lod1, int lod2, ale_pos depth,
			ale_pos depth_range) {

		d2::pixel comparison_color = al->get(f1)->get_image(lod1)->get_pixel(i, j);

		ale_pos best = -1;
		ale_pos best_depth = depth;

		for (ale_pos d = depth - depth_range; d < depth + depth_range; d += depth_range / 10) {

			if (!(d < 0))
				continue;
			
			point iw = _pt1.pw_unscaled(point(i, j, d));
			point is = _pt2.wp_unscaled(iw);

			if (!(is[2] < 0))
				continue;

			if (!al->get(f2)->get_image(lod2)->in_bounds(is.xy()))
				continue;

			ale_pos error = (comparison_color - al->get(f2)->get_image(lod2)->get_bl(is.xy())).norm();

			if (error < best || best == -1) {
				best = error;
				best_depth = d;
			}
		}

		return _pt1.pw_unscaled(point(i, j, best_depth));
	}

	/*
	 * Analyze space in a manner dependent on the score map.
	 */

	static void analyze_space_from_map(const char *d_out[], const char *v_out[],
				       std::map<const char *, pt> *d3_depth_pt,
				       std::map<const char *, pt> *d3_output_pt,
				       unsigned int f1, unsigned int f2, 
				       unsigned int i, unsigned int j, score_map _sm, int use_filler) {

		int accumulated_ambiguity = 0;
		int max_acc_amb = pairwise_ambiguity;

		pt _pt1 = al->get(f1)->get_t(0);
		pt _pt2 = al->get(f2)->get_t(0);

		if (_pt1.scale_2d() != 1)
			use_filler = 1;

		for(score_map::iterator smi = _sm.begin(); smi != _sm.end(); smi++) {

			point iw = smi->second.iw;
			point ip = smi->second.ip;
			// point is = smi->second.is;

			if (accumulated_ambiguity++ >= max_acc_amb)
				break;

			total_ambiguity++;

			ale_pos depth1 = _pt1.wc(iw)[2];
			ale_pos depth2 = _pt2.wc(iw)[2];

			ale_pos target_dim = calc_target_dim(iw, _pt1, d_out, v_out, d3_depth_pt, d3_output_pt);

			assert(target_dim > 0);

			int lod1 = calc_lod(depth1, _pt1, target_dim);
			int lod2 = calc_lod(depth2, _pt2, target_dim);

			while (lod1 < input_decimation_lower
			    || lod2 < input_decimation_lower) {
				target_dim *= 2;
				lod1 = calc_lod(depth1, _pt1, target_dim);
				lod2 = calc_lod(depth2, _pt2, target_dim);
			}


			if (lod1 >= (int) al->get(f1)->count()
			 || lod2 >= (int) al->get(f2)->count())
				continue;

			int multiplier = (unsigned int) floor(pow(2, primary_decimation_upper - lod1));

			ale_pos depth_range = calc_depth_range(iw, _pt1, _pt2);

			pt _pt1_lod = al->get(f1)->get_t(lod1);
			pt _pt2_lod = al->get(f2)->get_t(lod2);

			int im = i * multiplier;
			int jm = j * multiplier;

			for (int ii = 0; ii < multiplier; ii += 1)
			for (int jj = 0; jj < multiplier; jj += 1) {

				point refined_point = get_refined_point(_pt1_lod, _pt2_lod, im + ii, jm + jj, 
						f1, f2, lod1, lod2, depth1, depth_range);
				
				/*
				 * Re-evaluate target dimension.
				 */

				ale_pos target_dim_ = 
					calc_target_dim(refined_point, _pt1, d_out, v_out, d3_depth_pt, d3_output_pt);

				ale_pos depth1_ = _pt1.wc(refined_point)[2];
				ale_pos depth2_ = _pt2.wc(refined_point)[2];

				int lod1_ = calc_lod(depth1_, _pt1, target_dim_);
				int lod2_ = calc_lod(depth2_, _pt2, target_dim_);

				while (lod1_ < input_decimation_lower
				    || lod2_ < input_decimation_lower) {
					target_dim_ *= 2;
					lod1_ = calc_lod(depth1_, _pt1, target_dim_);
					lod2_ = calc_lod(depth2_, _pt2, target_dim_);
				}

				/*
				 * Attempt to refine space around the intersection point.
				 */

				space::traverse st = 
					refine_space(refined_point, target_dim_, use_filler || _pt1.scale_2d() != 1);

				ale_pos tc1 = al->get(f1)->get_t(0).trilinear_coordinate(st);
				ale_pos tc2 = al->get(f2)->get_t(0).trilinear_coordinate(st);


				assert(resolution_ok(al->get(f1)->get_t(0), tc1));
				assert(resolution_ok(al->get(f2)->get_t(0), tc2));
			}

		}
	}


	/*
	 * Initialize space and identify regions of interest for the adaptive
	 * subspace model.
	 */
	static void make_space(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		/*
		 * Variable indicating whether low-resolution filler space
		 * is desired to avoid aliased gaps in surfaces.
		 */

		int use_filler = d3_depth_pt->size() != 0
				   || d3_output_pt->size() != 0
				   || output_decimation_preferred > 0
				   || input_decimation_lower > 0
				   || !focus::is_trivial();

		fprintf(stderr, "[T=%lu]\n", (long unsigned) time(NULL));

		fprintf(stderr, "Subdividing 3D space");

		std::vector<pt> pt_outputs = make_pt_list(d_out, v_out, d3_depth_pt, d3_output_pt);

		/*
		 * Initialize root space.
		 */

		space::init_root();

		/*
		 * Special handling for experimental option 'subspace_traverse'.
		 */

		if (subspace_traverse) {
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

			fprintf(stderr, "Final spatial map size: %u\n", spatial_info_map.size());

			fprintf(stderr, ".\n");
			fprintf(stderr, "[T=%lu]\n", (long unsigned) time(NULL));

			return;
		}

		/*
		 * Subdivide space to resolve intensity matches between pairs
		 * of frames.
		 */

		for (unsigned int f1 = 0; f1 < d2::image_rw::count(); f1++)
		for (unsigned int f2 = 0; f2 < d2::image_rw::count(); f2++) {
			if (f1 == f2)
				continue;

			if (!d_out[f1] && !v_out[f1] && !d3_depth_pt->size()
			 && !d3_output_pt->size() && strcmp(pairwise_comparisons, "all"))
				continue;

			if (tc_multiplier == 0) {
				al->open(f1);
				al->open(f2);
			}

			/*
			 * Iterate over all points in the primary frame.
			 */

			for (unsigned int i = 0; i < al->get(f1)->get_image(primary_decimation_upper)->height(); i++)
			for (unsigned int j = 0; j < al->get(f1)->get_image(primary_decimation_upper)->width();  j++) {

				total_pixels++;

				/*
				 * Generate a map from scores to 3D points for
				 * various depths in f1.
				 */

				score_map _sm = p2f_score_map(f1, f2, i, j);

				/*
				 * Analyze space in a manner dependent on the score map.
				 */

				analyze_space_from_map(d_out, v_out, d3_depth_pt, d3_output_pt, 
						f1, f2, i, j, _sm, use_filler);

			}

			/*
			 * This ordering may encourage image f1 to be cached.
			 */

			if (tc_multiplier == 0) {
				al->close(f2);
				al->close(f1);
			}
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

		fprintf(stderr, "[T=%lu]\n", (long unsigned) time(NULL));
		/*
		 * Subspace model
		 */

		for (unsigned int i = 0; i < ou_iterations; i++)
			spatial_info_update();

		fprintf(stderr, "Final spatial map size: %u\n", spatial_info_map.size());
		fprintf(stderr, "[T=%lu]\n", (long unsigned) time(NULL));
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
