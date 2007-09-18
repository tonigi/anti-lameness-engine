// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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

/*
 * Help messages
 */

#define BETWEEN_SECTIONS "\n"
#define HEADER_SPACE ""

class help {
private:
	const char *package;
	const char *version;
	const char *invocation;
	FILE *help_stream;

	/*
	 * Stars
	 * 
	 * This function produces a line of stars for banner output.
	 */
	void stars(unsigned int n) {
		for (unsigned int i = 0; i < n; i++) {
			fprintf(help_stream, "*");
		}
		fprintf(help_stream, "\n");
	}


	/*
	 * Banner
	 * 
	 * This function aids in locating the start of help output.
	 */
	void banner(const char *name) {
		const char *package_banner = " Help Text, version ";
		const char *section_banner = "Section: ";

		int plen = strlen(package) + strlen(package_banner) + strlen(version);
		int slen = strlen(section_banner) + strlen(name);
		int len = (plen > slen) ? plen : slen;

		fprintf(help_stream, BETWEEN_SECTIONS);
		stars(len);
		fprintf(help_stream, "%s%s%s\n", package, package_banner, version);
		fprintf(help_stream, "%s%s\n", section_banner, name);
		stars(len);
	}

public:
	help(const char *package, const char *invocation, const char *version) {
		this->invocation = invocation;
		this->package = package;
		this->version = version;
		this->help_stream = stdout;
	}

	/*
	 * Describe how to use this program
	 */
	void usage() {
		banner("Usage");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Usage: %s [<options>] <original-frame> [<supplemental-frame> ...] <output-file>\n"
			"   or: %s [<help option> ...]\n"
			"   or: %s --version\n"
			BETWEEN_SECTIONS
			"Help options:\n"
			HEADER_SPACE
			"--hu              Usage (this message).\n"
			"--hf              Image files.\n"
			"--he              Exclusion regions.\n"
			"--ha              Alignment (not exposure-related).\n"
			"--hr              Rendering (not exposure-related).\n"
			"--hx              Exposure.\n"
			"--ht              Transformation data files.\n"
			"--hc              Control points.\n"
			"--hl              Filtering (PSFs, rendering chains).\n"
			"--hd              Devices.\n"
			"--hi              User Interfaces.\n"
			"--hp              Process details.\n"
			"--hs              Argument scope (Experimental).\n"
			"--hv              Video stream processing (Experimental).\n"
			"--h3              3D Modeling (Experimental).\n"
			"--hz              Undocumented options.\n"
			"--hA              Concatenate all help pages.\n"
			"\n",
			invocation, invocation, invocation);
	}

	void defaults() {
		banner("Defaults");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Default settings:\n"
			HEADER_SPACE
			"--q* options are no longer recognized.\n"
			"\n"
		       );
	}
	
	void file() {
		banner("File");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Cache options:\n"
			HEADER_SPACE
			"--cache <size>    Cache up to <size> megabytes of image data.   (256 is default)\n"
			BETWEEN_SECTIONS
			"Bit depth options:\n"
			HEADER_SPACE
			"--8bpc            Write 8 bit per channel output\n"
			"--16bpc           Write 16 bit per channel output [default]\n"
			BETWEEN_SECTIONS
			"Output format options:\n"
			HEADER_SPACE
#ifdef USE_MAGICK
			"--auto            Determine output file type automatically [default]\n"
			"--raw             Write raw PPM output\n"
			"--plain           Write plain PPM output\n"
#else
			"--raw             Write raw PPM output [default]\n"
			"--plain           Write plain PPM output\n"
#endif
			BETWEEN_SECTIONS
			"Incremental output:\n"
			HEADER_SPACE
			"--inc             Produce rough incremental output.\n"
			"--no-inc          Don't produce any incremental output.  [default]\n"

			BETWEEN_SECTIONS
			"Undefined values:\n"
			HEADER_SPACE
			"--def-nn <R>      Use nearest-neighbor defined values within\n"
			"                  radius <R>, zero outside.  Default radius is 0.\n"
			"\n"
		       );
	}
	void alignment() {
		banner("Alignment");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Alignment channel options:\n"
			HEADER_SPACE
			"--align-all       Align images using all color channels\n"
			"--align-green     Align images using the green channel\n"
			"--align-sum       Align images using a sum of channels [default]\n"
			BETWEEN_SECTIONS
			"Transformation options:\n"
			HEADER_SPACE
			"--translation     Only adjust the position of images\n"
			"--euclidean       Adjust the position and orientation of images [default]\n"
			"--projective      Use projective transformations.  Best quality, but slow.\n"
			BETWEEN_SECTIONS
			"Alignment following:\n"
			HEADER_SPACE
			"--follow          Frames align closely with their predecessor.  [default]\n"
			"--identity        Frames align closely with the original frame.\n"
			BETWEEN_SECTIONS
			"Alignment failure:\n"
			HEADER_SPACE
			"--fail-optimal    Frames beneath threshold are aligned optimally.  [default]\n"
			"--fail-default    Frames beneath threshold keep their default alignment.\n"
			BETWEEN_SECTIONS
			"Transformation file operations:\n"
			HEADER_SPACE
			"--trans-load <x>  Load initial transformation settings from file <x>\n"
			"--trans-save <x>  Save final transformation data in file <x>\n"
			BETWEEN_SECTIONS
			"Tunable parameters:\n"
			HEADER_SPACE
			"--metric=x        Set the alignment error metric exponent.       (2 is default)\n"
			"--threshold=x     Min. match threshold; a perfect match is 100.  (0 is default)\n"
			"--perturb-upper=x Perturbation upper bound pixels/arclength    (14%% is default)\n"
			"                     ('x%%' uses a fraction of the smallest image dimension.)\n"
			"--perturb-lower=x Perturbation lower bound pixels/arclength   (.125 is default)\n"
			"                     ('x%%' uses a fraction of the smallest image dimension.)\n"
			"--rot-upper=x     Rotation-specific upper bound in degrees    (32.0 is default)\n"
			"--bda-mult=x      Barrel distortion adjustment multiplier      (2.0 is default)\n"
			"--bda-rate=x      Barrel distortion rate of change maximum     (8.0 is default)\n"   
			"--lod-max=x       LOD scale factor is max(1, (2^floor(x))/perturb)  (1 is def.)\n"
			BETWEEN_SECTIONS
			"Certainty-weighted alignment:\n"
			HEADER_SPACE
			"--cw              Weight alignment error by certainty.\n"
			"--no-cw           Don't weight alignment error by certainty. [default]\n"
			BETWEEN_SECTIONS
			"Alignment weight maps:\n"
			HEADER_SPACE
			"--wm <f> <x> <y>  Use weight map image <f> at offset (<x>, <y>)\n"
			BETWEEN_SECTIONS
			"Frequency-weighted alignment:\n"
			HEADER_SPACE
			"--fl <h> <v> <a>  High-pass filters: horizontal <h>, vertical <v>, average <a>.\n"
			"                     Values should fall between 0 (pass all) and 1 (pass none).\n"
#ifndef USE_FFTW
			"\n"
			"                     NOTE: since this build of ALE does not link with FFTW,\n"
			"                           this option is not supported.  To use this option,\n"
			"                           first re-build with FFTW support.\n"
			"\n"
#endif
			"--flshow <o>      Write high-pass filtered data to file <o>.\n"
			BETWEEN_SECTIONS
			"Algorithmic alignment weighting:\n"
			HEADER_SPACE
			"--wmx <e> <r> <d> Write reference <r>, definition <d>, execute `<e> <f> <d>`,\n"
			"                  read weights <r> back.\n"
#ifndef USE_UNIX
			"\n"
			"                     NOTE: since this build was not configured with\n"
			"                           support for --wmx, this option is not supported.\n"
			"                           To use this option, re-build with support for --wmx.\n"
#endif
			BETWEEN_SECTIONS
			"Perturbation Type [experimental]:\n"
			HEADER_SPACE
			"--perturb-output  Apply perturbations in output image coordinates. [default]\n"
			"--perturb-source  Apply perturbations in source image coordinates.\n"
			BETWEEN_SECTIONS
			"Global searching:\n"
			HEADER_SPACE
			"--gs <type>       Set global search to <type>, one of:\n"
			"                     local     Local alignment only\n"
			"                     inner     Alignment reference image inner region\n"
			"                     outer     Alignment reference image outer region\n"
			"                     all       Union of inner and outer\n"
			"                     central   inner if below threshold or better; else, outer.\n"
			"                     defaults  'all' if default, 'local' otherwise.  [default]\n"
			"                     points    Align by control points.  Ignores gs-mo.\n"
			"--gs-mo <x>       Set <x> pixel min. overlap for global search.   (67%% default)\n"
			"                     ('x%%' uses a fraction of the total pixel count.)\n"
			BETWEEN_SECTIONS
			"Precise final alignment display:\n"
			HEADER_SPACE
			"--precise         Display precise final alignment.\n"
			"--no-precise      Do not display precise final alignment.  [default]\n"
			BETWEEN_SECTIONS
			"Multi-alignment [algorithms are not fully implemented]:\n"
			HEADER_SPACE
			"--ma-card <x>     Maximum number of multi-alignment elements.  (1 is default)\n"
			"--ma-cont <x>     Minimum contiguous area of each MA region.  (100 is default)\n"
			"\n"
		       );
	}

	void rendering() {
		banner("Rendering");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Mimicking ALE 0.6.0 merging and drizzling (see --hl for details):\n"
			HEADER_SPACE
			"   --dchain triangle:2     approximates merging.\n"
			"   --dchain fine:box:1     approximates drizzling.\n"
			BETWEEN_SECTIONS
			"Image extents:\n"
			HEADER_SPACE
			"--extend          Increase image extents to accommodate all pixel data.\n"
			"--no-extend       Don't increase extents; crop to original frame. [default]\n"
			BETWEEN_SECTIONS
			"Tunable parameters:\n"
			HEADER_SPACE
			"--scale=x         Scale images by the factor x, where x > 0.     (1 is default)\n"
			"--threshold=x     Min. match threshold; a perfect match is 100.  (0 is default)\n"
			BETWEEN_SECTIONS
			"Irani-Peleg iterative solver:\n"
			HEADER_SPACE
			"--ips <i>         Run <i> iterations.                            (1 is default)\n"
			"--ip-mean         Use the mean correction [default]\n"
			"--ip-median       Use the median correction\n"
			"--ip-wl <x>       Use weight limit <x>\n"
			"--ip-nowl         Use no weight limit [default]\n"
#if 0
			BETWEEN_SECTIONS
			"Unsharp Mask (was 'High-frequency Enhancement'):\n"
			HEADER_SPACE
			"--usm <m>         Apply an unsharp mask with multiplier <m>.\n"
			"                     (See also --device, --nlpsf, and --lpsf.)\n"
#endif
			BETWEEN_SECTIONS
			"Bayer pattern:\n"
			HEADER_SPACE
			"--bayer <b>       Set the Bayer pattern to <b>, one of:\n"
			"                    (clockwise from top left pixel)\n"
			"                      rgbg      Red-green-blue-green\n"
			"                      gbgr      Green-blue-green-red\n"
			"                      grgb      Green-red-green-blue\n"
			"                      bgrg      Blue-green-red-green\n"
			"                      none      RGB-RGB-RGB-RGB\n"
			"                  Default is none or device-specific.\n"
			BETWEEN_SECTIONS
			"Color adjustment:\n"
			HEADER_SPACE
			"--exp-mult=c,r,b   Adjust all channels by <c>, red by <r>, and blue by <b>.\n"
			"\n"
		);
	}
	void filtering() {
		banner("Filtering");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Point-spread functions (used with --ips; see --hr):\n"
			HEADER_SPACE
			"--lpsf <p>        Set linear colorspace point-spread function to <p>\n"
			"--nlpsf <p>       Set non-linear colorspace point-spread function to <p>\n"
			"                     Available point-spread functions:\n"
			"                        box=<diameter>\n"
			"                        circle=<diameter>\n"
			"                        gauss=<diameter of one std. deviation>\n"
			"                        stdin\n"
			"                        stdin_vg\n"
			"                        <p>+<p> (summation)\n"
			"                        <p>^<p> (convolution)\n"
			"                        <n>*<p> (multiplication by a scalar <n>)\n"
			"                     Default lpsf is either 'box=1.0' or device-specific.\n"
			"                     Default nlpsf is either disabled or device-specific.\n"
			BETWEEN_SECTIONS
			"Rendering chains:\n"
			HEADER_SPACE
			"--wt <w>          Set weight threshold for defined pixels (default is 0.0001).\n"
			"--dchain <g>      Use chain <g> to render the default output.\n"
			"--ochain <g> <o>  Use chain <g> to render output file <o>.\n"
			"--achain <g>      Use chain <g> to render the alignment reference image.\n"
			"--afilter <s>     Use SSF <s> to interpolate points in alignment.\n"
			"--3d-chain <g>    Use chain <g> by default to render 3d output.\n"
			"                  Example chains:\n"
			"                     triangle:2                   ALE 0.6.0 merging (roughly)\n"
			"                     fine:box:1                   Drizzling (roughly)\n"
			"                     fine:sinc*lanc:8             High-frequency preservation\n"                     
			"                     last:nex:sinc*lanc:8         Useful for video stabilization\n"
			"                     fine:box:1,triangle:2        Multi-resolution rendering\n"
			"                     fine:sinc*lanc:8,sinc*lanc:8 Multi-resolution rendering\n"
			"                     auto:box:1,fine:box:1,box:1  Range-resolution rendering.\n"
			"                  More precisely, chains <g> are one of:\n"
			"                     <g1>,<g2>      Chain <g1> where defined; <g2> elsewhere.\n"
			"                     <i>            Use rendering invariant <i>.\n"
			"                  Rendering invariants <i> are:\n"
			"                     avg:<e>        Avg. (mean) value using SSFE <e>.\n"
			"                     avgf:<x>:<e>   Avg. first up to weight <x> using SSFE <e>.\n"
			"                     first:<e>      First defined value using SSFE <e>.\n"
			"                     last:<e>       Last defined value using SSFE <e>.\n"
			"                     max:<e>        Maximum defined value using SSFE <e>.\n"
			"                     min:<e>        Minimum defined value using SSFE <e>.\n"
			"                     median:<e>     Median value using SSFE <e>.\n"
			"                     <e>            Same as avg:<e>.\n"
			"                  Scaled sampling filters with exclusion (SSFE) <e> are:\n"
			"                     ex:<s>         Use SSF <s>; honor exclusion regions.\n"
			"                     nex:<s>        Use SSF <s>; don't honor exclusion regions.\n"
			"                     <s>            Same as ex:<s>\n"
			"                  Scaled sampling filters (SSF) <s> are:\n"
			"                     auto:<f>       filter <f> from chain suffix resolutions.\n"
			"                     fine:<f>       filter <f> at output image resolution.\n"
			"                     coarse:<f>     filter <f> at resolution MIN(in, out).\n"
			"                     <f>            Same as coarse:<f>.\n"
			"                  Sampling filters <f> are:\n"
			"                     sinc           Sinc filter. (very large diameter)\n"
			"                     lanc:<x>       Lanczos, diameter <x>.\n"
			"                     triangle:<x>   Triangle, diameter <x>.\n"
			"                     box:<x>        Box, diameter <x>.\n"
			"                     gauss:<x>      Gaussian, standard deviation <x>.\n"
			"                     zero           Zero function.\n"
			"                     <f>*<f>        Pointwise multiplication (windowing).\n"
			"                  Defaults:\n"
			"                     dchain         auto:triangle:2,fine:triangle:2,triangle:2\n"
			"                     3d-chain       fine:triangle:2,fine:gauss:0.75,triangle:2\n"
			"                     achain         triangle:2\n"
			"                     afilter        internal (approximates triangle:2)\n"
			"\n"
			);
	}
	void device() {
		banner("Device");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Device (may set PSF, Bayer pattern, exposure, and view angle):\n"
			HEADER_SPACE
			"--device <d>      Set the capture device to <d>.\n"
			"                     Available devices (* expect linear inputs):\n"
			"                        canon_300d *\n"
			"                        canon_300d+50mm_1.4 *\n"
			"                        canon_300d+50mm_1.4@1.4 *\n"
			"                        canon_300d+50mm_1.8 *\n"
			"                        canon_300d+85mm_1.8 *\n"
			"                        nikon_d50 *\n"
			"                        ov7620 *\n"
			"                        xvp610_320x240\n"
			"                        xvp610_640x480\n"
			"\n"
		       );
	}
	void exclusion() {
		banner("Exclusion");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Rendering-coordinate exclusion regions:\n"
			HEADER_SPACE
			"--ex <args>       Exclude a specified volume.\n"
			"                     <args> are space-delimited horizontal,\n"
			"                     vertical, and frame limits:\n"
			"                         <xmin> <xmax> <ymin> <ymax> <fmin> <fmax>\n"
			"                     using unscaled rendering spatial coordinates.\n"
			"--crop <args>     Exclude the spatial complement of an area over a\n"
			"                  specified set of frames.  <args> are:\n"
			"                         <xmin> <xmax> <ymin> <ymax> <fmin> <fmax>\n"
			"                     using unscaled rendering spatial coordinates.\n"
			BETWEEN_SECTIONS
			"Frame-coordinate exclusion regions:\n"
			HEADER_SPACE
			"--fex <args>      Exclude a specified volume.\n"
			"                     <args> are space-delimited horizontal,\n"
			"                     vertical, and frame limits:\n"
			"                         <xmin> <xmax> <ymin> <ymax> <fmin> <fmax>\n"
			"                     using unscaled frame spatial coordinates.\n"
			"--fcrop <args>    Exclude the spatial complement of an area over a\n"
			"                  specified set of frames.  <args> are:\n"
			"                         <xmin> <xmax> <ymin> <ymax> <fmin> <fmax>\n"
			"                     using unscaled frame spatial coordinates.\n"
			"\n"
		       );
	}
	void exposure() {
		banner("Exposure");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Certainty-weighted rendering:\n"
			HEADER_SPACE
			"--cx <x>          Render with certainty exponent <x>. (default is 1)\n"
			"--no-cx           Render with uniform certainty.\n"
			BETWEEN_SECTIONS
			"Exposure registration:\n"
			HEADER_SPACE
			"--exp-register    Register exposure between frames.  [default]\n"
			"--exp-noregister  Assume uniform exposure across all frames.\n"
			"--exp-meta-only   Use only meta-data for registering exposure.\n"
			BETWEEN_SECTIONS
			"Range extension:\n"
			HEADER_SPACE
			"--exp-extend      Extend to include all calculated values.\n"
			"--exp-noextend    Restrict to the original frame's range.  [default]\n"
			BETWEEN_SECTIONS
			"Exposure value meta-data:\n"
			HEADER_SPACE
			"--ev <x>          Set ISO 100 equivalent EV <x>.  (default is 0)\n"
			"--black <x>       Set black level <x> as a fraction of saturation.  (default 0)\n"
			"\n"
		       );
	}
	void tdf() {
		banner("Transformation data files");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Version 2 syntax overview:\n"
			HEADER_SPACE
			"{version string}\n"
			"{supplemental frame 1 transformation}\n"
			"{supplemental frame 2 transformation}\n"
			"...\n"
			BETWEEN_SECTIONS
			"Version 3 syntax overview:\n"
			HEADER_SPACE
			"{version string}\n"
			"{original frame transformation}\n"
			"{supplemental frame 1 transformation}\n"
			"{supplemental frame 2 transformation}\n"
			"...\n"
			BETWEEN_SECTIONS
			"Version string:\n"
			HEADER_SPACE
			"V <x>             Transformation data file version <x>.\n"
			BETWEEN_SECTIONS
			"Transformation overview:\n"
			HEADER_SPACE
			"{barrel/pincushion distortion command (optional; version 3 only)}\n"
			"{projective, euclidean, or default command}\n"
			BETWEEN_SECTIONS
			"Barrel/pincushion distortion (version 3 only):\n"
			HEADER_SPACE
			"B <n> <coeff2> <coeff3> ... <coeff(n+1)>\n"
			BETWEEN_SECTIONS
			"Projective command:\n"
			HEADER_SPACE
			"P <xmax> <ymax> <tlx> <tly> <blx> <bly> <brx> <bry> <trx> <try>\n"
			BETWEEN_SECTIONS
			"Euclidean command:\n"
			HEADER_SPACE
			"E <xmax> <ymax> <xoffset> <yoffset> <angle>\n"
			"\n"
			BETWEEN_SECTIONS
			"Comments:\n"
			HEADER_SPACE
			"# Lines beginning with '#' are comments.\n"
			BETWEEN_SECTIONS
			"Example:\n"
			HEADER_SPACE
			"# Version 3 transformation data file.\n"
			"V 3\n"
			"# Original frame, with barrel/pincushion distortion correction\n"
			"B 3 0.1 0 -0.1\n"
			"D\n"
			"# Supplemental frame 1: shift a 640x480 image right by 100 pixels\n"
			"E 640 480 100 0 0\n"
			"\n"
			);
	}
	void visp() {
		banner("Video stream processing");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Video stream processing [Experimental]:\n"
			HEADER_SPACE
			"--visp <args>     Process a video sequence.\n"
			"                     <args> are:\n"
			"                         <chain> <stabilization-type> <prefix> <suffix>\n"
			"                     <chain> is a rendering chain.  (see --hl)\n"
			"                     <stabilization-type> is one of:\n"
			"                         ma:<x>     Moving average over 2*<x> + 1 frames\n"
			"                         sf:<x>     Stabilize to single frame number <x>\n"
			"                         identity   Same as ma:0\n"
			"                     <prefix> is an output file prefix\n"
			"                     <suffix> is an output file suffix\n"
			"--visp-scale=<x>  Use scale <x> for VISP output.  (default is 1.0)\n"
			"--exshow          For single-invariant chains, show --ex regions dimmed.\n"
			"\n");
	}
	void interface() {
		banner("User Interface");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"User Interfaces:\n"
			HEADER_SPACE
			"--ui=<type>       Set user interface to <type>, one of:\n"
			"                     quiet\n"
			"                     stream\n"
			"                     tty [default]\n"
			"                     log\n"
#ifndef USE_IOCTL
			"\n"
			"                     NOTE: since ALE was compiled without terminal size check,\n"
			"                     --ui=tty will behave identically to --ui=stream.\n"
			"                     For additional output, recompile with terminal size check.\n"
#endif
			"\n");
	}
	void cp() {
		banner("Control Points");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Control point files:\n"
			HEADER_SPACE
			"--cpf-load=<f>    Load control point data from file <f>\n"
			"\n");
	}
	void d3() {
		banner("3D Modeling [Experimental]");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Rendering:\n"
			HEADER_SPACE
			"--3dv <n> <o>     Render, to file <o>, colors as viewed from frame <n>.\n"
			"--3dd <n> <o>     Render, to file <o>, depths as viewed from frame <n>.\n"
			"--3dvp <args> <o> Render, to file <o>, colors viewed with projective <args>.\n"
			"--3ddp <args> <o> Render, to file <o>, depths viewed with projective <args>.\n"
			"                    <args> are W H V x y z P Y R:\n"
			"                         W    image width.\n"
			"                         H    image height.\n"
			"                         V    camera view angle.\n"
			"                         x    translation x component.\n"
			"                         y    translation y component.\n"
			"                         z    translation z component.\n"
			"                         P    rotation around x-axis.\n"
			"                         Y    rotation around y-axis.\n"
			"                         R    rotation around z-axis.\n"
			"--occ-norm        Normalize output with accumulated occupancy.        [default]\n"
			"--occ-nonorm      Don't normalize output with accumulated occupancy.\n"
			"--et <x>          Set encounter threshold <x> for defined pixels.[default is 0]\n"
			"--3dpx <args>     Exclude a specified spatial volume following full-scene\n"
			"                  reconstruction.  <args> are:\n"
			"                         <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>\n"
			"--3d-filter       Use filtering for 3D color output.                  [default]\n"
			"--3d-nofilter     Don't use filtering for 3D color output.\n"
			"--3d-dmr <x>      Set radius for filtering median depth to <x>      [default 0]\n"
			"--3d-fmr <x>      Set radius for filtering median diff to <x>       [default 0]\n"
			"--focus <ft> <op> Create focus region with type <ft> and options <op>:\n"
			"                     Focus type:\n"
			"                        d <d>           focus at distance <d>\n"
			"                        p <x> <y>       focus at point (x, y)\n"
			"                     Space-separated options may include zero or more of:\n"
			"                        ci=<ci>  camera index               [default 0]\n"
			"                        fr=<fr>  focal range                [default 0]\n"
			"                        ht=<ht>  horizontal tilt gradient   [default 0]\n"
			"                        vt=<vt>  vertical tilt gradient     [default 0]\n"
			"                        ap=<ap>  aperture diameter          [default 3]\n"
			"                        sc=<sc>  sample count               [default 3]\n"
			"                        sx=<sx>  start x coordinate      [default -Inf]\n"
			"                        ex=<ex>  end x coordinate         [default Inf]\n"
			"                        sy=<sy>  start y coordinate      [default -Inf]\n"
			"                        ey=<ey>  end y coordinate         [default Inf]\n"
			"                        sd=<sd>  start depth                [default 0]\n"
			"                        ed=<ed>  end depth                [default Inf]\n"
			"                        sr=<sr>  view sample randomization:\n"
			"                            aperture  per aperture            [default]\n"
			"                            pixel     per pixel\n"
			"                        fs=<fs>  focal statistic, one of:\n"
			"                            mean      mean of views           [default]\n"
			"                            median    median of views\n"
			"\n"                                 
			BETWEEN_SECTIONS
			"Camera parameters:\n"
			HEADER_SPACE
			"--view-angle <x>  Set the initial diagonal view angle to <x> degrees.\n"
			"                     (Default is 43.7 degrees or device-specific.)\n"
			"--cpp-upper=<x>   Set upper bound <x> for camera parameter perturbation,\n"
			"                  in pixels or degrees [default is 32]\n"
			"--cpp-lower=<x>   Set lower bound <x> for camera parameter perturbation,\n"
			"                  in pixels or degrees [default is 0.125]\n"
			"--cpp-err-mean    Use RMS error to determine camera parameters.       [default]\n"
			"--cpp-err-median  Use median error to determine camera parameters.               \n"
			"--va-upper=<x>    View-angle perturbation upper bound in degrees   [default 32]\n"
			"--st <x>          Set stereo threshold to <x> pixels.            [default is 4]\n"
			"--vp-adjust       Adjust the view point                               [default]\n"
			"--vp-noadjust     Do not adjust the view point\n"
			"--vo-adjust       Adjust the view orientation                         [default]\n"
			"--vo-noadjust     Do not adjust the view orientation\n"
			BETWEEN_SECTIONS
			"Transformation file operations:\n"
			HEADER_SPACE
			"--3d-trans-load=x Load initial transformation settings from file x\n"
			"--3d-trans-save=x Save final transformation data in file x\n"
			BETWEEN_SECTIONS
			"Model rules:\n"
			HEADER_SPACE
			"--di-upper <x>    Decimate primary input resolution by at most 2^x  [default 0]\n"
			"--di-lower <x>    Decimate input resolutions by at least 2^x     [default is 0]\n"
			"--do-try <x>      Decimate output resolution by 2^x if possible  [default is 0]\n"
			"--oc              Clip scene to output regions.\n"
			"--no-oc           Do not clip scene to output regions.                [default]\n"
			"--fc <x>          Set front-clip to <x> (0 < x < 1)              [default is 0]\n"
			"--rc <x>          Set rear-clip to <x> (1 < x < inf)           [default is inf]\n"
			"--fx <x>          Set falloff exponent to <x>                    [default is 1]\n"
			"--tcem <x>        Set third-camera error multiplier to <x>       [default is 0]\n"
			"--oui <x>         Set occupancy update iterations to <x>        [default is 10]\n"
			"--pa <x>          Set pairwise ambiguity to <x>                  [default is 3]\n"
			"--pc <type>       Set the type of pairwise comparisons:\n"
			"                     auto     Determine comparisons automatically.    [default]\n"
			"                     all      Perform all comparisons.\n"
			"\n");
	}
	void scope() {
		banner("Argument scope [experimental]");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"List of arguments admitting scoping:\n"
			HEADER_SPACE
			"--gs              (see --ha for more details)\n"
			"--gs-mo           (see --ha for more details)\n"
			"--threads         (see --hp for more details)\n"
			"--per-cpu         (see --hp for more details)\n"
			"--perturb-upper   (see --ha for more details)\n"
			"--ev              (see --hx for more details)\n"
			"--black           (see --hx for more details)\n"
			BETWEEN_SECTIONS
			"Implicit file scope (implementation may be buggy):\n"
			HEADER_SPACE
			"arg1 file1 arg2   <arg2> applies only to files later than <file1>.\n"
			"                    Example:\n"
			"\n"
			"                         --ev 10 file1 --ev 20 file2\n"
			"\n"
			"                    (file1 has EV 10; file2 has EV 20.)\n"
			"\n"
			BETWEEN_SECTIONS
			"Hidden scope (implementation may be buggy):\n"
			HEADER_SPACE
			"[ args ]          Hide the effects of <args> options within [ ... ] scope\n"
			"                    Example:\n"
			"\n"
			"                         file1 [ --ev 20 file2 ] file3\n"
			"\n"
			"                    (The EV argument is limited to file2.)\n"
			"\n"
			BETWEEN_SECTIONS
			"Exposed scope (implementation may be buggy):\n"
			HEADER_SPACE
			"[ arg1 < arg2 > ] Expose the effects of <arg2> outside of [ ... ] scope\n"
			BETWEEN_SECTIONS
			"Fully-exposed scope (implementation may be buggy):\n"
			HEADER_SPACE
			"{ arg1 }          Same as [ < arg1 > ].\n"
			"\n");
	}
	void process() {
		banner("Process details");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Profiling:\n"
			HEADER_SPACE
			"--profile         Output performance data\n"
			BETWEEN_SECTIONS
			"Thread count:\n"
			"\n"
			"      If the CPU count cannot be determined, the default thread count\n"
			"      is 4.  Otherwise, the default is one thread per CPU.\n"
			"\n"
			HEADER_SPACE
			"--threads <n>     Use <n> threads.\n"
#ifndef USE_PTHREAD
			"\n"
			"                     NOTE: since this build of ALE does not link with a\n"
			"                           threading library, this option is not supported.\n"
			"                           To use this option, first rebuild with support\n"
			"                           for threads.\n"
			"\n"
#endif
			"--per-cpu <n>     Use <n> threads for each detected CPU.\n"
#ifndef USE_PTHREAD
			"\n"
			"                     NOTE: since this build of ALE does not link with a\n"
			"                           threading library, this option is not supported.\n"
			"                           To use this option, first rebuild with support\n"
			"                           for threads.\n"
			"\n"
#endif
			"\n");
	}
	void undocumented() {
		banner("Undocumented");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Point-spread functions:\n"
			HEADER_SPACE
			"--psf-match <args> Can be used to evaluate PSFs. [details are undocumented]\n"
			"                      <args> are:\n"
			"                          <mr> <mg> <mb> <ar> <ag> <ab>\n"
			"                      where:\n"
			"                          r[calibrated] = r[input] * <mr> + <ar>\n"
			BETWEEN_SECTIONS
			"Projective transformation calculator:\n"
			HEADER_SPACE
			"--ptcalc\n"
			BETWEEN_SECTIONS
			"Traverse subspaces for 3D candidate selection:\n"
			HEADER_SPACE
			"--subspace-traverse\n"
			"\n"
			);
	}
};

#undef BETWEEN_SECTIONS
#undef HEADER_SPACE
