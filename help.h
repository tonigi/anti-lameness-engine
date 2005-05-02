// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * Help messages
 */

#define BETWEEN_SECTIONS "\n"
#define HEADER_SPACE ""

class help {
private:
	const char *version;
	const char *invocation;
	FILE *help_stream;

	/*
	 * Banner
	 * 
	 * This function aids in locating the start of help output.
	 */
	void banner(const char *name) {
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"***********************************\n"
			"ALE Help Text, version %s\n"
			"Section:  %s\n"
			"***********************************\n",
		       version, name);
	}

public:
	help(const char *invocation, const char *version) {
		this->invocation = invocation;
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
			"--hq              Default settings.\n"
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
			"--hv              Video stream processing (Experimental).\n"
			"--h3              3D Modeling (Very Experimental).\n"
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
			"--q0              Low quality, high speed. [default]\n"
			"--qn              Low noise, moderate speed.\n"
			"--q1              Moderate quality and speed.\n"
			"--q2              High quality, low speed.\n"
			"--qr              Range-extended high quality.\n"
			BETWEEN_SECTIONS
			"q0 defaults:\n"
			HEADER_SPACE
			"   --dchain fine:box:1,triangle:2\n"
			"   --achain triangle:2\n"
			"   --mc 30\n"
			"   --ips 0\n"
			"   --exp-noextend\n"
			"   --no-cx\n"
			BETWEEN_SECTIONS
			"Low noise defaults:\n"
			HEADER_SPACE
			"   --dchain sinc*lanc:6\n"
			"   --achain sinc*lanc:6\n"
			"   --mc 50\n"
			"   --ips 0\n"
			"   --exp-noextend\n"
			"   --no-cx\n"
			BETWEEN_SECTIONS
			"q1 defaults:\n"
			HEADER_SPACE
			"   --dchain fine:sinc*lanc:6,sinc*lanc:6\n"
			"   --achain sinc*lanc:6\n"
			"   --mc 50\n"
			"   --ips 0\n"
			"   --exp-noextend\n"
			"   --no-cx\n"
			BETWEEN_SECTIONS
			"q2 defaults:\n"
			HEADER_SPACE
			"   --dchain sinc*lanc:8\n"
			"   --achain sinc*lanc:8\n"
			"   --no-mc\n"
			"   --ips 4\n"
			"   --exp-noextend\n"
			"   --no-cx\n"
			BETWEEN_SECTIONS
			"Range-extended defaults:\n"
			HEADER_SPACE
			"   --dchain sinc*lanc:8\n"
			"   --achain sinc*lanc:8\n"
			"   --no-mc\n"
			"   --ips 6\n"
			"   --exp-extend\n"
			"   --cx 0.7\n"
			"\n"
		       );
	}
	
	void file() {
		banner("File");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Bit depth options:\n"
			HEADER_SPACE
			"--8bpc            Write 8 bit per channel output [default]\n"
			"--16bpc           Write 16 bit per channel output\n"
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
			"--inc             Produce incremental output.  [default]\n"
			"--no-inc          Don't produce incremental output.\n"
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
			"--identity        Frames align closely with the original frame.  [default]\n"
			"--follow          Frames align closely with their immediate predecessor.\n"
			BETWEEN_SECTIONS
			"Alignment failure:\n"
			HEADER_SPACE
			"--fail-optimal    Frames beneath threshold are aligned optimally.  [default]\n"
			"--fail-default    Frames beneath threshold keep their default alignment.\n"
			BETWEEN_SECTIONS
			"Transformation file operations:\n"
			HEADER_SPACE
			"--trans-load=x    Load initial transformation settings from file x\n"
			"--trans-save=x    Save final transformation data in file x\n"
			BETWEEN_SECTIONS
			"Monte Carlo alignment (see --hq for defaults):\n"
			HEADER_SPACE
			"--mc <x>          Align using, on average, x%% of available pixels (0 < x < 100)\n"
			"--no-mc           Align using all pixels.\n"
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
			"--bda-mult=x      Barrel distortion adjustment multiplier   (0.0001 is default)\n"
			"--bda-rate=x      Barrel distortion rate of change maximum  (0.0004 is default)\n"   
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
			"                     NOTE: since this build was not configured for POSIX,\n"
			"                           this option is not supported.  To use this option,\n"
			"                           first re-build with POSIX=1.\n"
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
			"                     local     Local alignment only [default]\n"
			"                     inner     Alignment reference image inner region\n"
			"                     outer     Alignment reference image outer region\n"
			"                     all       Union of inner and outer\n"
			"                     central   inner if below threshold or better; else, outer.\n"
			"                     points    Align by control points.  Ignores gs-mo.\n"
			"--gs-mo <x>       Set <x> pixel min. overlap for global search. (16 is default)\n"
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
			"Irani-Peleg iterative solver (see --hq for defaults):\n"
			HEADER_SPACE
			"--ips <i>         Run <i> iterations.  (see also --hx, --hl and --hd)\n"
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
			"--bayer <b>       Set the bayer pattern to <b>, one of:\n"
			"                    (clockwise from top left pixel)\n"
			"                      rgbg      Red-green-blue-green\n"
			"                      gbgr      Green-blue-green-red\n"
			"                      grgb      Green-red-green-blue\n"
			"                      bgrg      Blue-green-red-green\n"
			"                      none      RGB-RGB-RGB-RGB\n"
			"                  Default is none or device-specific.\n"
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
			"                        stdin\n"
			"                        stdin_vg\n"
			"                        <p>+<p> (summation)\n"
			"                        <p>^<p> (convolution)\n"
			"                        <n>*<p> (multiplication by a scalar <n>)\n"
			"                     Default lpsf is either 'box=1.0' or device-specific.\n"
			"                     Default nlpsf is either disabled or device-specific.\n"
			BETWEEN_SECTIONS
			"Incremental rendering chains:\n"
			HEADER_SPACE
			"--wt <w>          Set weight threshold for defined pixels    (default is 0.8).\n"
			"--dchain <g>      Use chain <g> to render the default output.\n"
			"--ochain <g> <o>  Use chain <g> to render output file <o>.\n"
			"--achain <g>      Use chain <g> to render the alignment reference image.\n"
			"--afilter <s>     Use SSF <s> to interpolate points in alignment.\n"
			"                  Example chains:\n"
			"                     triangle:2                   ALE 0.6.0 merging (roughly)\n"
			"                     fine:box:1                   Drizzling (roughly)\n"
			"                     fine:sinc*lanc:8             High-frequency preservation\n"                     
			"                     last:nex:sinc*lanc:8         Useful for video stabilization\n"
			"                     fine:box:1,triangle:2        Multi-resolution rendering\n"
			"                     fine:sinc*lanc:8,sinc*lanc:8 Multi-resolution rendering\n"
			"                  More precisely, chains <g> are one of:\n"
			"                     <g1>,<g2>      Chain <g1> where defined; <g2> elsewhere.\n"
			"                     <i>            Use rendering invariant <i>.\n"
			"                  Rendering invariants <i> are:\n"
			"                     avg:<e>        Avg. (mean) value using SSFE <e>.\n"
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
			"                     fine:<f>       filter <f> at output image resolution.\n"
			"                     coarse:<f>     filter <f> at resolution MIN(in, out).\n"
			"                     <f>            Same as coarse:<f>.\n"
			"                  Sampling filters <f> are:\n"
			"                     sinc           Sinc filter. (very large diameter)\n"
			"                     lanc:<x>       Lanczos, diameter <x>.\n"
			"                     triangle:<x>   Triangle, diameter <x>.\n"
			"                     box:<x>        Box, diameter <x>.\n"
			"                     zero           Zero function.\n"
			"                     <f>*<f>        Pointwise multiplication (windowing).\n"
			"                  Defaults:\n"
			"                     dchain         (see --hq)\n"
			"                     achain         (see --hq)\n"
			"                     afilter        internal (approximates triangle:2)\n"
			"\n"
			);
	}
	void device() {
		banner("Device");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Device (may set PSF, bayer pattern, exposure, and view angle):\n"
			HEADER_SPACE
			"--device <d>      Set the capture device to <d>.\n"
			"                     Available devices:\n"
			"                        xvp610_320x240\n"
			"                        xvp610_640x480\n"
			"                        ov7620_raw_linear\n"
			"                        canon_300d_raw_linear\n"
			"                        canon_300d_raw_linear+50mm_1.4\n"
			"                        canon_300d_raw_linear+50mm_1.4@1.4\n"
			"                        canon_300d_raw_linear+50mm_1.8\n"
			"                        canon_300d_raw_linear+85mm_1.8\n"
			"\n"
		       );
	}
	void exclusion() {
		banner("Exclusion");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Exclusion regions:\n"
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
			"\n"
		       );
	}
	void exposure() {
		banner("Exposure");
		fprintf(help_stream, 
			BETWEEN_SECTIONS
			"Certainty-weighted rendering (see --hq for defaults):\n"
			HEADER_SPACE
			"--cx <x>          Render with certainty exponent <x>.\n"
			"--no-cx           Render with uniform certainty.\n"
			BETWEEN_SECTIONS
			"Options that may be useful when using --cx:\n"
			"--ips <i>         Uses one-sided certainty.  (see --hr for details)\n"
			"--exp-extend      Extends the output range.  (see below for details)\n"
			BETWEEN_SECTIONS
			"Exposure registration:\n"
			HEADER_SPACE
			"--exp-register    Register exposure between frames.  [default]\n"
			"--exp-noregister  Assume uniform exposure across all frames.\n"
			"--exp-meta-only   Use only image metadata for registering exposure.\n"
			BETWEEN_SECTIONS
			"Range extension (see --hq for defaults):\n"
			HEADER_SPACE
			"--exp-extend      Extend range to include all calculated values.\n"
			"--exp-noextend    Restrict to the original frame's range.\n"
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
			"                     <chain> is an incremental rendering chain.  (see --hl)\n"
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
			"                     stream\n"
			"                     tty [default]\n"
#ifndef USE_IOCTL
			"\n"
			"                     NOTE: since ALE was compiled without ioctl support,\n"
			"                     --ui=tty will behave identically to --ui=stream.\n"
			"                     For additional output, recompile with IOCTL=1.\n"
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
		banner("3D Modeling (very experimental)");
		fprintf(help_stream,
			BETWEEN_SECTIONS
			"Rendering:\n"
			HEADER_SPACE
			"--3dv <n> <o>     Render, to file <o>, colors as viewed from frame <n>.\n"
			"--3dd <n> <o>     Render, to file <o>, depths as viewed from frame <n>.\n"
			BETWEEN_SECTIONS
			"Capture device view angle:\n"
			HEADER_SPACE
			"--view-angle <x>  Set the diagonal view angle to <x> degrees.\n"
			"                     (Default is 43.7 degrees or device-specific.)\n"
			BETWEEN_SECTIONS
			"Control point parameters:\n"
			HEADER_SPACE
			"--cpp-upper=<x>   Set perturbation upper bound <x> for control-point realignment,\n"
			"                  in pixels or degrees [default is 32]\n"
			"--cpp-lower=<x>   Set perturbation lower bound <x> for control-point realignment,\n"
			"                  in pixels or degrees [default is 0.125]\n"
			"--va-upper=<x>    View-angle perturbation upper bound in degrees  [default is 32]\n"
			"--st <x>          Set stereo threshold <x>.                        [default is 4]\n"
			BETWEEN_SECTIONS
			"Model costs:\n"
			HEADER_SPACE
			"--ecm <x>         Set edge length cost multiplier <x>.         [default is 0.001]\n"
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
			"Color adjustment:\n"
			HEADER_SPACE
			"--exp-mult=c,r,b   Adjust all channels by <c>, red by <r>, and blue by <b>.\n"
			"\n"
			);
	}
};

#undef BETWEEN_SECTIONS
#undef HEADER_SPACE
