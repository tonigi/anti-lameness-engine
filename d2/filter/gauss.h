// code by HJ Hornbeck, based on code copyright David Hilvert

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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
*/

#ifndef __d2_filter_gauss_h__
#define __d2_filter_gauss_h__

/*
 * A Gaussian filter.
 */

class gauss : public filter {
private:
    ale_real sigma;
    ale_real diameter;    /** measured in standard deviations **/

    /*
     * The heavy-lifting function:
     */
    ale_real _calc(ale_pos p) const {

        /** catch the trivial cases **/
        if ( p == 0 )                return 1;
/**         if ( p > (sigma*diameter) )         return 0;
             ^^ unnecessary, gauss is well-behaved
                outside its radius **/

        /** calculate the rest **/
        return exp( -p*p / (sigma*sigma) );

    }

    /** extend the heavy-lifting function to handle a coordinate **/
    ale_real _calc(point p) const {

        return _calc(  sqrt( (p[0]*p[0]) + (p[1]*p[1]) )  );

/** doesn't work as well:
        return _calc( p[0] ) * _calc( p[1] ); **/
    }

public:

    /*
     * Size of filter support, or how big a radius does this filter
     * touch?
     */
    ale_real support() const {
        return sigma*diameter;
    }

    /*
     * Response of filter at point p
     */
    virtual ale_real response(point p) const {
        return _calc(p);
    }

    /** only compare standard deviations, since this version fixes
         the diameter across all instances **/
    virtual int equals(const filter *f) const {
        if (typeid(*f) == typeid(*this))
            return ((gauss *)f)->sigma == sigma;
        return 0;
    }

    gauss(ale_real sigma) {
        this->sigma = sigma;
        this->diameter = 2;    /** fixed, standard for imaging **/
    }

};
#endif
