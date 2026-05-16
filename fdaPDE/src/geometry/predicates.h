// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __FDAPDE_ROBUST_PREDICATES_H
#define __FDAPDE_ROBUST_PREDICATES_H

#include "header_check.h"

namespace fdapde {

    // Robust adaptive floating-point geometric predicates (Shewchuk)
namespace robust {

    inline void two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bv = x - a;
        y = (a - (x - bv)) + (b - bv);
    }
    inline void fast_two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        y = b - (x - a);
    }
    inline void two_diff(double a, double b, double& x, double& y) {
        x = a - b;
        double bv = a - x;
        y = (a - (x + bv)) + (bv - b);
    }
    inline void two_diff_tail(double a, double b, double ab, double& tail) {
        double bvirt = a - ab;
        double avirt = ab + bvirt;
        double bround = bvirt - b;
        double around = a - avirt;
        tail = around + bround;
    }

    inline void split(double a, double& hi, double& lo) {
        constexpr double splitter = (static_cast<double>(1u << 27) + 1.0);
        double c = splitter * a;
        double abig = c - a;
        hi = c - abig;
        lo = a - hi;
    }
    inline void two_product(double a, double b, double& x, double& y) {
        x = a * b;
        double ah, al, bh, bl;
        split(a, ah, al);
        split(b, bh, bl);
        double err1 = x - (ah * bh);
        double err2 = err1 - (al * bh);
        double err3 = err2 - (ah * bl);
        y = (al * bl) - err3;
    }
    inline void square(double a, double& x, double& y) {
        two_product(a, a, x, y);
    }

    // helper: sum/diff of two double-double (h,l) ⊕ (H,L)
    inline void two_two_sum(double a1, double a0, double b1, double b0,
                            double& x3, double& x2, double& x1, double& x0) {
        double _i, _j, _0;
        two_sum(a0, b0, _i, x0);
        two_sum(a1, b1, _j, _0);
        two_sum(_i, _0, x2, x1);
        fast_two_sum(_j, x2, x3, x2);
    }
    inline void two_two_diff(double a1, double a0, double b1, double b0,
                            double& x3, double& x2, double& x1, double& x0) {
    double _i, _j, _0;
        two_diff(a0, b0, _i, x0);
        two_diff(a1, b1, _j, _0);
        two_sum(_i, _0, x2, x1);
        fast_two_sum(_j, x2, x3, x2);
    }

    inline int fast_expansion_sum_zeroelim(int elen, const double* e,
                                        int flen, const double* f,
                                        double* h) {
        double Q, Qnew, hh, bvirt;
        int eindex = 0, findex = 0, hindex = 0;
        double enow = e[0], fnow = f[0];
        if ((fnow > enow) == (fnow > -enow)) { 
            Q = enow; 
            eindex++; 
            enow = (eindex < elen) ? e[eindex] : 0.0; 
        }
        else  { 
            Q = fnow; 
            findex++; 
            fnow = (findex < flen) ? f[findex] : 0.0; 
        }
        while (eindex < elen && findex < flen) {
            if ((fnow > enow) == (fnow > -enow)) { 
                two_sum(Q, enow, Qnew, hh); 
                eindex++; 
                enow = (eindex < elen) ? e[eindex] : 0.0; 
            }
            else  { 
                two_sum(Q, fnow, Qnew, hh); 
                findex++; 
                fnow = (findex < flen) ? f[findex] : 0.0; 
            }
            if (hh != 0.0) { 
                h[hindex++] = hh; 
            }
            Q = Qnew;
        }
        while (eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            eindex++; 
            enow = (eindex < elen) ? e[eindex] : 0.0;
            if (hh != 0.0) { 
                h[hindex++] = hh; 
            }
            Q = Qnew;
        }
        while (findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            findex++; 
            fnow = (findex < flen) ? f[findex] : 0.0;
            if (hh != 0.0) { 
                h[hindex++] = hh; 
            }
            Q = Qnew;
        }
        if (Q != 0.0 || hindex == 0) { 
            h[hindex++] = Q; 
        }
        return hindex;
    }

    inline int scale_expansion_zeroelim(int elen, const double* e, double b, double* h) {
        double Q, sum, product1, product0;
        int eindex = 0, hindex = 0;
        double enow = e[0];
        two_product(enow, b, Q, h[0]); 
        hindex = 1;
        for (eindex = 1; eindex < elen; ++eindex) {
            enow = e[eindex];
            two_product(enow, b, product1, product0);
            two_sum(Q, product0, sum, h[hindex]); 
            if (h[hindex] != 0.0) ++hindex;
            fast_two_sum(product1, sum, Q, h[hindex]); 
            if (h[hindex] != 0.0) ++hindex;
        }
        if (Q != 0.0) 
            h[hindex++] = Q;
        return hindex;
    }

    inline double estimate(int elen, const double* e) {
        double Q = 0.0;
        for (int i = 0; i < elen; ++i) 
            Q += e[i];
        return Q;
    }

    inline double absolute(double x) { return x < 0.0 ? -x : x; }


    template <typename PointT>
        requires(internals::is_subscriptable<PointT, int>)
    constexpr double orient2d_fast(const PointT& A, const PointT& B, const PointT& C) {
          const double acx = A[0] - C[0];
          const double acy = A[1] - C[1];
          const double bcx = B[0] - C[0];
          const double bcy = B[1] - C[1];
          return acx * bcy - acy * bcx;
    }

    template <typename PointT>
        requires(internals::is_subscriptable<PointT, int>)
    constexpr int orient2d_sign(const PointT& A, const PointT& B, const PointT& C) {
          constexpr double eps  = std::numeric_limits<double>::epsilon();
          constexpr double errA = (3.0 + 16.0 * eps) * eps;
          constexpr double errB = (2.0 + 12.0 * eps) * eps;
          constexpr double errC = (9.0 + 64.0  * eps) * eps * eps;   
          constexpr double rErr = (3.0 + 8.0   * eps) * eps;

          // fast computation
          const double acx = A[0] - C[0];
          const double acy = A[1] - C[1];
          const double bcx = B[0] - C[0];
          const double bcy = B[1] - C[1];
          const double p1 = acx * bcy;
          const double p2 = acy * bcx;
          const double det = p1 - p2;
          const double detsum = absolute(acx * bcy) + absolute(acy * bcx);

          if (det >  errA * detsum) return +1;
          if (det < -errA * detsum) return -1;

          double acxtail, acytail, bcxtail, bcytail;
          two_diff_tail(A[0], C[0], acx, acxtail);
          two_diff_tail(A[1], C[1], acy, acytail);
          two_diff_tail(B[0], C[0], bcx, bcxtail);
          two_diff_tail(B[1], C[1], bcy, bcytail);

          if ((acxtail == 0.0) && (acytail == 0.0) && (bcxtail == 0.0) && (bcytail == 0.0)) 
               return (det > 0.0) ? +1 : (det < 0.0 ? -1 : 0);
          
          double det2 = det + (acx * bcytail + acytail * bcx) - (acy * bcxtail + acxtail * bcy);

          double scale = std::max({fabs(acx), fabs(acy), fabs(bcx), fabs(bcy)});
          if (std::fabs(det2 / (scale * scale)) < 1e-12) return 0;

          const double errbound = errB * detsum + rErr * absolute(det2);
          if (det2 >  errbound) return +1;
          if (det2 < -errbound) return -1; 

          // ---- exact fallback with expansions ----
          // (acx+acxtail)*(bcy+bcytail) - (acy+acytail)*(bcx+bcxtail)
          auto append_two_product = [](double x, double y, int curLen, double* cur) -> int {
                    double p1, p0, t[2];
                    two_product(x, y, p1, p0);  // p1=hi, p0=lo
                    t[0] = p0; t[1] = p1;
                    return fast_expansion_sum_zeroelim(curLen, cur, 2, t, cur);
          };

          double P[16]; int Plen = 0; 
          Plen = append_two_product(acx,      bcy,      Plen, P);
          if (bcytail != 0.0) 
                    Plen = append_two_product(acx, bcytail,  Plen, P);
          if (acxtail != 0.0) 
                    Plen = append_two_product(acxtail, bcy, Plen, P);
          if (acxtail != 0.0 && bcytail != 0.0)
                    Plen = append_two_product(acxtail, bcytail,  Plen, P);

          double Q[16]; int Qlen = 0; 
          Qlen = append_two_product(acy, bcx, Qlen, Q);
          if (bcxtail != 0.0) 
                    Qlen = append_two_product(acy, bcxtail, Qlen, Q);
          if (acytail != 0.0) 
                    Qlen = append_two_product(acytail, bcx, Qlen, Q);
          if (acytail != 0.0 && bcxtail != 0.0)
                    Qlen = append_two_product(acytail,  bcxtail,  Qlen, Q);

          // exact difference: P - Q
          double Qneg[16];
          for (int i = 0; i < Qlen; ++i) 
              Qneg[i] = -Q[i];

          double F[32]; 
          int Flen = fast_expansion_sum_zeroelim(Plen, P, Qlen, Qneg, F);
          const double det_exact = F[Flen - 1];

          return (det_exact > 0.0) ? +1 : (det_exact < 0.0 ? -1 : 0);
    }


    template <typename PointT>
        requires(internals::is_subscriptable<PointT, int>)
    constexpr double incircle_fast_det(const PointT& A, const PointT& B, const PointT& C, const PointT& D) {
        constexpr double eps  = std::numeric_limits<double>::epsilon();
        constexpr double errA = (10.0 + 96.0 * eps) * eps;

        const double adx = A[0] - D[0];
        const double ady = A[1] - D[1];
        const double bdx = B[0] - D[0];
        const double bdy = B[1] - D[1];
        const double cdx = C[0] - D[0];
        const double cdy = C[1] - D[1];

        const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
        const double cdxady = cdx * ady, adxcdy = adx * cdy;
        const double adxbdy = adx * bdy, bdxady = bdx * ady;

        const double alift = adx*adx + ady*ady;
        const double blift = bdx*bdx + bdy*bdy;
        const double clift = cdx*cdx + cdy*cdy;

        const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

        const double permanent = (absolute(bdxcdy) + absolute(cdxbdy)) * absolute(alift) + (absolute(cdxady) + absolute(adxcdy)) * absolute(blift) + (absolute(adxbdy) + absolute(bdxady)) * absolute(clift);

        const double errbound = errA * permanent;  
        if (det > errbound || -det > errbound) 
            return det;

        // force a fallback (incircle_adapt)
        return std::numeric_limits<double>::quiet_NaN();
    }


    inline double incircle_adapt(const double pa[2], const double pb[2],const double pc[2], const double pd[2],double permanent)
    {
        constexpr double eps  = std::numeric_limits<double>::epsilon();
        constexpr double iccerrboundB = (4.0  + 48.0  * eps) * eps;          
        constexpr double iccerrboundC = (44.0 + 576.0 * eps) * eps * eps;    
        constexpr double resulterrbound = (3.0  + 8.0   * eps) * eps;          

        double adx, bdx, cdx, ady, bdy, cdy;
        double det, errbound;

        double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
        double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
        double bc[4], ca[4], ab[4];
        double bc3, ca3, ab3;
        double axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
        int axbclen, axxbclen, aybclen, ayybclen, alen;
        double bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
        int bxcalen, bxxcalen, bycalen, byycalen, blen;
        double cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
        int cxablen, cxxablen, cyablen, cyyablen, clen;
        double abdet[64];
        int ablen;
        double fin1[1152], fin2[1152];
        double *finnow, *finother, *finswap;
        int finlength;

        double adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
        double adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
        double adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
        double aa[4], bb[4], cc[4];
        double aa3, bb3, cc3;
        double ti1, tj1;
        double ti0, tj0;
        double u[4], v[4];
        double u3, v3;
        double temp8[8], temp16a[16], temp16b[16], temp16c[16];
        double temp32a[32], temp32b[32], temp48[48], temp64[64];
        int temp8len, temp16alen, temp16blen, temp16clen;
        int temp32alen, temp32blen, temp48len, temp64len;
        double axtbb[8], axtcc[8], aytbb[8], aytcc[8];
        int axtbblen, axtcclen, aytbblen, aytcclen;
        double bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
        int bxtaalen, bxtcclen, bytaalen, bytcclen;
        double cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
        int cxtaalen, cxtbblen, cytaalen, cytbblen;
        double axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
        int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
        double axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
        int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
        double axtbctt[8], aytbctt[8], bxtcatt[8];
        double bytcatt[8], cxtabtt[8], cytabtt[8];
        int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
        double abt[8], bct[8], cat[8];
        int abtlen, bctlen, catlen;
        double abtt[4], bctt[4], catt[4];
        int abttlen, bcttlen, cattlen;
        double abtt3, bctt3, catt3;
        double negate;

        double bvirt;
        double avirt, bround, around;
        double c;
        double abig;
        double ahi, alo, bhi, blo;
        double _i, _j;
        double _0;

        adx = pa[0] - pd[0];
        bdx = pb[0] - pd[0];
        cdx = pc[0] - pd[0];
        ady = pa[1] - pd[1];
        bdy = pb[1] - pd[1];
        cdy = pc[1] - pd[1];

        two_product(bdx, cdy, bdxcdy1, bdxcdy0);
        two_product(cdx, bdy, cdxbdy1, cdxbdy0);
        two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
        bc[3] = bc3;
        axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
        axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
        aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
        ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
        alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

        two_product(cdx, ady, cdxady1, cdxady0);
        two_product(adx, cdy, adxcdy1, adxcdy0);
        two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
        ca[3] = ca3;
        bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
        bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
        bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
        byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
        blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

        two_product(adx, bdy, adxbdy1, adxbdy0);
        two_product(bdx, ady, bdxady1, bdxady0);
        two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
        ab[3] = ab3;
        cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
        cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
        cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
        cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
        clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

        ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
        finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

        det = estimate(finlength, fin1);
        errbound = iccerrboundB * permanent;
        if ((det >= errbound) || (-det >= errbound)) 
            return det;

        two_diff_tail(pa[0], pd[0], adx, adxtail);
        two_diff_tail(pa[1], pd[1], ady, adytail);
        two_diff_tail(pb[0], pd[0], bdx, bdxtail);
        two_diff_tail(pb[1], pd[1], bdy, bdytail);
        two_diff_tail(pc[0], pd[0], cdx, cdxtail);
        two_diff_tail(pc[1], pd[1], cdy, cdytail);
        if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) 
            return det;

        errbound = iccerrboundC * permanent + resulterrbound * absolute(det);
        det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
            + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
            + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
        if ((det >= errbound) || (-det >= errbound)) 
            return det;

        finnow = fin1;
        finother = fin2;

        if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0)) {
            square(adx, adxadx1, adxadx0);
            square(ady, adyady1, adyady0);
            two_two_sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
            aa[3] = aa3;
        }
        if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0)) {
            square(bdx, bdxbdx1, bdxbdx0);
            square(bdy, bdybdy1, bdybdy0);
            two_two_sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
            bb[3] = bb3;
        }
        if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0)) {
            square(cdx, cdxcdx1, cdxcdx0);
            square(cdy, cdycdy1, cdycdy0);
            two_two_sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
            cc[3] = cc3;
        }

        if (adxtail != 0.0) {
            axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
            temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx, temp16a);

            axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
            temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

            axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
            temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

            temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (adytail != 0.0) {
            aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
            temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady, temp16a);

            aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
            temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

            aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
            temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

            temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (bdxtail != 0.0) {
            bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
            temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx,temp16a);

            bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
            temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

            bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
            temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

            temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (bdytail != 0.0) {
            bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
            temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy, temp16a);

            bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
            temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

            bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
            temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

            temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (cdxtail != 0.0) {
            cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
            temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx, temp16a);

            cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
            temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

            cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
            temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

            temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (cdytail != 0.0) {
            cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
            temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy, temp16a);

            cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
            temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

            cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
            temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

            temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }

        if ((adxtail != 0.0) || (adytail != 0.0)) {
            if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0)) {
                two_product(bdxtail, cdy, ti1, ti0);
                two_product(bdx, cdytail, tj1, tj0);
                two_two_sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
                u[3] = u3;
                negate = -bdy;
                two_product(cdxtail, negate, ti1, ti0);
                negate = -bdytail;
                two_product(cdx, negate, tj1, tj0);
                two_two_sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
                v[3] = v3;
                bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

                two_product(bdxtail, cdytail, ti1, ti0);
                two_product(cdxtail, bdytail, tj1, tj0);
                two_two_diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
                bctt[3] = bctt3;
                bcttlen = 4;
            } else {
                bct[0] = 0.0;
                bctlen = 1;
                bctt[0] = 0.0;
                bcttlen = 1;
            }

            if (adxtail != 0.0) {
                temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
                axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
                temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx, temp32a);
                temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
                if (bdytail != 0.0) {
                    temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
                    temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, temp16a);
                    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
                if (cdytail != 0.0) {
                    temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
                    temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, temp16a);
                    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
                temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, temp32a);
                axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
                temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx, temp16a);
                temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, temp16b);
                temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
            }
            if (adytail != 0.0) {
                temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
                aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
                temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady, temp32a);
                temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;

                temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, temp32a);
                aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
                temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady, temp16a);
                temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, temp16b);
                temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
            }
        }

        if ((bdxtail != 0.0) || (bdytail != 0.0)) {
            if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0)) {
                two_product(cdxtail, ady, ti1, ti0);
                two_product(cdx, adytail, tj1, tj0);
                two_two_sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
                u[3] = u3;
                negate = -cdy;
                two_product(adxtail, negate, ti1, ti0);
                negate = -cdytail;
                two_product(adx, negate, tj1, tj0);
                two_two_sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
                v[3] = v3;
                catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

                two_product(cdxtail, adytail, ti1, ti0);
                two_product(adxtail, cdytail, tj1, tj0);
                two_two_diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
                catt[3] = catt3;
                cattlen = 4;
            } else {
                cat[0] = 0.0;
                catlen = 1;
                catt[0] = 0.0;
                cattlen = 1;
            }

            if (bdxtail != 0.0) {
                temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
                bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
                temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx, temp32a);
                temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
                if (cdytail != 0.0) {
                    temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
                    temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, temp16a);
                    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
                if (adytail != 0.0) {
                    temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
                    temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, temp16a);
                    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
                temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, temp32a);
                bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
                temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx, temp16a);
                temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, temp16b);
                temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
            }
            if (bdytail != 0.0) {
                temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
                bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
                temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy, temp32a);
                temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;

                temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, temp32a);
                bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
                temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy, temp16a);
                temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, temp16b);
                temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
            }
        }
        if ((cdxtail != 0.0) || (cdytail != 0.0)) {
            if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0)) {
                two_product(adxtail, bdy, ti1, ti0);
                two_product(adx, bdytail, tj1, tj0);
                two_two_sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
                u[3] = u3;
                negate = -ady;
                two_product(bdxtail, negate, ti1, ti0);
                negate = -adytail;
                two_product(bdx, negate, tj1, tj0);
                two_two_sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
                v[3] = v3;
                abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

                two_product(adxtail, bdytail, ti1, ti0);
                two_product(bdxtail, adytail, tj1, tj0);
                two_two_diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
                abtt[3] = abtt3;
                abttlen = 4;
            } else {
                abt[0] = 0.0;
                abtlen = 1;
                abtt[0] = 0.0;
                abttlen = 1;
            }

            if (cdxtail != 0.0) {
                temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
                cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
                temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx, temp32a);
                temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
                if (adytail != 0.0) {
                    temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
                    temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, temp16a);
                    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
                if (bdytail != 0.0) {
                    temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
                    temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, temp16a);
                    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }

                temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, temp32a);
                cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
                temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx, temp16a);
                temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, temp16b);
                temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
            }
            if (cdytail != 0.0) {
                temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
                cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
                temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy, temp32a);
                temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;

                temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, temp32a);
                cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
                temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy, temp16a);
                temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, temp16b);
                temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
                finswap = finnow; 
                finnow = finother; 
                finother = finswap;
            }
        }

        return finnow[finlength - 1];
    }


    template <typename PointT>
        requires(internals::is_subscriptable<PointT, int>)
    constexpr int incircle_sign_ccw(const PointT& A, const PointT& B, const PointT& C, const PointT& D) {
        double det = incircle_fast_det(A,B,C,D);
        if (std::isnan(det)) {
            const double pa[2]{ A[0], A[1] };
            const double pb[2]{ B[0], B[1] };
            const double pc[2]{ C[0], C[1] };
            const double pd[2]{ D[0], D[1] };

            const double adx = pa[0] - pd[0], ady = pa[1] - pd[1];
            const double bdx = pb[0] - pd[0], bdy = pb[1] - pd[1];
            const double cdx = pc[0] - pd[0], cdy = pc[1] - pd[1];
            const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
            const double cdxady = cdx * ady, adxcdy = adx * cdy;
            const double adxbdy = adx * bdy, bdxady = bdx * ady;
            const double alift = adx*adx + ady*ady;
            const double blift = bdx*bdx + bdy*bdy;
            const double clift = cdx*cdx + cdy*cdy;
            const auto absv = [](double x){ return x < 0 ? -x : x; };
            const double permanent = (absv(bdxcdy) + absv(cdxbdy)) * absv(alift) + (absv(cdxady) + absv(adxcdy)) * absv(blift) + (absv(adxbdy) + absv(bdxady)) * absv(clift);

            det = incircle_adapt(pa,pb,pc,pd,permanent);
        }
        return (det > 0.0) ? +1 : (det < 0.0 ? -1 : 0);
    }

    // assumes counterclockwise ordering of A, B, C
    template <typename PointT>
        requires(internals::is_subscriptable<PointT, int>)
    constexpr bool in_circle_ccw(const PointT& A, const PointT& B, const PointT& C, const PointT& D) {
        return incircle_sign_ccw(A,B,C,D) > 0;
    }


}  // namespace robust

}  // namespace fdapde

#endif // __FDAPDE_ROBUST_PREDICATES_H