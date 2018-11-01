#ifndef JACOBIAN_H_INCLUDED
#define JACOBIAN_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "fda-io.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "ekg-proc.h"
#include "lapacke.h"

using namespace std;
double rmin, metBC_coeff, dr;

//  for LAPACKE_dgbsv(): jac[ (kl + ku + 1) + (ldab - 1)*j + i ]  =  jac[ i, j ]
void set_jacCMabpclean(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		       int npts, int kl, int ku, int ldab);
void set_jacCMabpslow(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab);
void set_jacCMabpfast(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab);
void set_jacCM_ab_slow(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		       int npts, int kl, int ku, int ldab);
inline void set_jac_alphaCM(VD& jac, const VD& f_xi, const VD& f_pi,
			    const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			    int npts, int kl, int ku, int ldab);
inline void set_jac_betaCM(VD& jac, const VD& f_xi, const VD& f_pi,
			   const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			   int npts, int kl, int ku, int ldab);
inline void set_jac_psiCM(VD& jac, const VD& f_xi, const VD& f_pi,
			  const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			  int npts, int kl, int ku, int ldab);
inline int jac_ind(int j, int k) { return (4 + j + 6*k); }
// ***********************  JACOBIAN FUNCTIONS  ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_aa(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		  MAPID& r, int k)
{
  return r[NEG2INDRSQ] - 8*M_PI*sq(f_pi[k]) +
    2*pw4(f_ps[k])*sq(ddr_c(f_be,r,k) - f_be[k]*r[-k]) / (3*sq(f_al[k]));
}

inline dbl jac_aa_pm(const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r, int k, int p_m)
{
  return r[INDRSQ] + p_m*r[INDR]*((ddr_c(f_ps,r,k)/f_ps[k]) + r[-k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_bb(const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r, int k)
{
  return r[NEG2INDRSQ] - r[-k]*(2*r[-k] + 6*(ddr_c(f_ps,r,k)/f_ps[k]) - (ddr_c(f_al,r,k)/f_al[k]));
}

inline dbl jac_bb_pm(const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r, int k, int p_m)
{
  return r[INDRSQ] + p_m*r[IN2DR]*(2*r[-k] + 6*(ddr_c(f_ps,r,k)/f_ps[k]) - (ddr_c(f_al,r,k)/f_al[k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_pp(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		  MAPID& r, int k)
{
  return r[NEG2INDRSQ] + M_PI*(sq(f_xi[k]) + sq(f_pi[k])) +
    5*pw4(f_ps[k])*sq(ddr_c(f_be,r,k) - f_be[k]*r[-k]) / (12*sq(f_al[k]));
}

inline dbl jac_pp_pm(const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r, int k, int p_m)
{
  return r[INDRSQ] + p_m*r[INDR]*r[-k];
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// **********************************************************
// **********************************************************
//                    POPULATING JACOBIAN
// **********************************************************
// **********************************************************
////////////////////////////////////////////////////////////////////////////////////////////////

inline void set_jac_alphaCM(VD& jac, const VD& f_xi, const VD& f_pi,
			    const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			    int npts, int kl, int ku, int ldab)
{
  int k = kl + ku, one_past_last = npts - 3;
  int jm1, j = 0, jp1 = 1;
  // col 0
  jac[k] = -3*r[IN2DR];
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  // col 1
  jm1 = j; j = jp1; ++jp1;
  k += kl + 1;
  jac[++k] = 4*r[IN2DR];
  jac[++k] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  // col 2
  jm1 = j; j = jp1; ++jp1;
  k += kl;
  jac[++k] = -1*r[IN2DR];
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  jm1 = j; j = jp1; ++jp1;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(f_al)
    k += kl;
    jac[++k] = 0;
    jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
    jac[++k] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
    jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-3
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 1*r[IN2DR];
  jm1 = j; ++j;
  // col N-2
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = -4*r[IN2DR];
  // col N-1
  k += 1 + kl;
  jac[++k] = 0;
  jac[++k] = jac_aa_pm(f_al, f_be, f_ps, r, j, 1);
  jac[++k] = 3*r[IN2DR] + r[INRMAX];
  return;
}

inline void set_jac_betaCM(VD& jac, const VD& f_xi, const VD& f_pi,
			   const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			   int npts, int kl, int ku, int ldab)
{
  int k = kl + ku, one_past_last = npts - 3;
  int jm1, j = 0, jp1 = 1;
  // col 0
  jac[k] = 1;
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  // col 1
  jm1 = j; j = jp1; ++jp1;
  k += kl + 1;
  jac[++k] = 0;
  jac[++k] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  // col 2
  jm1 = j; j = jp1; ++jp1;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  jm1 = j; j = jp1; ++jp1;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(f_be)   
    k += kl;
    jac[++k] = 0;
    jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
    jac[++k] = jac_bb(f_al, f_be, f_ps, r, j);
    jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
    jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-3
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 1*r[IN2DR];
  jm1 = j; ++j;
  // col N-2
  
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++k] = -4*r[IN2DR];
  // col N-1
  k += 1 + kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(f_al, f_be, f_ps, r, j, 1);
  jac[++k] = 3*r[IN2DR] + r[INRMAX];
  return;
}


inline void set_jac_psiCM(VD& jac, const VD& f_xi, const VD& f_pi,
			  const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			  int npts, int kl, int ku, int ldab)
{
  int k = kl + ku, one_past_last = npts - 3;
  int jm1, j = 0, jp1 = 1;

  // col 0
  jac[k] = -3*r[IN2DR];
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  // col 1
  jm1 = j; j = jp1; ++jp1;
  k += kl + 1;
  jac[++k] = 4*r[IN2DR];
  jac[++k] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  // col 2
  jm1 = j; j = jp1; ++jp1;

  k += kl;
  jac[++k] = -1*r[IN2DR];
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 0;
  jm1 = j; ++jp1;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(f_ps)
    
    k += kl;
    jac[++k] = 0;
    jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
    jac[++k] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
    jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-3

  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++k] = 1*r[IN2DR];
  jm1 = j; ++j;
  // col N-2
  
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++k] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++k] = -4*r[IN2DR];

  // col N-1
  k += 1 + kl;
  jac[++k] = 0;
  jac[++k] = jac_pp_pm(f_al, f_be, f_ps, r, j, 1);
  jac[++k] = 3*r[IN2DR] + r[INRMAX];
  return;
}

void set_jacCMabpclean(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		       int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 3;
  int kal = 4; // kl+ku
  int kbe = npts*7 + 2; // npts*ldab + kl
  int kps = npts*14 + 2; // 2*npts*ldab + kl
  int jm1, j = 0, jp1 = 1;
  // EXTRA SPACE
  jac[kbe] = 0; jac[kps] = 0;
  jac[++kbe] = 0; jac[++kps] = 0; // skip al b/c no need to set these 2
  // ROW 0, COL 0
  jac[kal] = -3*r[IN2DR]; jac[++kbe] = 1; jac[++kps] = -3*r[IN2DR];
  // ROW 1, COL 0
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  // ROW 2, COL 0
  jac[++kal] = 0; jac[++kbe] = 0; jac[++kps] = 0;
  // ***********************************************************
  // NEXT COL: 1
  kal += 4; kbe += 3; kps += 3;
  jm1 = j; j = jp1; ++jp1;
  jac[kbe] = 0; jac[kps] = 0; // skip al b/c no need to set this 1
  // ROW 0, COL 1
  jac[kal] = 4*r[IN2DR]; jac[++kbe] = 0; jac[++kps] = 4*r[IN2DR];
  // ROW 1, COL 1
  jac[++kal] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++kbe] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++kps] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  // ROW 2, COL 1
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  // ROW 3, COL 1
  jac[++kal] = 0; jac[++kbe] = 0; jac[++kps] = 0;
  // ***********************************************************
  // NEXT COL: 2
  kal += 3; kbe += 3; kps += 3;
  jm1 = j; j = jp1; ++jp1;
  // ROW 0, COL 2
  jac[kal] = -1*r[IN2DR]; jac[kbe] = 0; jac[kps] = -1*r[IN2DR];
  // ROW 1, COL 2
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
  // ROW 2, COL 2
  jac[++kal] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++kbe] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++kps] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  // ROW 3, COL 2
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  // ROW 4, COL 2
  jac[++kal] = 0; jac[++kbe] = 0; jac[++kps] = 0;
  // ***********************************************************
  // NEXT COL: j
  kal += 3; kbe += 3; kps += 3;
  jm1 = j; ++jp1;
  for (j = 3; j < one_past_last; ++j) {
    // ROW j-2, COL j
    jac[kal] = 0; jac[kbe] = 0; jac[kps] = 0;
    // ROW j-1, COL j
    jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
    jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
    jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
    // ROW j, COL j
    jac[++kal] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    jac[++kbe] = jac_bb(f_al, f_be, f_ps, r, j);
    jac[++kps] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    // ROW j+1, COL j
    jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
    jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
    jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
    // ROW j+2, COL j
    jac[++kal] = 0; jac[++kbe] = 0; jac[++kps] = 0;
    // ***********************************************************
    // NEXT COL: j+1
    kal += 3; kbe += 3; kps += 3;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // ROW lastpt-4, COL lastpt-2
  jac[kal] = 0; jac[kbe] = 0; jac[kps] = 0;
  // ROW lastpt-3, COL lastpt-2
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
  // ROW lastpt-2, COL lastpt-2
  jac[++kal] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++kbe] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++kps] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  // ROW lastpt-1, COL lastpt-2
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jp1, -1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jp1, -1);
  // ROW lastpt, COL lastpt-2
  jac[++kal] = 1*r[IN2DR]; jac[++kbe] = 1*r[IN2DR]; jac[++kps] = 1*r[IN2DR];
  // ***********************************************************
  // NEXT COL: lastpt-1
  kal += 3; kbe += 3; kps += 3;
  jm1 = j; j = jp1;
  // ROW lastpt-3, COL lastpt-1
  jac[kal] = 0; jac[kbe] = 0; jac[kps] = 0;
  // ROW lastpt-2, COL lastpt-1
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, jm1, 1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, jm1, 1);
  // ROW lastpt-1, COL lastpt-1
  jac[++kal] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  jac[++kbe] = jac_bb(f_al, f_be, f_ps, r, j);
  jac[++kps] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
  // ROW lastpt, COL lastpt-1
  jac[++kal] = -4*r[IN2DR]; jac[++kbe] = -4*r[IN2DR]; jac[++kps] = -4*r[IN2DR];
  // ETRA SPACE
  jac[++kal] = 0; jac[++kbe] = 0; // no need to set f_ps here
  // ***********************************************************
  // NEXT COL: lastpt
  kal += 3; kbe += 3; kps += 4;
  // let j here be jm1 (save time w/o incrementing)
  dbl cR = 3*r[IN2DR] + r[INRMAX];
  // ROW lastpt-2, COL lastpt
  jac[kal] = 0; jac[kbe] = 0; jac[kps] = 0;
  // ROW lastpt-1, COL lastpt
  jac[++kal] = jac_aa_pm(f_al, f_be, f_ps, r, j, 1);
  jac[++kbe] = jac_bb_pm(f_al, f_be, f_ps, r, j, 1);
  jac[++kps] = jac_pp_pm(f_al, f_be, f_ps, r, j, 1);
  // ROW lastpt, COL lastpt-1
  jac[++kal] = cR; jac[++kbe] = cR; jac[++kps] = cR;
  // EXTRA SPACE
  jac[++kal] = 0; jac[++kbe] = 0;
  jac[++kal] = 0; jac[++kbe] = 0; // no need to set f_ps here
  // DONE!!!
  return;
}


void set_jacCMabpfast(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int k, kbe, kps, j = 0, jbe = npts, jps = 2*npts;
  // ROW 0, COL 0
  k = j; kbe = jbe; kps = jps;
  jac[jac_ind(j,k)] = r[JAC_N00];
  jac[jac_ind(jbe,kbe)] = 1;
  jac[jac_ind(jps,kps)] = r[JAC_N00];
  // ROW 0, COL 1
  k = j + 1; kbe = jbe + 1; kps = jps + 1;
  jac[jac_ind(j,k)] = r[JAC_N01];
  jac[jac_ind(jbe,kbe)] = 0;
  jac[jac_ind(jps,kps)] = r[JAC_N01];
  // ROW 0, COL 2
  k = j + 2; kbe = jbe + 2; kps = jps + 2;
  jac[jac_ind(j,k)] = r[JAC_N02];
  jac[jac_ind(jbe,kbe)] = 0;
  jac[jac_ind(jps,kps)] = r[JAC_N02];

  double drin_p, drin_m, dps_ps, dlogp6_a, p4db2_a2;
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    jps = jbe + npts;

    drin_p = r[INDRSQ] + r[-k]*r[INDR];
    drin_m = r[INDRSQ] - r[-k]*r[INDR];
    dps_ps = ddr_c(f_ps,r,k) / f_ps[k];
    dlogp6_a = 6*dps_ps - (ddr_c(f_al,r,k)/f_al[k]);
    p4db2_a2 = sq(f_ps[k]) * (ddr_c(f_be,r,k) - r[-k]*f_be[k]) / f_al[k];
    p4db2_a2 = sq(p4db2_a2);    
    
    // ROW j, COL j-1
    k = j-1; kbe = jbe-1; kps = jps-1;
    jac[jac_ind(j,k)] = drin_m - r[-k]*dps_ps;
    jac[jac_ind(jbe,kbe)] = drin_m - r[IN2DR]*dlogp6_a;
    jac[jac_ind(jps,kps)] = drin_m;
    // ROW j, COL j
    k = j; kbe = jbe; kps = jps;
    jac[jac_ind(j,k)] = r[NEG2INDRSQ] + r[TWO_THIRDS]*p4db2_a2 - r[EIGHT_PI]*sq(f_pi[k]);
    jac[jac_ind(jbe,kbe)] = r[NEG2INDRSQ] - r[-k]*(2*r[-k] + dlogp6_a);
    jac[jac_ind(jps,kps)] = r[NEG2INDRSQ] + r[FIVE_TWELFTHS]*p4db2_a2 + M_PI*(sq(f_xi[k]) + sq(f_pi[k]));
    // ROW j+1, COL j
    k = j + 1; kbe = jbe + 1; kps = jps + 1;
    jac[jac_ind(j,k)] = drin_p + r[-k]*dps_ps;
    jac[jac_ind(jbe,kbe)] = drin_p + r[IN2DR]*dlogp6_a;
    jac[jac_ind(jps,kps)] = drin_p;
  }
  j = one_past_last;
  jbe = j + npts;
  jps = jbe + npts;
  // ROW lastpt, COL lastpt-2
  k = j - 2; kbe = jbe - 2; kps = jps - 2;
  jac[jac_ind(j,k)] = r[JAC_RRM2];
  jac[jac_ind(jbe,kbe)] = r[JAC_RRM2];
  jac[jac_ind(jps,kps)] = r[JAC_RRM2];
  // ROW lastpt, COL lastpt-1
  k = j - 1; kbe = jbe - 1; kps = jps - 1;
  jac[jac_ind(j,k)] = r[JAC_RRM1];
  jac[jac_ind(jbe,kbe)] = r[JAC_RRM1];
  jac[jac_ind(jps,kps)] = r[JAC_RRM1];
  // ROW lastpt, COL lastpt
  k = j; kbe = jbe; kps = jps;
  jac[jac_ind(j,k)] = r[JAC_RR];
  jac[jac_ind(jbe,kbe)] = r[JAC_RR];
  jac[jac_ind(jps,kps)] = r[JAC_RR];
  return;
}

void set_jacCMabpslow(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		       int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int k, kbe, kps, j = 0, jbe = npts, jps = 2*npts;
  // ROW 0, COL 0
  k = j; kbe = jbe; kps = jps;
  jac[jac_ind(j,k)] = r[JAC_N00];
  jac[jac_ind(jbe,kbe)] = 1;
  jac[jac_ind(jps,kps)] = r[JAC_N00];
  // ROW 0, COL 1
  k = j + 1; kbe = jbe + 1; kps = jps + 1;
  jac[jac_ind(j,k)] = r[JAC_N01];
  jac[jac_ind(jbe,kbe)] = 0;
  jac[jac_ind(jps,kps)] = r[JAC_N01];
  // ROW 0, COL 2
  k = j + 2; kbe = jbe + 2; kps = jps + 2;
  jac[jac_ind(j,k)] = r[JAC_N02];
  jac[jac_ind(jbe,kbe)] = 0;
  jac[jac_ind(jps,kps)] = r[JAC_N02];
  
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    jps = jbe + npts;
    // ROW j, COL j-1
    k = j-1; kbe = jbe-1; kps = jps-1;
    jac[jac_ind(j,k)] = jac_aa_pm(f_al, f_be, f_ps, r, j, -1);
    jac[jac_ind(jbe,kbe)] = jac_bb_pm(f_al, f_be, f_ps, r, j, -1);
    jac[jac_ind(jps,kps)] = jac_pp_pm(f_al, f_be, f_ps, r, j, -1);
    // ROW j, COL j
    k = j; kbe = jbe; kps = jps;
    jac[jac_ind(j,k)] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    jac[jac_ind(jbe,kbe)] = jac_bb(f_al, f_be, f_ps, r, j);
    jac[jac_ind(jps,kps)] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    // ROW j+1, COL j
    k = j + 1; kbe = jbe + 1; kps = jps + 1;
    jac[jac_ind(j,k)] = jac_aa_pm(f_al, f_be, f_ps, r, j, 1);
    jac[jac_ind(jbe,kbe)] = jac_bb_pm(f_al, f_be, f_ps, r, j, 1);
    jac[jac_ind(jps,kps)] = jac_pp_pm(f_al, f_be, f_ps, r, j, 1);
  }
  j = one_past_last;
  jbe = j + npts;
  jps = jbe + npts;
  // ROW lastpt, COL lastpt-2
  k = j - 2; kbe = jbe - 2; kps = jps - 2;
  jac[jac_ind(j,k)] = r[JAC_RRM2];
  jac[jac_ind(jbe,kbe)] = r[JAC_RRM2];
  jac[jac_ind(jps,kps)] = r[JAC_RRM2];
  // ROW lastpt, COL lastpt-1
  k = j - 1; kbe = jbe - 1; kps = jps - 1;
  jac[jac_ind(j,k)] = r[JAC_RRM1];
  jac[jac_ind(jbe,kbe)] = r[JAC_RRM1];
  jac[jac_ind(jps,kps)] = r[JAC_RRM1];
  // ROW lastpt, COL lastpt
  k = j; kbe = jbe; kps = jps;
  jac[jac_ind(j,k)] = r[JAC_RR];
  jac[jac_ind(jbe,kbe)] = r[JAC_RR];
  jac[jac_ind(jps,kps)] = r[JAC_RR];
  return;
}

void set_jacCM_ab_slow(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		       int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int j = 0, jbe = npts;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = r[JAC_N00];
  jac[jac_ind(jbe,jbe)] = 1;
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = r[JAC_N01];
  jac[jac_ind(jbe,jbe + 1)] = 0;
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = r[JAC_N02];
  jac[jac_ind(jbe,jbe + 2)] = 0;
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = jac_aa_pm(f_al, f_be, f_ps, r, j, -1);
    jac[jac_ind(jbe,jbe - 1)] = jac_bb_pm(f_al, f_be, f_ps, r, j, -1);
    // ROW j, COL j
    jac[jac_ind(j,j)] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    jac[jac_ind(jbe,jbe)] = jac_bb(f_al, f_be, f_ps, r, j);
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = jac_aa_pm(f_al, f_be, f_ps, r, j, 1);
    jac[jac_ind(jbe,jbe + 1)] = jac_bb_pm(f_al, f_be, f_ps, r, j, 1);
  }
  j = one_past_last;
  jbe = j + npts;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = r[JAC_RRM2];
  jac[jac_ind(jbe,jbe - 2)] = r[JAC_RRM2];
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = r[JAC_RRM1];
  jac[jac_ind(jbe,jbe - 1)] = r[JAC_RRM1];
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = r[JAC_RR];
  jac[jac_ind(jbe,jbe)] = r[JAC_RR];
  return;
}




#endif
