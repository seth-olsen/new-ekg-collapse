#ifndef OLD_JACOBIAN_H_INCLUDED
#define OLD_JACOBIAN_H_INCLUDED

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
#include "jacobian.h"

using namespace std;

inline void set_jac_alphaCM(VD& jac, const VD& f_xi, const VD& f_pi,
			    const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			    int npts, int kl, int ku, int ldab);
inline void set_jac_betaCM(VD& jac, const VD& f_xi, const VD& f_pi,
			   const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			   int npts, int kl, int ku, int ldab);
inline void set_jac_psiCM(VD& jac, const VD& f_xi, const VD& f_pi,
			  const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
			  int npts, int kl, int ku, int ldab);
void set_jacCMabpclean(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		       int npts, int kl, int ku, int ldab);

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

#endif
