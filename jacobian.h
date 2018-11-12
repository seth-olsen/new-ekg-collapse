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

//  for LAPACKE_dgbsv(): jac[ (kl + ku + 1) + (ldab - 1)*j + i ]  =  jac[ i, j ]
void set_jacCMabpslow(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab);
void set_jacCMabpfast(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab);
void set_jacCM_ab_slow(VD& jac, const VD& f_xi, const VD& f_pi,
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
  return r[NEG2INDRSQ] - r[EIGHT_PI]*sq(f_pi[k]) +
    r[TWO_THIRDS]*pw4(f_ps[k])*sq(ddr_c(f_be,r,k) - f_be[k]*r[-k]) / sq(f_al[k]);
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
    r[FIVE_TWELFTHS]*pw4(f_ps[k])*sq(ddr_c(f_be,r,k) - f_be[k]*r[-k]) / sq(f_al[k]);
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

void set_jacCMabpfast(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int j = 0, jbe = npts, jps = 2*npts;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = r[JAC_N00];
  jac[jac_ind(jbe,jbe)] = 1;
  jac[jac_ind(jps,jps)] = r[JAC_N00];
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = r[JAC_N01];
  jac[jac_ind(jbe,jbe + 1)] = 0;
  jac[jac_ind(jps,jps + 1)] = r[JAC_N01];
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = r[JAC_N02];
  jac[jac_ind(jbe,jbe + 2)] = 0;
  jac[jac_ind(jps,jps + 2)] = r[JAC_N02];

  double drin_p, drin_m, dps_ps, dlogp6_a, p4db2_a2;
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    jps = jbe + npts;

    drin_p = r[INDRSQ] + r[-j]*r[INDR];
    drin_m = r[INDRSQ] - r[-j]*r[INDR];
    dps_ps = ddr_c(f_ps,r,j) / f_ps[j];
    dlogp6_a = 6*dps_ps - (ddr_c(f_al,r,j)/f_al[j]);
    p4db2_a2 = sq(f_ps[j]) * (ddr_c(f_be,r,j) - r[-j]*f_be[j]) / f_al[j];
    p4db2_a2 = sq(p4db2_a2);    
    
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = drin_m - r[-j]*dps_ps;
    jac[jac_ind(jbe,jbe - 1)] = drin_m - r[IN2DR]*dlogp6_a;
    jac[jac_ind(jps,jps - 1)] = drin_m;
    // ROW j, COL j
    jac[jac_ind(j,j)] = r[NEG2INDRSQ] + r[TWO_THIRDS]*p4db2_a2 - r[EIGHT_PI]*sq(f_pi[j]);
    jac[jac_ind(jbe,jbe)] = r[NEG2INDRSQ] - r[-j]*(2*r[-j] + dlogp6_a);
    jac[jac_ind(jps,jps)] = r[NEG2INDRSQ] + r[FIVE_TWELFTHS]*p4db2_a2 + M_PI*(sq(f_xi[j]) + sq(f_pi[j]));
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = drin_p + r[-j]*dps_ps;
    jac[jac_ind(jbe,jbe + 1)] = drin_p + r[IN2DR]*dlogp6_a;
    jac[jac_ind(jps,jps + 1)] = drin_p;
  }
  
  j = one_past_last;
  jbe = j + npts;
  jps = jbe + npts;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = r[JAC_RRM2];
  jac[jac_ind(jbe,jbe - 2)] = r[JAC_RRM2];
  jac[jac_ind(jps,jps - 2)] = r[JAC_RRM2];
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = r[JAC_RRM1];
  jac[jac_ind(jbe,jbe - 1)] = r[JAC_RRM1];
  jac[jac_ind(jps,jps - 1)] = r[JAC_RRM1];
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = r[JAC_RR];
  jac[jac_ind(jbe,jbe)] = r[JAC_RR];
  jac[jac_ind(jps,jps)] = r[JAC_RR];
  return;
}

void set_jacCMabpslow(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, MAPID& r,
		      int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int j = 0, jbe = npts, jps = 2*npts;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = r[JAC_N00];
  jac[jac_ind(jbe,jbe)] = 1;
  jac[jac_ind(jps,jps)] = r[JAC_N00];
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = r[JAC_N01];
  jac[jac_ind(jbe,jbe + 1)] = 0;
  jac[jac_ind(jps,jps + 1)] = r[JAC_N01];
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = r[JAC_N02];
  jac[jac_ind(jbe,jbe + 2)] = 0;
  jac[jac_ind(jps,jps + 2)] = r[JAC_N02];
  
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    jps = jbe + npts;
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = jac_aa_pm(f_al, f_be, f_ps, r, j, -1);
    jac[jac_ind(jbe,jbe - 1)] = jac_bb_pm(f_al, f_be, f_ps, r, j, -1);
    jac[jac_ind(jps,jps - 1)] = jac_pp_pm(f_al, f_be, f_ps, r, j, -1);
    // ROW j, COL j
    jac[jac_ind(j,j)] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    jac[jac_ind(jbe,jbe)] = jac_bb(f_al, f_be, f_ps, r, j);
    jac[jac_ind(jps,jps)] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, r, j);
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = jac_aa_pm(f_al, f_be, f_ps, r, j, 1);
    jac[jac_ind(jbe,jbe + 1)] = jac_bb_pm(f_al, f_be, f_ps, r, j, 1);
    jac[jac_ind(jps,jps + 1)] = jac_pp_pm(f_al, f_be, f_ps, r, j, 1);
  }
  j = one_past_last;
  jbe = j + npts;
  jps = jbe + npts;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = r[JAC_RRM2];
  jac[jac_ind(jbe,jbe - 2)] = r[JAC_RRM2];
  jac[jac_ind(jps,jps - 2)] = r[JAC_RRM2];
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = r[JAC_RRM1];
  jac[jac_ind(jbe,jbe - 1)] = r[JAC_RRM1];
  jac[jac_ind(jps,jps - 1)] = r[JAC_RRM1];
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = r[JAC_RR];
  jac[jac_ind(jbe,jbe)] = r[JAC_RR];
  jac[jac_ind(jps,jps)] = r[JAC_RR];
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
