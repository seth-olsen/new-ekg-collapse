#ifndef EKG_CLEAN_H_INCLUDED
#define EKG_CLEAN_H_INCLUDED

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
#include "jacobian.h"
#include "ekg-proc.h"
#include "lapacke.h"

using namespace std;

// perform gauss-seidel update on xi, pi (xpp version for including f_ps)
void hyp_solve_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt);
void hyp_solve_px(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt);
void hyp_solve_px_fast(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		       const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt);
void hyp_solve(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
	       const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt);
// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
// get inf-norm of residuals, residuals vector set
dbl full_res0clean(VD& residuals, const VD& old_xi, const VD& old_pi,
		   const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		   const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		   MAPID& r, int lastpt);
dbl get_hyp_res(VD& residuals, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		MAPID& r, int lastpt);
dbl get_hyp_res_fast(VD& residuals, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		     const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		     MAPID& r, int lastpt);
void get_ell_res_abpclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, MAPID& r, int lastpt);
void get_ell_res_abp_fast(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			  const VD& f_ps, MAPID& r, int lastpt);
void get_ell_res_pbaclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, MAPID& r, int lastpt);
double get_ell_res_ind(VD& res_al, VD& res_be, VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
		       const VD& f_ps, MAPID& r, int lastpt);
inline void get_resAl(VD& res_al, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      MAPID& r, int lastpt);
inline void get_resBe(VD& res_be, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      MAPID& r, int lastpt);
inline void get_resPs(VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      MAPID& r, int lastpt);

inline void apply_up_field(const VD& deltas, VD& field, int npts);
inline void apply_up_abp_join(const VD& deltas, VD& f_al, VD& f_be, VD& f_ps, int npts);
inline void apply_up_abp_join_cn(const VD& deltas, const VD& old_al, const VD& old_be, const VD& old_ps,
				 VD& f_al, VD& f_be, VD& f_ps, VD& cn_al, VD& cn_be, VD& cn_ps, int npts);
int ell_solve_abpclean(VD& jac, VD& res_ell, const VD& old_al, const VD& old_be, const VD& old_ps,
		       const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		       VD& cn_al, VD& cn_be, VD& cn_ps, MAPID& r, int npts,
		       int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb);
int ell_solve_abpclean_join(VD& jac, VD& res_ell, const VD& old_al, const VD& old_be, const VD& old_ps,
			    const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
			    VD& cn_al, VD& cn_be, VD& cn_ps, MAPID& r, int npts,
			    int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb);

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// perform gauss-seidel update (xi and pi only)
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt)
{
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_xi & cn_xi
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
    // update f_pi & cn_pi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_xi, f_xi, r, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  sommerfeld(old_pi, f_pi, r, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_px(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt)
{
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_pi & cn_pi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    // update f_xi & cn_xi
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_pi, f_pi, r, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  sommerfeld(old_xi, f_xi, r, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_px_fast(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		       const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt)
{
  double lam6_part;
  double ps2_m, ps2 = sq(cn_ps[0]), ps2_p = sq(cn_ps[1]);
  double al_ps2_m, al_ps2 = cn_al[0] / ps2, al_ps2_p = cn_al[1] / ps2_p;
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    lam6_part = r[LAM6VAL] * (d_c(cn_be,k) + cn_be[k]*(4*r[DRVAL]*r[-k] + 6*d_c(cn_ps,k)/cn_ps[k]));
    ps2_m = ps2; ps2 = ps2_p; ps2_p = sq(cn_ps[k+1]);
    al_ps2_m = al_ps2; al_ps2 = al_ps2_p; al_ps2_p = cn_al[k+1] / ps2_p;
    // update f_pi & cn_pi
    f_pi[k] = ( old_pi[k]*(1 - lam6_part) + (r[LAM2VAL] / sq(r[k]*ps2)) *
		(sq(r[k+1]*ps2_p)*(al_ps2_p*cn_xi[k+1] + cn_be[k+1]*cn_pi[k+1])
		 - sq(r[k-1]*ps2_m)*(al_ps2_m*cn_xi[k-1] + cn_be[k-1]*cn_pi[k-1])) )
      / (1 + lam6_part);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    // update f_xi & cn_xi
    f_xi[k] = old_xi[k] + r[LAM2VAL]*(al_ps2_p*cn_pi[k+1] + cn_be[k+1]*cn_xi[k+1]
				      - al_ps2_m*cn_pi[k-1] - cn_be[k-1]*cn_xi[k-1]);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_pi, f_pi, r, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  sommerfeld(old_xi, f_xi, r, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
	       const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int lastpt)
{
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_pi & f_xi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    // update cn_pi & cn_xi
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_pi, f_pi, r, lastpt);
  sommerfeld(old_xi, f_xi, r, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
////////////////////////////////////////////////////////////////////////////////////////////////
dbl full_res0clean(VD& residuals, const VD& old_xi, const VD& old_pi,
		   const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		   const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		   MAPID& r, int lastpt)
{
  int kpi = lastpt + 1, kal = 2*kpi, kbe = 3*kpi, kps = 4*kpi;
  residuals[0] = dirichlet0res(f_xi);
  residuals[kpi] = neumann0res(f_pi, r);
  residuals[kal] = neumann0res(f_al, r);
  residuals[kbe] = dirichlet0res(f_be);
  residuals[kps] = neumann0res(f_ps, r);
  for (int k = 1; k < lastpt; ++k) {
    residuals[k] = fda_resXi(old_xi, f_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    residuals[kpi+k] = fda_resPi(old_pi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    residuals[kal+k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    residuals[kbe+k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    residuals[kps+k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, r, k);
  }
  residuals[lastpt] = sommerfeldres(old_xi, f_xi, r, lastpt);
  residuals[kpi] = sommerfeldres(old_pi, f_pi, r, lastpt);
  residuals[kal] = fdaR_resAl(f_al, r, lastpt);
  residuals[kbe] = fdaR_resBe(f_be, r, lastpt);
  residuals[kps] = fdaR_resPs(f_ps, r, lastpt);
  return max(  *max_element(residuals.begin(), residuals.end()),
	     -(*min_element(residuals.begin(), residuals.end()))  );
}

dbl get_hyp_res(VD& res_hyp, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		MAPID& r, int lastpt)
{
  int kpi = lastpt + 1;
  res_hyp[0] = dirichlet0res(f_xi);
  res_hyp[kpi] = neumann0res(f_pi, r);
  for (int k = 1; k < lastpt; ++k) {
    res_hyp[k] = fda_resXi(old_xi, f_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
    res_hyp[kpi + k] = fda_resPi(old_pi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, k);
  }
  res_hyp[lastpt] = sommerfeldres(old_xi, f_xi, r, lastpt);
  res_hyp[kpi + lastpt] = sommerfeldres(old_pi, f_pi, r, lastpt);
  return max(  *max_element(res_hyp.begin(), res_hyp.end()),
	      -(*min_element(res_hyp.begin(), res_hyp.end()))  );
}

dbl get_hyp_res_fast(VD& res_hyp, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		     const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		     MAPID& r, int lastpt)
{
  double lam6_part;
  double ps2_m, ps2 = sq(cn_ps[0]), ps2_p = sq(cn_ps[1]);
  double al_ps2_m, al_ps2 = cn_al[0] / ps2, al_ps2_p = cn_al[1] / ps2_p;
  int kpi = lastpt + 1;
  res_hyp[0] = dirichlet0res(f_xi);
  res_hyp[kpi] = neumann0res(f_pi, r);
  for (int k = 1; k < lastpt; ++k) {
    lam6_part = r[LAM6VAL] * (d_c(cn_be,k) + cn_be[k]*(4*r[DRVAL]*r[-k] + 6*d_c(cn_ps,k)/cn_ps[k]));
    ps2_m = ps2; ps2 = ps2_p; ps2_p = sq(cn_ps[k+1]);
    al_ps2_m = al_ps2; al_ps2 = al_ps2_p; al_ps2_p = cn_al[k+1] / ps2_p;
    
    res_hyp[k] = r[INDT]*( f_xi[k] - old_xi[k] -
			   r[LAM2VAL]*(al_ps2_p*cn_pi[k+1] + cn_be[k+1]*cn_xi[k+1]
				       - al_ps2_m*cn_pi[k-1] - cn_be[k-1]*cn_xi[k-1]) );
    res_hyp[kpi + k] = r[INDT]*( f_pi[k]*(lam6_part + 1) + old_pi[k]*(lam6_part - 1)
				 - (r[LAM2VAL] / sq(r[k]*ps2)) *
				 (sq(r[k+1]*ps2_p)*(al_ps2_p*cn_xi[k+1] + cn_be[k+1]*cn_pi[k+1])
				  - sq(r[k-1]*ps2_m)*(al_ps2_m*cn_xi[k-1] + cn_be[k-1]*cn_pi[k-1])) );
  }
  res_hyp[lastpt] = sommerfeldres(old_xi, f_xi, r, lastpt);
  res_hyp[kpi + lastpt] = sommerfeldres(old_pi, f_pi, r, lastpt);
  return max(  *max_element(res_hyp.begin(), res_hyp.end()),
	      -(*min_element(res_hyp.begin(), res_hyp.end()))  );
}

void get_ell_res_abpclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, MAPID& r, int lastpt)
{
  int kbe = lastpt + 1;
  int kps = 2*kbe;
  res_ell[0] = neumann0res(f_al, r);
  res_ell[kbe] = dirichlet0res(f_be);
  res_ell[kps] = neumann0res(f_ps, r);
  for (int k = 1; k < lastpt; ++k) {
    res_ell[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    res_ell[kbe + k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    res_ell[kps + k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, r, k);
  }
  res_ell[lastpt] = fdaR_resAl(f_al, r, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, r, lastpt);
  res_ell[kps + lastpt] = fdaR_resPs(f_ps, r, lastpt);
  return;
}

void get_ell_res_abp_fast(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			  const VD& f_ps, MAPID& r, int lastpt)
{
  double drbe_r, dps_ps, alin, p4db2_a;
  int kbe = lastpt + 1;
  int kps = 2*kbe;
  res_ell[0] = neumann0res(f_al, r);
  res_ell[kbe] = dirichlet0res(f_be);
  res_ell[kps] = neumann0res(f_ps, r);
  for (int k = 1; k < lastpt; ++k) {
    drbe_r = ddr_c(f_be,r,k) - r[-k]*f_be[k];
    dps_ps = ddr_c(f_ps,r,k) / f_ps[k];
    alin = 1 / f_al[k];
    p4db2_a = sq( sq(f_ps[k]) * drbe_r ) * alin; 
    
    
    res_ell[k] = ddr2_c(f_al,r,k) + 2*ddr_c(f_al,r,k)*(r[-k] + dps_ps)
      - r[TWO_THIRDS]*p4db2_a - 8*M_PI*f_al[k]*sq(f_pi[k]);
    
    res_ell[kbe + k] = ddr2_c(f_be,r,k) + drbe_r*(2*r[-k] + 6*dps_ps - ddr_c(f_al,r,k)*alin)
      + 12*M_PI*f_al[k]*f_xi[k]*f_pi[k]/sq(f_ps[k]);
    
    res_ell[kps + k] = ddr2_c(f_ps,r,k) + f_ps[k]*( 2*r[-k]*dps_ps + r[TWELFTH]*p4db2_a*alin
						    + M_PI*(sq(f_xi[k]) + sq(f_pi[k])) );
  }
  res_ell[lastpt] = fdaR_resAl(f_al, r, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, r, lastpt);
  res_ell[kps + lastpt] = fdaR_resPs(f_ps, r, lastpt);
  return;
}

// ******** NOW THIS JOINS P->B->A
void get_ell_res_pbaclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, MAPID& r, int lastpt)
{
  int kbe = lastpt + 1;
  int kal = 2*kbe;
  res_ell[0] = neumann0res(f_ps, r);
  res_ell[kbe] = dirichlet0res(f_be);
  res_ell[kal] = neumann0res(f_al, r);
  for (int k = 1; k < lastpt; ++k) {
    res_ell[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    res_ell[kbe + k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    res_ell[kal + k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, r, k);
  }
  res_ell[lastpt] = fdaR_resPs(f_ps, r, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, r, lastpt);
  res_ell[kal + lastpt] = fdaR_resAl(f_al, r, lastpt);
  return;
}
double get_ell_res_ind(VD& res_al, VD& res_be, VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
		       const VD& f_ps, MAPID& r, int lastpt)
{
  res_ps[0] = neumann0res(f_ps, r);
  res_be[0] = dirichlet0res(f_be);
  res_al[0] = neumann0res(f_al, r);
  for (int k = 1; k < lastpt; ++k) {
    res_ps[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    res_be[k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, r, k);
    res_al[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, r, k);
  }
  res_ps[lastpt] = fdaR_resPs(f_ps, r, lastpt);
  res_be[lastpt] = fdaR_resBe(f_be, r, lastpt);
  res_al[lastpt] = fdaR_resAl(f_al, r, lastpt);
  return max(norm_inf(res_ps), max(norm_inf(res_be), norm_inf(res_al)));
}
inline void get_resAl(VD& res_al, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      MAPID& r, int lastpt)
{
  res_al[0] = neumann0res(f_al, r);
  for (int k = 0; k < lastpt; ++k) { res_al[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, r, k); }
  res_al[lastpt] = fdaR_resAl(f_al, r, lastpt);
  return;
}
inline void get_resBe(VD& res_be, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      MAPID& r, int lastpt)
{
  res_be[0] = dirichlet0res(f_be);
  for (int k = 0; k < lastpt; ++k) { res_be[k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, r, k); }
  res_be[lastpt] = fdaR_resBe(f_al, r, lastpt);
  return;
}
inline void get_resPs(VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		     MAPID& r, int lastpt)
{
  res_ps[0] = neumann0res(f_ps, r);
  for (int k = 0; k < lastpt; ++k) { res_ps[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, r, k); }
  res_ps[lastpt] = fdaR_resPs(f_ps, r, lastpt);
  return;
}

inline void apply_up_field(const VD& deltas, VD& field, int npts)
{
  for (int k = 0; k < npts; ++k) { field[k] -= deltas[k]; }
  return;
}

inline void apply_up_abp_join(const VD& deltas, VD& f_al, VD& f_be, VD& f_ps, int npts)
{
  int kbe = npts;
  int kps = 2*npts;
  for (int k = 0; k < npts; ++k) {
    f_al[k] -= deltas[k];
    f_be[k] -= deltas[kbe + k];
    f_ps[k] -= deltas[kps + k];
  }
  return;
}

inline void apply_up_abp_join_cn(const VD& deltas, const VD& old_al, const VD& old_be, const VD& old_ps,
				 VD& f_al, VD& f_be, VD& f_ps, VD& cn_al, VD& cn_be, VD& cn_ps, int npts)
{
  int kbe = npts;
  int kps = 2*npts;
  for (int k = 0; k < npts; ++k) {
    f_al[k] -= deltas[k];
    f_be[k] -= deltas[kbe];
    f_ps[k] -= deltas[kps];
    cn_al[k] = 0.5 * (old_al[k] + f_al[k]);
    cn_be[k] = 0.5 * (old_be[k] + f_be[k]);
    cn_ps[k] = 0.5 * (old_ps[k] + f_ps[k]);
    ++kbe; ++kps;
  }
  return;
}

int ell_solve_abpclean(VD& jac, VD& res_ell, const VD& old_al, const VD& old_be, const VD& old_ps,
		       const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		       VD& cn_al, VD& cn_be, VD& cn_ps, MAPID& r, int npts,
		       int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  get_ell_res_abpclean_join(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, npts-1);
  set_jacCMabpclean(jac, f_xi, f_pi, f_al, f_be, f_ps, r, npts, kl, ku, ldab);
  int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs, &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
  apply_up_abp_join_cn(res_ell, old_al, old_be, old_ps, f_al, f_be, f_ps, cn_al, cn_be, cn_ps, npts);
  return info;
}

int ell_solve_abpclean_join(VD& jac, VD& res_ell, const VD& old_al, const VD& old_be, const VD& old_ps,
			    const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
			    VD& cn_al, VD& cn_be, VD& cn_ps, MAPID& r, int npts,
			    int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  get_ell_res_abpclean_join(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, npts-1);
  set_jacCMabpclean(jac, f_xi, f_pi, f_al, f_be, f_ps, r, npts, kl, ku, ldab);
  int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs, &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
  apply_up_abp_join_cn(res_ell, old_al, old_be, old_ps, f_al, f_be, f_ps, cn_al, cn_be, cn_ps, npts);
  return info;
}

int ell_solve_indep(VD& jac, VD& res_ell, const VD& old_al, const VD& old_be, const VD& old_ps,
		    const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		    VD& cn_al, VD& cn_be, VD& cn_ps, MAPID& r, int npts,
		    int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  return 0;
}

#endif