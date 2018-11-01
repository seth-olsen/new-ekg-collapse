#ifndef SOLVERS_H_INCLUDED
#define SOLVERS_H_INCLUDED

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
#include "ekg-clean.h"

using namespace std;

typedef int (*SOLVER)(VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , const VD& , MAPID& ,
		      int , int , int , dbl , dbl , dbl ,
		      int , int , int , int , int , vector<int>& , int);

int solve_t0_fast(const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		  MAPID& r, int lastpt, int ell_maxit, dbl ell_tol)
{
  string error_response = "idc";
  lapack_int N_0 = 3*(lastpt + 1);
  lapack_int kl = 2;
  lapack_int ku = 2;
  lapack_int nrhs = 1;
  lapack_int ldab = 2*kl + ku + 1;
  lapack_int ldb_0 = N_0;
  vector<lapack_int> ipiv_0(N_0);
  VD res_0(ldb_0, 0);
  lapack_int info = 0;

  get_ell_res_abp_fast(res_0, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
  dbl res = norm_inf(res_0);
  int ell_itn = 0;
  while (res > ell_tol) {
    VD jac(ldab*N_0, 0);
    set_jacCMabpfast(jac, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt + 1, kl, ku, ldab);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N_0, kl, ku, nrhs,
			 &jac[0], ldab, &ipiv_0[0], &res_0[0], ldb_0);
    if (info != 0) { cout << ell_itn << "\nERROR: cannot solve initial elliptics\ninfo = " << info << endl; }
    
    apply_up_abp_join(res_0, f_al, f_be, f_ps, lastpt + 1);
    get_ell_res_abp_fast(res_0, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
    res = norm_inf(res_0);
    if (++ell_itn > ell_maxit) {
      cout << "\nSTUCK AT t=0 with res = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }
  }
  return ell_itn;
}

int solve_t0_slow(const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		  MAPID& r, int lastpt, int ell_maxit, dbl ell_tol)
{
  string error_response = "idc";
  lapack_int N_0 = 3*(lastpt + 1);
  lapack_int kl = 2;
  lapack_int ku = 2;
  lapack_int nrhs = 1;
  lapack_int ldab = 2*kl + ku + 1;
  lapack_int ldb_0 = N_0;
  vector<lapack_int> ipiv_0(N_0);
  VD res_0(ldb_0, 0);
  lapack_int info = 0;

  get_ell_res_abpclean_join(res_0, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
  dbl res = norm_inf(res_0);
  int ell_itn = 0;
  while (res > ell_tol) {
    VD jac(ldab*N_0, 0);
    set_jacCMabpslow(jac, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt + 1, kl, ku, ldab);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N_0, kl, ku, nrhs,
			 &jac[0], ldab, &ipiv_0[0], &res_0[0], ldb_0);
    if (info != 0) { cout << ell_itn << "\nERROR: cannot solve initial elliptics\ninfo = " << info << endl; }
    
    apply_up_abp_join(res_0, f_al, f_be, f_ps, lastpt + 1);
    get_ell_res_abpclean_join(res_0, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
    res = norm_inf(res_0);
    if (++ell_itn > ell_maxit) {
      cout << "\nSTUCK AT t=0 with res = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }
  }
  return ell_itn;
}


int solve_E_fast(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		 VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		 VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		 VD& res_hyp, VD& res_ell, const VD& jac_zero, MAPID& r,
		 int lastpt, int maxit, int i, dbl t, dbl tol, dbl ell_tol,
		 int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = tol + 1;
  while (res > tol) {
    hyp_itn = 0; ell_itn = 0;
    while (res > tol) {
      hyp_solve_px_fast(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			r, lastpt);
      res = get_hyp_res_fast(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			     r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << t << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    
    get_ell_res_abp_fast(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
    res = norm_inf(res_ell);
    while (res > ell_tol) {
      jac = jac_zero;
      set_jacCMabpfast(jac, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt + 1, kl, ku, ldab);
      //set_jacCMabpclean(jac, f_xi, f_pi, f_al, f_be, f_ps, r, npts, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << t << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_abp_join(res_ell, f_al, f_be, f_ps, lastpt + 1);
      get_ell_res_abp_fast(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << t << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set3_cn(old_ps, old_be, old_al, f_ps, f_be, f_al, cn_ps, cn_be, cn_al, lastpt + 1);
    res = get_hyp_res_fast(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			   r, lastpt);
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << t << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

int solve_E_slow(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		 VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		 VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		 VD& res_hyp, VD& res_ell, const VD& jac_zero, MAPID& r,
		 int lastpt, int maxit, int i, dbl t, dbl tol, dbl ell_tol,
		 int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = tol + 1;
  while (res > tol) {
    hyp_itn = 0; ell_itn = 0;
    while (res > tol) {
      hyp_solve_px(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
		   r, lastpt);
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << t << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    
    get_ell_res_abpclean_join(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
    res = norm_inf(res_ell);
    while (res > ell_tol) {
      jac = jac_zero;
      set_jacCMabpslow(jac, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt + 1, kl, ku, ldab);
      //set_jacCMabpclean(jac, f_xi, f_pi, f_al, f_be, f_ps, r, npts, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << t << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_abp_join(res_ell, f_al, f_be, f_ps, lastpt + 1);
      get_ell_res_abpclean_join(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << t << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set3_cn(old_ps, old_be, old_al, f_ps, f_be, f_al, cn_ps, cn_be, cn_al, lastpt + 1);
    res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
		      r, lastpt);
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << t << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

int solve_H_slow(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		 VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		 VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		 VD& res_hyp, VD& res_ell, const VD& jac_zero, MAPID& r,
		 int lastpt, int maxit, int i, dbl t, dbl tol, dbl ell_tol,
		 int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = tol + 1;
  while (res > tol) {
    hyp_itn = 0; ell_itn = 0;
    while (res > tol) {
      hyp_solve_ppx_slow(old_xi, old_pi, old_ps, f_xi, f_pi, f_ps,
			 cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, lastpt);
      res = get_hyp_res_ppx(res_hyp, old_xi, old_pi, old_ps, f_xi, f_pi, f_ps,
			    cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << t << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    
    get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
    res = norm_inf(res_ell);
    while (res > ell_tol) {
      jac = jac_zero;
      set_jacCM_ab_slow(jac, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt + 1, kl, ku, ldab);
      //set_jacCMabpclean(jac, f_xi, f_pi, f_al, f_be, f_ps, r, npts, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << t << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_ab(res_ell, f_al, f_be, lastpt + 1);
      get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << t << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set2_cn(old_be, old_al, f_be, f_al, cn_be, cn_al, lastpt + 1);
    res = get_hyp_res_ppx(res_hyp, old_xi, old_pi, old_ps, f_xi, f_pi, f_ps,
			  cn_xi, cn_pi, cn_al, cn_be, cn_ps, r, lastpt);
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << t << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

#endif
