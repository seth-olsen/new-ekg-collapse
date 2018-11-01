#ifndef NEG2INDRSQ
#define INDRSQ 32500
#define IN2DR 32501
#define INDR 32502
#define INDT 32503
#define RMIN 0
#define RMAX 32505
#define DRVAL 32506
#define DTVAL 32507
#define LAMVAL 32508
#define LAM2VAL 32509
#define LAM6VAL 32510
#define CSOMM 32511
#define INRMAX 32512
#define NEG2INDRSQ 32513
#define WRITEDR 32514
#define TWO_THIRDS 32515
#define TWELFTH 32516
#define FOUR_THIRDS 32517
#define FIVE_TWELFTHS 32518
#define EIGHT_PI 32519
#define TWELVE_PI 32520
#define JAC_RR 32521
#define JAC_RRM1 32522
#define JAC_RRM2 32523
#define JAC_N00 32524
#define JAC_N01 32525
#define JAC_N02 32526
#define DT_TWELVE 32527
#define CPSI_RHS 32528
#define CSOMM_RHS 32529
#define CSOMM_OLD 32530
#endif

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <cmath> // for ICs
#include <vector> // for everything
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "lapacke.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "jacobian.h"
#include "ekg-clean.h"
#include "ekg-proc.h"
#include "solvers.h"

int main(int argc, char **argv)
{
  // **********************************************************
  // ******************** PARAMETERS **************************
  // **********************************************************
  
  // user-set parameters
  str outfile = "ekg";
  int lastpt = 1000; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 2000; // time steps
  int save_step = 4; // write only every (save_step)th time step
  dbl lam = 0.25; // dt/dr
  dbl r2m = 0.0;
  dbl rmin = 0.0;
  dbl rmax = 100.0;
  dbl dspn = 0.5; // dissipation coefficient
  dbl tol = 0.00000001; // iterative method tolerance
  dbl ell_tol = 0.001*tol;
  // for computing residual: 0 = inf-norm, 1 = 1-norm, 2 = 2-norm
  int resnorm_type = 0;
  int maxit = 25; // max iterations for debugging
  dbl ic_Dsq = 4.0; // gaussian width
  dbl ic_r0 = 25.0; // gaussian center
  dbl ic_Amp = 0.004; // gaussian amplitude
  int check_step = 10; // for monitoring invariant mass
  // note: set bools in command line with integers 1=true or 0=false
  bool psi_hyp = false; // update psi with hyperbolic evolution eqn after IC?
  bool zero_pi = false; // zero initial time derivative?
  bool somm_cond = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool dspn_psi = false; // dissipate psi (only activated if psi_hyp=true)?
  bool dr3_up = false; // update pi with d/dr^3 scheme?
  bool static_metric = false; // ignore scalar field's effect on metric?
  bool clean_hyp = false; // use clean hyperbolic update functions (slower)?
  bool clean_ell = false; // use clean hyperbolic update functions (slower)?
  bool write_res = false; // write residuals?
  bool write_ricci = false; // write ricci?
  bool write_itn = false; // write itn counts?
  bool write_mtot = false; // write total mass?
  bool write_maspect = true; // write mass aspect?
  bool write_outnull = true; // write outgoing null expansion?
  bool write_xp = true; // write xi and pi?
  bool write_abp = true; // write metric fields (alpha, beta, psi)?
  bool write_ires_xp = false; // write ires for xi and pi?
  bool write_ires_abp = false; // write ires for metric variables?
  // variable to hold constant across resolutions
  str hold_const = "lambda"; // "lambda", "dt", or "dr"
  int nresn = 1; // 1, 2, or 3
  int resn0 = 1, resn1 = 2, resn2 = 4; // in order of priority
  int *resns[3] = {&resn0, &resn1, &resn2};
  
  map<str, str *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<str, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt",&save_pt},
      {"-nsteps",&nsteps}, {"-save_step",&save_step},
      {"-resnorm_type",&resnorm_type}, {"-maxit",&maxit},
      {"-check_step",&check_step}, {"-nresn",&nresn},
      {"-resn0",&resn0}, {"-resn1",&resn1}, {"-resn2",&resn2}};
  map<str, dbl *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ell_tol",&ell_tol},
      {"-ic_Dsq",&ic_Dsq}, {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<str, bool *> p_bool { {"-psi_hyp",&psi_hyp}, {"-zero_pi",&zero_pi},
      {"-somm_cond",&somm_cond}, {"-dspn_bound",&dspn_bound}, {"-dr3_up",&dr3_up}, {"-dspn_psi",&dspn_psi},
      {"-write_res",&write_res},{"-write_ricci",&write_ricci}, {"-write_itn",&write_itn},
      {"-write_mtot",&write_mtot},{"-write_maspect",&write_maspect}, {"-write_outnull",&write_outnull},
      {"-write_xp",&write_xp}, {"-write_abp",&write_abp},
      {"-write_ires_xp",&write_ires_xp}, {"-write_ires_abp",&write_ires_abp},
      {"-clean_hyp",&clean_hyp}, {"-clean_ell",&clean_ell}};
  map<str, str> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);

  // OBTAIN RESOLUTION FACTORS
  vector<int> resolutions(nresn);
  for (int k = 0; k < nresn; ++k) {
    resolutions[k] = *resns[k];
  }

  // check that grid size (lastpt = npts-1) is divisible by save_pt 
  if (lastpt % save_pt != 0) {
    cout << "ERROR: save_pt = " << save_pt << " entered for grid size " << lastpt << endl;
    save_pt -= lastpt % save_pt;
    cout << "--> corrected: save_pt = " << save_pt << endl;
  }
  // check that resnorm_type is valid
  if (resnorm_type < 0 || resnorm_type > 2) {
    cout << "ERROR: resnorm_type=" << save_pt << " not valid -> using inf-norm" << endl;
    resnorm_type = 0;
  }

  // SET SOLVER
  SOLVER ekg_solver = solve_E_fast;
  if (psi_hyp) { ekg_solver = solve_H_slow; }
  else if (clean_hyp && clean_ell) { ekg_solver = solve_E_slow; }

  // HERE psi_hyp = True means do NOT dissipate the point next to r = 0
  void (*apply_dissipation)(const VD&, const VD&, VD&, VD&, int, dbl);
  if (dspn_bound) { apply_dissipation = dissipationB_xp; }
  else { apply_dissipation = dissipationNB_xp; }

  int n_ell = 3;
  if (psi_hyp) { n_ell = 2; }
  int n_hyp = 5 - n_ell;

// **********************************************************
//                      WRITE START
// **********************************************************

  // bbhutil parameters for writing data to sdf
  int lastwr = lastpt/save_pt;
  int wr_shape = lastwr + 1;
  int *bbh_shape = &wr_shape;
  int bbh_rank = 1;
  dbl coord_lims[2] = {rmin, rmax};
  dbl *coords = &coord_lims[0];
  dbl wr_dr = (rmax - rmin) / ((dbl) lastwr);

  int vlen = ((write_xp) ? wr_shape : 1);
  VD wr_xi(vlen, 0.0), wr_pi(vlen, 0.0);
  vlen = ((write_abp) ? wr_shape : 1);
  VD wr_al(vlen, 0.0), wr_be(vlen, 0.0), wr_ps(vlen, 0.0);

  VD ricci(((write_ricci) ? wr_shape : 1), 0.0);
  VD maspect(((write_maspect) ? wr_shape : 1), 0.0);
  VD outnull(((write_outnull) ? wr_shape : 1), 0.0);

  vlen = ((write_ires_xp) ? wr_shape : 1);
  VD ires_xi(vlen, 0.0), ires_pi(vlen, 0.0);
  vlen = ((write_ires_abp) ? wr_shape : 1);
  VD ires_al(vlen, 0.0), ires_be(vlen, 0.0), ires_ps(vlen, 0.0);

  // SET WRITE FUNCTION & POPULATE WRITE POINTERS
  vector<dbl *> wr_ptrs({&wr_xi[0], &wr_pi[0],
	&wr_al[0], &wr_be[0], &wr_ps[0], &ricci[0], &maspect[0], &outnull[0],
	&ires_xi[0], &ires_pi[0], &ires_al[0], &ires_be[0], &ires_ps[0]});
  WR_FN wr_fn = get_wr_fn(wr_ptrs, write_xp, write_abp, write_ricci, write_maspect,
			  write_outnull, write_ires_xp, write_ires_abp);  
  int num_wr = wr_ptrs.size();

  // SAME FOR WRITE RESIDUAL WRITING
  vector<dbl *> wr_res_ptrs({&wr_xi[0], &wr_pi[0],
	&wr_al[0], &wr_be[0], &wr_ps[0]});
  WR_RES_FN wr_res_fn = get_wr_res_fn(wr_res_ptrs, write_res,
				      write_xp, write_abp);
  int num_wr_res = wr_res_ptrs.size();
// **********************************************************
//                       WRITE END
// **********************************************************
  
// **************************************************************
// **************************************************************
//                 LOOP PROGRAM OVER RESOLUTIONS
// **************************************************************
// **************************************************************

  str outfile0 = outfile;
  int lastpt0 = lastpt; 
  int save_pt0 = save_pt;
  int nsteps0 = nsteps;
  int save_step0 = save_step;
  dbl lam0 = lam;
  // will be used repeatedly:
  dbl lam2 = 0.5 * lam;
  dbl one_third = 1 / 3.0;
  dbl lam6 = one_third * lam2;
  dbl eight_M_PI = 8 * M_PI;
  
  for (int factor : resolutions) {
    if (hold_const == "lambda") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      outfile = to_string(factor) + "-" + outfile0;
    }
    else if (hold_const == "dt") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      lam = lam0 * factor;
      outfile = to_string(factor) + "dr-" + outfile0;
    }
    else if (hold_const == "dr") {
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      lam = lam0 / ((dbl) factor);
      outfile = to_string(factor) + "dt-" + outfile0;
    }
    else { cout << "ERROR: hold_const must be 'lambda' or 'dt' or 'dr'" << endl; }

// **********************************************************
// **********************************************************
//                      PARAMETER MAP
// **********************************************************
// **********************************************************
    
    // derived parameters
    int npts = lastpt + 1;
    dbl norm_factor = 1 / ((dbl) 5*npts);
    if (resnorm_type == 1) { norm_factor = sqrt(norm_factor); }
    dbl dr = (rmax - rmin) / ((dbl) lastpt);
    dbl dt = lam * dr;
    MAPID r { {DRVAL,dr}, {DTVAL,dt}, {RMIN,rmin}, {RMAX,rmax},
	      {LAMVAL,lam}, {LAM2VAL,lam2}, {LAM6VAL,lam6}, {WRITEDR,wr_dr} };

    r[INDR] = 1 / dr;
    r[IN2DR] = 0.5*r[INDR];
    r[INDRSQ] = sq(r[INDR]);
    r[INDT] = 1 / dt;
    r[INRMAX] = 1 / rmax;
    r[NEG2INDRSQ] = -2 * r[INDRSQ];
    r[CSOMM] = 0.75*lam + 0.5*dt*r[INRMAX]; // for field's outer bc
    r[TWO_THIRDS] = 2 * one_third;
    r[FOUR_THIRDS] = 2 * r[TWO_THIRDS];
    r[TWELFTH] = 0.25 * one_third;
    r[FIVE_TWELFTHS] = 5 * r[TWELFTH];
    r[EIGHT_PI] = eight_M_PI;
    r[TWELVE_PI] = 1.5 * r[EIGHT_PI];
    r[JAC_RR] = (3 * r[IN2DR]) + r[INRMAX];
    r[JAC_RRM1] = -4 * r[IN2DR];
    r[JAC_RRM2] = r[IN2DR];
    r[JAC_N00] = -3 * r[IN2DR];
    r[JAC_N01] = 4 * r[IN2DR];
    r[JAC_N02] = -1 * r[IN2DR];
    r[DT_TWELVE] = r[TWELFTH] * dt;
    r[CPSI_RHS] = 1 / r[JAC_RR];
    r[CSOMM_RHS] = 1 / (1 + r[CSOMM]);
    r[CSOMM_OLD] = 1 - r[CSOMM];
    
// **********************************************************
//                      WRITE START
// **********************************************************
    
    // OUTPUT parameter data
    param_print(outfile,lastpt,save_pt,nsteps,save_step,lam,r2m,rmin,rmax,
		dspn,tol,maxit,ic_Dsq,ic_r0,ic_Amp,check_step,dr,dt,
		psi_hyp,zero_pi,somm_cond,dspn_bound,clean_hyp,clean_ell);
    // FILE NAMES for FIELDS, MASPECT, METRIC VARIABLES, IRES
    vector<str> filenames = get_filenames(outfile, write_xp, write_abp, write_ricci, write_maspect,
					  write_outnull, write_ires_xp, write_ires_abp);
    vector<char *> wr_names(num_wr);
    for (int k = 0; k < num_wr; ++k) { wr_names[k] = &filenames[k][0]; }
    // FILE NAMES for RESIDUALS
    vector<str> resfilenames = get_resfilenames(outfile, write_res, write_xp, write_abp);
    vector<char *> wr_res_names(num_wr_res);
    for (int k = 0; k < num_wr_res; ++k) { wr_res_names[k] = &resfilenames[k][0]; }
    // FILE NAMES for ITN, MTOT
    str itn_file = "itns-" + outfile + ".csv";
    ofstream ofs_itn;
    if (write_itn) {
      ofs_itn.open(itn_file, ofstream::out);
      ofs_itn << "step,time,hyp_itn,ell_itn" << endl;
    }
    str mass_file = "mass-" + outfile + ".csv";
    ofstream ofs_mass;
    if (write_mtot) {
      ofs_mass.open(mass_file, ofstream::out);
      ofs_mass << "save=" << save_step <<","<< "check=" << check_step << endl;
    }
// **********************************************************
//                       WRITE END
// **********************************************************


    // **********************************************************
    // ***************** OBJECT DECLARATIONS ********************
    // **********************************************************
    
    // fields and residuals
    VD old_xi(npts, 0), old_pi(npts, 0);
    VD old_al(npts, 1), old_be(npts, 0), old_ps(npts, 1);
    VD f_xi(npts, 0), f_pi(npts, 0);
    // initial const (flat space) guess for metric vars
    VD f_al(npts, 1), f_be(npts, 0), f_ps(npts, 1);
    VD cn_xi(npts, 0), cn_pi(npts, 0);
    VD cn_al(npts, 1), cn_be(npts, 0), cn_ps(npts, 1);
    VD res_hyp(n_hyp*npts, 0);
    vlen = ((write_ires_xp) ? npts : 1);
    VD older_xi(vlen, 0), older_pi(vlen, 0);
    vlen = ((write_ires_abp) ? npts : 1);
    VD older_ps(vlen, 1);
    vlen = ((write_res) ? 5*npts : 1);
    VD residuals(vlen, 0);
  
    // lapack object declaration
    lapack_int N = n_ell*npts;
    lapack_int kl = 2;
    lapack_int ku = 2;
    lapack_int nrhs = 1;
    lapack_int ldab = 2*kl + ku + 1;
    lapack_int ldb = N;
    vector<lapack_int> ipiv(N);
    VD res_ell(ldb, 0);
    VD jac_zero(ldab*N, 0);

    time_t start_time = time(NULL); // time for rough performance measure
    
// **********************************************************
// ******************* INITIAL DATA ************************
// **********************************************************

    dbl t = 0; // declare position and time variables
    int j, itn = 0;
    for (j = 1; j < npts; ++j) {
      r[j] = rmin + j*dr;
      r[-j] = 1 / r[j];
      f_xi[j] = ic_xi(r[j], ic_Amp, ic_Dsq, ic_r0);
      if (!zero_pi) { f_pi[j] = ic_pi(r[j], ic_Amp, ic_Dsq, ic_r0); }
    }
    if (r[lastpt] != rmax) { cout << "\n**\n****ERROR: r[J] != rmax****\n**\n" << endl; }
    // ASSUMING rmin = 0
    dirichlet0(f_xi);
    neumann0(f_pi);

    // SOLVE ELLIPTIC EQUATIONS FOR t=0  
    if (!static_metric) {
      itn = solve_t0_fast(f_xi, f_pi, f_al, f_be, f_ps, r, lastpt, maxit, 0.1*ell_tol);
    }
    
    // set old_f = f = cn_f for writing initial step
    old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
    cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;
    if (write_ires_xp) {
      older_xi = old_xi;
      older_pi = old_pi;
    }
    if (write_ires_abp) { older_ps = old_ps; }
// **********************************************************
// ******************* TIME STEPPING ************************
// *******************   & WRITING   ************************
// **********************************************************
  
    gft_set_multi(); // start bbhutil file i/o
    for (int i = 0; i < nsteps; ++i) {
      t = i * dt;
// **********************************************************
//                      WRITE START
// **********************************************************
      if (i % save_step == 0) {
	// get coarsened arrays
	wr_fn(older_xi, older_pi, older_ps, old_xi, old_pi, old_al, old_be, old_ps,
	      f_xi, f_pi, f_al, f_be, f_ps, maspect, ires_xi,
	      ires_pi, ires_al, ires_be, ires_ps, wr_xi, wr_pi, wr_al, wr_be, wr_ps,
	      r, save_pt, lastwr);
	if (write_ricci) {
	  ricci[0] = sRicci(f_xi, f_pi, f_ps, 0);//outgoing_null0(f_al, f_be, f_ps, inv2dr, r);
	  int s = save_pt;
	  for (int k = 1; k < lastwr; ++k) {
	    ricci[k] = sRicci(f_xi, f_pi, f_ps, s); // outgoing_null(f_al, f_be, f_ps, s, inv2dr, r);
	    s += save_pt;
	  }
	  ricci[lastwr] = sRicci(f_xi, f_pi, f_ps, lastpt);//outgoing_nullR(f_al, f_be, f_ps, s, inv2dr, r);
	}
	if (write_outnull) {
	  outnull[0] = outgoing_null_f(f_al, f_be, f_ps, r, 0);
	  int s = save_pt;
	  for (int k = 1; k < lastwr; ++k) {
	    outnull[k] = outgoing_null(f_al, f_be, f_ps, r, s);
	    s += save_pt;
	  }
	  outnull[lastwr] = outgoing_null_b(f_al, f_be, f_ps, r, lastpt);
	}
	  
	// write field, maspect, ires outfiles
	wr_step(num_wr, wr_names, t, bbh_shape, bbh_rank, coords, wr_ptrs);
	// write residuals
	if (write_res) {
	  wr_res_fn(residuals, wr_xi, wr_pi, wr_al, wr_be, wr_ps, save_pt, wr_shape);
	  wr_step(num_wr_res, wr_res_names, t, bbh_shape, bbh_rank, coords, wr_res_ptrs);
	}
      }
// **********************************************************
//                       WRITE END
// **********************************************************

// ******************************************************************
// ******************************************************************
//         SOLVE HYPERBOLIC & ELLIPTIC EQUATIONS ITERATIVELY
// ******************************************************************
// ******************************************************************
      // save time n-1 for ires with different d/dt scheme
      if (write_ires_xp) {
	older_xi = old_xi;
	older_pi = old_pi;
      }
      if (write_ires_abp) { older_ps = old_ps; }

      itn = (*ekg_solver)(old_xi, old_pi, old_al, old_be, old_ps, f_xi, f_pi, f_al, f_be, f_ps,
			  cn_xi, cn_pi, cn_al, cn_be, cn_ps, res_hyp, res_ell, jac_zero,
			  r, lastpt, maxit, i, t, tol, ell_tol, N, kl, ku, nrhs, ldab, ipiv, ldb);

      // record itn count of this sweep
      if (write_itn) { ofs_itn << i <<","<< t <<","<< itn << endl; }
// *****************************************************************
// ****************** ITERATIVE SOLUTION COMPLETE ******************
// *****************************************************************
    
    // **************************************************************************
    // *********************** kreiss-oliger DISSIPATION ************************
    // **************************************************************************
      (*apply_dissipation)(old_xi, old_pi, f_xi, f_pi, lastpt-1, dspn);
    }
    
// ***********************************************************************
// ***********************************************************************
// *********************** FULL SOLUTION COMPLETE ************************
// ***********************************************************************
// ***********************************************************************

    // ******************** DONE TIME STEPPING *********************
    // WRITE final time step
    if (nsteps % save_step == 0) {
      wr_fn(older_xi, older_pi, older_ps, old_xi, old_pi, old_al, old_be, old_ps,
	    f_xi, f_pi, f_al, f_be, f_ps, maspect, ires_xi,
	    ires_pi, ires_al, ires_be, ires_ps, wr_xi, wr_pi, wr_al, wr_be, wr_ps,
	    r, save_pt, lastwr);
      if (write_ricci) {
	ricci[0] = sRicci(f_xi, f_pi, f_ps, 0);
	int s = save_pt;
	for (int k = 1; k < lastwr; ++k) {
	  ricci[k] = sRicci(f_xi, f_pi, f_ps, s);
	  s += save_pt;
	}
	ricci[lastwr] = sRicci(f_xi, f_pi, f_ps, lastpt);
      }
      if (write_outnull) {
	outnull[0] = outgoing_null_f(f_al, f_be, f_ps, r, 0);
	int s = save_pt;
	for (int k = 1; k < lastwr; ++k) {
	  outnull[k] = outgoing_null(f_al, f_be, f_ps, r, s);
	  s += save_pt;
	}
	outnull[lastwr] = outgoing_null_b(f_al, f_be, f_ps, r, lastpt);
      }
      wr_step(num_wr, wr_names, t, bbh_shape, bbh_rank, coords, wr_ptrs);
      if (write_res) {
	wr_res_fn(residuals, wr_xi, wr_pi, wr_al, wr_be, wr_ps, save_pt, wr_shape);
	wr_step(num_wr_res, wr_res_names, t, bbh_shape, bbh_rank, coords, wr_res_ptrs);
      }
    }
    // CLOSE outfiles
    gft_close_all();
    if (write_itn) { ofs_itn.close(); }
    if (write_mtot) { ofs_mass.close(); }
    // PRINT resolution runtime and number of steps reaching maxit
    cout << factor << "-"+outfile+" written in "
	 << difftime(time(NULL),start_time) << " seconds" << endl;
  }
  // ******************** DONE LOOPING OVER RESOLUTIONS *********************
  return 0;
}
