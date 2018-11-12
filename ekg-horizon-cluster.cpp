#ifndef NEG2INDRSQ
#define NORM_FACTOR 32504
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
#define DSPN_WEIGHT 32531
#define EUP_WEIGHT 32532
#define HYPTOL 32533
#define ELLTOL 32534
#endif

#ifndef RESN_FACTOR
#define LASTPOINT 32535
#define SAVEPOINT 32536
#define LASTSTEP 32537
#define SAVESTEP 32538
#define NUM_POINTS 32539
#define NUM_ELL 32540
#define NUM_HYP 32541
#define RESN_FACTOR 32542
#define LASTWRITE 32543
#define WRITESHAPE 32544
#define LP_N 32545
#define LP_KL 32546
#define LP_KU 32547
#define LP_NRHS 32548
#define LP_LDAB 32549
#define LP_LDB 32550
#define MAX_ITN 32551
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
  int lastpt = 500; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 1000; // time steps
  int save_step = 2; // write only every (save_step)th time step
  dbl lam = 0.25; // dt/dr
  dbl r2m = 0.0;
  dbl rmin = 0.0;
  dbl rmax = 50.0;
  dbl dspn = 0.5; // dissipation coefficient
  dbl tol = 0.00000001; // iterative method tolerance
  dbl ell_tol = 0.001*tol;
  dbl ell_up_weight = 0.5;
  // for computing residual: 0 = inf-norm, 1 = 1-norm, 2 = 2-norm
  int resnorm_type = 0;
  int maxit = 50; // max iterations for debugging
  dbl ic_Dsq = 25.0; // gaussian width
  dbl ic_r0 = 15.0; // gaussian center
  dbl ic_Amp = 0.004; // gaussian amplitude
  int check_step = 10; // for monitoring invariant mass
  // note: set bools in command line with integers 1=true or 0=false
  bool psi_hyp = true; // update psi with hyperbolic evolution eqn after IC?
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
  bool write_maspect = false; // write mass aspect?
  bool write_outnull = false; // write outgoing null expansion?
  bool write_xp = false; // write xi and pi?
  bool write_abp = false; // write metric fields (alpha, beta, psi)?
  bool write_ires_xp = false; // write ires for xi and pi?
  bool write_ires_abp = false; // write ires for metric variables?
  bool horizon_search = true; // search for apparent horizon after each step?
  // variable to hold constant across resolutions
  str hold_const = "lambda"; // "lambda", "dt", or "dr"
  int nresn = 1; // 1, 2, or 3
  int resn0 = 16, resn1 = 4, resn2 = 1; // in order of priority
  int *resns[3] = {&resn0, &resn1, &resn2};
  
  map<str, str *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<str, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt",&save_pt},
      {"-nsteps",&nsteps}, {"-save_step",&save_step},
      {"-resnorm_type",&resnorm_type}, {"-maxit",&maxit},
      {"-check_step",&check_step}, {"-nresn",&nresn},
      {"-resn0",&resn0}, {"-resn1",&resn1}, {"-resn2",&resn2}};
  map<str, dbl *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin}, {"-rmax",&rmax},
      {"-dspn",&dspn}, {"-tol",&tol}, {"-ell_tol",&ell_tol}, {"-ell_up_weight",&ell_up_weight},
      {"-ic_Dsq",&ic_Dsq}, {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<str, bool *> p_bool { {"-psi_hyp",&psi_hyp}, {"-zero_pi",&zero_pi},
      {"-somm_cond",&somm_cond}, {"-dspn_bound",&dspn_bound}, {"-dr3_up",&dr3_up}, {"-dspn_psi",&dspn_psi},
      {"-write_res",&write_res},{"-write_ricci",&write_ricci}, {"-write_itn",&write_itn},
      {"-write_mtot",&write_mtot},{"-write_maspect",&write_maspect}, {"-write_outnull",&write_outnull},
      {"-write_xp",&write_xp}, {"-write_abp",&write_abp}, {"-static_metric",&static_metric},
      {"-write_ires_xp",&write_ires_xp}, {"-write_ires_abp",&write_ires_abp},
      {"-clean_hyp",&clean_hyp}, {"-clean_ell",&clean_ell}, {"-horizon_search",&horizon_search} };
  map<str, str> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);

  // OBTAIN RESOLUTION FACTORS
  vector<int> resolutions(nresn);
  for (int k = 0; k < nresn; ++k) {
    resolutions[k] = *resns[k];
  }

  // HERE psi_hyp = True means do NOT dissipate the point next to r = 0
  //void (*apply_dissipation)(const VD&, const VD&, VD&, VD&, int, dbl);
  //if (dspn_bound) { apply_dissipation = dissipationB_xp; }
  //else { apply_dissipation = dissipationNB_xp; }

  int n_ell = 3;
  if (psi_hyp) { n_ell = 2; }
  int n_hyp = 5 - n_ell;
  
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
    MAPID r { /*{WRITEDR,wr_dr},*/ {HYPTOL,tol}, {ELLTOL,ell_tol},
	      {EUP_WEIGHT,ell_up_weight}, {DSPN_WEIGHT,dspn} };
    set_rmap(r, lastpt, dr, dt, lam, rmin, rmax);

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

    VD ricci_vec(npts, 0);
    str ricci_file = "maxRicci-" + outfile + ".csv";
    ofstream ofs_ricci;
    ofs_ricci.open(ricci_file, ofstream::out);
    ofs_ricci << "i,t,max_Ricci" << endl;

    time_t start_time = time(NULL); // time for rough performance measure
    
// **********************************************************
// ******************* INITIAL DATA ************************
// **********************************************************

    dbl t = 0; // declare position and time variables
    for (int j = 1; j < npts; ++j) {
      r[j] = rmin + j*dr;
      r[-j] = 1 / r[j];
      f_xi[j] = ic_xi(r[j], ic_Amp, ic_Dsq, ic_r0);
      if (!zero_pi) { f_pi[j] = ic_pi(r[j], ic_Amp, ic_Dsq, ic_r0); }
    }
    if (r[lastpt] != rmax) { cout << "\n**\n****ERROR: r[J] != rmax****\n**\n" << endl; }
    // ASSUMING rmin = 0
    dirichlet0(f_xi);
    neumann0(f_pi);
    cout << outfile << endl;

    // SOLVE ELLIPTIC EQUATIONS FOR t=0  
    int itn = solve_t0_slow(f_xi, f_pi, f_al, f_be, f_ps, r, lastpt, 10*maxit);
    if (itn < 0) {
      if (itn > -npts) {
	record_horizon(r, outfile, lastpt, save_pt, nsteps, save_step, maxit,
		       ic_Dsq, ic_r0, ic_Amp, psi_hyp, zero_pi, somm_cond,
		       dspn_bound, clean_hyp, clean_ell, f_ps, -itn, 0, 0, 0);
      }
      else if (itn == -npts) {
	record_horizon(r, outfile, lastpt, save_pt, nsteps, save_step, maxit,
		       ic_Dsq, ic_r0, ic_Amp, psi_hyp, zero_pi, somm_cond,
		       dspn_bound, clean_hyp, clean_ell, f_ps, 0, 0, 0, 0);
      }
      else { cout << "\nUNDETERMINED ERROR IN T = 0 ELLIPTIC CONSTRAINTS" << endl; }
      return itn;
    }
    
// **********************************************************
// ******************* TIME STEPPING ************************
// *******************   & WRITING   ************************
// **********************************************************

    for (int i = 1; i < nsteps; ++i) {
      t = i * dt;
// ******************************************************************
// ******************************************************************
//         SOLVE HYPERBOLIC & ELLIPTIC EQUATIONS ITERATIVELY
// ******************************************************************
// ******************************************************************

      itn = solve_Hsearch(old_xi, old_pi, old_al, old_be, old_ps, f_xi, f_pi, f_al, f_be, f_ps,
			  cn_xi, cn_pi, cn_al, cn_be, cn_ps, res_hyp, res_ell, jac_zero,
			  r, lastpt, maxit, i, N, kl, ku, nrhs, ldab, ipiv, ldb);
      
// *****************************************************************
// ****************** ITERATIVE SOLUTION COMPLETE ******************
// *****************************************************************

    // **************************************************************************
    // *********************** kreiss-oliger DISSIPATION ************************
    // **************************************************************************
      dissipationNB_xp(old_xi, old_pi, f_xi, f_pi, lastpt-1, dspn);

    // **************************************************************************
    // ************************ APPARENT HORIZON SEARCH *************************
    // **************************************************************************
      if (itn < 0 || itn == maxit) {
	if (outgoing_null_b(f_al, f_be, f_ps, r, lastpt) <= 0) {
	  record_horizon(r, outfile, lastpt, save_pt, nsteps, save_step, maxit,
			 ic_Dsq, ic_r0, ic_Amp, psi_hyp, zero_pi, somm_cond,
			 dspn_bound, clean_hyp, clean_ell, f_ps, lastpt, itn, i, t);
	  return -lastpt;
	}
	int k = lastpt;
	while (--k > 0) {
	  if (outgoing_null(f_al, f_be, f_ps, r, k) <= 0) {
	    record_horizon(r, outfile, lastpt, save_pt, nsteps, save_step, maxit,
			   ic_Dsq, ic_r0, ic_Amp, psi_hyp, zero_pi, somm_cond,
			   dspn_bound, clean_hyp, clean_ell, f_ps, k, itn, i, t);
	    return -k;
	  }
	}
	if (outgoing_null_f(f_al, f_be, f_ps, r, 0) <= 0) {
	  record_horizon(r, outfile, lastpt, save_pt, nsteps, save_step, maxit,
			 ic_Dsq, ic_r0, ic_Amp, psi_hyp, zero_pi, somm_cond,
			 dspn_bound, clean_hyp, clean_ell, f_ps, 0, itn, i, t);
	  return -npts;
	}
      }

      // **************************************************************************
      // ************************ SAVE MAX VALUE OF RICCI *************************
      // **************************************************************************
      if (i % save_step == 0) {
	for (int j = 0; j < npts; ++j) {
	  ricci_vec[j] = sRicci(f_xi, f_pi, f_ps, j);
	}
	ofs_ricci << i <<","<< t <<","<< r[EIGHT_PI]*norm_inf(ricci_vec) << endl;
      }


      
    }
    
// ***********************************************************************
// ***********************************************************************
// *********************** FULL SOLUTION COMPLETE ************************
// ***********************************************************************
// ***********************************************************************

    ofs_ricci.close();
    // PRINT resolution runtime and number of steps reaching maxit
    cout << outfile+" written in " << difftime(time(NULL),start_time) << " seconds" << endl;
  }
  // ******************** DONE LOOPING OVER RESOLUTIONS *********************
  return 0;
}
