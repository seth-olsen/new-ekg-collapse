#ifndef FDA_IO_H_INCLUDED
#define FDA_IO_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "fda-fns.h"
#include "ekg-fns.h"
#include "ekg-proc.h"

using namespace std;

typedef string str;

typedef void (*WR_FN)(const VD&, const VD&, const VD&, const VD&, const VD&,
		      const VD&, const VD&, const VD&, const VD&, const VD&,
		      const VD&, const VD&, const VD&, VD&,
		      VD&, VD&, VD&, VD&, VD&, VD&, VD&, VD&, VD&, VD&,
		      MAPID&, int, int);

typedef void (*WR_RES_FN)(const VD&, VD&, VD&, VD&, VD&, VD&, int, int);

void param_collect(char **source, int num, map<str, str>& dest) {
  for (int arg = 1; arg < num; ++arg) {
    if (source[arg][0] == '-') {
      dest[source[arg]] = source[arg+1];
    }
  }
}

void param_set(map<str, str>& p_all, map<str, str *>& p_str,
	       map<str, int *>& p_int, map<str, dbl *>& p_dbl,
	       map<str, bool *>& p_bool) {
  for (pair<str, str> p : p_all) {
    if (p_str.count(p.first)) { *p_str[p.first] = p.second; }
    else if (p_int.count(p.first)) { *p_int[p.first] = atoi(&p.second[0]); }
    else if (p_dbl.count(p.first)) { *p_dbl[p.first] = atof(&p.second[0]); }
    else if (p_bool.count(p.first)) { *p_bool[p.first] = (bool) atoi(&p.second[0]); }
  }
}

void param_print(str outfile, int lastpt, int save_pt, int nsteps,
		 int save_step, dbl lam, dbl r2m, dbl rmin, dbl rmax,
		 dbl dspn, dbl tol, int maxit, dbl ic_Dsq, dbl ic_r0,
		 dbl ic_Amp, int check_step, dbl dr, dbl dt, bool psi_hyp,
		 bool zero_pi, bool somm_cond, bool dspn_bound,
		 bool clean_hyp, bool clean_ell)
{
  str p_h = (psi_hyp) ? "true" : "false";
  str z_p = (zero_pi) ? "true" : "false";
  str s_c = (somm_cond) ? "true" : "false";
  str d_b = (dspn_bound) ? "true" : "false";
  str c_h = (clean_hyp) ? "true" : "false";
  str c_e = (clean_ell) ? "true" : "false";
  str param_str = "\noutfile name = " + outfile + "\ngrid size = " +
    to_string(lastpt) + " (" + to_string(save_pt) + "/write)\ntime steps = " +
    to_string(nsteps) + " (" + to_string(save_step) + "/write)\nlambda = " +
    to_string(lam) + "\nr2m = " + to_string(r2m) + "\nrmin = " + to_string(rmin)
    + "\nrmax = " + to_string(rmax) + "\ndissipation = " + to_string(dspn) +
    "\niterative tolerance = " + to_string(tol) + "\nmaximum iterations = " +
    to_string(maxit) + "\nic_Dsq = " + to_string(ic_Dsq) + "\nic_r0 = " +
    to_string(ic_r0) + "\nic_Amp = " + to_string(ic_Amp) + "\nmass check step = "
    + to_string(check_step) + "\nmaximum evolution time = " + to_string(nsteps*dt)
    + "\ndr = " + to_string(dr) + "\ndt = " + to_string(dt) +
    "\noptions:\nhyperbolic psi evolution = " + p_h + "\nzero pi_0 = " + z_p +
    "\nsommerfeld bc = " + s_c + "\ndissipation at bound = " + d_b
    + "\nclean hyperbolic update functions = " + c_h
    + "\nclean elliptic update functions = " + c_e + "\n";
  cout << param_str;
  ofstream specs;
  str specs_name = outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << param_str;
  specs.close();
  return;
}

// write fields using bbhutil
void wr_step(int nfields, const vector<char *>& files,
	     dbl time, int *shape, int rank, dbl *coordinates,
	     const vector<dbl *>& fields) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}

// read fields using bbhutil
void read_step(const vector<char *>& files, int times[], const vector<dbl *>& fields, int nfields) {
  for (int k = 0; k < nfields; ++k) {
    gft_read_brief(files[k], times[k], fields[k]);
  }
  return;
}

// GET COARSENED ARRAY FOR WRITING
void get_wr_f(const VD& f, VD& wr, int wr_shape, int save_pt)
{
  int k, s = 0;
  for (k = 0; k < wr_shape; ++k) {
    wr[k] = f[s];
    s += save_pt;
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////

void wr_fn0(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  return;
}
void wr_fn1(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k, s = 0;
  for (k = 0; k < lastwr+1; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    s += save_pt;
  }
  return;
}
void wr_fn2(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k, s = 0;
  for (k = 0; k < lastwr+1; ++k) {
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    s += save_pt;
  }
  return;
}
void wr_fn3(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k, s = 0;
  for (k = 0; k < lastwr+1; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    s += save_pt;
  }
  return;
}
void wr_fn4(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k, s = save_pt;
  
  for (k = 1; k < lastwr; ++k) {
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    s += save_pt;
    
  }
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn5(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k = 0, s = 0;
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  
  s += save_pt;
  for (k = 1; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn6(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k = 0, s = 0;
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  
  s += save_pt;
  for (k = 1; k < lastwr; ++k) {
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    
    s += save_pt;
  }
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn7(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  int k = 0, s = 0;
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  
  s += save_pt;
  for (k = 1; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn8(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  dbl al_ps2, lam6_part; // ADD
  for (k = 1; k < lastwr; ++k) {
    
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
  }
  return;
}
void wr_fn9(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	    const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	    const VD& f_xi, const VD& f_pi, const VD& f_al,
	    const VD& f_be, const VD& f_ps, VD& maspect,
	    VD& ires_xi, VD& ires_pi, VD& ires_al,
	    VD& ires_be, VD& ires_ps, VD& wr_xi,
	    VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	    int save_pt, int lastwr)
{
  
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  return;
}
void wr_fn10(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
  }
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  return;
}
void wr_fn11(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  return;
}
void wr_fn12(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
    
  }
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn13(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
    
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn14(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
    
  }
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn15(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  int k, s = save_pt;
  
  dbl al_ps2, lam6_part;
  for (k = 1; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    s += save_pt;
    
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn16(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  for (k = 2; k < lastwr; ++k) {
    
    s += save_pt;
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  }
  return;
}
void wr_fn17(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  for (k = 2; k < lastwr; ++k) {
    
    s += save_pt;
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  }
  s += save_pt;
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  return;
}
void wr_fn18(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  for (k = 2; k < lastwr; ++k) {
    
    s += save_pt;
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  }
  s += save_pt;
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  return;
}
void wr_fn19(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  for (k = 2; k < lastwr; ++k) {
    
    s += save_pt;
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  }
  s += save_pt;
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  return;
}
void wr_fn20(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn21(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn22(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn23(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn24(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  for (k = 2; k < lastwr; ++k) {
    s += save_pt;
    
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  }
  return;
}
void wr_fn25(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  return;
}
void wr_fn26(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    s += save_pt;
  }
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  return;
}
void wr_fn27(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  return;
}
void wr_fn28(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn29(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn30(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}
void wr_fn31(const VD& older_xi, const VD& older_pi, const VD& older_ps,
	     const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be, const VD& old_ps,
	     const VD& f_xi, const VD& f_pi, const VD& f_al,
	     const VD& f_be, const VD& f_ps, VD& maspect,
	     VD& ires_xi, VD& ires_pi, VD& ires_al,
	     VD& ires_be, VD& ires_ps, VD& wr_xi,
	     VD& wr_pi, VD& wr_al, VD& wr_be, VD& wr_ps, MAPID& r,
	     int save_pt, int lastwr)
{
  wr_xi[0] = f_xi[0];
  wr_pi[0] = f_pi[0];
  wr_al[0] = f_al[0];
  wr_be[0] = f_be[0];
  wr_ps[0] = f_ps[0];
  ires_xi[0] = f_xi[0];
  ires_pi[0] = -3*f_pi[0] + 4*f_pi[1] - f_pi[2];
  ires_al[0] = -3*f_al[0] + 4*f_al[1] - f_al[2];
  ires_be[0] = f_be[0];
  ires_ps[0] = -3*f_ps[0] + 4*f_ps[1] - f_ps[2];
  int k = 1, s = save_pt;
  
  dbl al_ps2, lam6_part;
  wr_xi[k] = f_xi[s];
  wr_pi[k] = f_pi[s];
  wr_al[k] = f_al[s];
  wr_be[k] = f_be[s];
  wr_ps[k] = f_ps[s];
  maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
  al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
  lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
  ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
  ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			r[DRVAL], r[s], al_ps2, lam6_part);
  ires_al[k] = iresalpha_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_be[k] = iresbeta_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  ires_ps[k] = irespsi_f(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
  
  s += save_pt;
  for (k = 2; k < lastwr; ++k) {
    wr_xi[k] = f_xi[s];
    wr_pi[k] = f_pi[s];
    wr_al[k] = f_al[s];
    wr_be[k] = f_be[s];
    wr_ps[k] = f_ps[s];
    maspect[k] = mass_aspect(f_al, f_be, f_ps, r, s);
    al_ps2 = old_al[s] / sq(old_ps[s]); // ADD
    lam6_part = r[LAM6VAL]*(d_c(old_be,s) + old_be[s]*(4*r[DRVAL]*r[-s] + 6*log(old_ps[s+1]/old_ps[s-1])));
    ires_xi[k] = iresxi_c(older_xi, old_xi, old_pi, old_al, old_be, old_ps, f_xi, s, r[LAMVAL], al_ps2);
    ires_pi[k] = irespi_c(older_pi, old_xi, old_pi, old_al, old_be, old_ps, f_pi, s, r[LAMVAL],
			  r[DRVAL], r[s], al_ps2, lam6_part);
    ires_al[k] = iresalpha_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_be[k] = iresbeta_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    ires_ps[k] = irespsi_c(f_xi, f_pi, f_al, f_be, f_ps, s, r[LAMVAL], r[DRVAL], r[s]);
    
    s += save_pt;
  }
  wr_xi[lastwr] = f_xi[s];
  wr_pi[lastwr] = f_pi[s];
  wr_al[lastwr] = f_al[s];
  wr_be[lastwr] = f_be[s];
  wr_ps[lastwr] = f_ps[s];
  maspect[lastwr] = mass_aspectR(f_al, f_be, f_ps, r, s);
  return;
}

WR_FN get_wr_fn(vector<dbl *>& wr_ptrs, bool write_xp, bool write_abp,
		bool write_ricci, bool write_maspect, bool write_outnull,
		bool write_ires_xp, bool write_ires_abp)
{
  int wr_fn_code = 0, ind = 0;
  auto it = wr_ptrs.begin();
  if (write_xp) {
    wr_fn_code += 1;
    ind += 2;
  }
  else { wr_ptrs.erase(it+ind, it+ind+2); }
  
  if (write_abp) {
    wr_fn_code += 2;
    ind += 3;
  }
  else { wr_ptrs.erase(it+ind, it+ind+3); }
  
  if (write_ricci) { ind += 1; }
  else { wr_ptrs.erase(it+ind); }
  
  if (write_maspect) {
    wr_fn_code += 4;
    ind += 1;
  }
  else { wr_ptrs.erase(it+ind); }

  if (write_outnull) { ind += 1; }
  else { wr_ptrs.erase(it+ind); }
  
  if (write_ires_xp) {
    wr_fn_code += 8;
    ind += 2;
  }
  else { wr_ptrs.erase(it+ind, it+ind+2); }
  
  if (write_ires_abp) { wr_fn_code += 16; }
  else { wr_ptrs.erase(it+ind, it+ind+3); }
  
  if (wr_fn_code == 0) { return &wr_fn0; }
  else if (wr_fn_code == 1) { return &wr_fn1; }
  else if (wr_fn_code == 2) { return &wr_fn2; }
  else if (wr_fn_code == 3) { return &wr_fn3; }
  else if (wr_fn_code == 4) { return &wr_fn4; }
  else if (wr_fn_code == 5) { return &wr_fn5; }
  else if (wr_fn_code == 6) { return &wr_fn6; }
  else if (wr_fn_code == 7) { return &wr_fn7; }
  else if (wr_fn_code == 8) { return &wr_fn8; }
  else if (wr_fn_code == 9) { return &wr_fn9; }
  else if (wr_fn_code == 10) { return &wr_fn10; }
  else if (wr_fn_code == 11) { return &wr_fn11; }
  else if (wr_fn_code == 12) { return &wr_fn12; }
  else if (wr_fn_code == 13) { return &wr_fn13; }
  else if (wr_fn_code == 14) { return &wr_fn14; }
  else if (wr_fn_code == 15) { return &wr_fn15; }
  else if (wr_fn_code == 16) { return &wr_fn16; }
  else if (wr_fn_code == 17) { return &wr_fn17; }
  else if (wr_fn_code == 18) { return &wr_fn18; }
  else if (wr_fn_code == 19) { return &wr_fn19; }
  else if (wr_fn_code == 20) { return &wr_fn20; }
  else if (wr_fn_code == 21) { return &wr_fn21; }
  else if (wr_fn_code == 22) { return &wr_fn22; }
  else if (wr_fn_code == 23) { return &wr_fn23; }
  else if (wr_fn_code == 24) { return &wr_fn24; }
  else if (wr_fn_code == 25) { return &wr_fn25; }
  else if (wr_fn_code == 26) { return &wr_fn26; }
  else if (wr_fn_code == 27) { return &wr_fn27; }
  else if (wr_fn_code == 28) { return &wr_fn28; }
  else if (wr_fn_code == 29) { return &wr_fn29; }
  else if (wr_fn_code == 30) { return &wr_fn30; }
  else if (wr_fn_code == 31) { return &wr_fn31; }
  else {
    cout << "ERROR: invalid write function code -> WRITING XI & PI" << endl;
    return &wr_fn1;
  }
}

void wr_res_fn0(const VD& residuals, VD& wr_xi, VD& wr_pi,
		VD& wr_al, VD& wr_be, VD& wr_ps,
		int save_pt, int wr_shape)
{
  return;
}
void wr_res_fn1(const VD& residuals, VD& wr_xi, VD& wr_pi,
		VD& wr_al, VD& wr_be, VD& wr_ps,
		int save_pt, int wr_shape)
{
  int k, s = 0, next = 5*save_pt;
  for (k = 0; k < wr_shape; ++k) {
    wr_xi[k] = residuals[s];
    wr_pi[k] = residuals[s+1];
    s += next;
  }
  return;
}
void wr_res_fn2(const VD& residuals, VD& wr_xi, VD& wr_pi,
		VD& wr_al, VD& wr_be, VD& wr_ps,
		int save_pt, int wr_shape)
{
  int k, s = 0, next = 5*save_pt;
  for (k = 0; k < wr_shape; ++k) {
    wr_al[k] = residuals[s+2];
    wr_be[k] = residuals[s+3];
    wr_ps[k] = residuals[s+4];
    s += next;
  }
  return;
}
void wr_res_fn3(const VD& residuals, VD& wr_xi, VD& wr_pi,
		VD& wr_al, VD& wr_be, VD& wr_ps,
		int save_pt, int wr_shape)
{
  int k, s = 0, next = 5*save_pt;
  for (k = 0; k < wr_shape; ++k) {
    wr_xi[k] = residuals[s];
    wr_pi[k] = residuals[s+1];
    wr_al[k] = residuals[s+2];
    wr_be[k] = residuals[s+3];
    wr_ps[k] = residuals[s+4];
    s += next;
  }
  return;
}

WR_RES_FN get_wr_res_fn(vector<dbl *>& wr_res_ptrs, bool write_res,
			bool write_xp, bool write_abp)
{
  int wr_fn_code = 0, ind = 0;
  auto it = wr_res_ptrs.begin();
  if (write_res) {
    if (write_xp) {
      wr_fn_code += 1;
      ind += 2;
    }
    else { wr_res_ptrs.erase(it+ind, it+ind+2); }
    
    if (write_abp) { wr_fn_code += 2; }
    else { wr_res_ptrs.erase(it+ind, it+ind+3); }
  }
  else { wr_res_ptrs.erase(wr_res_ptrs.begin(), wr_res_ptrs.end()); }
  
  if (wr_fn_code == 0) { return &wr_res_fn0; }
  else if (wr_fn_code == 1) { return &wr_res_fn1; }
  else if (wr_fn_code == 2) { return &wr_res_fn2; }
  else if (wr_fn_code == 3) { return &wr_res_fn3; }
  else {
    cout << "ERROR: invalid write residual function code -> NOT writing res" << endl;
    return &wr_res_fn0;
  }
}

vector<str> get_filenames(str outfile, bool write_xp, bool write_abp, bool write_ricci,
			  bool write_maspect, bool write_outnull, bool write_ires_xp, bool write_ires_abp)
{
  vector<str> filenames({});
  if (write_xp) {
    filenames.push_back(("xi-" + outfile + ".sdf"));
    filenames.push_back(("pi-" + outfile + ".sdf"));
  }
  if (write_abp) {
    filenames.push_back(("al-" + outfile + ".sdf"));
    filenames.push_back(("be-" + outfile + ".sdf"));
    filenames.push_back(("ps-" + outfile + ".sdf"));
  }
  if (write_ricci) { filenames.push_back(("ricci-" + outfile + ".sdf")); }
  if (write_maspect) { filenames.push_back(("maspect-" + outfile + ".sdf")); }
  if (write_outnull) { filenames.push_back(("outnull-" + outfile + ".sdf")); }
  if (write_ires_xp) {
    filenames.push_back(("ires_xi-" + outfile + ".sdf"));
    filenames.push_back(("ires_pi-" + outfile + ".sdf"));
  }
  if (write_ires_abp) {
    filenames.push_back(("ires_al-" + outfile + ".sdf"));
    filenames.push_back(("ires_be-" + outfile + ".sdf"));
    filenames.push_back(("ires_ps-" + outfile + ".sdf"));
  }
  return filenames;
}

vector<str> get_resfilenames(str outfile, bool write_res, bool write_xp, bool write_abp)
{
  vector<str> resfilenames({});
  if (write_res) {
    if (write_xp) {
      resfilenames.push_back(("res_xi-" + outfile + ".sdf"));
      resfilenames.push_back(("res_pi-" + outfile + ".sdf"));
    }
    if (write_abp) {
      resfilenames.push_back(("res_al-" + outfile + ".sdf"));
      resfilenames.push_back(("res_be-" + outfile + ".sdf"));
      resfilenames.push_back(("res_ps-" + outfile + ".sdf"));
    }
  }
  return resfilenames;
}



#endif


