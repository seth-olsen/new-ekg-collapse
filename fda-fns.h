#ifndef FDA_FNS_H_INCLUDED
#define FDA_FNS_H_INCLUDED

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
#include <vector> // for everything
#include <cmath> // for ICs
#include <map>
#include <string>

using namespace std;

typedef double dbl;

typedef vector<dbl> VD;

typedef map<int, double> MAPID;

// x^n and 1/x^n
inline dbl sq(dbl x) { return (x*x); }
inline dbl pw3(dbl x) { return (x*x*x); }
inline dbl pw4(dbl x) { return (x*x*x*x); }
inline dbl pw5(dbl x) { return (x*x*x*x*x); }
inline dbl pw6(dbl x) { return (x*x*x*x*x*x); }

inline dbl norm_inf(const VD& vec)
{
  return max( *max_element(vec.begin(), vec.end()),
	     -(*min_element(vec.begin(), vec.end())) );
}

// initialize r with r{ {WRITEDR, wr_dr} } then pass to set
void set_rmap(MAPID& r, int lastpt,
	      dbl dr, dbl dt, dbl lam, dbl rmin, dbl rmax)
{
  string error_response = "idc";
  if ((lastpt*dr) != (rmax - rmin)) {
    cout << "\nERROR: dr != (rmax-rmin)/lastpt\n" << endl;
    cout << endl << "continue? " << endl;
    cin >> error_response;
  }
  if (dt != lam*dr) {
    cout << "\nERROR: lam != dt/dr\n" << endl;
    cout << endl << "continue? " << endl;
    cin >> error_response;
  }
  dbl one_third = 1 / 3.0;
  r[RMIN] = rmin; r[RMAX] = rmax;
  r[DRVAL] = dr; r[DTVAL] = dt; r[LAMVAL] = lam;
  r[LAM2VAL] = 0.5 * lam;
  r[LAM6VAL] = one_third * r[LAM2VAL];
  r[INDR] = 1 / dr;
  r[IN2DR] = 0.5*r[INDR];
  r[INDRSQ] = sq(r[INDR]);
  r[INDT] = 1 / dt;
  r[INRMAX] = 1 / rmax;
  r[NEG2INDRSQ] = -2 * r[INDRSQ];
  r[CSOMM] = 0.75*lam + 0.5*dt*r[INRMAX];
  r[TWO_THIRDS] = 2 * one_third;
  r[FOUR_THIRDS] = 2 * r[TWO_THIRDS];
  r[TWELFTH] = 0.25 * one_third;
  r[FIVE_TWELFTHS] = 5 * r[TWELFTH];
  r[EIGHT_PI] = 8 * M_PI;
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
  return;
}

// ******** DIFFERENCES ********

inline dbl d_c(const VD& u, int ind)
{ return u[ind+1] - u[ind-1]; }

inline dbl d_f(const VD& u, int ind)
{ return -3*u[ind] + 4*u[ind+1] - u[ind+2]; }

inline dbl d_b(const VD& u, int ind)
{ return 3*u[ind] - 4*u[ind-1] + u[ind-2]; }

inline dbl d2_c(const VD& u, int ind)
{ return u[ind-1] - 2*u[ind] + u[ind+1]; }

inline dbl d2_f(const VD& u, int ind)
{ return 2*u[ind] - 5*u[ind+1] + 4*u[ind+2] - u[ind+3]; }

inline dbl d2_b(const VD& u, int ind)
{ return 2*u[ind] - 5*u[ind-1] + 4*u[ind-2] - u[ind-3]; }

inline dbl dlog_c(const VD& u, int ind)
{ return log(u[ind+1] / u[ind-1]); }
inline dbl dlog_f(const VD& u, int ind)
{ return -3*log(u[ind]) + 4*log(u[ind+1]) - log(u[ind+2]); }
inline dbl dlog_b(const VD& u, int ind)
{ return 3*log(u[ind]) - 4*log(u[ind+1]) + log(u[ind+2]); }

// DERIVATIVES

inline dbl ddr_c(const VD& u, MAPID& r, int ind)
{ return r[IN2DR]*(u[ind+1] - u[ind-1]); }

inline dbl ddr_f(const VD& u, MAPID& r, int ind)
{ return r[IN2DR]*(-3*u[ind] + 4*u[ind+1] - u[ind+2]); }

inline dbl ddr_b(const VD& u, MAPID& r, int ind)
{ return r[IN2DR]*(3*u[ind] - 4*u[ind-1] + u[ind-2]); }

inline dbl ddr2_c(const VD& u, MAPID& r, int ind)
{ return r[INDRSQ]*(u[ind-1] - 2*u[ind] + u[ind+1]); }

inline dbl ddr2_f(const VD& u, MAPID& r, int ind)
{ return r[INDRSQ]*(2*u[ind] - 5*u[ind+1] + 4*u[ind+2] - u[ind+3]); }

inline dbl ddr2_b(const VD& u, MAPID& r, int ind)
{ return r[INDRSQ]*(2*u[ind] - 5*u[ind-1] + 4*u[ind-2] - u[ind-3]); }

inline dbl ddrlog_c(const VD& u, MAPID& r, int ind)
{ return r[IN2DR]*log(u[ind+1] / u[ind-1]); }
inline dbl ddrlog_f(const VD& u, MAPID& r, int ind)
{ return r[IN2DR]*(-3*log(u[ind]) + 4*log(u[ind+1]) - log(u[ind+2])); }
inline dbl ddrlog_b(const VD& u, MAPID& r, int ind)
{ return r[IN2DR]*(3*log(u[ind]) - 4*log(u[ind+1]) + log(u[ind+2])); }

// d(r*u)/dr * 2*dr
inline dbl d_ru_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]*(r+dr) - u[ind-1]*(r-dr); }

inline dbl d_ru_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]*r + 4*u[ind+1]*(r+dr) - u[ind+2]*(r+2*dr); }

inline dbl d_ru_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]*r - 4*u[ind-1]*(r-dr) + u[ind-2]*(r-2*dr); }

// d(u/r)/dr * 2*dr
inline dbl d_urinv_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]/(r+dr) - u[ind-1]/(r-dr); }

inline dbl d_urinv_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]/r + 4*u[ind+1]/(r+dr) - u[ind+2]/(r+2*dr); }

inline dbl d_urinv_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]/r - 4*u[ind-1]/(r-dr) + u[ind-2]/(r-2*dr); }

// d(r^2*u)/dr * 2*dr
inline dbl d_r2u_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]*sq(r+dr) - u[ind-1]*sq(r-dr); }

inline dbl d_r2u_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]*sq(r) + 4*u[ind+1]*sq(r+dr) - u[ind+2]*sq(r+2*dr); }

inline dbl d_r2u_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]*sq(r) - 4*u[ind-1]*sq(r-dr) + u[ind-2]*sq(r-2*dr); }

// d(u/r^2)/dr * 2*dr
inline dbl d_ur2inv_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]/sq(r+dr) - u[ind-1]/sq(r-dr); }

inline dbl d_ur2inv_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]/sq(r) + 4*u[ind+1]/sq(r+dr) - u[ind+2]/sq(r+2*dr); }

inline dbl d_ur2inv_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]/sq(r) - 4*u[ind-1]/sq(r-dr) + u[ind-2]/sq(r-2*dr); }

// d(r^2*u^4)/dr * 2*dr
inline dbl d_r2u4_c(const VD& u, int ind, dbl dr, dbl r)
{ return ( sq(r+dr)*pw4(u[ind+1]) - sq(r-dr)*pw4(u[ind-1]) ); }

inline dbl d_r2u4_f(const VD& u, int ind, dbl dr, dbl r)
{ return ( -3*sq(r)*pw4(u[ind]) + 4*sq(r+dr)*pw4(u[ind+1]) - sq(r+2*dr)*pw4(u[ind+2]) ); }

inline dbl d_r2u4_b(const VD& u, int ind, dbl dr, dbl r)
{ return ( 3*sq(r)*pw4(u[ind]) - 4*sq(r-dr)*pw4(u[ind-1]) + sq(r-2*dr)*pw4(u[ind-2]) ); }

// d(u/v^2)/dr * 2*dr
inline dbl d_uv2inv_c(const VD& u, const VD& v, int ind)
{ return ( u[ind+1]/sq(v[ind+1]) - u[ind-1]/sq(v[ind-1]) ); }

inline dbl d_uv2inv_f(const VD& u, const VD& v, int ind)
{ return ( -3*u[ind]/sq(v[ind]) + 4*u[ind+1]/sq(v[ind+1]) - u[ind+2]/sq(v[ind+2]) ); }

inline dbl d_uv2inv_b(const VD& u, const VD& v, int ind)
{ return ( 3*u[ind]/sq(v[ind]) - 4*u[ind-1]/sq(v[ind-1]) + u[ind-2]/sq(v[ind-2]) ); }


// CRANK-NICHOLSON
inline dbl cn(const VD& old_f, const VD& f, int ind)
{
  return 0.5*(old_f[ind] + f[ind]);
}

inline void set2_cn(const VD& old_f1, const VD& old_f2,
		    const VD& f1, const VD& f2,
		    VD& cn_f1, VD& cn_f2, int npts)
{
  for (int k = 0; k < npts; ++k) {
      cn_f1[k] = 0.5*(old_f1[k] + f1[k]);
      cn_f2[k] = 0.5*(old_f2[k] + f2[k]);
  }
  return;
}
inline void set3_cn(const VD& old_f1, const VD& old_f2, const VD& old_f3,
		    const VD& f1, const VD& f2, const VD& f3,
		    VD& cn_f1, VD& cn_f2, VD& cn_f3, int npts)
{
  for (int k = 0; k < npts; ++k) {
      cn_f1[k] = 0.5*(old_f1[k] + f1[k]);
      cn_f2[k] = 0.5*(old_f2[k] + f2[k]);
      cn_f3[k] = 0.5*(old_f3[k] + f3[k]);
  }
  return;
}

// *********** BOUNDARY CONDITIONS ****************


inline void dirichlet0(VD& field)
{ field[0] = 0; }

inline void neumann0(VD& field)
{ field[0] = (4*field[1] - field[2]) / 3.0; }

inline void sommerfeld(const VD& oldfield, VD& field, MAPID& r,
		       int ind)
{
  field[ind] = r[CSOMM_RHS]*( r[LAMVAL]*(field[ind-1] + oldfield[ind-1])
			      - 0.25*r[LAMVAL]*(field[ind-2] + oldfield[ind-2])
			      + r[CSOMM_OLD]*oldfield[ind] );
  return;
}

inline dbl dirichlet0res(const VD& field)
{ return field[0]; }

inline dbl neumann0res(const VD& field, MAPID& r)
{ return ddr_f(field, r, 0); }

inline dbl sommerfeldres(const VD& oldfield, const VD& field,
			 MAPID& r, int ind)
{
  return 0.5*(field[ind] + oldfield[ind]) +
    r[ind]*( r[INDT]*(field[ind] - oldfield[ind])
	    + 0.5*(ddr_b(field,r,ind) + ddr_b(oldfield,r,ind)) );
}


// ********** DISSIPATION FUNCTIONS **************


// kreiss-oliger dissipation (p.23 choptuik notes)
inline dbl dissipate(dbl eps, const VD& u, int ind)
{ return -0.0625 * eps * ( u[ind-2] - 4*u[ind-1] + 6*u[ind]
			       - 4*u[ind+1] + u[ind+2] ); }

inline dbl symdiss1(dbl eps, const VD& u)
{ return -0.0625 * eps * ( u[3] - 4*u[2] + 7*u[1] - 4*u[0] ); }

inline dbl antidiss1(dbl eps, const VD& u)
{ return -0.0625 * eps * ( u[3] - 4*u[2] + 5*u[1] - 4*u[0] ); }

inline dbl symdiss0(dbl eps, const VD& u)
{ return -0.0625 * eps * ( 2*u[2] - 8*u[1] + 6*u[0] ); }

inline dbl antidiss0(dbl eps, const VD& u)
{ return -0.0625 * eps * ( 6*u[0] ); }



#endif

