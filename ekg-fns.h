#ifndef EKG_FNS_H_INCLUDED
#define EKG_FNS_H_INCLUDED

#include "fda-fns.h"
#include <vector> // for everything
#include <cmath> // for ICs


////////////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************
// **********************   FDA IMPLEMENTATIONS   **************************
// *************************************************************************
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_xi(const VD& old_xi, const VD& cn_xi, const VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int k)
{
  return old_xi[k] + r[LAM2VAL]*(cn_al[k+1]*cn_pi[k+1]/sq(cn_ps[k+1]) + cn_be[k+1]*cn_xi[k+1]
				 - cn_al[k-1]*cn_pi[k-1]/sq(cn_ps[k-1]) - cn_be[k-1]*cn_xi[k-1]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_pi(const VD& old_pi, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
		  const VD& cn_be, const VD& cn_ps, MAPID& r, int k)
{
  double lam6_part = r[LAM6VAL] * (d_c(cn_be,k) + cn_be[k]*(4*r[DRVAL]*r[-k] + 6*d_c(cn_ps,k)/cn_ps[k]));
  return ( old_pi[k]*(1 - lam6_part) + (r[LAM2VAL] / sq(r[k]*sq(cn_ps[k]))) *
    (sq(r[k+1]*sq(cn_ps[k+1]))*(cn_al[k+1]*cn_xi[k+1]/sq(cn_ps[k+1]) + cn_be[k+1]*cn_pi[k+1])
     - sq(r[k-1]*sq(cn_ps[k-1]))*(cn_al[k-1]*cn_xi[k-1]/sq(cn_ps[k-1]) + cn_be[k-1]*cn_pi[k-1])) )
    / (1 + lam6_part);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_hyp_ps(const VD& old_ps, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
		      const VD& cn_be, const VD& cn_ps, MAPID& r, int k)
{
  double dt12_part = r[DT_TWELVE] * (ddr_c(cn_be,r,k) + 2*cn_be[k]*r[-k]);
  return ( old_ps[k]*(1 + dt12_part) + r[DTVAL]*cn_be[k]*ddr_c(cn_ps,r,k) ) / (1 - dt12_part);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_hyp_ps(const VD& f_ps, MAPID& r, int k)
{
  return r[CPSI_RHS]*( 1 - r[CPSI_RRM1]*f_ps[k-1] - r[CPSI_RRM2]*f_ps[k-2] );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_hyp_resPs(const VD& old_ps, const VD& f_ps, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
			 const VD& cn_be, const VD& cn_ps, MAPID& r, int k)
{
  return r[INDT]*(f_ps[k] - old_ps[k]) - cn_be[k]*ddr_c(cn_ps,r,k)
    - r[DT_TWELVE]*(ddr_c(cn_be,r,k) + 2*cn_be[k]*r[-k])*(f_ps[k] + old_ps[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_hyp_resPs(const VD& f_ps, MAPID& r, int k)
{
  return r[INRMAX]*( r[CPSI_RR]*f_ps[k] + r[CPSI_RRM1]*f_ps[k-1] + r[CPSI_RRM2]*f_ps[k-2] - 1 );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resXi(const VD& old_xi, const VD& f_xi, const VD& cn_xi, const VD& cn_pi,
		     const VD& cn_al, const VD& cn_be, const VD& cn_ps, MAPID& r, int k)
{
  return r[INDT]*( f_xi[k] - old_xi[k]
		   - r[LAM2VAL]*(cn_al[k+1]*cn_pi[k+1]/sq(cn_ps[k+1]) + cn_be[k+1]*cn_xi[k+1]
				 - cn_al[k-1]*cn_pi[k-1]/sq(cn_ps[k-1]) - cn_be[k-1]*cn_xi[k-1]) );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resPi(const VD& old_pi, const VD& f_pi, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
		     const VD& cn_be, const VD& cn_ps, MAPID& r, int k)
{
  return r[INDT]*( f_pi[k] - old_pi[k]
		   - ( (r[LAM2VAL] / sq(r[k]*sq(cn_ps[k]))) *
		       (sq(r[k+1]*sq(cn_ps[k+1]))*(cn_al[k+1]*cn_xi[k+1]/sq(cn_ps[k+1]) + cn_be[k+1]*cn_pi[k+1])
			- sq(r[k-1]*sq(cn_ps[k-1]))*(cn_al[k-1]*cn_xi[k-1]/sq(cn_ps[k-1]) + cn_be[k-1]*cn_pi[k-1])) )
		   + ( (f_pi[k] + old_pi[k]) * r[LAM6VAL] *
		       (d_c(cn_be,k) + cn_be[k]*(4*r[DRVAL]*r[-k] + 6*d_c(cn_ps,k)/cn_ps[k])) )  );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resPs(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		     MAPID& r, int k)
{
  return ddr2_c(f_ps,r,k) + M_PI * f_ps[k] * (sq(f_xi[k]) + sq(f_pi[k])) +
    2*r[-k]*ddr_c(f_ps,r,k) +
    pw5(f_ps[k]) * sq(ddr_c(f_be,r,k) - r[-k]*f_be[k]) / (12*sq(f_al[k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resBe(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		     MAPID& r, int k)
{
  return ddr2_c(f_be,r,k)
    + 12 * M_PI * f_al[k] * f_xi[k] * f_pi[k] / sq(f_ps[k])
    + ( (ddr_c(f_be,r,k) - r[-k]*f_be[k]) *
	(2*r[-k] + 6*(ddr_c(f_ps,r,k)/f_ps[k]) - (ddr_c(f_al,r,k)/f_al[k])) );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resAl(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		     MAPID& r, int k)
{
  return ddr2_c(f_al,r,k)
    - 8 * M_PI * f_al[k] * sq(f_pi[k])
    + 2 * ddr_c(f_al,r,k) * (r[-k] + (ddr_c(f_ps,r,k)/f_ps[k]))
    - (2 * pw4(f_ps[k]) * sq(ddr_c(f_be,r,k) - r[-k]*f_be[k]) / (3*f_al[k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_resPs(const VD& f_ps, MAPID& r, int k)
{ return ddr_b(f_ps,r,k) + r[INRMAX]*(f_ps[k] - 1); }
inline dbl RxfdaR_resPs(const VD& f_ps, MAPID& r, int k)
{ return r[k]*ddr_b(f_ps,r,k) + f_ps[k] - 1; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_resBe(const VD& f_be, MAPID& r, int k)
{ return ddr_b(f_be,r,k) + r[INRMAX]*f_be[k]; }
inline dbl RxfdaR_resBe(const VD& f_be, MAPID& r, int k)
{ return r[k]*ddr_b(f_be,r,k) + f_be[k]; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_resAl(const VD& f_al, MAPID& r, int k)
{ return ddr_b(f_al,r,k) + r[INRMAX]*(f_al[k] - 1); }
inline dbl RxfdaR_resAl(const VD& f_al, MAPID& r, int k)
{ return r[k]*ddr_b(f_al,r,k) + f_al[k] - 1; }
////////////////////////////////////////////////////////////////////////////////////////////////
// **********************************************************
// **********************************************************
//             INITIAL AND BOUNDARY CONDITIONS
// **********************************************************
// **********************************************************
////////////////////////////////////////////////////////////////////////////////////////////////
// for gaussian field or sin(coeff*r)*cos(coeff*t)/(coeff*r)
inline dbl ic_sol(dbl r, dbl amp, dbl dsq, dbl r0)
{ return amp * exp(-(r - r0)*(r - r0)/dsq); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl ic_xi(dbl r, dbl amp, dbl dsq, dbl r0)
{ return -2 * (r - r0) * amp * exp(-(r - r0)*(r - r0)/dsq) / dsq; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl ic_pi(dbl r, dbl amp, dbl dsq, dbl r0)
{ return ic_xi(r, amp, dsq, r0) + ic_sol(r, amp, dsq, r0)/r; }
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// ******** IRES FUNCTIONS *************************

// *****IMPORTANT NOTE:************
// anything with sq(psi) or /sq() may be wrong
// had to replace sqin and did so in a rush and not sure now if the replacements were right


// centered ires_f1 = f1 - ires_c(f1, f2, k, c1, d1, c2, d2)
//                   - oldf1 - ires_c(oldf1, oldf2, k, c1, d1, c2, d2)
inline dbl iresxi_c(const VD& older_xi, const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be,
		    const VD& old_ps, const VD& f_xi, int k, dbl lam, dbl al_ps2)
{ 
  return f_xi[k] - older_xi[k]
    - lam*( old_pi[k]*d_c(old_al,k)/sq(old_ps[k]) + al_ps2*d_c(old_pi,k)
	    - 2*al_ps2*old_pi[k]*(log(old_ps[k+1]/old_ps[k-1])) 
	    + old_xi[k]*d_c(old_be,k) + old_be[k]*d_c(old_xi,k) );
}

inline dbl irespi_c(const VD& older_pi, const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be,
		    const VD& old_ps, const VD& f_pi, int k, dbl lam, dbl dr, dbl r, dbl al_ps2, dbl lam6_part)
{
  return f_pi[k] - older_pi[k] + 4*old_pi[k]*lam6_part
    - lam*( (old_be[k]*old_pi[k] + al_ps2*old_xi[k])*(2/r + 4*(log(old_ps[k+1]/old_ps[k-1])))
	    + old_xi[k]*d_c(old_al,k)/sq(old_ps[k]) + al_ps2*d_c(old_xi,k)
	    - 2*al_ps2*old_xi[k]*(log(old_ps[k+1]/old_ps[k-1])) 
	    + old_pi[k]*d_c(old_be,k) + old_be[k]*d_c(old_pi,k) );
}

inline dbl irespsihyp_c(const VD& older_ps, const VD& old_ps, const VD& f_ps, int k, dbl lam6_part)
{ return f_ps[k] - older_ps[k] - old_ps[k]*lam6_part; }

inline dbl irespsi_c(const VD& xi, const VD& pi,
		     const VD& alpha, const VD& beta,
		     const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return d2_c(psi,k) + sq(dr)*psi[k]*M_PI*(sq(xi[k]) + sq(pi[k])) + dr*d_c(psi,k)/r
    + pw5(psi[k])*sq(r*d_urinv_c(beta,k,dr,r)) / (48*sq(alpha[k]));
}

inline dbl iresbeta_c(const VD& xi, const VD& pi,
		      const VD& alpha, const VD& beta,
		      const VD& psi, int k, dbl lam, dbl dr, dbl r)
{ return d2_c(beta,k) + sq(dr)*12*M_PI*alpha[k]*xi[k]*pi[k] / sq(psi[k]) +
    0.25*r*d_urinv_c(beta,k,dr,r)*(6*d_c(psi,k)/psi[k] - d_c(alpha,k)/alpha[k] + 4*dr/r); }

inline dbl iresalpha_c(const VD& xi, const VD& pi,
		       const VD& alpha, const VD& beta,
		       const VD& psi, int k, dbl lam, dbl dr, dbl r)
{ return d2_c(alpha,k) - sq(dr)*8*M_PI*alpha[k]*sq(pi[k])
    + 0.25*d_c(alpha,k)*(d_c(psi,k)/psi[k] + 4*dr/r)
    - pw4(psi[k])*sq(r*d_urinv_c(beta,k,dr,r)) / (6*alpha[k]); }


inline double irespsi_f(const vector<double>& xi, const vector<double>& pi,
			const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, double lam, double dr, double r) {
  return d2_f(psi,ind) + sq(dr)*psi[ind]*M_PI*(sq(xi[ind]) + sq(pi[ind])) + dr*d_f(psi,ind)/r
    + pw5(psi[ind])*sq(r*d_urinv_f(beta,ind,dr,r)) / (48*sq(alpha[ind])); }

inline double iresbeta_f(const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double lam, double dr, double r) {
  return d2_f(beta,ind) + sq(dr)*12*M_PI*alpha[ind]*xi[ind]*pi[ind] / sq(psi[ind]) +
    0.25*r*d_urinv_f(beta,ind,dr,r)*(6*d_f(psi,ind)/psi[ind] - d_f(alpha,ind)/alpha[ind] + 4*dr/r); }

inline double iresalpha_f(const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double lam, double dr, double r) {
  return d2_f(alpha,ind) - sq(dr)*8*M_PI*alpha[ind]*sq(pi[ind])
    + 0.25*d_f(alpha,ind)*(d_f(psi,ind)/psi[ind] + 4*dr/r)
    - pw4(psi[ind])*sq(r*d_urinv_f(beta,ind,dr,r)) / (6*alpha[ind]); }

// *********************************************
// **************  DIAGNOSTICS  ****************
// *********************************************
inline dbl mass_aspect(const VD& alpha, const VD& beta, const VD& psi, MAPID& r, int k)
{
  return (0.5*sq(r[k])*r[INDRSQ])*( (r[k]*pw6(psi[k]) / (18*sq(alpha[k]))) * sq(0.5*d_c(beta,k) - r[DRVAL]*beta[k]*r[-k])
			      - d_c(psi,k) * d_ru_c(psi,k,r[DRVAL],r[k]) );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl mass_aspectR(const VD& alpha, const VD& beta, const VD& psi, MAPID& r, int k)
{
  return (0.5*sq(r[k])*r[INDRSQ])*( r[k]*pw6(psi[k]) * sq(0.5*d_b(beta,k) - r[DRVAL]*beta[k]*r[-k]) / (18*sq(alpha[k]))
			      - d_b(psi,k) * d_ru_b(psi,k,r[DRVAL],r[k]) );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null(const VD& alpha, const VD& beta,
			 const VD& psi, MAPID& r, int k)
{
  return (r[k]*ddr_c(beta,r,k) - beta[k])/(3*alpha[k]) + 1/sq(psi[k]) + 2*r[k]*ddr_c(psi,r,k)/pw3(psi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null_f(const VD& alpha, const VD& beta,
			   const VD& psi, MAPID& r, int k)
{
  return (r[k]*ddr_f(beta,r,k) - beta[k])/(3*alpha[k]) + 1/sq(psi[k]) + 2*r[k]*ddr_f(psi,r,k)/pw3(psi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null_b(const VD& alpha, const VD& beta,
			   const VD& psi, MAPID& r, int k)
{
  return (r[k]*ddr_b(beta,r,k) - beta[k])/(3*alpha[k]) + 1/sq(psi[k]) + 2*r[k]*ddr_b(psi,r,k)/pw3(psi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl sRicci(const VD& f_xi, const VD& f_pi, const VD& f_ps, int k)
{
  return (sq(f_xi[k]) - sq(f_pi[k])) / pw4(f_ps[k]);
}

#endif
