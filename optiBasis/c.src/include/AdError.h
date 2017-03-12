#pragma once 

#include "SimulationRegion.h"
#include <omp.h>

#include "IkError.h"

template <typename HatComputer>
class AdError 
{
public:
  typedef HatComputer HatType;
  AdError (const double & beta,
	   const vector<int> & KK,
	   const vector<HatComputer> & hat_comput,
	   const int & l_cut,
	   const int b_style,
	   const int numb_threads);
  double estimate (const double & q2,
		   const int & natoms,
		   const SimulationRegion<double> & region) ;
  double estimate (vector<double > & result,
		   const double & q2,
		   const int & natoms,
		   const SimulationRegion<double> & region) ;
private:
  double beta;
  vector<int > KK;
  const vector<HatComputer> & hat_comput;
  int l_cut;
  int b_style;
  int func_numb_threads;
  double sum_o1;
  vector<double > sum_o1_deriv;
  void prepare_sum (const SimulationRegion<double> & region);
  void prepare_sum_deriv (const SimulationRegion<double> & region) {};
  void prepare_sum_pppm (const SimulationRegion<double> & region) {};
  void prepare_sum_pppm_deriv (const SimulationRegion<double> & region) {};
  void compute_G (double * gm, const int ip, const int jp, const int kp,
		  const SimulationRegion<double> & region) const;
  void rect_vec (double * mm, const int ip, const int jp, const int kp,
		 const SimulationRegion<double> & region) const;
};

template<typename HatComputer> 
AdError<HatComputer>::
AdError (const double & beta_,
	 const vector<int> & KK_,
	 const vector<HatComputer> & hat_comput_,
	 const int & l_cut_,
	 const int b_style_,
	 const int numb_threads)
    : beta (beta_),
      KK (KK_),
      hat_comput (hat_comput_),
      l_cut (l_cut_),
      b_style (b_style_),
      func_numb_threads (numb_threads),
      sum_o1 (0)
{
  
}

template<typename HatComputer> 
void
AdError<HatComputer>::
rect_vec (double * mm,
	  const int ip,
	  const int jp,
	  const int kp,
	  const SimulationRegion<double> & region) const
{
  const double * rec_boxt = region.getRecBoxTensor();
  mm[0] = MathUtilities::dot<double> (ip, jp, kp, rec_boxt[0], rec_boxt[3], rec_boxt[6]);
  mm[1] = MathUtilities::dot<double> (ip, jp, kp, rec_boxt[1], rec_boxt[4], rec_boxt[7]);
  mm[2] = MathUtilities::dot<double> (ip, jp, kp, rec_boxt[2], rec_boxt[5], rec_boxt[8]);
}

template<typename HatComputer> 
void
AdError<HatComputer>::
compute_G (double * gm,
	   const int ip,
	   const int jp,
	   const int kp,
	   const SimulationRegion<double> & region) const
{
  gm[0] = gm[1] = gm[2] = 0;
  if (abs(ip) + abs(jp) + abs(kp) == 0) return;

  double mm[3] = {0};
  rect_vec (mm, ip, jp, kp, region);

  double mm2 = MathUtilities::dot<double> (mm, mm);
  double expp = exp (- M_PI*M_PI/beta/beta*mm2) / (mm2);

  for (int dd = 0; dd < 3; ++dd) {
    gm[dd] = -4. * M_PI * mm[dd] * expp;
  }
}

template<typename HatComputer> 
void
AdError<HatComputer>::
prepare_sum (const SimulationRegion<double> & region)
{
  double V = region.getVolume ();
  double tmp_sum_o1 = 0;
  double scale = M_PI * M_PI / beta / beta;
  const double * rec_boxt = region.getRecBoxTensor();

#pragma omp parallel for reduction (+:tmp_sum_o1) num_threads(func_numb_threads)
  for (int m0 = -int(KK[0]/2); m0 <= KK[0]/2; ++m0){
    for (int m1 = -int(KK[1]/2); m1 <= KK[1]/2; ++m1){
      for (int m2 = -int(KK[2]/2); m2 <= KK[2]/2; ++m2){
	if (abs(m0) + abs(m1) + abs(m2) == 0) continue;
	int idxmm[3];
	idxmm[0] = m0;
	idxmm[1] = m1;
	idxmm[2] = m2;

	double mm[3] = {0};
	rect_vec (mm, idxmm[0], idxmm[1], idxmm[2], region);
	double mm2 = MathUtilities::dot<double> (mm, mm);
	double expp = exp (- scale * mm2) / (mm2);
	double gm[3];
	for (int dd = 0; dd < 3; ++dd) {
	  gm[dd] = -4. * M_PI * mm[dd] * expp;
	}
	double gm2 = MathUtilities::dot<double > (gm, gm);
	
	double o1e = 0;

	for (int ll = -l_cut; ll <= l_cut; ll += 1){
	  if (ll == 0) continue;

	  double zl[3];
	  zl[0] = hat_comput[0].value(ll * KK[0] + idxmm[0]) / hat_comput[0].value(idxmm[0]);
	  zl[1] = hat_comput[1].value(ll * KK[1] + idxmm[1]) / hat_comput[1].value(idxmm[1]);
	  zl[2] = hat_comput[2].value(ll * KK[2] + idxmm[2]) / hat_comput[2].value(idxmm[2]);

	  double fm_scale = expp;
	  for (int alpha = 0; alpha < 3; ++alpha){
	    double tmp_scale = zl[alpha] * fm_scale * -4. * M_PI * ll * KK[alpha];
	    double fm[3];
	    fm[0] = tmp_scale * rec_boxt[alpha*3+0];
	    fm[1] = tmp_scale * rec_boxt[alpha*3+1];
	    fm[2] = tmp_scale * rec_boxt[alpha*3+2];
	    fm[0] += gm[0] * zl[alpha];
	    fm[1] += gm[1] * zl[alpha];
	    fm[2] += gm[2] * zl[alpha];
	    o1e += 1. * gm2 * zl[alpha] * zl[alpha];
	    o1e += MathUtilities::dot<double > (fm, fm);
	  }
	}
	tmp_sum_o1 += o1e;
      }
    }
  }
  
  sum_o1 = tmp_sum_o1;
  sum_o1 /= (2. * M_PI * V * 2. * M_PI * V);
}

template<typename HatComputer> 
double
AdError<HatComputer>::
estimate (const double & q2,
	  const int & natoms,
	  const SimulationRegion<double> & region) 
{
  if (b_style == B_STYLE::SPME){
    prepare_sum (region);
  }
  else if (b_style == B_STYLE::PPPM){
    prepare_sum_pppm (region);
  }
  else {
    cerr << "unknow b style, exit" <<endl;
    exit(1);
  }

  return sqrt (sum_o1 * q2 * q2 / double(natoms)) * ElectrostaticConvertion;
}

template<typename HatComputer> 
double
AdError<HatComputer>::
estimate (vector<double > & result,
	  const double & q2,
	  const int & natoms,
	  const SimulationRegion<double> & region) 
{
  if (b_style == B_STYLE::SPME){
    prepare_sum_deriv (region);
  }
  else if (b_style == B_STYLE::PPPM){
    prepare_sum_pppm_deriv (region);
  }
  else {
    cerr << "unknow b style, exit" <<endl;
    exit(1);
  }

  unsigned basis_size = hat_comput[0].basis_size();
  result.resize (basis_size);
  
  for (unsigned pp = 0; pp < basis_size; ++pp){
    result[pp] = sum_o1_deriv[pp] * 0.5 / sqrt (sum_o1) * sqrt(q2 * q2 / double(natoms)) * ElectrostaticConvertion ;
  }

  return sqrt (sum_o1 * q2 * q2 / double(natoms)) * ElectrostaticConvertion;
}
