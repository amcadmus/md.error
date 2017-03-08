#pragma once

#include "SimulationRegion.h"

const static double ElectrostaticConvertion = 138.93545756169981341199;

template <typename HatComputer>
class IkError 
{
public:
  typedef HatComputer HatType;
  IkError (const double & beta,
	   const vector<int> & KK,
	   const vector<HatComputer> & hat_comput,
	   const int & l_cut);
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
  double sum_o1;
  vector<double > sum_o1_deriv;
  void prepare_sum (const SimulationRegion<double> & region);
  void prepare_sum_deriv (const SimulationRegion<double> & region);
  void compute_G (double * gm, const int ip, const int jp, const int kp,
		  const SimulationRegion<double> & region) const;
  void rect_vec (double * mm, const int ip, const int jp, const int kp,
		 const SimulationRegion<double> & region) const;
};



template<typename HatComputer> 
IkError<HatComputer>::
IkError (const double & beta_,
	 const vector<int> & KK_,
	 const vector<HatComputer> & hat_comput_,
	 const int & l_cut_)
    : beta (beta_),
      KK (KK_),
      hat_comput (hat_comput_),
      l_cut (l_cut_)
{
  
}

template<typename HatComputer> 
void
IkError<HatComputer>::
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
IkError<HatComputer>::
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
IkError<HatComputer>::
prepare_sum (const SimulationRegion<double> & region)
{
  double V = region.getVolume ();
  sum_o1 = 0;

  for (int m0 = -int(KK[0]/2); m0 <= KK[0]/2; ++m0){
    double ip = m0;
    for (int m1 = -int(KK[1]/2); m1 <= KK[1]/2; ++m1){
      double jp = m1;
      for (int m2 = -int(KK[2]/2); m2 <= KK[2]/2; ++m2){
	double kp = m2;
	if (abs(m0) + abs(m1) + abs(m2) == 0) continue;
	int idxmm[3];
	idxmm[0] = m0;
	idxmm[1] = m1;
	idxmm[2] = m2;

	double gm[3];
	compute_G (gm, ip, jp, kp, region);
	double gm2 = MathUtilities::dot<double > (gm, gm);	
	
	double o1e = 0;
	
	for (int dd = 0; dd < 3; ++dd){
	  double hat_phi0 = hat_comput[dd].value (idxmm[dd]);
	  double hat_phi02i = 1./(hat_phi0 * hat_phi0);
	  for (int ll = -l_cut; ll <= l_cut; ll += 1){
	    if (ll == 0) continue;
	    double hat_phil = hat_comput[dd].value(idxmm[dd] + ll * KK[dd]);
	    o1e += hat_phil * hat_phil * hat_phi02i;
	    // cout << hat_phil * hat_phil * hat_phi02i << endl;
	  }
	}
	sum_o1 += 2. * gm2 * o1e;
      }
    }
  }
  sum_o1 /= (2. * M_PI * V * 2. * M_PI * V);
}


template<typename HatComputer> 
void
IkError<HatComputer>::
prepare_sum_deriv (const SimulationRegion<double> & region)
{
  double V = region.getVolume ();
  unsigned basis_size = hat_comput[0].basis_size();
  sum_o1 = 0;
  sum_o1_deriv.resize (basis_size, 0);
  fill (sum_o1_deriv.begin(), sum_o1_deriv.end(), 0.);

  for (int m0 = -int(KK[0]/2); m0 <= KK[0]/2; ++m0){
    double ip = m0;
    for (int m1 = -int(KK[1]/2); m1 <= KK[1]/2; ++m1){
      double jp = m1;
      for (int m2 = -int(KK[2]/2); m2 <= KK[2]/2; ++m2){
	double kp = m2;
	if (abs(m0) + abs(m1) + abs(m2) == 0) continue;
	int idxmm[3];
	idxmm[0] = m0;
	idxmm[1] = m1;
	idxmm[2] = m2;

	double gm[3];
	compute_G (gm, ip, jp, kp, region);
	double gm2 = MathUtilities::dot<double > (gm, gm);	
	
	double o1e = 0;
	vector<double> o1e_deriv (basis_size, 0);
	
	for (int dd = 0; dd < 3; ++dd){
	  const double &	hat_phi0   = hat_comput[dd].value (idxmm[dd]);
	  const vector<double>&	hat_basis0 = hat_comput[dd].basis_value (idxmm[dd]);
	  double		hat_phi02i = 1./(hat_phi0 * hat_phi0);
	  double		hat_phi03i = 1./(hat_phi0 * hat_phi0 * hat_phi0);
	  for (int ll = -l_cut; ll <= l_cut; ll += 1){
	    if (ll == 0) continue;
	    const double &		hat_phil   = hat_comput[dd].value(idxmm[dd] + ll * KK[dd]);
	    const vector<double> &	hat_basisl = hat_comput[dd].basis_value(idxmm[dd] + ll * KK[dd]);
	    o1e += hat_phil * hat_phil * hat_phi02i;
	    for (unsigned pp = 0; pp < basis_size; ++pp){
	      o1e_deriv[pp] += (2 * hat_phil * hat_phi02i * hat_basisl[pp] - 2 * hat_phil * hat_phil * hat_phi03i * hat_basis0[pp] );
	    }
	  }
	}
	sum_o1 += 2. * gm2 * o1e;
	for (unsigned pp = 0; pp < basis_size; ++pp){
	  sum_o1_deriv[pp] += 2. * gm2 * o1e_deriv[pp];
	}
      }
    }
  }

  sum_o1 /= (2. * M_PI * V * 2. * M_PI * V);
  for (unsigned pp = 0; pp < basis_size; ++pp){
    sum_o1_deriv[pp] = sum_o1_deriv[pp] / (2. * M_PI * V * 2. * M_PI * V);
  }
}

template<typename HatComputer> 
double
IkError<HatComputer>::
estimate (const double & q2,
	  const int & natoms,
	  const SimulationRegion<double> & region) 
{
  prepare_sum (region);
  return sqrt (sum_o1 * q2 * q2 / double(natoms)) * ElectrostaticConvertion;
}

template<typename HatComputer> 
double
IkError<HatComputer>::
estimate (vector<double > & result,
	  const double & q2,
	  const int & natoms,
	  const SimulationRegion<double> & region) 
{
  prepare_sum_deriv (region);

  unsigned basis_size = hat_comput[0].basis_size();
  result.resize (basis_size);
  
  for (unsigned pp = 0; pp < basis_size; ++pp){
    result[pp] = sum_o1_deriv[pp] * 0.5 / sqrt (sum_o1) * sqrt(q2 * q2 / double(natoms)) * ElectrostaticConvertion ;
  }

  return sqrt (sum_o1 * q2 * q2 / double(natoms)) * ElectrostaticConvertion;
}

