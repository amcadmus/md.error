#pragma once

#include <vector>
#include <cassert>
#include <fftw3.h>
#include "SymmHermiteBase.h"

using namespace std;

class HatSymmHermiteBase 
{
public:
  HatSymmHermiteBase (const int & CC,
		      const int & nbins,
		      const int & KK,
		      const int & over_sample);
  void set_value (const vector<double > & vi,
		  const vector<double > & di);
  const double &	 value		(const int & mm) const	{return result_buff[abs(mm)];}
  unsigned		 basis_size	() const		{return 2 * (nbins-1);}
  const vector<double> & basis_value	(const int & mm) const	{return fhvdt[abs(mm)];}
private:
  int CC;
  int nbins;
  int KK;
  int over_sample;
  double hh;
  int MM;
  vector<double > fhv0;			// size MM
  vector<vector<double > > fhv;		// size nbins-1 x MM
  vector<vector<double > > fhd;		// size nbins-1 x MM
  vector<vector<double > > fhvd;	// size 2(nbins - 1) x MM
  vector<vector<double > > fhvdt;	// size MM x 2(nbins-1)
  vector<double > result_buff;		// size MM
};


HatSymmHermiteBase::
HatSymmHermiteBase(const int & CC_,
		   const int & nbins_,
		   const int & KK_,
		   const int & over_sample_)
    : CC (CC_), nbins (nbins_), KK (KK_), over_sample (over_sample_)
{
  hh = CC / double (nbins);
  MM = KK * over_sample;
  double mh = KK / double(MM);
  vector<vector<double> >  ihv(nbins), ihd(nbins);
  SymmHermiteBaseV sbv;
  SymmHermiteBaseD sbd;

  for (int ii = 0; ii < nbins; ++ii){
    ihv[ii].resize (MM, 0);
    ihd[ii].resize (MM, 0);
  }
  // assemble input vectors
  for (int ii = 0; ii < nbins; ++ii){
    for (int jj = 0; jj < MM; ++jj){
      double xx = jj * mh;
      if (xx > KK/2) xx -= KK;
      ihv[ii][jj] = sbv.value (xx, ii, hh);
      ihd[ii][jj] = sbd.value (xx, ii, hh);
    }
    // char file_name[1024];
    // sprintf (file_name, "out.%03d", ii);
    // FILE * fp = fopen (file_name, "w");
    // for (int jj = 0; jj < MM; ++jj){
    //   fprintf (fp, "%f %f %f\n", jj * mh, ihv[ii][jj], ihd[ii][jj]);
    // }
    // fclose (fp);
  }

  fhv0.resize (MM);
  fhv.resize(nbins-1);
  fhd.resize(nbins-1);
  for (int ii = 0; ii < nbins-1; ++ii){
    fhv[ii].resize (MM, 0);
    fhd[ii].resize (MM, 0);
  }

  fftw_complex * in, *out;
  fftw_plan fplan;
  in  = (fftw_complex * ) fftw_malloc (sizeof(fftw_complex) * MM);
  out = (fftw_complex * ) fftw_malloc (sizeof(fftw_complex) * MM);
  for (int ii = 0; ii < MM; ++ii) in[ii][0] = in[ii][1] = 0;

  fplan = fftw_plan_dft_1d (MM, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  for (int ii = 0; ii < nbins; ++ii){
    if (ii == 0){
      for (int jj = 0; jj < MM; ++jj){
	in[jj][0] = ihv[0][jj];
      }
      fftw_execute (fplan);
      for (int jj = 0; jj < MM; ++jj){
	fhv0[jj] = out[jj][0] / double (MM);
      }
    }
    else {
      for (int jj = 0; jj < MM; ++jj){
	in[jj][0] = ihv[ii][jj];
      }
      fftw_execute (fplan);
      for (int jj = 0; jj < MM; ++jj){
	fhv[ii-1][jj] = out[jj][0] / double (MM);
      }
      for (int jj = 0; jj < MM; ++jj){
	in[jj][0] = ihd[ii][jj];
      }
      fftw_execute (fplan);
      for (int jj = 0; jj < MM; ++jj){
	fhd[ii-1][jj] = out[jj][0] / double (MM);
      }
    }
  }

  fftw_destroy_plan (fplan);
  fftw_free (in);
  fftw_free (out);

  fhvd = fhv;
  fhvd.insert (fhvd.end(), fhd.begin(), fhd.end());

  fhvdt.resize (MM);
  for (int jj = 0; jj < MM; ++jj){
    fhvdt[jj].resize(basis_size());
    for (unsigned ii = 0; ii < basis_size(); ++ii){
      fhvdt[jj][ii] = fhvd[ii][jj];
    }
  }
}

void 
HatSymmHermiteBase::
set_value (const vector<double > & vi,
	   const vector<double > & di)
{
  result_buff = fhv0;
  
  assert (vi.size() == fhv.size());
  assert (di.size() == fhd.size());

  for (int jj = 0; jj < MM; ++jj){
    for (int ii = 0; ii < int(vi.size()); ++ii){
      result_buff[jj] += vi[ii] * fhv[ii][jj] + di[ii] * fhd[ii][jj];
    }
  }
}


