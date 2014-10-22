#ifndef __nb_kernels_h_wanghan__
#define __nb_kernels_h_wanghan__

#include "nb_interaction_acc_128s.h"
#include "nb_interaction_lj.h"
// #include "nb_interactions.h"

template <int MDDIM,
	  typename ValueType,
	  typename Geometric,
	  typename Acceleration,
	  typename ElectrostaticType,
	  typename VanderWaalsType,
	  typename ComputeMode> 
void
nb_pair_force (const ValueType *		eleParam_,
	       const ValueType *		vdwParam_,
	       const ValueType *		cutoff,
	       const int			ni,
	       const int *			neighbor_index_data,
	       const int *			neighbor_index,
	       const ValueType *		dof_,
	       const ValueType *		charge,
	       const int *			vdwtype,
	       const int			nvdw,
	       ValueType *			force,
	       ValueType *			force_shift,
	       ValueType *			energy)
{
  // nb_pair_force_impl <MDDIM, ValueType,
  // 		      ElectrostaticType,
  // 		      VanderWaalsType,
  // 		      ComputeMode>
  //     (eleParam, vdwParam, cutoff, idx0, idx1, dof, charge, force, energy, fshift, Geometric(), Acceleration());
  VanderWaals_Interaction <Acceleration, VanderWaalsType> vdw_inter;
  
  const typename Acceleration::ValueType * dof,*vdwParam;
  dof = static_cast <const typename Acceleration::ValueType *> (dof_);
  // eleParam = static_cast <const typename Acceleration::ValueType *> (eleParam_);
  vdwParam = static_cast <const typename Acceleration::ValueType *> (vdwParam_);
  
  typename Acceleration::DataType cix, ciy, ciz;
  typename Acceleration::DataType fix, fiy, fiz;

  for (int iindex = 0; iindex < ni; ++iindex){
    load_data_s3_a1<Acceleration> (dof, iindex, &cix, &ciy, &ciz);

    cix = sub<Acceleration> (cix, init<Acceleration>(2.0));
    ciy = sub<Acceleration> (ciy, init<Acceleration>(2.0));
    ciz = sub<Acceleration> (ciz, init<Acceleration>(2.0));
    fix = init<Acceleration>();
    fiy = init<Acceleration>();
    fiz = init<Acceleration>();

    int vdw_ishift = vdwtype[iindex] * nvdw * 2;

    for (int jj = neighbor_index_data[iindex*2]; jj < neighbor_index_data[iindex*2+1]; jj += Acceleration::index_stride){
      typename Acceleration::IndexType jindex;
      copy_index<Acceleration> (jj, neighbor_index, &jindex);

      typename Acceleration::DataType cjx, cjy, cjz;
      typename Acceleration::DataType diffx, diffy, diffz;

      load_data_s3_afull<Acceleration> (dof, jindex, &cjx, &cjy, &cjz);

      diffx = sub<Acceleration> (cix, cjx);
      diffy = sub<Acceleration> (ciy, cjy);
      diffz = sub<Acceleration> (ciz, cjz);

      typename Acceleration::DataType r2, rinv, rinv2;

      r2 = rsq<Acceleration>(diffx, diffy, diffz);
      rinv = invsqrt<Acceleration>(r2);
      rinv2 = mul<Acceleration>(rinv, rinv);

      typename Acceleration::DataType c6, c12;
      typename Acceleration::IndexType param_index;
      index_table_trans<Acceleration> (&jindex, vdwtype, &param_index);
      index_increase<Acceleration> (&param_index, vdw_ishift);
    
      load_data_s2_afull<Acceleration> (vdwParam, param_index, &c6, &c12);
    
      typename Acceleration::DataType fscale;
      fscale = vdw_inter.fscale (rinv2, c6, c12);
    
      typename Acceleration::DataType tfx, tfy, tfz;
      tfx = mul<Acceleration> (fscale, diffx);
      tfy = mul<Acceleration> (fscale, diffy);
      tfz = mul<Acceleration> (fscale, diffz);

      // typename Acceleration::ValueType ofx[4], ofy[4], ofz[4];
      // get<Acceleration> (tfx, ofx);
      // get<Acceleration> (tfy, ofy);
      // get<Acceleration> (tfz, ofz);
      // for (int kk = 0; kk < 4; ++kk){
      //   printf ("%f %f %f   %f\n", ofx[kk], ofy[kk], ofz[kk], fscale[kk]);
      // }

      fix = add<Acceleration> (fix, tfx);
      fiy = add<Acceleration> (fiy, tfy);
      fiz = add<Acceleration> (fiz, tfz);
      
      decrease_data_s3_afull<Acceleration> (force, jindex, tfx, tfy, tfz);
    }

    increase_data_s3_a1<Acceleration> (force, force_shift, iindex, fix, fiy, fiz);
    
  }
}

#endif

