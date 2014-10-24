#ifndef __nb_kernels_h_wanghan__
#define __nb_kernels_h_wanghan__

#include "global_defines.h"

#include "nb_interaction_acc_128s.h"
#include "nb_interaction_vdw_cutoff.h"
// #include "nb_interactions.h"

template <int MDDIM,
	  typename ValueType,
	  typename Geometric,
	  typename Acceleration,
	  typename ElectrostaticType,
	  typename VanderWaalsType,
	  typename ComputeMode> 
void
nb_pair_force (const int					ni,
	       const int *					neighbor_index_data,
	       const int *					neighbor_index,
	       const VanderWaals_Parameter<Acceleration> 	vdw_param,
	       const Electrostatic_Control<Acceleration> 	ele_ctrl,
	       const VanderWaals_Control<Acceleration>		vdw_ctrl,
	       const typename Acceleration::ValueType *		dof,
	       const typename Acceleration::ValueType *		charge,
	       const int *					vdwtype,
	       typename Acceleration::ValueType *		force,
	       typename Acceleration::ValueType *		force_shift,
	       typename Acceleration::ValueType *		total_energy)
{
  VanderWaals_Interaction <Acceleration, VanderWaalsType, ComputeMode> vdw_inter (vdw_ctrl);
  
  // init internal vdw parameters
  const typename Acceleration::ValueType * vdwParam (vdw_param.data);
  const int nvdw (vdw_param.nvdw);

  typename Acceleration::DataType fscale, energy;
  energy = fscale = init<Acceleration>();

  for (int iindex = 0; iindex < ni; ++iindex){
    typename Acceleration::DataType cix, ciy, ciz;
    load_data_s3_a1<Acceleration> (dof, iindex, &cix, &ciy, &ciz);
    cix = sub<Acceleration> (cix, init<Acceleration>(2.0));
    ciy = sub<Acceleration> (ciy, init<Acceleration>(2.0));
    ciz = sub<Acceleration> (ciz, init<Acceleration>(2.0));

    typename Acceleration::DataType fix, fiy, fiz;
    fix = init<Acceleration>();
    fiy = init<Acceleration>();
    fiz = init<Acceleration>();

    int vdw_ishift = vdwtype[iindex] * nvdw;

    int jj;
    for (jj = neighbor_index_data[iindex*2];
	 jj < neighbor_index_data[iindex*2+1] && neighbor_index[jj+Acceleration::index_stride-1] >= 0;	// should be valid neighbor
	 jj += Acceleration::index_stride){
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

      {
	typename Acceleration::DataType c6, c12;
	typename Acceleration::IndexType param_index;
	index_table_trans<Acceleration> (jindex, vdwtype, &param_index);
	index_increase<Acceleration> (&param_index, vdw_ishift);
    
	load_data_s2_afull<Acceleration> (vdwParam, param_index, &c6, &c12);
    
	vdw_inter.force_energy (rinv2, c6, c12, &fscale, &energy);
      }
    
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
    
    if (jj < neighbor_index_data[iindex*2+1]){				// deal with phatom neighbor
      typename Acceleration::IndexType jindex;
      copy_index_capped<Acceleration> (jj, neighbor_index, &jindex);
      typename Acceleration::DataType mask;
      mask = cal_mask<Acceleration> (jj, neighbor_index);
      
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

      {
	typename Acceleration::DataType c6, c12;
	typename Acceleration::IndexType param_index;
	index_table_trans<Acceleration> (jindex, vdwtype, &param_index);
	index_increase<Acceleration> (&param_index, vdw_ishift);
    
	load_data_s2_afull<Acceleration> (vdwParam, param_index, &c6, &c12);
    
	vdw_inter.force_energy (rinv2, c6, c12, &fscale, &energy);
	fscale = apply_mask<Acceleration> (mask,fscale);
	energy = apply_mask<Acceleration> (mask,energy);
      }
    
      typename Acceleration::DataType tfx, tfy, tfz;
      tfx = mul<Acceleration> (fscale, diffx);
      tfy = mul<Acceleration> (fscale, diffy);
      tfz = mul<Acceleration> (fscale, diffz);

      fix = add<Acceleration> (fix, tfx);
      fiy = add<Acceleration> (fiy, tfy);
      fiz = add<Acceleration> (fiz, tfz);
      
      decrease_data_s3_afull<Acceleration> (force, jindex, tfx, tfy, tfz);
    }
    
    increase_data_s3_a1<Acceleration> (force, force_shift, iindex, fix, fiy, fiz);
  }
}

#endif

