#ifndef __nb_interaction_acc_h_wanghan__
#define __nb_interaction_acc_h_wanghan__

#include "nb_interaction_operators.h"

struct nb_interaction_accelteration_none_s_tag	;
struct nb_interaction_accelteration_none_d_tag	{} ;
struct nb_interaction_accelteration_128s_tag	;
struct nb_interaction_accelteration_128d_tag	{} ;
struct nb_interaction_accelteration_256s_tag	{} ;
struct nb_interaction_accelteration_256d_tag	{} ;

#include "nb_interaction_acc_none_s.h"
#include "nb_interaction_acc_128s.h"

struct nb_interaction_accelteration_sse2_s_tag
    :  public nb_interaction_accelteration_128s_tag {} ;
struct nb_interaction_accelteration_sse4_s_tag
    :  public nb_interaction_accelteration_128s_tag {} ;
struct nb_interaction_accelteration_avx_128_s_tag
    :  public nb_interaction_accelteration_128s_tag {} ;
struct nb_interaction_accelteration_avx_256_s_tag {} ;




#endif
