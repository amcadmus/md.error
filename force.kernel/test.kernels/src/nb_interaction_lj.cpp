#include "nb_interaction_lj.h"
#include "nb_interaction_acc.h"

template <>
typename nb_interaction_accelteration_128s_tag::DataType
LennardJones_6_12<nb_interaction_accelteration_128s_tag>::one_over_12 =
    init<nb_interaction_accelteration_128s_tag> (nb_interaction_accelteration_128s_tag::ValueType(1./12.));

template <>
typename nb_interaction_accelteration_128s_tag::DataType
LennardJones_6_12<nb_interaction_accelteration_128s_tag>::one_over_6 =
    init<nb_interaction_accelteration_128s_tag> (nb_interaction_accelteration_128s_tag::ValueType(1./6.));

