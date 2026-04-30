#ifndef Ros_sums
#define Ros_sums
#include "Ros_class.hpp"

namespace rosenbrock{

/*-----------------------Begin: sum_ak---------------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::sum_ak(const unsigned int& stage){
    // this function stores sum_{j}^{stage-1}a_{stage,j}\vec{k}_j in ak, so we first need to make all elements zero, 
    // and then take the sum for each component
    for (unsigned int eq = 0; eq <N_eqs ; eq++){
        ak[eq]=0.; 
        for (unsigned int j = 0; j < stage; j++){ ak[eq]+=RK_method::a[stage][j]*k[eq][j];  }
    }
    
}
/*-----------------------End: sum_ak---------------------------------*/


/*-----------------------Begin: sum_gk---------------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::sum_gk(const unsigned int& stage){
    // this function stores sum_{j}^{stage-1}g_{stage,j}\vec{k}_j in ak, so we first need to make all elements zero, and then take the sum for each component
    for (unsigned int eq = 0; eq <N_eqs ; eq++){
        gk[eq]=0.;  
        for (unsigned int j = 0; j < stage; j++){ gk[eq]+=RK_method::g[stage][j]*k[eq][j];  }
    }
    
}
/*-----------------------End: sum_ak---------------------------------*/

/*-----------------------Begin: sum_bk---------------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::sum_bk(){
    // this function stores sum_{i}^{s}b_{i}\vec{k}_i in bk and sum_{i}^{s}b_{i}^{\star}\vec{k}_i in bstark  
    for (unsigned int eq = 0; eq <N_eqs ; eq++){
        bk[eq]=0.;
        bstark[eq]=0.; 
        
        for (unsigned int i = 0; i < RK_method::s; i++){ 
            bk[eq]+=RK_method::b[i]*k[eq][i];  
            bstark[eq]+=RK_method::bstar[i]*k[eq][i];  
            }
        // std::cout<<bk[eq]<<"   "<<bstark[eq]<<"   "<< h<<"\n";std::cin.get();
    }

}
/*-----------------------End: sum_bk---------------------------------*/

}

#endif