#ifndef Ros_LU
#define Ros_LU
#include "Ros_class.hpp"

namespace rosenbrock{

/*--------Calculate the LU decomposition of (1-h*gamma*J) for this step--------------------------*/
template<class LD, class RK_method, step_controllers step_controller> 
void Solver<LD, RK_method, step_controller>::LU(){
    //initialize coefficient to 0
    std::vector<std::vector<LD>> coeff(N_eqs,std::vector<LD>(N_eqs,0));
   
    // Jac(J,dfdt,yprev,tn);

    // Find the LUP decomposition of   (I-\gamma*h*J)     
    for(unsigned int i=0; i<N_eqs ; i++){
        coeff[i][i]=1;
        for(unsigned int j=0; j<N_eqs ; j++){
            coeff[i][j]+=-h_trial*RK_method::gamma*J[i][j];
        }
    }
    LUP<LD>(coeff,L,U,P,_tiny<LD>); //LU decomposition of (1-h*gamma*J)
    Inverse_LU<LD>(L,U,P,_inv); // the inverse of (1-h*gamma*J)
}
/*---------------------------------------------------------------------------------------------*/

}

#endif
