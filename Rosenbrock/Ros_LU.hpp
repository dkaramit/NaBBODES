#ifndef Ros_LU
#define Ros_LU
#include "Ros_class.hpp"



/*--------Calculate the LU decomposition of (1-h*gamma*J) for this step--------------------------*/
Ros_Template
void Ros_Namespace::LU(){
    // initialize the lhs coefficient to unity
    for(int i=0; i<N_eqs ; i++){ for(int j=0; j<N_eqs ; j++){ if(i==j){ _coeff[i][j]=1; }else{_coeff[i][j]=0;}  }}
    // Find the LUP decomposition of   (I-\gamma*h*J)     
    Jac(J,dfdt,tmp_sol,tn);

    for(int i=0; i<N_eqs ; i++){
        for(int j=0; j<N_eqs ; j++){
            _coeff[i][j]+=-h*method.gamma*J[i][j];
        }
    }
    LUP<N_eqs,LD>(_coeff,L,U,P); //LU decomposition of (1-h*gamma*J)
    Inverse_LU<N_eqs,LD>(L,U,P,_inv); // the inverse of (1-h*gamma*J)
}
/*---------------------------------------------------------------------------------------------*/








#endif