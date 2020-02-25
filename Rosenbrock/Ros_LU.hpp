#ifndef Ros_LU
#define Ros_LU
#include "Ros_class.hpp"
#include "Jacobian/Jacobian.hpp"
#include "LU/LU.hpp"



/*--------Calculate the LU decomposition of (1-h*gamma*J) for this step--------------------------*/
_Ros_template_
_Ros_Func_::LU(){
    // initialize the lhs coefficient to unity
    for(int i=0; i<N_eqs ; i++){ for(int j=0; j<N_eqs ; j++){  if(i==j){ _coeff[i][j]=1; }else{_coeff[i][j]=0;}  }}

    // Find the LUP decomposition of   (I-\gamma*h*J)     
    for(int i=0; i<N_eqs ; i++){yn[i]=tmp_sol[i] ; }
    Jac(J,dfdt,yn,tn);

    for(int i=0; i<N_eqs ; i++){
        for(int j=0; j<N_eqs ; j++){
            _coeff[i][j]+=-h0*method.gamma*J[i][j];
    

        }
    }

    LUP<N_eqs,LD>(_coeff,L,U,P); 

}
/*---------------------------------------------------------------------------------------------*/








#endif