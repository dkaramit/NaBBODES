#ifndef Ros_calc_k
#define Ros_calc_k
#include "Ros_class.hpp"




/*-----------------------Begin: calc_Jk---------------------------------*/
_Ros_template_
_Ros_Func_::calc_Jk(){
    /*
    Calculate product of matrix with vector. 
    We use it to get J*gk.
    */
    for(int i=0; i<N_eqs ; i++){
        Jk[i]=0;
        for(int j=0; j<N_eqs ; j++){ Jk[i]+=J[i][j]*gk[j]; }

        
    }         
}
/*-----------------------End: calc_Jk---------------------------------*/

/*-----------------------Begin: calc_k---------------------------------*/
_Ros_template_
_Ros_Func_::calc_k(){
    // since LU decomposition is not updated as you try to find suitable step, put it ouside the step control loop (should be faster)!
    // LU();
    // calculate k for the other stages
    for(int stage = 0; stage < method.s; stage++){

        
        sum_ak(stage);
        sum_gk(stage);
        
        // setup the argument for dydt (it is evaluated at y_n+\sum a*k )
        for(int eq=0; eq<N_eqs ; eq++){yn[eq]=tmp_sol[eq] +ak[eq]  ; }
        /*--- Get the rhs terms ---*/
        // first term
        dydt(rhs1, yn  , tn+method.c[stage]*h0 );
        // second term
        for(int eq=0; eq<N_eqs ; eq++) { rhs2[eq] =  (h0*h0)*(method.gamma+sum_gamma[stage])*dfdt[eq]   ;  }
        // third term
        calc_Jk();
        // then the rhs becomes
        for(int eq=0; eq<N_eqs ; eq++) { rhs[eq]= rhs1[eq]*h0 + rhs2[eq] + Jk[eq]*h0  ;  }
        /*------------------------*/

        Solve_LU<N_eqs>(L,U,P,rhs, lu_sol);
        for( int eq = 0; eq < N_eqs; eq++ ){ k[eq][stage]=lu_sol[eq]; }
    }







    
          
}
/*-----------------------End: calc_k---------------------------------*/
#endif