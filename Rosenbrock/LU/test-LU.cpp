#include<iostream>
#include<array>
#include<cmath>
#include "LU.hpp"

using std::cout;
using std::endl;
using namespace rosenbrock;

#define LD double

// #define indmax // run ind_max test

// #define swap //run index_swap test

// #define perm //run permutation test

// #define lup //run LUP test

#define _rand //run random tests of Solve_LU

// #define inv_test // random tests of Inverse_LU




#ifdef _rand
    #include <cstdlib>
#endif

int main(){


    #ifdef indmax
        std::vector<LD> x(5,0);
        x={2,8,-1,2,50};
        cout<<ind_max<LD>(x,3)<<endl;
    #endif


    
    #ifdef swap
        std::vector<LD> x(5,0);
        x={2,1,-1,2,50};

        index_swap<LD>(x,4,1);
        for( LD i : x ){ cout<<i<<endl;}
    #endif
    

    #ifdef perm
        const unsigned int N=5;
        std::vector<LD> A(N,0); 
        A={1,2,5,8,3};
        std::vector<int>  P(N,0);
        P={2,4,0,3,1};

        std::vector<LD> Ap(N,0);

        apply_permutations_vector<LD>(A,P,Ap);
        for( int i =0 ; i<5 ; i++){ cout<<A[i]<<" ";}
        cout<<endl;
        for( int i =0 ; i<5 ; i++){ cout<<Ap[i]<<" ";}
        cout<<endl;
    #endif




    #ifdef lup
    const unsigned int N=5;
    std::vector<std::vector<LD>> M(N,std::vector<LD>(N,0)),L(N,std::vector<LD>(N,0)),U(N,std::vector<LD>(N,0));
    M={{
    { 0,  2,  2 , 3 , 5},
    {-3, -1,  1 , 5 , 9},
    { 1, -1,  1 , 4 , 7},
    { 1, -1,  1 , 0 , 2},
    { 1, -1,  1 , 0 , 3}
    }};


    std::vector<int> P(N,0);

    LUP<LD>(M,L,U,P);

    for( LD i : P ){ cout<< i<<' ';}
    cout<<endl;
    cout<<endl;
    cout<<endl;

    for (int i = 0; i < N; i++) {
    
    for (int j = 0; j < N; j++) {
        cout<< U[i][j]<<"\t";

    }
    cout<<endl;
    }
    cout<<endl;
    cout<<endl;

    for (int i = 0; i < N; i++) {
    
    for (int j = 0; j < N; j++) {
        cout<< L[i][j]<<"\t";

    }
    cout<<endl;
    }
    #endif


    #ifdef _rand
        /* 
            Random tests for Solve_LU. Basically run "runs" tests of Solve_LU with N number 
            of equations (with random coefficients of magnitude up to 100), and if 
            (M*x-b)/b > 1e-11, print it.
        */    
        const unsigned int runs=100000;
        std::vector<LD> err(runs,0);


        const unsigned int N=10;
        std::vector<std::vector<LD>> M(N,std::vector<LD>(N,0)),U(N,std::vector<LD>(N,0)),L(N,std::vector<LD>(N,0));
        std::vector<LD> b(N,0),x(N,0);
        std::vector<int> P(N,0);
        
        std::vector<LD> tmp(N,0);
        
        LD tmpE;

        for(unsigned int _r=0; _r<runs ; _r++){
            for (unsigned int i = 0; i < N; i++) { for (int j = 0; j < N; j++)  {
                M[i][j]= ( rand()/ ((LD) RAND_MAX ) -0.5 ) *100 ;  } 
                b[i]= (rand()/ ((LD) RAND_MAX ) -0.5 ) *100 ;  
            } 
        
            LUP<LD>(M,L,U,P);
            Solve_LU<LD>(L,U,P,b,x);

            err[_r]=0;
            for (int i = 0; i < N; i++){
                dot<LD>(M,x,tmp);
                tmpE= std::abs((tmp[i] - b[i])/b[i]) ;
                if(tmpE>err[_r] ) {err[_r] = tmpE ;}

            }
        
        
        }

        for(LD _err: err){  
            if(_err>1e-8){ cout<<_err<<endl;} 
        }
        

    #endif



    #ifdef inv_test 
    /*
    run tests for Inverse_LU. Basically run "runs" random tests of Inverse_LU of N*N matrix 
    (with random entries of magnitude up to 100), and if 
    (M*M^{-1}-1)>1e-12, print it.
    */
    
    
    const unsigned int N=10;
    std::vector<std::vector<LD>> M(N,std::vector<LD>(N,0)),L(N,std::vector<LD>(N,0)),U(N,std::vector<LD>(N,0)),invM(N,std::vector<LD>(N,0)),R(N,std::vector<LD>(N,0));

    std::vector<std::vector<LD>> Unit(N,std::vector<LD>(N,0));
    for (int i = 0; i < N; i++){Unit[i][i]=1;} 


    std::vector<int> P(N,0);
    LD tmp;
    
    for(int _run=0 ; _run<1000; ++_run)
    {
        for (int i = 0; i < N; i++) 
        { 
            // Unit is initialized as zero matrix. So put 1 at Unit[i][i].
            Unit[i][i]=1;
            for (int j = 0; j < N; j++)  
            {        
                // fil M with random numbers
                M[i][j]= ( rand()/ ((LD) RAND_MAX ) -0.5 ) *100 ;  

            }
        } 

        // LUP decomposition of M
        LUP<LD>(M,L,U,P);
        // Given LUP decomposition you can calculate the inverse. 
        Inverse_LU<LD>(L,U,P,invM);

        // calculate M*M^{-1}
        dot<LD>(M,invM,R);


        // print if M*M^{-1} - 1 > 10^{-10}
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                tmp=fabs(R[i][j]-Unit[i][j]);
                if(tmp>1e-8){ cout<< tmp << endl;}
            }
        }
    }
    #endif



    return 0;
}