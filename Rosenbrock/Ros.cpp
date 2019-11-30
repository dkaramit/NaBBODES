// This is how you run RK. 
#include<iostream>
#include<fstream>
#include<cmath>
#include"Ros.hpp"
#include "Jacobian/Jacobian.hpp"//this is treated as user input, since one may have an analytic form.

/*--------------------------------------------------------------------------------------------------------*/
class ROS3w{
    public:
        const int s=3;
        const int p=2;
        double c[3]={0,2/3.,4/3.};
        double b[3]={0.25,0.25,0.5 };
        double bstar[3]={ 0.746704703274 , 0.1144064078371 , 0.1388888888888 };
        double gamma=0.4358665215084;
        
        double a[3][3] , g[3][3];

        ~ROS3w(){};
        ROS3w(){
            
            for (int i = 0; i < s; i++){ for (int j = 0; j < s; j++){a[i][j]=0; g[i][j]=0; } }
            
            a[1][0]=2/3.;
            a[2][0]=2/3.;
            a[2][1]=2/3.;
            
            g[1][0]=0.3635068368900;
            g[2][0]=-0.8996866791992;
            g[2][1]=-0.1537997822626;
        };
};
/*--------------------------------------------------------------------------------------------------------*/
class ROS34PW2{
    public:
        const int s=4;
        const int p=3;
        double c[4];// c[i]=sum_{j} a[i][j]
        double b[4]={2.4212380706095346e-1 , -1.2232505839045147 , 1.5452602553351020 , 4.3586652150845900e-1};
        double bstar[4]={ 3.7810903145819369e-1 , -9.6042292212423178e-2 , 0.5 , 2.1793326075422950e-1};
        double gamma=4.3586652150845900e-1; 
        
        double a[4][4] , g[4][4];

        ~ROS34PW2(){};
        ROS34PW2(){
            
            for (int i = 0; i < s; i++){c[i]=0; for (int j = 0; j < s; j++){  a[i][j]=0; g[i][j]=0; } }
            
            a[1][0]=8.7173304301691801e-1;
            
            a[2][0]=8.4457060015369423e-1;
            a[2][1]=-1.1299064236484185e-1;
           
            a[3][2]=1.;

            for (int i = 0; i < s; i++){for (int j = 0; j < s; j++){  c[i]+=a[i][j]; } }
           

            g[1][0]=-8.7173304301691801e-1;
            
            g[2][0]=-9.0338057013044082e-1;
            g[2][1]=5.4180672388095326e-2;
            
            g[3][0]=2.4212380706095346e-1;
            g[3][1]=-1.2232505839045147;
            g[3][2]=5.4526025533510214e-1;
           
        };
};
/*--------------------------------------------------------------------------------------------------------*/
class RODASP{
    public:
        const int s=6;
        const int p=4;
        double c[6];// c[i]=sum_{j} a[i][j]
        double b[6]={ -8.0368370789e-2 , -5.6490613592e-2, 4.8828563004e-1 , 5.0571621148e-1 , -1.0714285714e-1 , 0.25};
        double bstar[6]={ -1.7644376488 , -4.7475655721e-1 , 2.3696918469 , 6.1950235906e-1 , 0.25 , 0  };
        double gamma=0.25; 
        
        double a[6][6] , g[6][6];

        ~RODASP(){};
        RODASP(){
            
            for (int i = 0; i < s; i++){c[i]=0; for (int j = 0; j < s; j++){  a[i][j]=0; g[i][j]=0; } }
            a[1][0]=0.75;
            
            a[2][0]= 8.6120400814152190e-2; 
            a[2][1]=0.1238795991858478;
            
            
            a[3][0]= 0.7749345355073236;    
            a[3][1] = 0.1492651549508680; 
            a[3][2] = -0.2941996904581916;
            
            a[4][0]= 5.308746682646142;     
            a[4][1] = 1.330892140037269;  
            a[4][2] = -5.374137811655562; 
            a[4][3]= -0.2655010110278497;

            a[5][0]= -1.764437648774483;    
            a[5][1] =-0.4747565572063027; 
            a[4][2] = 2.369691846915802;
            a[5][3]= 0.6195023590649829;    
            a[5][4] = 0.25;

            for (int i = 0; i < s; i++){for (int j = 0; j < s; j++){  c[i]+=a[i][j]; } }
           

            g[1][0]=-0.75;
            
            g[2][0]=-1.3551200000e-1;
            g[2][1]=-1.3799200000e-1;
            
            g[3][0]=-1.2560800000;
            g[3][1]=-2.5014500000e-1;
            g[3][2]=1.2209300000;
            
            g[4][0]=-7.0731800000;
            g[4][1]=-1.8056500000;
            g[4][2]=7.7438300000;
            g[4][3]=8.8500300000e-1;
            
            g[5][0]=1.6840700000;
            g[5][1]=4.1826600000e-1;
            g[5][2]=-1.8814100000;
            g[5][3]=-1.1378600000e-1;
            g[5][4]=-3.5714300000e-1;
           
           
        };
};
/*--------------------------------------------------------------------------------------------------------*/












// this is how the diffeq should look like
#define n_eqs 3 //number of equations
typedef double Array[n_eqs];//define an array type of length n_eqs
typedef void (*diffeq)(Array &lhs, Array &y  , double t );
//-------------------------------------------------------------------------//

void sys( Array &lhs, Array &y  , double t )
        {
            //lhs is an array that gets the return value (the left hand side of the equation)
            //y is an array with values of y
            //t is the value of the variable t
            
            lhs[0]=-20*y[0]*pow(t,3.);
            lhs[1]=5*y[0]*pow(t,2)+2*(-pow( y[1],2  )+pow( y[2],2 ) )*t;
            lhs[2]=15*y[0]*pow(t,2)+2*(pow( y[1],2  )-pow( y[2],2 ) )*t;

        };

#define initial_step_size 1e-4 
#define minimum_step_size 1e-11 
#define maximum_step_size 1e-2
#define maximum_No_steps 1000000
#define absolute_tolerance 1e-12
#define relative_tolerance 1e-12
#define beta 0.85
#define fac_max 3

#define METHOD ROS34PW2
// #define METHOD ROS3w // this is not very good...
// #define METHOD RODASP // this is slow. I'll have to look at the paper to see what happens

int main(int argc, const char** argv) {
    
    Array lhs;
    Array y0;
    y0[0]=5;
    y0[1]=10;
    y0[2]=0;
    Jacobian<diffeq,n_eqs> jac(sys);

    // this is not looing good.. Just pass the number of equations as argument to get rid of Array and Matrix...
    Ros<diffeq,n_eqs, METHOD ,Jacobian<diffeq,n_eqs> > System(sys,y0, 
     initial_step_size,  minimum_step_size,  maximum_step_size, maximum_No_steps, 
     absolute_tolerance, relative_tolerance, beta, fac_max);
    // System.next_step();
    System.solve();
    // 
    

    #if 0
    for (int i = 0; i < System.max_N; i++)
    {
       System.next_step();

       std::cout<<System.tn<<"\n";

       if(System.tn==1){break;} 
    }
    #endif

    
    #if 0
    for (int i = 0; i < System.current_step; i++)
    {

        std::cout<<i<<"  "<<System.steps[i]<<"\t\t\t\t";
        for(int eq =0; eq<n_eqs; eq++){ std::cout<< System.solution[eq][i]<< "\t\t\t"; }
        std::cout<<"\n";
        
        if(System.steps[i]==1){break;}
    }
    #endif


    std::ofstream f1,f2,f3,t,err;
    f1.open ("./0-test/y1.dat");
    f2.open ("./0-test/y2.dat");
    f3.open ("./0-test/y3.dat");
    t.open ("./0-test/t.dat");
    err.open ("./0-test/err.dat");
    

   
    for (int i = 0; i < System.current_step; i++)
    {
        f1 << System.solution[0][i] ;
        f1 << "\n";
        f2 << System.solution[1][i] ;
        f2 << "\n";
        f3 << System.solution[2][i] ;
        f3 << "\n";
        t << System.steps[i] ;
        t << "\n";
        err << System.err[i] ;
        err << "\n";
        
        // std::cout<<System.steps[i]<<"\t"<< System.solution[0][i] << "\t"<< System.solution[1][i] << "\t"<< System.solution[2][i] <<"\n";
        if(System.steps[i]==1){break;}
    }

    f1.close();
    f2.close();
    f3.close();
    t.close();

    // std::cout<<System.current_step<<std::endl;

    return 0;
 }