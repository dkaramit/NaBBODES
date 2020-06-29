#ifndef RKF_METHOD
#define RKF_METHOD


template<class LD>
class DormandPrince{
    public:
        const int s=7;
        int p=5;
        LD c[7]={0,1/5.,3/10.,4/5.,8/9.,1.,1.};
        LD b[7]={5179/57600.,0,7571/16695.,393/640.,-92097/339200.,187/2100.,1/40.};
        LD bstar[7]={ 35/384.,0.,500/1113.,125/192.,-2187/6784.,11/84.,0 };
        LD a[7][7];

        ~DormandPrince(){};
        DormandPrince(){
            
            for (int i = 0; i < s; i++){ for (int j = 0; j < s; j++){a[i][j]=0;} }
            
            a[1][0]=1/5.;
            
            a[2][0]=3/40.;
            a[2][1]=9/40.;
            
            a[3][0]=44/45.;
            a[3][1]=-56/15.;
            a[3][2]=32/9.;
            
            a[4][0]=19372/6561.;
            a[4][1]=-25360/2187.;
            a[4][2]=64448/6561.;
            a[4][3]=-212/729. ;       
            

            a[5][0]=9017/3168.;
            a[5][1]=-355/33.;
            a[5][2]=46732/5247.;
            a[5][3]=49/176.;
            a[5][4]=-5103/18656.;
            
            
            a[6][0]=35/384.;
            a[6][1]=0;
            a[6][2]=500/1113.;
            a[6][3]=125/192.;
            a[6][4]=-2187/6784.;
            a[6][5]=11/84.;
        };
};



template<class LD>
class CashKarp{
    public:
        const int s=6;
        int p=5;
        LD c[6]={0,1/5.,3/10.,3/5.,1.,7/8.};
        LD b[6]={37/378.,0,250/621.,125/594.,0,512/1771.};
        LD bstar[6]={2825/27648.,0.,18575/48384.,13525/55296,277/14336.,0.25};
        LD a[6][6];

        ~CashKarp(){};
        CashKarp(){
            
            for (int i = 0; i < s; i++){ for (int j = 0; j < s; j++){a[i][j]=0;} }
            
            a[1][0]=0.2;
            
            a[2][0]=3/40.;
            a[2][1]=9/40.;
            
            a[3][0]=0.3;
            a[3][1]=-0.9;
            a[3][2]=6/5.;
            
            a[4][0]=-11/54.;
            a[4][1]=2.5;
            a[4][2]=-70/27;
            a[4][3]=35/27. ;       
            

            a[5][0]=1631/55296.;
            a[5][1]=175/512.;
            a[5][2]=575/13824.;
            a[5][3]=44275/110592.;
            a[5][4]=253/4096.;
            
        };
};



template<class LD>
class RKF45{
    public:
        const int s=6;
        int p=5;
        LD c[6]={0,1/4.,3/8.,12/13.,1.,0.5};
        LD b[6]={16/135.,0,6656/12825.,28561/56430.,-9/50.,2/55.};
        LD bstar[6]={25/216.,0.,1408/2565.,2197/4104,-0.2,0};
        LD a[6][6];

        ~RKF45(){};
        RKF45(){
            
            for (int i = 0; i < s; i++){ for (int j = 0; j < s; j++){a[i][j]=0;} }
            
            a[1][0]=0.25;
            
            a[2][0]=3/32.;
            a[2][1]=9/32.;
            
            a[3][0]=1932/2197.;
            a[3][1]=-7200/2197.;
            a[3][2]=7296/2197.;
            
            a[4][0]=439/216.;
            a[4][1]=-8.;
            a[4][2]=3680/513.;
            a[4][3]=-845/4104.;       
            

            a[5][0]=-8/27.;
            a[5][1]=2.;
            a[5][2]=-3544/2565.;
            a[5][3]=1859/4104.;
            a[5][4]=-11/40.;
            
        };
};
#endif