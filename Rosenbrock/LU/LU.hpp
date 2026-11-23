#ifndef LU_Include
#define LU_Include
#include<vector>
#include<cmath>

namespace rosenbrock{

/*---------------------Functions I need for LU decomposition-------------------------------------------------*/
template<class LD>
size_t ind_max(std::vector<LD> &row, size_t len_col){
    /*   
    Find the index of the maximum of a list (row) of lentgth N up to len_col.
    */
    size_t _in=0;
    LD _max = std::abs(row[0]);

    for (size_t  i = 1; i < len_col; i++)
    {
        if(std::abs(row[i])>_max){_max=std::abs(row[i]); _in=i; }
    }
    

    return _in;
}


template<class LD>
void index_swap(std::vector<LD> &A, size_t index_1, size_t index_2){

    /* 
        index_swap takes an array and interchanges 
         A[index_1] with A[index_2].
    */
    LD tmp=A[index_1];
    A[index_1]=A[index_2];
    A[index_2]=tmp;


}

// template<class LD>
// void index_swap(std::vector<int> &A, size_t index_1, size_t index_2){

//     /* 
//         index_swap takes an array and interchanges 
//          A[index_1] with A[index_2].
//     */
//     size_t tmp=A[index_1];
//     A[index_1]=A[index_2];
//     A[index_2]=tmp;


// }

template<class LD>
void apply_permutations_vector(const std::vector<LD> &A, const std::vector<size_t>P, std::vector<LD> &Ap){
    /*
    Applies the permutations given by P from LUP
    to a list A of length N, and stores the result to Ap.
    
    Example:
    If we do this:

    LD A[]={1,2,5,8,3};
    int P[]={2,4,0,3,1};

    LD Ap[5];
    apply_permutations_vector(A,P,5,Ap)

    we get Ap={5,3,1,8,2}
    */
    size_t N=A.size();
    for (size_t i = 0; i < N; i++){Ap[i]=A[ P[i] ];}

}

template<class LD>
void Map( LD (*F)(LD), const std::vector<LD> &L, std::vector<LD> &FL){
    size_t N=L.size();
    FL.resize(N);
    for (size_t i = 0; i < N; i++){ FL[i] = F(L[i]); }
}

/*--------------------------------------------------------------------------------------------------------------*/

/*-----------------------------LUP decompositioning--------------------------------------------------------*/

template<class LD>
void LUP(const std::vector<std::vector<LD>> &M, std::vector<std::vector<LD>> &L ,std::vector<std::vector<LD>> &U, std::vector<size_t> &P, LD _tiny=0){
    size_t N=M.size();
    
    L.assign(N, std::vector<LD>(N, 0));
    U.assign(N, std::vector<LD>(N, 0));
    P.resize(N);

    // Initialize LU
    for (size_t  i = 0; i < N; i++){
        P[i]=i;
        for (size_t  j = 0; j < N; j++){
            if(i==j){L[i][j]=1;}
            if(i!=j){L[i][j]=0;}
            U[i][j]=M[i][j];
        }
    }
    std::vector<LD> _col(N),tmpU(N),tmpL(N);
    size_t len_col,pivot;
    

    for (size_t  k = 1; k < N; k++){
    for (size_t _r=k-1 ; _r<N ; _r++ ) { _col[_r-(k-1)]=std::abs(U[_r][k-1]); }//we need to convert the index of _col because we start the loop from k-1


    len_col=N-(k-1);
    pivot=ind_max<LD>( _col,len_col) + k - 1;
    // std::cout<<pivot<<std::endl;

    if (U[pivot][k-1] == 0 || std::abs(U[pivot][k-1]) < _tiny)  {continue;}

    if (pivot != k-1){ 
            
        index_swap(P,k-1,pivot);
        
        for (size_t _r=k-1 ; _r<N ; _r++ ) { tmpU[_r-(k-1)]= U[k-1][_r] ; }//we need to convert the index of tmpU because we start the loop from k-1

        for (size_t _r=k-1 ; _r<N ; _r++ ) { U[k-1][_r]=U[pivot][_r] ; }
        
        for (size_t _r=k-1 ; _r<N ; _r++ ) { U[pivot][_r]=tmpU[_r-(k-1)] ; }//we need to convert the index of tmpU because we start the loop from k-1

        for (size_t _r=0 ; _r<k-1 ; _r++ ) {tmpL[_r]= L[k-1][_r] ; }
        
        for (size_t _r=0 ; _r<k-1 ; _r++ ) {L[k-1][_r]=L[pivot][_r] ; }
        
        for (size_t _r=0 ; _r<k-1 ; _r++ ) {L[pivot][_r]=tmpL[_r] ; }
    }

    for ( size_t  i = k; i < N; i++){
    L[i][k-1]=U[i][k-1]/U[k-1][k-1];

    for (size_t j=k-1 ; j<N ; j++ ) {  U[i][j]=U[i][j]-L[i][k-1]*U[k-1][j] ; }

    }

    }
    
}
/*-------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------Solve-LU----------------------------------------------------------*/
template<class LD>
void Solve_LU(const std::vector<std::vector<LD>> &L, const std::vector<std::vector<LD>> &U, const std::vector<size_t> &P, const std::vector<LD> &b , std::vector<LD> &x){
    /*
    This solves M*x=b
    Input:
    L,U,P= LUP decomposition of M, which is the output of the function LUP.

    b=the right hand side of the equation
    N=the number of equations

    x=an array to store the solution of M*x=b
    */    

    size_t N=L.size();
    x.resize(N);

    LD tmps=0;

    // b and x must be different vectors. x first holds the permuted b, then d (L*d=b), then the solution.
    for (size_t i = 0; i < N; i++){ x[i]=b[ P[i] ]; }

    for(size_t i=1; i<N  ; i++){
        tmps=0;
        for (size_t j = 0; j < i; j++){ tmps +=L[i][j]*x[j]; }

        x[i]=x[i]-tmps;
    }


    x[N-1]  = x[N-1]/U[N-1][N-1] ;
    // if(N>=1){
        for (size_t i = N-1; i >= 1; i--){
            // std::cout<<v<<"\t"<<N<<"\n";
            tmps=0;
            for (size_t j = i; j < N; j++){ tmps += U[i-1][j]*x[j];}
            x[i-1]=(x[i-1]-tmps )/U[i-1][i-1];
        // }
    }

}
/*-------------------------------------------------------------------------------------------------------------------*/



/*------------------------------------------------Solve-LU----------------------------------------------------------*/
template<class LD>
void Inverse_LU(const std::vector<std::vector<LD>> &L, const std::vector<std::vector<LD>> &U, const std::vector<size_t> &P, std::vector<std::vector<LD>> &invM){
    /*
    Finds the Inverse matrix given its LU decomposition.
    Basically this solves M*M^{-1}=1

    Input:
    L,U,P= LUP decomposition of M, which is the output of the function LUP.

    N=dimension of the matrix (N*N)

    invM=an array to store the solution inverse matrix.
    */    
    size_t N=L.size();
    invM.assign(N, std::vector<LD>(N, 0));
    
    //     
    std::vector<LD> e(N,0);
    // for(size_t i=0 ; i< N ; ++i){ e[i]=0; } 
    std::vector<LD> x(N,0);

    for(size_t i=0 ; i< N ; ++i){
        e[i]=1;
        Solve_LU<LD>(L,U,P,e,x);

        for(size_t j=0 ; j<N ; ++j){
            invM[j][i]=x[j];
        }

        e[i]=0;
    }
}
/*-------------------------------------------------------------------------------------------------------------------*/
 
/*------------------------------------------------Product of two matrices----------------------------------------------------------*/
template<class LD>
void dot(const std::vector<std::vector<LD>> &A, const std::vector<std::vector<LD>> &B, std::vector<std::vector<LD>> &R){
    /*
    Calculates the product of two matrices.
    R=A*B
    */
    size_t N=A.size();
    R.assign(N, std::vector<LD>(N, 0));

    for(size_t i=0; i<N; ++i){
        for(size_t j=0; j<N; ++j){
            R[i][j]=0;
            for(size_t l=0; l<N; ++l){
                R[i][j] += A[i][l]*B[l][j];
            }
        }
    }
}
/*-------------------------------------------------------------------------------------------------------------------*/


/*------------------------------------------------Product of matrix with vector----------------------------------------------------------*/
template<class LD>
void dot(const std::vector<std::vector<LD>> &A, const std::vector<LD> &x, std::vector<LD> &b){
    /*
    Calculates the product of  matrix with vector.
    b=A*x
    */
    size_t N=A.size();
    b.assign(N,0);

    for(size_t i=0; i<N; ++i){
        b[i]=0;
        for(size_t j=0; j<N; ++j){
            b[i] += A[i][j]*x[j];
        }
    }
}
/*-------------------------------------------------------------------------------------------------------------------*/

}

#endif