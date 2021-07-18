for f in $(find . -regextype egrep -regex ".*\.hpp") ;do

        perl -pe 's/^#define.*R.*<.*$//g' $f |\
    perl -pe 's/Ros_Template/template<class diffeq, int N_eqs, class RK_method, class jacobian, class LD> 
/g'  |\
    perl -pe 's/Ros_Namespace/Ros<diffeq, N_eqs, RK_method,  jacobian, LD>/g' > tmp 
    
    cat tmp >$f

done