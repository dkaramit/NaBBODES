#set PATH variable to working directory


# METHOD=ROS34PW2
METHOD=RODASPR2
# METHOD=ROS3w


# LONG=  #leave empty to use doubles
LONG=long #set to long , in order to use long doubles

 

Path=$(PWD)
CC=g++
FLG=  -std=c++17  -I "$(Path)" -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3


all: Ros.run

Ros.run: makefile METHOD.hpp  Ros_calc_k.hpp  Ros_class.hpp  Ros_costructor.hpp  Ros.hpp  Ros_step_control-PI.hpp Ros_step_control-simple.hpp  Ros_steps.hpp  Ros_LU.hpp Ros_sums.hpp  Ros.cpp LU/LU.hpp Jacobian/Jacobian.hpp
	$(CC) -o "$(Path)/Ros.run" "$(Path)/Ros.cpp" $(FLG) 


clean:
	@[ -f "$(Path)/Ros.run" ] && rm "$(Path)/Ros.run" || true
	 