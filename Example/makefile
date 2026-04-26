METHOD=RODAS5
# METHOD=ROS34PW2
# METHOD=RODASPR2
# METHOD=ROS3w


# LONG=  #leave empty to use doubles
LONG=long #set to long , in order to use long doubles

 

Path=$(CURDIR)
PathRos=$(abspath $(Path)/..)


CC=g++
FLG=  -std=c++20  -I "$(PathRos)/Rosenbrock/Jacobian"  -I "$(Path)" -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3 -Wpedantic


all: TwoState.run

TwoState.run: makefile TwoState.cpp
	$(CC) -o "$(Path)/TwoState.run" "$(Path)/TwoState.cpp" $(FLG) 


clean:
	@[ -f "$(Path)/TwoState.run" ] && rm "$(Path)/TwoState.run" || true
	 