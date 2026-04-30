METHOD=RODAS5
# METHOD=RODASPR2
# METHOD=ROS34PW2
# METHOD=GRK4A
# METHOD=GRK4T
# METHOD=ROS3w

# LONG=  #leave empty to use doubles
LONG=long #set to long , in order to use long doubles


rootDir="../"
IncludeDir=$(rootDir)
IncludeHPP = $(shell find $(rootDir) -type f -name '*.hpp')


Path=$(CURDIR)
CC=g++

FLG=  -std=c++20  -I "$(rootDir)"  -lm  -DMETHOD=$(METHOD) -DLONG=$(LONG)  -O3 -Wpedantic

all: Ros.run Robertson.run Example.run

Ros.run: $(IncludeHPP) Ros.cpp makefile
	$(CC) -o "$(Path)/Ros.run" "$(Path)/Ros.cpp" $(FLG) 

Robertson.run: $(IncludeHPP) Robertson.cpp makefile
	$(CC) -o "$(Path)/Robertson.run" "$(Path)/Robertson.cpp" $(FLG) 


Example.run: $(IncludeHPP) Example.cpp makefile
	$(CC) -o "$(Path)/Example.run" "$(Path)/Example.cpp" $(FLG) 

clean:
	rm -f *.run	 