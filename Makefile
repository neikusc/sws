CCsingle = g++
CCFLAGS  = -O3 -Wall -ftemplate-depth-150 
#CCFLAGS  = -O0 -g -Wall
#CCFLAGS  = -pg  -Wall
 
LDFLAGS     = -L/usr/lib
IDFLAGS     = -I. -I/usr/include
LIBS_LAPACK = -lm -llapack -lblas -lg2c

all        : sws 
 
swt        : sws.cpp
	$(CCsingle) $(CCFLAGS) -o $@  swt.cpp $(LDFLAGS) $(IDFLAGS) $(LIBS_LAPACK)

clean:
	rm -f -r swt
