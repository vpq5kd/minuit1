# This Makefile builds two stand alone programs that have access to 
# the ROOT libraries.  To use this with your own code just substitute
# the name of your program below.

# here we access the root configuration, include files, and libraries
ROOTCFLAGS=$(shell root-config --cflags)
ROOTINC=$(shell root-config --incdir)
ROOTLIBDIR=$(shell root-config --libdir)
ROOTLIBS=$(shell root-config --libs) -lMinuit
ROOTLDFLAGS=$(shell root-config --ldflags)

ROOTC=$(ROOTCFLAGS) 
#-I$(ROOTINC)
ROOTLINK=-L$(ROOTLIBDIR) $(ROOTLIBS) $(ROOTLDFLAGS)

CPP=g++

default: expFit ex1 rootExample ex3

expFit: expFit.cpp
	$(CPP) -O -Wall $(ROOTC) -o expFit expFit.cpp $(ROOTLINK) 
ex1: ex1.cpp
	$(CPP) -O -Wall $(ROOTC) -o ex1 ex1.cpp $(ROOTLINK) 

ex3: ex3.cpp
	$(CPP) -O -Wall $(ROOTC) -o ex3 ex3.cpp $(ROOTLINK) 


rootExample: rootExample.cpp
	$(CPP) -O -Wall $(ROOTC) -o rootExample rootExample.cpp $(ROOTLINK) 
# note: just replace the -O flag with -g to build a debug version

clean: 
	rm -f expFit rootExample *~ *.d *.so
