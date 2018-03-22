INCDIR = -I ../../GXBERT/include 
INCDIR += -I ../../clhep/include
INCDIR += -I ../../external

CXX = g++
CXXFLAGS += -Wall -std=c++11 -g $(INCDIR)

.SUFFIXES: .cc .o
.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all
all: gxbert 

gxbert: main.o 
	$(CXX) $(CXXFLAGS) -o $@ $? -L../../lib/Linux-g++ -l G4hadronic_bert_cascade -l G4clhep 


.PHONY: lib
lib: libG4hadronic_bert_cascade.a

libG4hadronic_bert_cascade.a:
	cd ..; make;

.PHONY: clean
clean:
	rm -r gxbert main.o
