#
INCDIR += -I /home/tkoi/BERT/standalone/rev0/geant4-install/include/Geant4

CXX = g++
CXXFLAGS += -Wall -std=c++11 -g $(INCDIR)

.SUFFIXES: .cc .o
.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all
all: g4bert 

g4bert: main.o 
	$(CXX) $(CXXFLAGS) -o $@ $? -L/home/tkoi/BERT/standalone/rev0/geant4-install/lib64 -l G4processes -l G4particles -l G4global -l G4clhep 

.PHONY: clean
clean:
	rm -r g4bert main.o
