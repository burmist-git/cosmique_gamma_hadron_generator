ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX  = g++
CXX += -I./	
CXX += -I./obj/

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)

OUTLIB = ./obj/

#----------------------------------------------------#

all: cosmique_gamma_hadron_generator

.PHONY: printmakehelp_and_reminder
printmakehelp_and_reminder: cosmique_gamma_hadron_generator.C Makefileterzina
	$(info  /******************************************************************************/)
	$(info  * task --> printmakehelp_and_reminder: plots_effective_area.C Makefileterzina *)
	$(info  * $$@ ----> $@                                         *)
	$(info  * $$< --------------------------------> $<                 *)
	$(info  * $$^ --------------------------------> $^ *)
	$(info  /******************************************************************************/)

cosmique_gamma_hadron_generator: cosmique_gamma_hadron_generator.C
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS)

mergeDATfiles: mergeDATfiles.C
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS)

clean:
	rm -rf cosmique_gamma_hadron_generator
	rm -rf mergeDATfiles
	rm -rf *~
	rm -rf src/*~
	rm -rf $(OUTLIB)*.o
