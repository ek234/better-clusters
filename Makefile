srcfile := 1_true_clusters.cpp

CC=g++
objects = $(patsubst %.c,%.o,$(wildcard *.c))

DELPHES_DIR := /opt/madgraph/installation/MG5_aMC_v3_5_0/Delphes
FASTJET_DIR := /opt/madgraph/installation/MG5_aMC_v3_5_0/HEPTools/fastjet
PYTHIA8_DIR := /opt/madgraph/installation/MG5_aMC_v3_5_0/HEPTools/pythia8

PLATFORM_TYPE := $(shell uname -s)
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTINC := $(shell root-config --incdir)

INC_FJ := `$(FASTJET_DIR)/bin/fastjet-config --cxxflags`
LIB_FJ := `$(FASTJET_DIR)/bin/fastjet-config --libs`

INC_DP := -I$(DELPHES_DIR) -I$(DELPHES_DIR)/classes -I$(DELPHES_DIR)/external
LIB_DP := -L$(DELPHES_DIR) -lDelphes -Wl,-rpath=$(DELPHES_DIR)

all: Test_project

Test_project: $(objects) $(srcfile)
	$(CXX) -O3 -g $(objects) $(srcfile) $(ROOTFLAGS) $(INC_DP) $(INC_FJ) $(ROOTLIBS) $(LIB_DP) $(LIB_FJ) -lfastjet -lNsubjettiness -Wl,-rpath=$(DELPHES_DIR) -o $(srcfile).ex

%.o: %.C %.h
	rm  lib/*; gcc -fPIC -O2 -c $(INC_FJ) $< -o $@

clean:
	rm $(FILES) Test_project lib/*
