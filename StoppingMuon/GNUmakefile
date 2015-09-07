#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = MichelReco_StoppingMuonCint
SHLIB = libMichelReco_StoppingMuon.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS  = -I../
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVIDR)
INCFLAGS += -I$(LARLITE_USERDEVIDR)/MichelReco/MichelCluster/Fmwk
INCFLAGS += -I$(LARLITE_USERDEVIDR)/MichelReco/MichelCluster/Algo

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

LDFLAGS += -L$(LARLITE_LIBDIR) -lMichelReco_MichelClusterFmwk
LDFLAGS += -L$(LARLITE_LIBDIR) -lMichelReco_MichelClusterAlgo

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
LDFLAGS += $(shell larlite-config --libs)
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE


