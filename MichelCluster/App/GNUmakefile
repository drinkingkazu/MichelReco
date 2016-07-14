#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = MichelReco_MichelClusterAppCint
SHLIB = libMichelReco_MichelClusterApp.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += -I../
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell recotool-config --includes)
INCFLAGS += $(shell basictool-config --includes)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += -L$(LARLITE_LIBDIR) -lMichelReco_MichelClusterFmwk
LDFLAGS += -L$(LARLITE_LIBDIR) -lRecoTool_MCComp
LDFLAGS += -L$(LARLITE_LIBDIR) -lBasicTool_GeoAlgo
LDFLAGS += $(shell root-config --libs)
#LDFLAGS += -lSpectrum
# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE

