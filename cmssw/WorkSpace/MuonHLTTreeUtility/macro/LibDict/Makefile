ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)

CXX           = g++
CXXFLAGS      = -O -Wall -fPIC $(ROOTCFLAGS)
LD            = g++
LDFLAGS       = -O $(ROOTLDFLAGS)
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

libTreeLib.so:  MyLibLinkDef.h
		@rootcint -f MyLibDict.cc -c map.h vector.h MyLibLinkDef.h
		$(CXX) $(LDFLAGS) $(SOFLAGS) $(CXXFLAGS) MyLibDict.cc -o $@ 

