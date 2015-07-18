# makefile for tagging performance histogram filler
# Author: Dan Guest (dguest@cern.ch)
# created: Thu Dec 19 17:34:32 EST 2013

# --- set dirs
BUILD          := build
SRC          := src
INC          := include
DICT         := dict
OUTPUT       := bin
DELPHES      := delphes

# --- make sure we have delphes
LIB_DELPHES := $(DELPHES)/libDelphes.so
ifeq (,$(wildcard $(LIB_DELPHES) ))
$(error "couldn't find '$(DELPHES)' folder, make sure you've linked it")
endif
DELPHES_ABSOLUTE := $(abspath $(DELPHES))

# --- HACKS ----
# CXXFLAG_HACKS := -Wno-literal-suffix #hdf5 header sets this off

#  set search path
vpath %.cxx  $(SRC)
vpath %.hh   $(INC)
vpath %.h    $(INC)
vpath %Dict.h $(DICT)
vpath %Dict.cxx $(DICT)
# vpath %.h $(DELPHES)
# vpath %.h $(DELPHES)/external

# --- hdf and ndhist
HDF_INFO := $(shell h5c++ -showconfig | grep 'Installation point:')
HDF_PATH := $(strip $(shell echo ${HDF_INFO} | cut -d ':' -f 2 ))
ifndef HDF_PATH
$(error "couldn't find HDF, quitting")
endif

# ND_HIST_DIR      := $(CURDIR)/ndhist
# ND_HIST_INC      := $(ND_HIST_DIR)/include
# ND_HIST_LIB      := $(ND_HIST_DIR)/lib

# --- load in root config
ROOTCFLAGS    := $(shell root-config --cflags)
# would be nice to avoid linking everything, but that will probably cause
# problems...
ROOTLIBS      := $(shell root-config --libs)
# ROOTLIBS      := -L$(shell root-config --libdir)
ROOTLIBS      += -lCore -lTree -lRIO
ROOTLIBS      += -lCint		# don't know why we need this...
ROOTLDFLAGS   := $(shell root-config --ldflags)

# --- set compiler and flags (roll c options and include paths together)
CXX          ?= g++
CXXFLAGS     := -O2 -Wall -fPIC -I$(INC) -I$(ND_HIST_INC) -g -std=c++11
CXXFLAGS     += $(CXXFLAG_HACKS)
CXXFLAGS     += -I$(DELPHES) -I$(DELPHES)/external
# LIBS         := -L$(ND_HIST_LIB) -Wl,-rpath,$(ND_HIST_LIB) -lndhist
LIBS         := $(shell ndhist-config --libs)
CXXFLAGS     += $(shell ndhist-config --cflags)

LDFLAGS      := #-Wl,--no-undefined

CXXFLAGS     += -I$(HDF_PATH)/include
LIBS         += -L$(HDF_PATH)/lib -Wl,-rpath,$(HDF_PATH)/lib
LIBS         += -L$(DELPHES_ABSOLUTE) -Wl,-rpath,$(DELPHES_ABSOLUTE)
LIBS         += -lDelphes

# --- HDF5 needed for hist saving
LIBS         += -lhdf5_cpp -lhdf5

# --- rootstuff
CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS         += $(ROOTLIBS)

# ---- define objects
GEN_OBJ     := ExRootTreeReader.o misc_func.o truth_tools.o AllPlanes.o
TOP_OBJ     += delphes-tracking-plots-build.o
TOP_OBJ     += delphes-vertex-plots-build.o

# stuff used for the c++ executable
ALL_EXE    := delphes-tracking-plots-build delphes-vertex-plots-build

GEN_OBJ_PATHS := $(GEN_OBJ:%=$(BUILD)/%)
TOP_OBJ_PATHS := $(TOP_OBJ:%=$(BUILD)/%)
ALL_EXE_PATHS := $(ALL_EXE:%=$(OUTPUT)/%)

all: $(ALL_EXE_PATHS)

$(OUTPUT)/delphes-%: $(GEN_OBJ_PATHS) $(BUILD)/delphes-%.o
	@mkdir -p $(OUTPUT)
	@echo "linking $^ --> $@"
	@$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)


# --------------------------------------------------

# compile rule
$(BUILD)/%.o: %.cxx
	@echo compiling $<
	@mkdir -p $(BUILD)
	@$(CXX) -c $(CXXFLAGS) $< -o $@

# use auto dependency generation
ALLOBJ       := $(GEN_OBJ)
DEP = $(BUILD)

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),rmdep)
include  $(ALLOBJ:%.o=$(DEP)/%.d)
endif
endif

DEPTARGSTR = -MT $(BUILD)/$*.o -MT $(DEP)/$*.d
$(DEP)/%.d: %.cxx
	@echo making dependencies for $<
	@mkdir -p $(DEP)
	@$(CXX) -MM -MP $(DEPTARGSTR) $(CXXFLAGS) $(PY_FLAGS) $< -o $@

# clean
.PHONY : clean rmdep all
CLEANLIST     = *~ *.o *.o~ *.d core
clean:
	rm -fr $(CLEANLIST) $(CLEANLIST:%=$(BUILD)/%) $(CLEANLIST:%=$(DEP)/%)
	rm -fr $(BUILD) $(DICT)

rmdep:
	rm -f $(DEP)/*.d
