# Compiler and flags
CC= g++

# NETCDF library directories
NETCDF_INCLUDEDIR=/opt/homebrew/include
NETCDF_LIBDIR=/opt/homebrew/lib

# Library files to include
LDFILES= -lnetcdf

# Local files
FILES= SQuadGen.cpp \
	   GridElements.cpp \
	   CSRefinementMap.cpp \
       CubedSphereGrid.cpp \
	   IcosahedralFlagGrid.cpp \
	   RefineGrid.cpp \
       RefinementTemplateCUBIT.cpp \
	   RefinementTemplateLOWCONN.cpp \
	   RefinementTemplateLOWCONNOLD.cpp \
       SpringDynamics.cpp \
	   MeshUtilities.cpp \
       Tessellate.cpp\
       lodepng.cpp \
	   netcdf.cpp \
	   ncvalues.cpp

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR) -fsanitize=address
LDFLAGS+= -L$(NETCDF_LIBDIR) -fsanitize=address

include Make.defs

##
## Build instructions
##
SQuadGen: $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)
	mv SQuadGen ..

##
## Clean
##
clean:
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

