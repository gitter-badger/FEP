BDIR = ./bin
SDIR = ./src
IDIR = ./inc
LDIR = ./lib
ODIR = ./obj

export OMPI_CXX=clang++

#===============================
#	  XERCES C++ XML PARSING 
#===============================

XERCES_DIR  = $(LDIR)/xerces-c-3.1.2/build
XERCES_INC = -I $(XERCES_DIR)/include/

XERCES_LIBS += -L $(XERCES_DIR)/lib

#===============================
#			PUMI 
#===============================

PUMI_DIR = ~/Documents/GitHub/core
PUMI_INC = -I $(PUMI_DIR)/apf -I $(PUMI_DIR)/gmi -I $(PUMI_DIR)/mds \
			-I $(PUMI_DIR)/pcu

PUMI_LIBS += -L $(PUMI_DIR)/build/apf -L $(PUMI_DIR)/build/gmi -L \
		$(PUMI_DIR)/build/mds -L $(PUMI_DIR)/build/pcu -L $(LDIR)

#===============================
#		GOOGLE TEST 
#===============================

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR =  ~/Documents/GitHub/googletest

# Where to find user code.
USER_DIR = $(GTEST_DIR)/samples

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include -isystem $(GTEST_DIR)

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread --pedantic-errors

MPIFLAGS += -Wall -g -I $(IDIR) --pedantic-errors

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = $(BDIR)/a4

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
				$(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

all : $(TESTS)

#there are other things in lib, so only remove libraries explictly by name
clean :
	rm -vf $(TESTS) $(BDIR)/* $(LDIR)/gtest.a $(LDIR)/gtest_main.a $(ODIR)/*.o

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(ODIR)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
			$(GTEST_DIR)/src/gtest-all.cc -o $@

$(ODIR)/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
			$(GTEST_DIR)/src/gtest_main.cc -o $@

$(LDIR)/gtest.a : $(ODIR)/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(LDIR)/gtest_main.a : $(ODIR)/gtest-all.o $(ODIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.
.PHONY : PumiBuild

PumiBuild :
	$(MAKE) -C $(PUMI_DIR)/build


$(ODIR)/MeshAdjReorder.o : $(SDIR)/MeshAdjReorder.cc \
					  $(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/MeshBuilder_unittest.o : $(SDIR)/MeshBuilder_unittest.cc \
					  $(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/MeshBuilder.o : $(SDIR)/MeshBuilder.cc $(IDIR)/MeshBuilder.h \
					  $(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/SimpleRectMesh_unittest.o : $(SDIR)/SimpleRectMesh_unittest.cc \
					  $(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

#use the phony rule to ensure that the latest pumi library is built
$(ODIR)/a4.o : $(SDIR)/a4.cc $(GTEST_HEADERS) .PHONY
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@


$(BDIR)/a4 : $(ODIR)/a4.o $(ODIR)/MeshBuilder.o $(ODIR)/MeshBuilder_unittest.o \
						   $(ODIR)/SimpleRectMesh_unittest.o $(LDIR)/gtest.a \
						   $(ODIR)/MeshAdjReorder.o
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(PUMI_LIBS) $^ -lmds -lapf -lgmi -lpcu -lpthread  -o $@



Goofy : $(BDIR)/test_xml_reader.out

$(ODIR)/test_xml_reader.o : $(SDIR)/test_xml_reader.cc
	mpicxx $(CPPFLAGS) $(CXXFLAGS) -I $(IDIR) $(PUMI_INC) $(XERCES_INC) -c $< -o $@

$(ODIR)/vtkSAX2Handler.o : $(SDIR)/vtkSAX2Handler.cc
	mpicxx $(CPPFLAGS) $(CXXFLAGS) -I $(IDIR) $(PUMI_INC) $(XERCES_INC) -c $< -o $@

$(BDIR)/test_xml_reader.out : $(ODIR)/test_xml_reader.o $(ODIR)/vtkSAX2Handler.o
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(XERCES_LIBS) $^\
			-lpthread -lxerces-c -o $@

#parser debugging compile
Factfinder : $(BDIR)/parse.out

$(ODIR)/parse.o : $(SDIR)/parseVtk.cc
	mpicxx -v $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) $(XERCES_INC) -c $< -o $@

$(BDIR)/parse.out : $(ODIR)/parse.o
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(PUMI_LIBS) $(XERCES_LIBS) $^ \
			-lmds -lapf -lgmi -lpcu -lpthread -lxerces-c -o $@

