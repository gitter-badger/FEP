BDIR = ./bin
SDIR = ./src
IDIR = ./inc
LDIR = ./lib
ODIR = ./src/obj


GTEST_DIR = ~/Documents/GitHub/googletest

GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

CPP_FLAGS += -isystem $(GTEST_DIR)/include

PUMI_DIR = ~/Documents/GitHub/core
PUMI_INC = -I $(PUMI_DIR)/apf -I $(PUMI_DIR)/gmi -I $(PUMI_DIR)/mds \
			-I $(PUMI_DIR)/pcu

#--verbose adds a lot of information to output and expands on include search path
CXXFLAGS = -Wall -g -I $(IDIR) $(PUMI_INC) --pedantic-errors -lmds -lapf -lgmi -lpcu


CC = mpicxx

LIBS = -L $(PUMI_DIR)/build/apf -L $(PUMI_DIR)/build/gmi -L \
		$(PUMI_DIR)/build/mds -L $(PUMI_DIR)/build/pcu


#use DEPS to require a rebuild whenever a header file changes
_DEPS = TestMeshBuilder.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

#hold all the objects other than server
#this is so that whenever something changes the server is recompiled
_OBJ = TestMeshBuilder.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

all: ./bin/a4

$(BDIR)/a4 : ./src/obj/a4.o $(OBJ)
	$(CC) $< -o $@ $(LIBS) $(CFLAGS) 

$(ODIR)/a4.o : ./src/a4.cc 
	$(CC) -c $(CFLAGS) $< -o $@ 

$(ODIR)/TestMeshBuilder.o : ./src/TestMeshBuilder.cc 
	$(CC) -c $(CFLAGS) $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/server $(ODIR)/*.o
