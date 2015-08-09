BDIR = ./bin
SDIR = ./src
IDIR = ./inc
LDIR = ./lib
ODIR = ./src/obj

#===============================
#			PUMI 
#===============================

PUMI_DIR = ~/Documents/GitHub/core
PUMI_INC = -I $(PUMI_DIR)/apf -I $(PUMI_DIR)/gmi -I $(PUMI_DIR)/mds \
			-I $(PUMI_DIR)/pcu

PUMI_LIBS += -L $(PUMI_DIR)/build/apf -L $(PUMI_DIR)/build/gmi -L \
		$(PUMI_DIR)/build/mds -L $(PUMI_DIR)/build/pcu -L $(LDIR)

#--verbose adds a lot of information to output and expands on include search path
#-L$(Lgmi)/libgmi  -L$(Lapf)/libapf  -L$(Lpcu)/libpcu  -L$(Lmds)/libmds


CFLAGS = -Wall -g -I $(IDIR) --pedantic-errors
CXX = mpicxx
#LIBS =  $(Lapf)/libapf.a $(Lgmi)/libgmi.a $(Lmds)/libmds.a $(Lpcu)/libpcu.a

#use DEPS to require a rebuild whenever a header file changes
_DEPS = 
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

#hold all the objects other than server
#this is so that whenever something changes the server is recompiled
_OBJ = 
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

#`pkg-config --cflags --libs libmds`

all: ./bin/a1tet ./bin/a1quad ./bin/a1mix

$(BDIR)/a1tet : ./src/obj/a1tet.o
	$(CXX) $(LIBS) $(CFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a1tet.o : ./src/a1tet.cc 
	$(CXX) $(CFLAGS) $(PUMI_INC) -c $< -o $@ 

$(BDIR)/a1quad : ./src/obj/a1quad.o 
	$(CXX) $(LIBS) $(CFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a1quad.o : ./src/a1quad.cc 
	$(CXX) $(CFLAGS) $(PUMI_INC) -c $< -o $@ 

$(BDIR)/a1mix : ./src/obj/a1mix.o
	$(CXX) $(LIBS) $(CFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a1mix.o : ./src/a1mix.cc 
	$(CXX) $(CFLAGS) $(PUMI_INC) -c $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/server $(ODIR)/*.o
