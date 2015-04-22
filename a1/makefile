BDIR = ./bin
SDIR = ./src
IDIR = ./inc
LDIR = ./lib
ODIR = ./src/obj

apf = /home/shivaebola/Documents/GitHub/core/apf/
gmi = /home/shivaebola/Documents/GitHub/core/gmi/
mds = /home/shivaebola/Documents/GitHub/core/mds/
pcu = /home/shivaebola/Documents/GitHub/core/pcu/
export CPLUS_INCLUDE_PATH = $(apf):$(gmi):$(mds):$(pcu)

Lapf = ~/Documents/GitHub/core/build/apf
Lgmi = ~/Documents/GitHub/core/build/gmi
Lmds = ~/Documents/GitHub/core/build/mds
Lpcu = ~/Documents/GitHub/core/build/pcu

export LD_LIBRARY_PATH = "$(Lapf):$(Lgmi):$(Lmds):$(Lpcu)"

#--verbose adds a lot of information to output and expands on include search path
#-L$(Lgmi)/libgmi  -L$(Lapf)/libapf  -L$(Lpcu)/libpcu  -L$(Lmds)/libmds


CFLAGS = -Wall -g -I $(IDIR) --pedantic-errors -lmds -lapf -lgmi -lpcu  
CC = mpicxx
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

$(BDIR)/a1tet : ./src/obj/a1tet.o $(OBJ)
	$(CC) $< -o $@ $(LIBS) $(CFLAGS) 

$(ODIR)/a1tet.o : ./src/a1tet.cc 
	$(CC) -c $(CFLAGS) $< -o $@ 

$(BDIR)/a1quad : ./src/obj/a1quad.o $(OBJ)
	$(CC) $< -o $@ $(LIBS) $(CFLAGS) 

$(ODIR)/a1quad.o : ./src/a1quad.cc 
	$(CC) -c $(CFLAGS) $< -o $@ 

$(BDIR)/a1mix : ./src/obj/a1mix.o $(OBJ)
	$(CC) $< -o $@ $(LIBS) $(CFLAGS) 

$(ODIR)/a1mix.o : ./src/a1mix.cc 
	$(CC) -c $(CFLAGS) $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/server $(ODIR)/*.o