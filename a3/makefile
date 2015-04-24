BDIR = ./bin
SDIR = ./src
IDIR = ./inc
LDIR = ./lib
ODIR = ./src/obj

PUMI_DIR = ~/Documents/GitHub/core
PUMI_INC = -I $(PUMI_DIR)/apf -I $(PUMI_DIR)/gmi -I $(PUMI_DIR)/mds \
			-I $(PUMI_DIR)/pcu

#--verbose adds a lot of information to output and expands on include search path
CFLAGS = -Wall -g -I $(IDIR) $(PUMI_INC) --pedantic-errors -lmds -lapf -lgmi -lpcu
CC = mpicxx

LIBS = -L $(PUMI_DIR)/build/apf -L $(PUMI_DIR)/build/gmi -L \
		$(PUMI_DIR)/build/mds -L $(PUMI_DIR)/build/pcu


#use DEPS to require a rebuild whenever a header file changes
_DEPS = 
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

#hold all the objects other than main target
_OBJ = 
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

#`pkg-config --cflags --libs libmds`

all: ./bin/a3

$(BDIR)/a3 : ./src/obj/a3.o $(OBJ)
	$(CC) $< -o $@ $(LIBS) $(CFLAGS) 

$(ODIR)/a3.o : ./src/a3.cc 
	$(CC) -c $(CFLAGS) $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/server $(ODIR)/*.o