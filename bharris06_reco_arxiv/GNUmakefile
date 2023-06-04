
# COMPILER/LINKER
CXX=g++

# C-FLAG helper routines include include directories
CFLAGS=`root-config --cflags`
CFLAGS+=`geant4-config --cflags-without-gui`
CFLAGS+= -std=c++14

# OTHER INCLUDES
# CENNS Stuff
INCLUDES=-I$(CENNS10_BASEDIR)

# LIBRARIES TO LINK
LDLIBS=`root-config --libs`
LDLIBS+=`geant4-config --libs`
LDLIBS+= -L$(CENNS10_BASEDIR)/lib -lMyLib

# PROGRAM NAME(S)
EXE_LIST=reco_vertex plot_results plot_hit_totals plot_hit_positions corr_vertex chaos

# get source
EXE_SRCS=$(addsuffix .cxx , $(EXE_LIST))

# get other sources
SOURCE=$(filter-out $(EXE_SRCS), $(wildcard *.cxx))
# name of object files for these sources
OBJS=$(SOURCE:.cxx=.o)

# program binary
EXE_BIN=$(addprefix bin/, $(EXE_LIST))

# default target are the executables
all: $(EXE_BIN)

# rule to make executables
bin/%: %.cxx
	@echo "<<< compile and link $@ >>>"
	$(CXX) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDLIBS)

# rule to make objects
%.o: %.cxx %.h
	@echo "<<< compiling $< >>>"
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< $@

clean:
	@echo "<<< Cleaning >>>"
	rm -f $(EXE_BIN) *.o
