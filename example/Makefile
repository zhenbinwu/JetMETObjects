DEL       = $(shell echo $(DELPHES))
CXX       = g++
CXXFLAGS  = -Wno-deprecated $(shell root-config --cflags) -std=c++0x
CXXFLAGS += -I$(DEL)/external/ -I$(DEL) -DSTANDALONE
LDFLAGS   =
LD        = g++

LIBS      = $(shell root-config --glibs)
LIBS     += -L$(DEL) -lDelphes -L../lib/ -lJetMETObjects


PROGRAM   = DelHATS_JetMET
OBJS      = $(PROGRAM).o

all:            $(PROGRAM)


%.o : %.C $(HEADERS)


%.o : %.cc $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<


$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@$(LD) $(OBJS) $(LIBS) -o $(PROGRAM)
	@echo "done"

clean:;         @rm -f $(OBJS) $(PROGRAM) core 

