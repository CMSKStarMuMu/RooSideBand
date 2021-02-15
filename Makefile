#object_files (testEffi.o)
OBJECTS := $(wildcard *.o)

#root_stuff (root libraries and needed root options)
ROOTLIBS := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lFoam -lRooFit -lRooFitCore
ROOTCINT := $(shell which rootcint)

#exe_files
EXECUTABLEFIT      := testSidebandFit
EXECUTABLEREAD     := readSideband
LIBBERN            := libRooBernsteinSideband.so
CLASS		   := RooBernsteinSideband
CLASSDICT	   := $(CLASS)Dictionary.cxx
CLASSCB 	   := RooDoubleCBFast
CLASSCBDICT	   := $(CLASSCB)Dictionary.cxx

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(CLASS).cxx  $(CLASSDICT) $(CLASSCB).cxx  $(CLASSCBDICT) 

	
all: $(CLASSDICT)   $(CLASSCBDICT) $(EXECUTABLEFIT)  $(EXECUTABLEREAD)

dict: $(CLASSDICT)   $(CLASSCBDICT) 


$(CLASSDICT): $(CLASS).h $(CLASS)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^

$(CLASSCBDICT): $(CLASSCB).h $(CLASSCB)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^


$(EXECUTABLEFIT): $(EXECUTABLEFIT).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.


$(EXECUTABLEREAD): $(EXECUTABLEREAD).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.


#$(EXMAKEHIST): $(EXMAKEHIST).cc 
#	$(CXX) $(CXXFLAGS)  -o $@  $^ $(ROOTLIBS) $(ROOTFLAGS) -I.

$(LIBBERN): $(CLASS).cxx
	$(CXX) $(CXXFLAGS) -fPIC -shared -o $@ $^ $(ROOTLIBS) $(ROOTFLAGS) -I.


#cleaning options
.PHONY: clean cleanall
clean:
	rm -f $(OBJECTS) && rm -f $(EXECUTABLEFIT) $(EXECUTABLEREAD)  $(CLASSDICT)

