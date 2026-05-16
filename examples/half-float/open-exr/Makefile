#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

HALF_DIR = ../../../thirdparty/imath/v3.2.2
CXXFLAGS =  -std=c++17
LIBS =  -lhdf5  -lz -ldl -lm
INCLUDES = -I/usr/local/include -I/usr/include -I./ -I$(HALF_DIR)

test: half-float
	./half-float

half.o : $(HALF_DIR)/half.cpp
	$(CXX) $(INCLUDES) -o half.o  $(CXXFLAGS) -c $(HALF_DIR)/half.cpp

half-float.o : half-float.cpp
	$(CXX) $(INCLUDES) -o half-float.o  $(CXXFLAGS) -c half-float.cpp

half-float: half-float.o half.o
	$(CXX) $^ $(LIBS) -o $@

clean:
	@$(RM) *.o *.h5 half-float

dist-clean: clean
.PHONY: test

