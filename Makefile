INCS=`pkg-config --cflags polyclipping`
LIBS=`pkg-config --libs polyclipping`
CXXFLAGS=-std=c++11 -O2 -W -Wall

svg2nc : svg2nc.cpp svg.cpp path_util.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ $^ $(LIBS)
