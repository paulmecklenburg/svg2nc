CXXFLAGS=-std=c++11 -W -Wall `PKG_CONFIG_PATH=/opt/netsurf/lib/pkgconfig pkg-config --cflags libsvgtiny polyclipping`

svg2nc : svg2nc.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ `PKG_CONFIG_PATH=/opt/netsurf/lib/pkgconfig pkg-config --libs libsvgtiny libwapcaplet polyclipping`
