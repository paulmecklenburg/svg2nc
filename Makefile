INCS=`PKG_CONFIG_PATH=/opt/netsurf/lib/pkgconfig pkg-config --cflags libsvgtiny polyclipping`
LIBS=`PKG_CONFIG_PATH=/opt/netsurf/lib/pkgconfig pkg-config --libs libsvgtiny libwapcaplet polyclipping`
CXXFLAGS=-std=c++11 -O2 -W -Wall

svg2nc : svg2nc.cpp uniq_seg.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ $^ $(LIBS) -lann

.PHONY: test
test: uniq_seg_test

uniq_seg_test : uniq_seg.cpp uniq_seg.h
	$(CXX) $(CXXFLAGS) -DTEST_ONLY $(INCS) -o $@ $< $(LIBS) && ./$@ && echo PASS
