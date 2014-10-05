svg2nc
======

## About svg2nc

Command line tool to covert scalable vector graphics (svg) files to g-code nc files for a 3-axis nc mill. svg2nc assumes a cylindrical bit. All units are inches.


## Using svg2nc

svg2nc generates imperial unit nc files for a 3-axis mill from an svg file. It's easiest to understand with an example:

<pre>
./svg2nc --verbose             \
  --color-elevation=FF0000:.25 \
  --color-elevation=00FF00:.5  \
  --feed-rate=12               \
  --diameter=.125              \
  --clearance-space=.2         \
  --material-thickness=.5      \
  --ps-file=test-debug.ps      \
  --nc-file=test.nc            \
  --through-elevation=-.003    \
  --max-pass-depth=.3          \
  test.svg
</pre>

This command reads in test.svg and treats red (FF0000) as .25" and green (00FF00) as .5" elevations. The material thickness is explicitly set to .5". No cut will remove more than .3" inches at a time. Deeper cuts are done in passes. While moving, the tool will be lifted .2" above the surface. The results are written to two files: test-debug.ps and test.nc. The former is a PostScript file which is useful for sanity checking the results.

If no material thickness is specified, svg assumes the material is as thick as the highest elevation given.


## Known Issues and Limitations

svg2nc uses libsvgtiny to parse svg files. libsvgtiny does not correctly handle a relative move after a path close. Inkscape generates files that contain these. The fix is trivial, but is not yet in mainline libsvgtiny.

svg2nc relies on libsvgtiny to handle svg input. Paths and basic shapes are supported. Clones are ignored. Layers that are marked hidden are not actually ignored.

svg2nc breaks the design into parts and cuts each separately. It will not necessarily order the cuts of nested parts correctly. Do not place parts inside of pockets of other parts. Cutting a set of concentric rings will not work.

svg2nc does not support color gradients. The elevation of each path must be uniform.


## Building svg2nc

svg2nc depends on
* Clipper (libpolyclipping-dev on debian).
* Boost (libboost-dev on debian)
* ANN (libann-dev on debian)
* libsvgtiny

The following instructions assume a debian-based system and were tested on Ubuntu.

### git, g++, ANN, Clipper, Boost, Expat, and gperf
ANN, Clipper and Boost are used directly by svg2nc. Expat and gperf are requried by libsvgtiny.
<pre>
$ sudo apt-get install --assume-yes libann-dev libpolyclipping-dev libboost-dev libexpat1-dev gperf git g++
</pre>

### libsvgtiny
<pre>
$ mkdir tmp; cd tmp

$ git clone git://git.netsurf-browser.org/buildsystem
$ cd buildsystem && sudo make install && cd ..

$ git clone git://git.netsurf-browser.org/libwapcaplet
$ cd libwapcaplet && make && sudo make install && cd ..

$ git clone git://git.netsurf-browser.org/libparserutils
$ cd libparserutils && make && sudo make install && cd ..

$ git clone git://git.netsurf-browser.org/libhubbub.git
$ cd libhubbub && make && sudo make install && cd ..

$ git clone git://git.netsurf-browser.org/libdom
$ cd libdom && make && sudo make install && cd ..

$ git clone git://git.netsurf-browser.org/libsvgtiny.git
$ cd libsvgtiny && make && sudo make install && cd ..
</pre>

### svg2nc
<pre>
$ git clone https://github.com/paulmecklenburg/svg2nc.git
$ cd svg2nc && make
</pre>
