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
  --svg-file=test-plan.svg     \
  --nc-file=test.nc            \
  --through-elevation=-.003    \
  --max-pass-depth=.3          \
  test.svg
</pre>

This command reads in test.svg and treats red (FF0000) as .25" and green (00FF00) as .5" elevations. The material thickness is explicitly set to .5". No cut will remove more than .3" inches at a time. Deeper cuts are done in passes. While moving, the tool will be lifted .2" above the surface. The results are written to two files: test-plan.svg and test.nc. The former is a svg file depicting the plan which is useful for sanity checking the results.

If no material thickness is specified, svg assumes the material is as thick as the highest elevation given.


## Known Issues and Limitations

svg2nc relies on nanosvg to handle svg input.

Parts placed within holes of other parts are cut first.

svg2nc does not support color gradients. The elevation of each path must be uniform.


## Building svg2nc

svg2nc depends on
* Clipper
* Boost (libboost-dev on debian)

The following instructions assume a debian-based system and were tested on Ubuntu 14.04.

### CMake, git, g++, and Boost
Boost is used directly by svg2nc. CMake is required to build clipper.
<pre>
$ sudo apt-get install --assume-yes libboost-dev git g++ cmake
</pre>

### Clipper
Clipper is available as a Debian package, but as of November 2015, it is version v6.1.3a which has bugs that make svg2nc unusable. Instead download the [latest version](http://www.angusj.com/delphi/clipper.php).
<pre>
$ makedir clipper_ver6.2.1; cd clipper_ver6.2.1
$ unzip ../clipper_ver6.2.1.zip
$ cd cpp
$ cmake .
# Open clipper.hpp in an editor and uncomment the line "#define use_lines".
$ make
$ sudo make install
$ sudo ldconfig
</pre>

### svg2nc and nanosvg
<pre>
$ git clone https://github.com/paulmecklenburg/svg2nc.git
$ cd svg2nc && wget https://raw.githubusercontent.com/memononen/nanosvg/master/src/nanosvg.h && make
</pre>
