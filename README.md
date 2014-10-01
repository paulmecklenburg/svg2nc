svg2nc
======

## About svg2nc

Command line tool to covert scalable vector graphics (svg) files to g-code nc files for a 3-axis nc mill. svg2nc assumes a cylindrical bit and all units are inches.


## Using svg2nc

svg2nc reliease on libsvgtiny to handle svg input. Paths and basic shapes are supported. Clones are ignored. Layers that are marked hidden are not actually ignored.

libsvgtiny does not correctly handle a relative move after a path close. Inkscape generates files that contain this sequence. I have a fix and am working on getting permission from my employer to share it. Until this is done, expect bizarre results.


## Building svg2nc

Dependencies (on debian based systems):
  * clipper (libpolyclipping-dev).
