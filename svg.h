#ifndef _SVG_H
#define _SVG_H

#include <map>

#include "common.h"
#include "polyclipping/clipper.hpp"

bool SvgToPolygons(const char *file_name,
                   const std::map<uint32_t, double> &color_to_elevation,
                   const bool as_drawn,
                   std::map<double, ClipperLib::Paths> *layers,
                   ClipperLib::cInt *width,
                   ClipperLib::cInt *height) MUST_USE_RESULT;

#endif  // _SVG_H
