#ifndef _SVG_H
#define _SVG_H

#include <map>
#include <set>
#include <vector>

#include "common.h"
#include "polyclipping/clipper.hpp"

const double kSvgUnitsPerInch = 90;

bool SvgToPolygons(const char *file_name,
                   const std::map<uint32_t, double> &color_to_elevation,
                   const std::set<uint32_t> &ignore_colors,
                   const bool as_drawn,
                   std::map<double, ClipperLib::Paths> *layers,
                   std::vector<ClipperLib::IntPoint> *delay_points,
                   ClipperLib::cInt *width,
                   ClipperLib::cInt *height) MUST_USE_RESULT;

#endif  // _SVG_H
