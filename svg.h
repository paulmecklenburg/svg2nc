#ifndef _SVG_H
#define _SVG_H

#include <map>

#include "common.h"
#include "polyclipping/clipper.hpp"
extern "C" {
#include "svgtiny.h"
}  // extern "C"

const double kSvgUnitsPerInch = 90;

inline double SvgToInches(double x) {
  return x / kSvgUnitsPerInch;
}

inline ClipperLib::cInt SvgToQuanta(double x) {
  return InchesToQuanta(SvgToInches(x));
}

struct svgtiny_diagram *LoadSvg(const char *) MUST_USE_RESULT;

bool SvgToPolygons(const struct svgtiny_diagram &,
                   const std::map<uint32_t, double> &color_to_elevation,
                   const bool as_drawn,
                   std::map<double, ClipperLib::Paths> *) MUST_USE_RESULT;

#endif  // _SVG_H
