#ifndef _PATH_UTIL_H
#define _PATH_UTIL_H

#include "common.h"
#include "polyclipping/clipper.hpp"

bool SubtractFrom(const ClipperLib::Paths &,
                  ClipperLib::Paths *) MUST_USE_RESULT;
bool UnionInto(const ClipperLib::Paths &,
               ClipperLib::Paths *) MUST_USE_RESULT;
bool CopyAndForceOrientation(const ClipperLib::Paths &paths,
                             bool orientation,
                             ClipperLib::Paths *out) MUST_USE_RESULT;
bool ComputeOffset(const ClipperLib::Paths &,
                   double,
                   ClipperLib::Paths *) MUST_USE_RESULT;

#endif  // _PATH_UTIL_H
