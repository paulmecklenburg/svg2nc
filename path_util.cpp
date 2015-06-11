#include "path_util.h"

#include <algorithm>

#include "polyclipping/clipper.hpp"

using ClipperLib::Clipper;
using ClipperLib::ClipperOffset;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::PolyNode;
using ClipperLib::PolyTree;
using ClipperLib::ctDifference;
using ClipperLib::ctUnion;
using ClipperLib::etClosedPolygon;
using ClipperLib::jtRound;
using ClipperLib::ptSubject;
using ClipperLib::ptClip;

// TODO: Convert most functions to exit(1) instead of returning bool.

bool SubtractFrom(const Paths &a, Paths *b) {
  if (a.empty() || b->empty())
    return true;
  Clipper c;
  return c.AddPaths(*b, ptSubject, true) &&
    c.AddPaths(a, ptClip, true) &&
    c.Execute(ctDifference, *b);
}

bool UnionInto(const Paths &paths, Paths *result) {
  if (paths.empty())
    return true;

  if (result->empty()) {
    *result = paths;
    return true;
  }

  Clipper c;
  if (!c.AddPaths(*result, ptSubject, true) ||
      !c.AddPaths(paths, ptClip, true)) {
    return false;
  }
  Paths solution;
  if (!c.Execute(ctUnion, solution)) {
    return false;
  }
  result->swap(solution);
  return true;
}

namespace {
  void CopyAndForceOrientation(const PolyNode &node, bool orientation, Paths *out) {
    for (const auto &child : node.Childs) {
      out->push_back(child->Contour);
      if (Orientation(out->back()) != orientation) {
        Path *b = &out->back();
        std::reverse(b->begin(), b->end());
      }
      CopyAndForceOrientation(*child, !orientation, out);
    }
  }
}  // namespace

bool CopyAndForceOrientation(const Paths &paths, bool orientation, Paths *out) {
  // Use Clipper to convert paths to a PolyTree.
  Clipper c;
  PolyTree solution;
  if (!c.AddPaths(paths, ptSubject, true) ||
      !c.Execute(ctUnion, solution))
    return false;
  out->clear();
  CopyAndForceOrientation(solution, orientation, out);
  return true;
}

bool ComputeOffset(const Paths &paths, double amount, Paths *result) {
  // Previous operations can leave small artifacts (e.g. self-intersecting
  // polygons) which ClipperOffset cannot handle. CleanPolygons fixes at least
  // some of these cases.
  Paths cleaned(paths.size());  // CleanPolygons does not resize 'cleaned'.
  CleanPolygons(paths, cleaned);
  Paths tmp_paths;
  if (!CopyAndForceOrientation(cleaned, true, &tmp_paths))
    return false;

  ClipperOffset co;
  co.ArcTolerance = kQuantaPerInch / 1000;
  co.AddPaths(tmp_paths, jtRound, etClosedPolygon);
  co.Execute(*result, InchesToQuanta(amount));
  return true;
}
