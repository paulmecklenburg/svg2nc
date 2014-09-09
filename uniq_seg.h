#ifndef _UNIQ_SEG_H
#define _UNIQ_SEG_H

#include <map>
#include <tuple>

#include <boost/icl/interval_set.hpp>
#include "polyclipping/clipper.hpp"

// Tool to remove duplicate spans of segments.
class UniqSeg {
public:
  // Given a segment [p0, p1], emit any portions of the span not already seen
  // to output. Also record this segment for future calls RemoveRedundant.
  void RemoveRedundant(const ClipperLib::IntPoint &p0, const ClipperLib::IntPoint &p1,
                       ClipperLib::Paths *output);

  struct Line
    : public std::tuple<ClipperLib::cInt, ClipperLib::cInt, ClipperLib::cInt> {
    // a*x + b*y + c = 0.
    Line(ClipperLib::cInt a, ClipperLib::cInt b, ClipperLib::cInt c)
      : tuple(a, b, c) {}
    ClipperLib::cInt &a() { return std::get<0>(*this); }
    ClipperLib::cInt &b() { return std::get<1>(*this); }
    ClipperLib::cInt &c() { return std::get<2>(*this); }
    const ClipperLib::cInt &a() const { return std::get<0>(*this); }
    const ClipperLib::cInt &b() const { return std::get<1>(*this); }
    const ClipperLib::cInt &c() const { return std::get<2>(*this); }
  };

private:
  std::map<Line, boost::icl::interval_set<ClipperLib::cInt>> seen_;
};

#endif  // _UNIQ_SEG_H
