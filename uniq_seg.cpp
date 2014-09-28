#include "uniq_seg.h"

#include <assert.h>
#include <utility>

using ClipperLib::cInt;
using ClipperLib::IntPoint;
using ClipperLib::Path;
using ClipperLib::Paths;
using boost::icl::interval;
using boost::icl::interval_set;

namespace {
  void Rot90(cInt x, cInt y, cInt *ox, cInt *oy) {
    *ox = -y;
    *oy = x;
  }

  cInt Abs(cInt a) {
    return a > 0 ? a : -a;
  }

  cInt Gcd(cInt a, cInt b) {
    if (b == 0)
      return a;
    return Gcd(b, a % b);
  }  

  UniqSeg::Line NormalizedLineFromPoints(const IntPoint &p0, const IntPoint &p1) {
    assert(p0 != p1);
    // Compute a,b, and c such that a*p0.X + b*p0.Y + c = 0.
    cInt a, b;
    Rot90(p0.X - p1.X, p0.Y - p1.Y, &a, &b);
    // Normalize.
    if (a == 0) {
      b = b / Abs(b);
    } else if (b == 0) {
      a = a / Abs(a);
    } else {
      const cInt coef = Gcd(Abs(a), Abs(b));
      a /= coef;
      b /= coef;
    }
    cInt c = -(a * p0.X + b * p0.Y);
    if (c < 0 ||
        (c == 0 && b < 0) ||
        (c == 0 && b == 0 && a < 0)) {
      a = -a;
      b = -b;
      c = -c;
    }
    return UniqSeg::Line(a, b, c);
  }

  cInt CoordOnLine(const IntPoint &p, const UniqSeg::Line &line) {
    cInt a, b;
    Rot90(line.a(), line.b(), &a, &b);
    return -(a * p.X + b * p.Y);
  }

  IntPoint PointFromLineAndCoord(const UniqSeg::Line &line, cInt c) {
    cInt a, b;
    Rot90(line.a(), line.b(), &a, &b);
    // a * X + b * Y + c == 0.
    // line.a() * X + line.b() * Y + line.c() = 0.
    // [ a     b ] [X] = [ -c ]
    // [ a() b() ] [Y]   [-c()]
    const cInt denom = a * line.b() - b * line.a();
    const cInt x = (-c * line.b() + -line.c() * -b) / denom;
    const cInt y = (-c * -line.a() + -line.c() * a) / denom;
    return IntPoint(x, y);
  }
}  // namespace

void UniqSeg::RemoveRedundant(
    const IntPoint &p0, const IntPoint &p1, Paths *output) {
  output->clear();
  if (p0 == p1)
    return;
  const UniqSeg::Line line = NormalizedLineFromPoints(p0, p1);
  cInt q0 = CoordOnLine(p0, line);
  cInt q1 = CoordOnLine(p1, line);
  bool reverse = false;
  if (q0 > q1) {
    std::swap(q0, q1);
    reverse = true;
  }
  interval<cInt>::type full_span(q0, q1);
  interval_set<cInt> remainder;
  remainder += full_span;
  auto &seen_spans = seen_[line];
  remainder -= seen_spans;
  for (const interval<cInt>::type &seg : remainder) {
    IntPoint r0 = PointFromLineAndCoord(line, seg.lower());
    IntPoint r1 = PointFromLineAndCoord(line, seg.upper());
    output->push_back(Path());
    Path *o = &output->back();
    if (reverse)
      std::swap(r0, r1);
    o->push_back(r0);
    o->push_back(r1);
  }
  if (reverse) {
    std::reverse(output->begin(), output->end());
  }
  seen_spans.insert(full_span);
}

#ifdef TEST_ONLY
#include <iostream>

bool ExpectEq(const IntPoint &a, const IntPoint &b) {
  if (a == b)
    return true;
  std::cout << "! " << a << " != " << b << std::endl;
  return false;
}

int main(int, char *[]) {
  static const IntPoint pts[] = {
    {5, 3}, {10, 10}, {15, 17}, {20, 24}, {25, 31},
  };

  {
    const UniqSeg::Line line = NormalizedLineFromPoints(pts[0], pts[1]);
    const cInt q = CoordOnLine(pts[2], line);
    assert(PointFromLineAndCoord(line, q) == pts[2]);
  }

  {
    UniqSeg uniq;
    Paths output;

    uniq.RemoveRedundant(pts[1], pts[2], &output);
    if (output.size() != 1 || output[0].size() != 2 ||
        !ExpectEq(output[0][0], pts[1]) ||
        !ExpectEq(output[0][1], pts[2]))
      return 1;

    uniq.RemoveRedundant(pts[1], pts[3], &output);
    if (output.size() != 1 || output[0].size() != 2 ||
        !ExpectEq(output[0][0], pts[2]) ||
        !ExpectEq(output[0][1], pts[3]))
      return 1;

    uniq.RemoveRedundant(pts[4], pts[0], &output);
    if (output.size() != 2 ||
        output[0].size() != 2 ||
        !ExpectEq(output[0][0], pts[4]) ||
        !ExpectEq(output[0][1], pts[3]) ||
        output[1].size() != 2 ||
        !ExpectEq(output[1][0], pts[1]) ||
        !ExpectEq(output[1][1], pts[0]))
      return 1;

    uniq.RemoveRedundant(pts[1], pts[3], &output);
    if (output.size() != 0)
      return 1;

    uniq.RemoveRedundant(IntPoint(0, 0), IntPoint(5, 0), &output);
    if (output.empty())
      return 1;

    uniq.RemoveRedundant(IntPoint(4, 0), IntPoint(0, 0), &output);
    if (!output.empty())
      return 1;

    uniq.RemoveRedundant(IntPoint(0, 0), IntPoint(0, 5), &output);
    if (output.empty())
      return 1;

    uniq.RemoveRedundant(IntPoint(0, 4), IntPoint(0, 0), &output);
    if (!output.empty())
      return 1;
  }

  return 0;
}
#endif  // TEST_ONLY
