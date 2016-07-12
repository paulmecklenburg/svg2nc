#include "svg.h"

#include <math.h>
#include <memory>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"
#include "path_util.h"

#include <string.h>
#include <math.h>
#define NANOSVG_IMPLEMENTATION  // Expands implementation
#include "nanosvg.h"

using ClipperLib::IntPoint;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::cInt;

namespace {
  double SvgToInches(double x) {
    return x / kSvgUnitsPerInch;
  }

  cInt SvgToQuanta(double x) {
    return InchesToQuanta(SvgToInches(x));
  }

  void BezierPoint(double x0, double y0,
                   double x1, double y1,
                   double x2, double y2,
                   double x3, double y3,
                   double t,
                   double *x, double *y) {
    const double omt = 1. - t;
    const double a = omt * omt * omt;
    const double b = 3. * omt * omt * t;
    const double c = 3. * omt * t * t;
    const double d = t * t *t;
    *x = a * x0 + b * x1 + c * x2 + d * x3;
    *y = a * y0 + b * y1 + c * y2 + d * y3;
  }

  void AddBezierPointsToPath(double x0, double y0,
                             double x1, double y1,
                             double x2, double y2,
                             double x3, double y3,
                             double t0, double t2,
                             Path *output) {
    double a0, b0, a1, b1, a2, b2;
    const double t1 = (t0 + t2) * .5;
    BezierPoint(x0, y0, x1, y1, x2, y2, x3, y3, t0, &a0, &b0);
    BezierPoint(x0, y0, x1, y1, x2, y2, x3, y3, t1, &a1, &b1);
    BezierPoint(x0, y0, x1, y1, x2, y2, x3, y3, t2, &a2, &b2);
    const double a_err = (a0 + a2) * .5 - a1;
    const double b_err = (b0 + b2) * .5 - b1;
    const double err = sqrt(a_err * a_err + b_err * b_err);
    if (err > .001) {
      AddBezierPointsToPath(x0, y0, x1, y1, x2, y2, x3, y3, t0, t1, output);
      AddBezierPointsToPath(x0, y0, x1, y1, x2, y2, x3, y3, t1, t2, output);
    } else {
      output->push_back(IntPoint(InchesToQuanta(a2), InchesToQuanta(b2)));
    }
  }

  void ConvertPath(const NSVGpath *input, double svg_height, Paths *output) {
    for (; input; input = input->next) {
      if (!input->closed) {
        fprintf(stderr, "warning: skipping open input path.\n");
        continue;
      }
      output->push_back(Path());
      Path *path = &output->back();
      for (int i = 0; i < input->npts-1; i += 3) {
        const float* p = &input->pts[i*2];
        AddBezierPointsToPath(SvgToInches(p[0]),
                              SvgToInches(svg_height - p[1]),
                              SvgToInches(p[2]),
                              SvgToInches(svg_height - p[3]),
                              SvgToInches(p[4]),
                              SvgToInches(svg_height - p[5]),
                              SvgToInches(p[6]),
                              SvgToInches(svg_height - p[7]),
                              0., 1., path);
      }
    }
  }

  bool BadPath(const Path &p) {
    return p.size() < 3;
  }

  bool BadPaths(const Paths &a) {
    for (const Path &p : a) {
      if (BadPath(p))
        return true;
    }
    return false;
  }
}  // namespace

bool SvgToPolygons(const char *file_name,
                   const std::map<uint32_t, double> &color_to_elevation,
                   const std::set<uint32_t> &ignore_colors,
                   const bool as_drawn,
                   std::map<double, Paths> *layers,
                   std::vector<IntPoint> *delay_points,
                   cInt *width, cInt *height) {
  // TODO: Switch to use inches instead of Inkscape 90 pixels per inch.
  std::unique_ptr<struct NSVGimage, decltype(&nsvgDelete)> diagram(
      nsvgParseFromFile(file_name, "px", kSvgUnitsPerInch), nsvgDelete);
  if (!diagram)
    return false;
  if (width)
    *width = SvgToQuanta(diagram->width);
  if (height)
    *height = SvgToQuanta(diagram->height);

  for (auto shape = diagram->shapes; shape; shape = shape->next) {
    if (shape->paths) {
      Paths paths;
      ConvertPath(shape->paths, diagram->height, &paths);
      if (!paths.empty() && !paths[0].empty() && shape->strokeDashCount > 0) {
        delay_points->push_back(paths[0][0]);
        continue;
      }
      if (BadPaths(paths)) {
        fprintf(stderr,
                "Input contains a path with fewer than 3 points. This may be a "
                "bug in nanosvg.\n");
        return false;
      }
      if (shape->fill.type != NSVG_PAINT_COLOR) {
        fprintf(stderr, "warning: ignoring shaped with non-solid fill color.\n");
        continue;
      }
      const unsigned int fill =
        (shape->fill.color & 0xFF) << 16 |
        (shape->fill.color & 0xFF00) |
        (shape->fill.color & 0xFF0000) >> 16;
      if (ignore_colors.count(fill)) {
        continue;
      }
      const auto elevation_itr = color_to_elevation.find(fill);
      if (elevation_itr == color_to_elevation.end()) {
        fprintf(stderr, "no elevation mapped for: %06X\n", fill);
        return false;
      }
      const double elevation = elevation_itr->second;
      if (as_drawn) {
        for (auto above_iter = layers->upper_bound(elevation);
             above_iter != layers->end(); ++above_iter) {
          if (!SubtractFrom(paths, &above_iter->second)) {
            fprintf(stderr, "unable to subtract paths from above layer.\n");
            return false;
          }
        }
      }
      if (elevation > 0) {
        Paths &layer_paths = (*layers)[elevation];
        if (!UnionInto(paths, &layer_paths)) {
          fprintf(stderr, "unable to merge paths into a common layer.\n");
          return false;
        }
      }
    }
  }
  return true;
}
