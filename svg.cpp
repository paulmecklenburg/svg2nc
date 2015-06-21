#include "svg.h"

#include <math.h>
#include <memory>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"
#include "path_util.h"
extern "C" {
#include "svgtiny.h"
}  // extern "C"

using ClipperLib::IntPoint;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::cInt;

namespace {

  bool ReadFileToString(const char *, std::string *) MUST_USE_RESULT;
  bool ReadFileToString(const char *path, std::string *result) {
    struct stat sb;
    if (stat(path, &sb)) {
      perror(path);
      return false;
    }

    FILE *fd = fopen(path, "rb");
    if (!fd) {
      perror(path);
      return false;
    }

    result->reserve(sb.st_size);
    result->resize(sb.st_size);
    const size_t n = fread(&((*result)[0]), 1, sb.st_size, fd);
    if (n != (size_t)sb.st_size) {
      perror(path);
      return false;
    }

    fclose(fd);
    return true;
  }

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

  bool ConvertPath(const float *, unsigned int, double, Paths *) MUST_USE_RESULT;
  bool ConvertPath(const float *input, unsigned int length, double svg_height, Paths *output) {
    // svg_height is used to invert the y axis.
    Path *path = nullptr;
    unsigned int j;
    for (j = 0; j < length; ) {
      switch ((int) input[j]) {
      case svgtiny_PATH_MOVE:
        if (path)
          return false;
        output->push_back(Path());
        path = &output->back();
      case svgtiny_PATH_LINE:
        if (!path || (j + 3) > length)
          return false;
        path->push_back(IntPoint(SvgToQuanta(input[j+1]),
                                 SvgToQuanta(svg_height - input[j+2])));
        j += 3;
        break;
      case svgtiny_PATH_CLOSE:
        if (!path || (j + 1) > length)
          return false;
        if (path->front() == path->back())
          path->pop_back();  // Sometimes the final point is duplicated.
        path = nullptr;
        j += 1;
        break;
      case svgtiny_PATH_BEZIER:
        {
          if (!path || path->empty() || (j + 7) > length)
            return false;
          const double x0 = QuantaToInches(path->back().X);
          const double y0 = QuantaToInches(path->back().Y);
          AddBezierPointsToPath(x0, y0,
                                SvgToInches(input[j+1]),
                                SvgToInches(svg_height - input[j+2]),
                                SvgToInches(input[j+3]),
                                SvgToInches(svg_height - input[j+4]),
                                SvgToInches(input[j+5]),
                                SvgToInches(svg_height - input[j+6]),
                                0., 1., path);
          j += 7;
          break;
        }
      default:
        return false;
      }
    }
    return true;
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


  struct svgtiny_diagram *LoadSvg(const char *path) {
    std::string contents;
    if (!ReadFileToString(path, &contents)) {
      return nullptr;
    }
    struct svgtiny_diagram *diagram = svgtiny_create();
    const svgtiny_code code = svgtiny_parse(diagram,
                                            contents.data(), contents.size(),
                                            path, 100000, 100000);
    switch (code) {
    case svgtiny_OK:
      return diagram;
    case svgtiny_OUT_OF_MEMORY:
      fprintf(stderr, "svgtiny_OUT_OF_MEMORY\n");
      break;
    case svgtiny_LIBDOM_ERROR:
      fprintf(stderr, "svgtiny_LIBDOM_ERROR\n");
      break;
    case svgtiny_NOT_SVG:
      fprintf(stderr, "svgtiny_NOT_SVG\n");
      break;
    case svgtiny_SVG_ERROR:
      fprintf(stderr, "svgtiny_SVG_ERROR: line %i: %s\n",
              diagram->error_line,
              diagram->error_message);
      break;
    }
    svgtiny_free(diagram);
    return nullptr;
  }

}  // namespace

bool SvgToPolygons(const char *file_name,
                   const std::map<uint32_t, double> &color_to_elevation,
                   const bool as_drawn,
                   std::map<double, Paths> *layers,
                   cInt *width, cInt *height) {
  std::unique_ptr<struct svgtiny_diagram, decltype(&svgtiny_free)> diagram(
      LoadSvg(file_name), svgtiny_free);
  if (!diagram)
    return false;
  if (width)
    *width = SvgToQuanta(diagram->width);
  if (height)
    *height = SvgToQuanta(diagram->height);

  for (unsigned i = 0; i < diagram->shape_count; i++) {
    const auto &shape = diagram->shape[i];
    if (shape.path) {
      Paths paths;
      if (!ConvertPath(shape.path, shape.path_length, diagram->height, &paths)) {
        fprintf(stderr, "unable to process path.\n");
        return false;
      }
      if (BadPaths(paths)) {
        fprintf(stderr,
                "Input contains a bad path. This may be a bug in "
                "libtinysvg. For example, arc commands are converted to a "
                "single line segment.\n");
        return false;
      }
      const auto elevation_itr = color_to_elevation.find(shape.fill);
      if (elevation_itr == color_to_elevation.end()) {
        fprintf(stderr, "no elevation mapped for: %06X\n", shape.fill);
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
