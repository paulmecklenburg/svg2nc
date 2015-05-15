#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>

#include <assert.h>
#include <ANN/ANN.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "polyclipping/clipper.hpp"
#include "uniq_seg.h"
extern "C" {
#include "svgtiny.h"
}  // extern "C"

using ClipperLib::Area;
using ClipperLib::cInt;
using ClipperLib::ctUnion;
using ClipperLib::Clipper;
using ClipperLib::ClipperOffset;
using ClipperLib::ctDifference;
using ClipperLib::ctIntersection;
using ClipperLib::etClosedLine;
using ClipperLib::etClosedPolygon;
using ClipperLib::etOpenButt;
using ClipperLib::IntPoint;
using ClipperLib::Orientation;
using ClipperLib::jtRound;
using ClipperLib::jtSquare;
using ClipperLib::ptSubject;
using ClipperLib::ptClip;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::PolyNode;
using ClipperLib::PolyTree;

// TODO: Improve the error reporting. In some cases the program should just
// abort instead of returning false up multple layers of function calls.

#define MUST_USE_RESULT __attribute__((warn_unused_result))

namespace {
  struct Config {
    std::map<uint32_t, double> color_to_elevation;
    double clearance_space = .25;
    double diameter = -1.;
    double feed_rate = 12.;
    double material_thickness = -1.;
    double mill_overlap = .1;
    double max_pass_depth = .25;
    double spindle_speed = 1000.;
    double through_elevation = 0;
    bool verbose = false;
    int part = -1;
    std::string svg_path;
    std::string output_ps_path;
    std::string output_nc_path;
  };

  void PrintUsage(FILE *stream, const char *program_name, int exit_code) {
    fprintf(stream, "Usage:  %s [OPTION]... [SVG FILE]\n", program_name);
    fprintf(
        stream,
        "  -h  --help                 Display this usage information.\n"
        "  -c --color-elevation=<hex color>:<inches>\n"
        "     Specify the elevation for a color.\n"
        "  -f --feed-rate=<inches per minute>\n"
        "     Specify the feed rate.\n"
        "  -d --diameter=<inches>     Set the tool DIAMETER.\n"
        "  -l --clearance-space=<inches>\n"
        "     Safety space to use above the part while moving.\n"
        "  -m --material-thickness=<inches>\n"
        "     Specify the material thickness.\n"
        "  -n --nc-file=<path>        Write the cut plan to a .nc file.\n"
        "  -o --mill-overlap=<fraction>\n"
        "     Fraction of a cut to overlap with adjacent cuts.\n"
        "  -p --ps-file=<path>        Write the cut plan to a .ps file.\n"
        "  -r --part=<number>         Only process the specified number part.\n"
        "  -s --spindle_speed=<rpm>   Set the spindle speed.\n"
        "  -t --through-elevatation=<inches>\n"
        "     The elevation to use while cutting holes and outlines.\n"
        "  -v --verbose               Enable verbose output.\n"
        "  -x --max-pass-depth=<depth>\n"
        "     The maximum depth to cut in a single pass.\n"
    );
    exit(exit_code);
  }

  void ParseColorElevationPair(const char *optarg,
                            std::map<uint32_t, double> *color_to_elevation) {
    uint32_t color;
    double elevation;
    if (2 == sscanf(optarg, "%x:%lg", &color, &elevation)) {
      if (elevation <= 0) {
        fprintf(stderr, "Elevation values must be greater than zero.\n");
        exit(EXIT_FAILURE);
      }
      (*color_to_elevation)[color] = elevation;
    } else {
      fprintf(stderr,
              "'%s' does not match the expected format "
              "<hex color>:<decimal height>\n",
              optarg);
      exit(EXIT_FAILURE);
    }
  }

  Config ParseArgs(int argc, char *argv[]) {
    const char *program_name = argv[0];
    Config config;
    int c;
    static const struct option long_options[] = {
      {"clearance-space", required_argument, nullptr, 'l'},
      {"color-elevation", required_argument, nullptr, 'c'},
      {"diameter", required_argument, nullptr, 'd'},
      {"feed-rate", required_argument, nullptr, 'f'},
      {"help", no_argument, nullptr, 'h'},
      {"material-thickness", required_argument, nullptr, 'm'},
      {"max-pass-depth", required_argument, nullptr, 'x'},
      {"mill-overlap", required_argument, nullptr, 'o'},
      {"nc-file", required_argument, nullptr, 'n'},
      {"part", required_argument, nullptr, 'r'},
      {"ps-file", required_argument, nullptr, 'p'},
      {"spindle-speed", required_argument, nullptr, 's'},
      {"through-elevation", required_argument, nullptr, 't'},
      {"verbose", no_argument, nullptr, 'v'},
      {nullptr, 0, nullptr, 0},
    };
    while ((c = getopt_long(
                argc, argv, "c:d:f:hl:m:n:o:p:r:s:t:vx:",
                long_options, nullptr)) != -1) {
      switch (c) {
      case '?':
        PrintUsage(stderr, program_name, EXIT_FAILURE);
        break;
      case 'c':
        ParseColorElevationPair(optarg, &config.color_to_elevation);
        break;
      case 'd':
        config.diameter = atof(optarg);
        break;
      case 'f':
        config.feed_rate = atof(optarg);
        break;
      case 'h':
        PrintUsage(stderr, program_name, EXIT_SUCCESS);
        break;
      case 'l':
        config.clearance_space = atof(optarg);
        break;
      case 'm':
        config.material_thickness = atof(optarg);
        break;
      case 'n':
        config.output_nc_path = optarg;
        break;
      case 'o':
        config.mill_overlap = atof(optarg);
        break;
      case 'p':
        config.output_ps_path = optarg;
        break;
      case 'r':
        config.part = atoi(optarg);
        break;
      case 's':
        config.spindle_speed = atof(optarg);
        break;
      case 't':
        config.through_elevation = atof(optarg);
        break;
      case 'v':
        config.verbose = true;
        break;
      case 'x':
        config.max_pass_depth = atof(optarg);
        break;
      }
    }

    if (optind != argc - 1) {
      PrintUsage(stderr, program_name, EXIT_FAILURE);
    }
    config.svg_path = argv[optind];
    if (config.material_thickness <= 0) {
      for (const auto &iter : config.color_to_elevation) {
        config.material_thickness = std::max(config.material_thickness,
                                             iter.second);
      }
      if (config.material_thickness <= 0) {
        fprintf(stderr, "Material thickness must be greater than zero.\n");
        exit(EXIT_FAILURE);
      }
    }
    if (config.diameter <= 0) {
      fprintf(stderr, "Bit diameter must be greater than zero.\n");
      exit(EXIT_FAILURE);
    }
    if (config.through_elevation > 0) {
      fprintf(stderr, "Through elevation cannot be greater than zero.\n");
      exit(EXIT_FAILURE);
    }
    if (config.mill_overlap < 0 || config.mill_overlap >= 1.) {
      fprintf(stderr,
              "Mill overlap must be greater than or equal to zero and less than"
              " one.\n");
      exit(EXIT_FAILURE);
    }
    if (config.max_pass_depth <= 0) {
      fprintf(stderr, "Max pass depth must be greater than zero.\n");
      exit(EXIT_FAILURE);
    }
    return config;
  }

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

  struct svgtiny_diagram *LoadSvg(const char *) MUST_USE_RESULT;
  struct svgtiny_diagram *LoadSvg(const char *path) {
    std::string svg_contents;
    if (!ReadFileToString(path, &svg_contents)) {
      return nullptr;
    }
    struct svgtiny_diagram *diagram = svgtiny_create();
    const svgtiny_code code = svgtiny_parse(diagram,
                                            svg_contents.data(), svg_contents.size(),
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

  const cInt kQuantaPerInch = 100000;
  const double kSvgUnitsPerInch = 90;

  inline double SvgToInches(double x) {
    return x / kSvgUnitsPerInch;
  }

  inline cInt InchesToQuanta(double x) {
    return x * kQuantaPerInch;
  }

  inline double QuantaToInches(cInt x) {
    return double(x) / kQuantaPerInch;
  }

  inline cInt SvgToQuanta(double x) {
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

  bool UnionInto(const Paths &, Paths *) MUST_USE_RESULT;
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
  
  bool SvgToPolygons(const struct svgtiny_diagram &,
                     const Config &,
                     std::map<double, Paths> *) MUST_USE_RESULT;
  bool SvgToPolygons(const struct svgtiny_diagram &diagram,
                     const Config &config,
                     std::map<double, Paths> *layers) {
    for (unsigned i = 0; i < diagram.shape_count; i++) {
      const auto &shape = diagram.shape[i];
      if (shape.path) {
        Paths paths;
        if (!ConvertPath(shape.path, shape.path_length, diagram.height, &paths)) {
          fprintf(stderr, "unable to process path.\n");
          return false;
        }
        const auto elevation_itr = config.color_to_elevation.find(shape.fill);
        if (elevation_itr == config.color_to_elevation.end()) {
          fprintf(stderr, "no elevation mapped for: %06X\n", shape.fill);
          return false;
        }
        Paths &layer_paths = (*layers)[elevation_itr->second];
        if (!UnionInto(paths, &layer_paths)) {
            fprintf(stderr, "unable to merge paths into a common layer.\n");
            return false;
        }
      }
    }
    return true;
  }

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

  bool ComputeOffset(const Paths &, double, Paths *) MUST_USE_RESULT;
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

  void CutToPolygon(const Paths &cut, double radius, Paths *result) {
    // Offset a path to a polygon representing the removed material.
    ClipperOffset co;
    co.ArcTolerance = kQuantaPerInch / 1000;
    co.AddPaths(cut, jtRound, etClosedLine);
    co.Execute(*result, InchesToQuanta(radius));
  }

  bool SubtractFrom(const Paths &, Paths *) MUST_USE_RESULT;
  bool SubtractFrom(const Paths &a, Paths *b) {
    if (a.empty() || b->empty())
      return true;
    Clipper c;
    return c.AddPaths(*b, ptSubject, true) &&
      c.AddPaths(a, ptClip, true) &&
      c.Execute(ctDifference, *b);
  }

  bool EnforceStandardMilling(Paths *) MUST_USE_RESULT;
  bool EnforceStandardMilling(Paths *edge_cuts) {
    return CopyAndForceOrientation(*edge_cuts, true, edge_cuts);
  }

  bool MillEdges(const Config &,
                 const Paths &,
                 Paths *) MUST_USE_RESULT;
  bool MillEdges(const Config &config,
                 const Paths &higher_layers,
                 Paths *cuts) {
    Paths edge_cuts;
    if (!ComputeOffset(higher_layers, config.diameter * .5, &edge_cuts))
      return false;
    if (!EnforceStandardMilling(&edge_cuts))
      return false;
    cuts->insert(cuts->end(), edge_cuts.begin(), edge_cuts.end());
    return true;
  }

  bool CanSlide(const IntPoint &p0, const IntPoint &p1,
                double radius, const Paths &layers_above) {
    if (layers_above.empty())
      return true;

    Paths swept_path;
    {  // Create a swept path for the tool. Ideal would use 'jtRound' with
      // exactly radius. Instead etOpenButt is used. We assume that the start
      // and end are already safe and make the path a wider (kEpsilon) to give
      // a bit of extra space. Meatspace is not as precise as math.
      ClipperOffset co;
      const Path line{p0, p1};
      co.AddPath(line, jtSquare /* doesn't actually matter */, etOpenButt);
      const double kEpsilon = 0.005;
      co.Execute(swept_path, InchesToQuanta(radius + kEpsilon));
    }

    Clipper c;
    if (!c.AddPaths(layers_above, ptSubject, true) ||
        !c.AddPaths(swept_path, ptClip, true))
      return false;  // TODO: report error.

    Paths intersection;
    return c.Execute(ctIntersection, intersection) && intersection.empty();
  }

  void MergeCut(const Path &cut,
                const Paths &layers_above,
                double diameter,
                Paths *all_cuts) {
    if (all_cuts->empty()) {
      all_cuts->push_back(cut);
      return;
    }

    std::vector<std::array<ANNcoord, 2>> point_data;
    for (const auto &p : *all_cuts) {
      for (const auto &pt : p) {
        point_data.push_back(std::array<ANNcoord, 2>{
            QuantaToInches(pt.X), QuantaToInches(pt.Y)});
      }
    }

    std::vector<ANNpoint> points;
    for (auto &pd : point_data) {
      points.push_back(pd.data());
    }

    ANNkd_tree tree(points.data(), points.size(), 2);

    ANNidx nearest_idx = 0;
    ANNdist nearest_dist = std::numeric_limits<ANNdist>::max();
    size_t cut_pos = 0;
    for (size_t i = 0; i < cut.size(); ++i) {
      const auto &pt = cut[i];
      std::array<ANNcoord, 2> q{QuantaToInches(pt.X), QuantaToInches(pt.Y)};
      ANNidx idx;
      ANNdist dist;
      tree.annkSearch(q.data(), 1, &idx, &dist);
      if (dist < nearest_dist) {
        nearest_dist = dist;
        nearest_idx = idx;
        cut_pos = i;
      }
    }

    Path *dst = nullptr;
    for (auto &p : *all_cuts) {
      if (nearest_idx < ANNidx(p.size())) {
        dst = &p;
        break;
      }
      nearest_idx -= p.size();
    }

    const int kMaxSlideDiameters = 2;
    if (sqrt(nearest_dist) <= diameter * kMaxSlideDiameters &&
        CanSlide(cut[cut_pos], (*dst)[nearest_idx], diameter * .5, layers_above)) {
      Path to_insert;
      to_insert.push_back((*dst)[nearest_idx]);
      to_insert.insert(to_insert.end(), cut.begin() + cut_pos, cut.end());
      to_insert.insert(to_insert.end(), cut.begin(), cut.begin() + cut_pos + 1);
      dst->insert(dst->begin() + nearest_idx, to_insert.begin(), to_insert.end());
    } else {
      all_cuts->push_back(cut);
    }
  }

  bool MillSurface(const Config &,
                   const Paths &,
                   const Paths &,
                   Paths *) MUST_USE_RESULT;
  bool MillSurface(const Config &config,
                   const Paths &surface,
                   const Paths &layers_above,
                   Paths *cuts) {
    Paths remaining(surface);
    while (!remaining.empty()) {
      // Shrink by less than radius.
      Paths cut;
      if (!ComputeOffset(remaining,
                         -(config.diameter * .5 * (1. - config.mill_overlap)),
                         &cut)) {
        return false;
      }

      // In some cases an area smaller than the bit is left and cut is empty.
      // To make sure this material is removed we cut the outer edge. This is
      // less efficient than is optimal, but it's a rare, small case.
      if (cut.empty())
        cut.swap(remaining);

      for (const auto &c : cut) {
        MergeCut(c, layers_above, config.diameter, cuts);
      }

      if (!remaining.empty()) {
        Paths mask;
        CutToPolygon(cut, config.diameter * .5, &mask);
        if (!SubtractFrom(mask, &remaining))
          return false;
      }
    }
    return true;
  }

  struct LayerCuts {
    // Direction of edge cuts is important.
    Paths edges;

    // Direction of surface cuts is not important.
    Paths surface;
  };

  bool ComputeCuts(const Config &,
                   const std::map<double, Paths> &,
                   std::map<double, LayerCuts> *) MUST_USE_RESULT;
  bool ComputeCuts(const Config &config,
                   const std::map<double, Paths> &layers,
                   std::map<double, LayerCuts> *elevation_to_cuts) {
    Paths higher_layer_union;
    for (auto iter = layers.rbegin(); iter != layers.rend(); ++iter) {
      const double elevation = iter->first;
      const Paths &polygons = iter->second;
      if (elevation < config.material_thickness) {
        Paths mill_area = polygons;
        LayerCuts *cuts = &(*elevation_to_cuts)[elevation];
        if (!higher_layer_union.empty()) {
          // Mill the edges of the layer above.
          if (!MillEdges(config, higher_layer_union, &cuts->edges))
            return false;

          Paths tmp;
          if (!ComputeOffset(higher_layer_union, config.diameter, &tmp))
            return false;
          if (!SubtractFrom(tmp, &mill_area))
            return false;
        }
        if (!MillSurface(config, mill_area, higher_layer_union, &cuts->surface))
          return false;
      }
      if (!UnionInto(polygons, &higher_layer_union))
        return false;
    }
    return MillEdges(config,
                     higher_layer_union,
                     &(*elevation_to_cuts)[config.through_elevation].edges);
  }

  void ClosedToOpen(Paths *paths) {
    for (auto &p : *paths) {
      if (p.back() != p.front()) {
        p.push_back(p.front());
      }
    }
  }
  
  bool OnOrInside(const Path &, const Path &) MUST_USE_RESULT;
  bool OnOrInside(const Path &inner, const Path &outer) {
    for (const auto &pt : inner) {
      // PointInPolygon uses doubles internally and as a result occasionally
      // lies about the status of a point.
      if (PointInPolygon(pt, outer) == 0)
        return false;
    }
    return true;
  }

  bool TrimCutsToBoundingBox(cInt, cInt, Paths *) MUST_USE_RESULT;
  bool TrimCutsToBoundingBox(cInt width, cInt height, Paths *cuts) {
    const Path bbox{{0, 0}, {width, 0}, {width, height}, {0, height}};
    Paths saved;
    saved.reserve(cuts->size());
    Clipper c;
    for (auto &cut : *cuts) {
      // This test is done to work around a clipper bug where segments are
      // dropped from open paths with the final point the same as the starting
      // point in the *output*. This only happens with paths that aren't clipped
      // at all.
      if (OnOrInside(cut, bbox)) {
        saved.push_back(Path());
        saved.back().swap(cut);
      } else {
        if (!c.AddPath(cut, ptSubject, false))
          return false;
      }
    }
    cuts->clear();
    PolyTree solution;
    if (!c.AddPath(bbox, ptClip, true) ||
        !c.Execute(ctIntersection, solution))
      return false;
    OpenPathsFromPolyTree(solution, *cuts);
    cuts->insert(cuts->end(), saved.begin(), saved.end());
    return true;    
  }

  struct CutPath {
    double elevation;
    Path path;
  };

  struct Part {
    std::vector<CutPath> interior;
    std::vector<CutPath> perimeter;
  };

  Part *ResizeGet(std::vector<Part> *t, size_t i) {
    if (t->size() <= i)
      t->resize(i + 1);
    return &(*t)[i];
  }

  bool AddCutPathToPart(
      const CutPath &, const Paths &, std::vector<Part> *) MUST_USE_RESULT;
  bool AddCutPathToPart(
      const CutPath &cp, const Paths &perimeters, std::vector<Part> *parts) {
    for (size_t i = 0; i < perimeters.size(); ++i) {
      const Path &perimeter = perimeters[i];
      if (cp.path.front() == cp.path.back()) {
        const Path closed_cp(cp.path.begin(), cp.path.end() - 1);
        const double cp_area = fabs(Area(closed_cp));
        const double perimeter_area = fabs(Area(perimeter));
        Paths unn{perimeter};
        if (!UnionInto(Paths{closed_cp}, &unn))
          return false;
        if (unn.size() > 1)
          continue;  // Disjoint
        const double union_area = fabs(Area(unn[0]));
        if (union_area < perimeter_area + cp_area * .5) {
          bool is_interior = perimeter_area - cp_area > perimeter_area * .00001;
          if (cp.elevation > 0. || is_interior) {
            ResizeGet(parts, i)->interior.push_back(cp);
          } else {
            ResizeGet(parts, i)->perimeter.push_back(cp);
          }
          return true;
        }
      } else {
        // PointInPolygon isn't 100% accurate, so a heuristic is used. Assume
        // parts are reasonably spaced. The first intersection is accepted.
        for (const auto &pt : cp.path) {
          if (PointInPolygon(pt, perimeter) != 0) {
            ResizeGet(parts, i)->interior.push_back(cp);
            return true;
          }
        }
      }
    }
    ResizeGet(parts, perimeters.size())->interior.push_back(cp);
    return true;
  }

  bool FindDisjointParts(const std::vector<CutPath> &, std::vector<Part> *) MUST_USE_RESULT;
  bool FindDisjointParts(const std::vector<CutPath> &all, std::vector<Part> *parts) {
    parts->clear();
    Paths unn;
    for (const auto &cp : all) {
      if (cp.elevation <= 0. && cp.path.front() == cp.path.back()) {
        Paths paths{cp.path};
        paths.front().pop_back();
        if (!UnionInto(paths, &unn))
          return false;        
      }
    }
    for (const auto &cp : all) {
      if (!AddCutPathToPart(cp, unn, parts))
        return false;
    }
    return true;
  }

  std::vector<double> PassElevations(const Config &config, double elevation) {
    const double delta = config.material_thickness - elevation;
    const int passes = ceil(delta / config.max_pass_depth);
    std::vector<double> result;
    for (int i = passes - 1; i >= 0; --i) {
      result.push_back(i * (delta / passes) + elevation);
    }
    return result;
  }

  cInt DistanceSquared(const IntPoint &a, const IntPoint &b) {
    const cInt dx = a.X - b.X;
    const cInt dy = a.Y - b.Y;
    return dx*dx + dy*dy;
  }

  IntPoint LastPosition(std::vector<CutPath> &history) {
    if (history.empty()) {
      return IntPoint{0, 0};
    } else {
      return history.back().path.back();
    }
  }

  cInt ClosestPointInPath(
      const IntPoint &pt, const Path &path, size_t *index = nullptr) {
    const bool is_loop = path.front() == path.back();
    const size_t end = is_loop ? path.size() : 1;
    cInt min_d2 = std::numeric_limits<cInt>::max();
    for (size_t i = 0; i < end; ++i) {
      const cInt d2 = DistanceSquared(pt, path[i]);
      if (d2 < min_d2) {
        if (index)
          *index = i;
        min_d2 = d2;
      }
    }
    return min_d2;
  }

  cInt DistanceSquaredToCutPaths(const IntPoint &last_pos,
                                 const std::vector<CutPath> &paths) {
    cInt min_d2 = std::numeric_limits<cInt>::max();
    for (const auto &cp : paths) {
      min_d2 = std::min(min_d2, ClosestPointInPath(last_pos, cp.path));
    }
    return min_d2;
  }

  Part RemoveClosestPart(const IntPoint &last_pos, std::vector<Part> *parts) {
    size_t closest = 0;
    cInt min_d2 = std::numeric_limits<cInt>::max();
    for (size_t i = 0; i < parts->size(); ++i) {
      const Part &part = (*parts)[i];
      const cInt d2 = std::min(
          DistanceSquaredToCutPaths(last_pos, part.interior),
          DistanceSquaredToCutPaths(last_pos, part.perimeter));
      if (d2 < min_d2) {
        min_d2 = d2;
        closest = i;
      }
    }
    Part res = (*parts)[closest];
    if (closest + 1 != parts->size()) {
      std::swap((*parts)[closest], parts->back());
    }
    parts->pop_back();
    return res;
  }

  void MoveNearestCutPathAndLayer(const Config &config,
                                  std::vector<CutPath> *input,
                                  std::vector<CutPath> *output) {
    const IntPoint last_pos = LastPosition(*output);
    // Find the point in the remaining cuts closest to the current position.
    // For non-loop cuts, only consider starting at the beginning of the cut.
    CutPath *cp = nullptr;
    size_t ind;
    cInt min_d2 = std::numeric_limits<cInt>::max();
    for (auto &cpt : *input) {
      size_t index;
      const cInt d2 = ClosestPointInPath(last_pos, cpt.path, &index);
      if (d2 < min_d2) {
        ind = index;
        cp = &cpt;
        min_d2 = d2;
      }
    }

    CutPath new_cp;
    new_cp.elevation = cp->elevation;
    if (ind == 0) {
      new_cp.path = cp->path;
    } else {
      new_cp.path.insert(new_cp.path.end(),
                         cp->path.begin() + ind, cp->path.end());
      new_cp.path.insert(new_cp.path.end(),
                         cp->path.begin() + 1, cp->path.begin() + ind + 1);
    }

    // Remove cp from input.
    if (cp != &input->back()) {
      std::swap(*cp, input->back());
    }
    input->pop_back();

    for (const double elevation : PassElevations(config, new_cp.elevation)) {
      new_cp.elevation = elevation;
      output->push_back(new_cp);
    }
  }

  double PathLength(const Path &cut) {
    double total = 0;
    for (size_t i = 0; i + 1 < cut.size(); ++i) {
      total += sqrt(DistanceSquared(cut[i], cut[i+1]));
    }
    return QuantaToInches(total);
  }

  void AddOpenPathsToPs(const Path &cut, FILE *fp) {
    for (const auto &pt : cut) {
      fprintf(fp, "%f %f %s\n",
              QuantaToInches(pt.X) * 72,
              QuantaToInches(pt.Y) * 72,
              &pt == &cut.front() ? "newpath moveto" : "lineto");
    }
    fprintf(fp, "stroke\n");
  }

  void WriteCutsToPs(const std::string &path,
                     const std::vector<CutPath> &ordered_cuts) {
    FILE *fp = fopen(path.c_str(), "w");
    if (!fp) {
      perror(path.c_str());
      exit(EXIT_FAILURE);
    }
    fprintf(fp, "0.2 setlinewidth\n");
    IntPoint last(0, 0);
    bool green = true;
    for (const auto &cp : ordered_cuts) {
      // Move
      Path move{last, cp.path.front()};
      fprintf(fp, "1 0 0 setrgbcolor\n");
      AddOpenPathsToPs(move, fp);
      // Cut
      fprintf(fp, "%s setrgbcolor\n", green ? "0 1 0" : "0 0 1");
      green = !green;
      AddOpenPathsToPs(cp.path, fp);
      last = cp.path.back();
    }

    fclose(fp);
  }

  void WriteCutsToNc(const Config &config,
                     const std::vector<CutPath> &ordered_cuts) {
    FILE *fp = fopen(config.output_nc_path.c_str(), "w");
    if (!fp) {
      perror(config.output_nc_path.c_str());
      exit(EXIT_FAILURE);
    }

    fprintf(fp,
            ";; %s\n"
            "G40 (disable tool radius compensation)\n"
            "G49 (disable tool length compensation)\n"
            "G80 (cancel modal motion)\n"
            "G54 (select coordinate system 1)\n"
            "G90 (disable incremental moves)\n"
            "G20 (imperial)\n"
            "G61 (exact path mode)\n"
            "F%.3f (feed rate)\n"
            "S%.1f (spindle speed)\n"
            "M3 (start spindle)\n"
            "G04 p1 (wait for 1 second)\n",
            config.svg_path.c_str(), config.feed_rate, config.spindle_speed);

    const double safe_elevation =
      config.material_thickness + config.clearance_space;
    IntPoint last{0, 0};
    for (const auto &cp : ordered_cuts) {
      const auto &first = cp.path.front();
      if (first != last) {
        fprintf(fp, "G0 Z%f\n", safe_elevation);
        fprintf(fp, " X%f Y%f\n",
                QuantaToInches(first.X),
                QuantaToInches(first.Y));
      }
      fprintf(fp, "G1 Z%f\n", cp.elevation);
      for (size_t i = 1; i < cp.path.size(); ++i) {
        const auto &cur = cp.path[i];
        fprintf(fp, " X%f Y%f\n", QuantaToInches(cur.X), QuantaToInches(cur.Y));
      }
      last = cp.path.back();
    }

    fprintf(fp,
            "G0 Z%f\n"
            "M5 (stop spindle)\n"
            "M2 (end program)\n",
            safe_elevation);
    fclose(fp);
  }
}  // namespace

int main(int argc, char *argv[]) {
  const Config config = ParseArgs(argc, argv);
  struct svgtiny_diagram *diagram = LoadSvg(config.svg_path.c_str());
  if (!diagram)
    return EXIT_FAILURE;

  std::map<double, Paths> layers;
  const cInt width = SvgToQuanta(diagram->width);
  const cInt height = SvgToQuanta(diagram->height);
  {
    const bool r = SvgToPolygons(*diagram, config, &layers);
    svgtiny_free(diagram);
    if (!r) {
      fprintf(stderr, "svg processing failed.\n");
      return EXIT_FAILURE;
    }
  }

  // Compute layer by layer cuts as closed loops.
  std::map<double, LayerCuts> elevation_to_cuts;
  if (!ComputeCuts(config, layers, &elevation_to_cuts)) {
    fprintf(stderr, "failed to compute cut paths.\n");
    return EXIT_FAILURE;
  }

  for (auto &iter : elevation_to_cuts) {
    // Convert the cuts to open seqment sequences.
    ClosedToOpen(&iter.second.edges);
    ClosedToOpen(&iter.second.surface);

    // Clip to the size of the input drawing.
    if (!TrimCutsToBoundingBox(width, height, &iter.second.edges) ||
        !TrimCutsToBoundingBox(width, height, &iter.second.surface)) {
      fprintf(stderr, "failed trimming cuts to bounding box.\n");
      return EXIT_FAILURE;
    }
  }

  {  // When multiple cuts are on top of each other, only keep the bottom cut.
    UniqSeg uniq_seg;
    for (auto &iter : elevation_to_cuts) {
      Paths *edges = &iter.second.edges;
      Paths new_edges;
      for (const Path &path : *edges) {
        Paths new_path;
        for (size_t i = 0; i + 1 < path.size(); ++i) {
          Paths tmp;
          uniq_seg.RemoveRedundant(path[i], path[i + 1], &tmp);
          for (const Path &p : tmp) {
            assert(p.size() == 2);
            if (!new_path.empty() && new_path.back().back() == p.front()) {
              new_path.back().push_back(p.back());
            } else {
              new_path.push_back(p);
            }
          }
        }
        if (new_path.size() > 1 &&
            new_path.front().front() == new_path.back().back()) {
          new_path.back().insert(new_path.back().end(),
                                 new_path.front().begin(),
                                 new_path.front().end());
          new_path.front().swap(new_path.back());
          new_path.pop_back();
        }
        new_edges.insert(new_edges.end(), new_path.begin(), new_path.end());
      }
      edges->swap(new_edges);
    }
  }

  // Separate cuts into separate pieces that are surrounded by full depth cuts.
  std::vector<CutPath> all_ordered_cuts;
  {
    // Fill unordered_cuts with all cuts from all layers.
    std::vector<CutPath> unordered_cuts;
    for (const auto &layer : elevation_to_cuts) {
      const double elevation = layer.first;
      for (const auto &p : layer.second.edges)
        unordered_cuts.push_back(CutPath{elevation, p});
      for (const auto &p : layer.second.surface)
        unordered_cuts.push_back(CutPath{elevation, p});
    }

    std::vector<Part> parts;
    if (!FindDisjointParts(unordered_cuts, &parts)) {
      fprintf(stderr, "failed computing disjoint parts.\n");
      return EXIT_FAILURE;      
    }

    if (config.verbose) {
      printf("Number of distinct parts: %lu\n", parts.size());
    }

    if (config.part >= 0 && config.part < (signed)parts.size()) {
      parts = {parts[config.part]};
    }

    while (!parts.empty()) {
      Part part = RemoveClosestPart(
          LastPosition(all_ordered_cuts), &parts);
      while (!part.interior.empty()) {
        MoveNearestCutPathAndLayer(
            config, &part.interior, &all_ordered_cuts);
      }
      while (!part.perimeter.empty()) {
        MoveNearestCutPathAndLayer(
            config, &part.perimeter, &all_ordered_cuts);
      }
    }
  }

  if (config.verbose) {
    double cut_distance = 0;
    for (const auto &cp : all_ordered_cuts) {
      cut_distance += PathLength(cp.path);
    }
    printf("Total cut distance: %f\n", cut_distance);
    unsigned seconds = cut_distance / config.feed_rate * 60;
    const unsigned minutes = seconds / 60;
    seconds = seconds % 60;
    printf("Total cut time excluding moves: %u:%02u\n", minutes, seconds);
  }

  if (!config.output_ps_path.empty()) {
    // TODO: Replace this with output to .svg so that the drawing size can be correct.
    WriteCutsToPs(config.output_ps_path, all_ordered_cuts);
  }

  if (!config.output_nc_path.empty()) {
    WriteCutsToNc(config, all_ordered_cuts);
  }

  return EXIT_SUCCESS;
}
