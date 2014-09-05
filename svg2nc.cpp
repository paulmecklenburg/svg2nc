#include <map>
#include <set>
#include <string>

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "polyclipping/clipper.hpp"
extern "C" {
#include "svgtiny.h"
}

using ClipperLib::cInt;
using ClipperLib::ctUnion;
using ClipperLib::Clipper;
using ClipperLib::ClipperOffset;
using ClipperLib::ctDifference;
using ClipperLib::etClosedLine;
using ClipperLib::etClosedPolygon;
using ClipperLib::IntPoint;
using ClipperLib::jtRound;
using ClipperLib::ptSubject;
using ClipperLib::ptClip;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::SimplifyPolygons;

// TODO: Rename height to thickness and elevation.
// TODO: All cuts should be masked by the boundary.
// TOOD: Remove redundant segments (one above the other).

// By convention all paths are type EvenOdd.

#define MUST_USE_RESULT __attribute__((warn_unused_result))

namespace {
  struct Config {
    std::map<uint32_t, double> color_to_height;
    double diameter = -1.;
    double material_height = -1.;
    double non_overlap = .9;  // Needs a better name
    std::string svg_path;
  };

  void PrintUsage(FILE *stream, const char *program_name, int exit_code) {
    fprintf(stream, "Usage:  %s [OPTION]... [FILE]\n", program_name);
    fprintf(stream,
            "  -h  --help             Display this usage information.\n"
            "  -c --color-height=<hex color>:<height>\n"
            "  -d --diameter=DIAMETER Set the tool DIAMETER.\n"
            "  -m --material-height=H Material thickness.\n");
    exit(exit_code);
  }

  void ParseColorHeightPair(const char *optarg,
                            std::map<uint32_t, double> *color_to_height) {
    uint32_t color;
    double height;
    if (2 == sscanf(optarg, "%x:%lg", &color, &height)) {
      if (height <= 0) {
        fprintf(stderr, "Height values must be greater than zero.\n");
        exit(EXIT_FAILURE);
      }
      (*color_to_height)[color] = height;
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
      {"help", no_argument, NULL, 'h'},
      {"color-height", required_argument, NULL, 'c'},
      {"diameter", required_argument, NULL, 'd'},
      {"material-height", required_argument, NULL, 'm'},
      {NULL, 0, NULL, 0},
    };
    while ((c = getopt_long(
                argc, argv, "hc:d:m:", long_options, nullptr)) != -1) {
      switch (c) {
      case 'h':        
        PrintUsage(stderr, program_name, EXIT_SUCCESS);
        break;
      case '?':
        PrintUsage(stderr, program_name, EXIT_FAILURE);
        break;
      case 'c':
        ParseColorHeightPair(optarg, &config.color_to_height);
        break;
      case 'd':
        config.diameter = atof(optarg);
        break;
      case 'm':
        config.material_height = atof(optarg);
        break;
      }
    }

    if (optind != argc - 1) {
      PrintUsage(stderr, program_name, EXIT_FAILURE);
    }
    config.svg_path = argv[optind];
    if (config.material_height <= 0) {
      for (const auto &iter : config.color_to_height) {
        config.material_height = std::max(config.material_height, iter.second);
      }
      if (config.material_height <= 0) {
        fprintf(stderr,
                "Material height (thickness) must be greater than zero.\n");
        exit(EXIT_FAILURE);
      }
    }
    if (config.diameter <= 0) {
      fprintf(stderr, "Bit diameter must be specified and greater than zero.\n");
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

  bool ConvertPath(const float *, unsigned int, double, Paths *) MUST_USE_RESULT;
  bool ConvertPath(const float *input, unsigned int length, double svg_height, Paths *output) {
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
        // Invert y axis.
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
        // TOOD: Add support for bezier.
        fprintf(stderr, "(skipping) bezier..\n");
        if (!path || (j + 7) > length)
          return false;
        // http://cairographics.org/manual/cairo-Paths.html#cairo-curve-to
        // http://en.wikipedia.org/wiki/B%C3%A9zier_curve
        // cairo_curve_to(cr,
        //                scale * path->path[j + 1],
        //                scale * path->path[j + 2],
        //                scale * path->path[j + 3],
        //                scale * path->path[j + 4],
        //                scale * path->path[j + 5],
        //                scale * path->path[j + 6]);
        j += 7;
        break;
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
        const auto height_itr = config.color_to_height.find(shape.fill);
        if (height_itr == config.color_to_height.end()) {
          fprintf(stderr, "no height mapped for: %06X\n", shape.fill);
          return false;
        }
        Paths &layer_paths = (*layers)[height_itr->second];
        if (!UnionInto(paths, &layer_paths)) {
            fprintf(stderr, "unable merge paths in a common layer.\n");
            return false;
        }
      }
    }
    return true;
  }

  bool SameOrInside(const Path &, const Path &) MUST_USE_RESULT;
  bool SameOrInside(const Path &inner, const Path &outer) {
    for (const auto &pt : inner) {
      if (PointInPolygon(pt, outer) == 0)
        return false;  // point outside
    }
    return true;  // All points on or inside.
  }

  // Assuming all is one or more polygons (which may have holes). Find one
  // disjoint polygon and copy it and any corresponding holes to island.
  // bool FindIsland(const Paths &all, Paths *island) {
  //   island->clear();
  //   Paths unn;
  //   for (const Path &path : all) {
  //     Paths paths{path};
  //     if (!UnionInto(paths, &unn))
  //       return false;
  //   }
  //   if (unn.empty())
  //     return false;
  //   // unn is now a set of polygons without any holes.
  //   for (const Path &path : all) {
  //     if (SameOrInside(path, unn[0])) {
  //       island->push_back(path);
  //     }
  //   }
  //   return true;
  // }

  void ComputeOffset(const Paths &paths, double amount, Paths *result) {
    ClipperOffset co;
    co.ArcTolerance = kQuantaPerInch / 1000;
    // TODO: Do we need to check orientations here?
    co.AddPaths(paths, jtRound, etClosedPolygon);
    co.Execute(*result, InchesToQuanta(amount));
  }

  // Offset a path to a polygon.
  void CutToPolygon(const Paths &cut, double radius, Paths *result) {
    ClipperOffset co;
    co.ArcTolerance = kQuantaPerInch / 1000;
    // TODO: Do we need to check orientations here?
    co.AddPaths(cut, jtRound, etClosedLine);
    co.Execute(*result, InchesToQuanta(radius));
  }

  bool SubtractFrom(const Paths &, Paths *) MUST_USE_RESULT;
  bool SubtractFrom(const Paths &a, Paths *b) {
    Clipper c;
    return c.AddPaths(*b, ptSubject, true) &&
      c.AddPaths(a, ptClip, true) &&
      c.Execute(ctDifference, *b);
  }

  void MillEdges(const Config &config,
                 const Paths &higher_layers,
                 Paths *cuts) {
    Paths edge_cuts;
    ComputeOffset(higher_layers, config.diameter * .5, &edge_cuts);
    cuts->insert(cuts->end(), edge_cuts.begin(), edge_cuts.end());
  }

  bool MillSurface(const Config &,
                   const Paths &,
                   Paths *) MUST_USE_RESULT;
  bool MillSurface(const Config &config,
                   const Paths &surface,
                   Paths *cuts) {
    Paths remaining(surface);
    while (!remaining.empty()) {
      // Shrink by less than radius.
      Paths cut;
      ComputeOffset(remaining, -(config.diameter * .5 * config.non_overlap), &cut);

      // In some cases an area smaller than the bit is left and cut is empty.
      // To make sure this material is removed we cut the outer edge. This is
      // less efficient than is optimal, but it's a rare, small case.
      if (cut.empty())
        cut.swap(remaining);

      cuts->insert(cuts->end(), cut.begin(), cut.end());

      Paths mask;
      CutToPolygon(cut, config.diameter * .5, &mask);
      if (!SubtractFrom(mask, &remaining))
        return false;
    }
    return true;
  }

  bool ComputeCuts(const Config &,
                   const std::map<double, Paths> &,
                   std::map<double, Paths> *) MUST_USE_RESULT;
  bool ComputeCuts(const Config &config,
                   const std::map<double, Paths> &layers,
                   std::map<double, Paths> *height_to_cuts) {
    Paths higher_layer_union;
    for (auto iter = layers.rbegin(); iter != layers.rend(); ++iter) {
      const double height = iter->first;
      const Paths &polygons = iter->second;
      if (height < config.material_height) {
        Paths mill_area = polygons;
        Paths *cuts = &(*height_to_cuts)[height];
        if (!higher_layer_union.empty()) {
          // Mill the edges of the layer above.
          MillEdges(config, higher_layer_union, cuts);

          Paths tmp;
          ComputeOffset(higher_layer_union, config.diameter, &tmp);
          if (!SubtractFrom(tmp, &mill_area))
            return false;
        }
        if (!MillSurface(config, mill_area, cuts))
          return false;
      }
      if (!UnionInto(polygons, &higher_layer_union))
        return false;
    }
    MillEdges(config, higher_layer_union, &(*height_to_cuts)[0]);
    return true;
  }

  void ToDebugPs(const std::map<double, Paths> &height_to_cuts) {
    FILE *fp = fopen("debug.ps", "w");
    if (!fp)
      return;

    for (const auto &layer : height_to_cuts) {
      const double color = 1. - layer.first;
      for (const auto &cut : layer.second) {
        fprintf(fp, "0.2 setlinewidth\n");
        for (const auto &pt : cut) {
          fprintf(fp, "%f %f %s\n",
                  QuantaToInches(pt.X) * 72,
                  QuantaToInches(pt.Y) * 72,
                  &pt == &*cut.begin() ? "newpath moveto" : "lineto");
        }
        fprintf(fp, "closepath\n");
        fprintf(fp, "%f setgray\n", color);
        fprintf(fp, "stroke\n");
      }
    }

    fclose(fp);
  }
}  // namespace

int main(int argc, char *argv[]) {
  const Config config = ParseArgs(argc, argv);
  struct svgtiny_diagram *diagram = LoadSvg(config.svg_path.c_str());
  if (!diagram)
    return EXIT_FAILURE;

  std::map<double, Paths> layers;
  bool r = SvgToPolygons(*diagram, config, &layers);
  svgtiny_free(diagram);

  for (auto iter : layers) {
    SimplifyPolygons(iter.second);
  }
  std::map<double, Paths> height_to_cuts;
  r = r && ComputeCuts(config, layers, &height_to_cuts);
  ToDebugPs(height_to_cuts);
  return r ? EXIT_SUCCESS : EXIT_FAILURE;
}
