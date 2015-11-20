#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include "polyclipping/clipper.hpp"
#include "common.h"
#include "path_util.h"
#include "svg.h"
#include "uniq_seg.h"

using ClipperLib::Area;
using ClipperLib::cInt;
using ClipperLib::Clipper;
using ClipperLib::ClipperOffset;
using ClipperLib::ctIntersection;
using ClipperLib::etClosedLine;
using ClipperLib::etOpenButt;
using ClipperLib::IntPoint;
using ClipperLib::jtRound;
using ClipperLib::jtSquare;
using ClipperLib::ptSubject;
using ClipperLib::ptClip;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::PolyTree;

namespace {
  struct Config {
    bool as_drawn = false;
    std::map<uint32_t, double> color_to_elevation;
    std::set<uint32_t> ignore_colors;
    double clearance_space = .25;
    double diameter = -1.;
    double feed_rate = 12.;
    double material_thickness = -1.;
    double mill_overlap = .1;
    double max_pass_depth = .25;
    double spindle_speed = 1000.;
    double through_elevation = 0;
    double blend_tolerance = 0;
    bool verbose = false;
    std::string svg_path;
    std::string output_svg_path;
    std::string output_nc_path;
  };

  void PrintUsage(FILE *stream, const char *program_name, int exit_code) {
    fprintf(stream, "Usage:  %s [OPTION]... [SVG FILE]\n", program_name);
    fprintf(
        stream,
        "  -h  --help                 Display this usage information.\n"
        "  -a --as-drawn              Allow holes and pockets to be drawn on.\n"
        "  -b --blend-tolerance=<inches> Allow paths to deviate.\n"
        "  -c --color-elevation=<hex color>:<inches>\n"
        "     Specify the elevation for a color.\n"
        "  -d --diameter=<inches>     Set the tool DIAMETER.\n"
        "  -e --hole-color=<hex color> Color for drawn-on through holes.\n"
        "  -f --feed-rate=<inches per minute>\n"
        "     Specify the feed rate.\n"
        "  -g --svg-file=<path>        Write the cut plan to a .svg file.\n"
        "  -i --ignore-color=<hex color> Ignore this color.\n"
        "  -l --clearance-space=<inches>\n"
        "     Safety space to use above the part while moving.\n"
        "  -m --material-thickness=<inches>\n"
        "     Specify the material thickness.\n"
        "  -n --nc-file=<path>        Write the cut plan to a .nc file.\n"
        "  -o --mill-overlap=<fraction>\n"
        "     Fraction of a cut to overlap with adjacent cuts.\n"
        "  -s --spindle_speed=<rpm>   Set the spindle speed.\n"
        "  -t --through-elevatation=<inches>\n"
        "     The elevation to use while cutting holes and outlines.\n"
        "  -v --verbose               Enable verbose output.\n"
        "  -x --max-pass-depth=<depth>\n"
        "     The maximum depth to cut in a single pass.\n"
    );
    exit(exit_code);
  }

  uint32_t ParseColor(const char *arg) {
    uint32_t color;
    if (1 == sscanf(arg, "%x", &color)) {
      return color;
    } else {
      fprintf(stderr,
              "'%s' does not match the expected format <hex color>\n",
              arg);
      exit(EXIT_FAILURE);
    }
  }

  void ParseColorElevationPair(const char *arg,
                               std::map<uint32_t, double> *color_to_elevation) {
    uint32_t color;
    double elevation;
    if (2 == sscanf(arg, "%x:%lg", &color, &elevation)) {
      if (elevation <= 0) {
        fprintf(stderr, "Elevation values must be greater than zero.\n");
        exit(EXIT_FAILURE);
      }
      (*color_to_elevation)[color] = elevation;
    } else {
      fprintf(stderr,
              "'%s' does not match the expected format "
              "<hex color>:<decimal height>\n",
              arg);
      exit(EXIT_FAILURE);
    }
  }

  Config ParseArgs(int argc, char *argv[]) {
    const char *program_name = argv[0];
    Config config;
    int c;
    std::string hole_color;
    static const struct option long_options[] = {
      {"as-drawn", no_argument, nullptr, 'a'},
      {"blend-tolerance", required_argument, nullptr, 'b'},
      {"clearance-space", required_argument, nullptr, 'l'},
      {"color-elevation", required_argument, nullptr, 'c'},
      {"diameter", required_argument, nullptr, 'd'},
      {"feed-rate", required_argument, nullptr, 'f'},
      {"help", no_argument, nullptr, 'h'},
      {"hole-color", required_argument, nullptr, 'e'},
      {"ignore-color", required_argument, nullptr, 'i'},
      {"material-thickness", required_argument, nullptr, 'm'},
      {"max-pass-depth", required_argument, nullptr, 'x'},
      {"mill-overlap", required_argument, nullptr, 'o'},
      {"nc-file", required_argument, nullptr, 'n'},
      {"spindle-speed", required_argument, nullptr, 's'},
      {"svg-file", required_argument, nullptr, 'g'},
      {"through-elevation", required_argument, nullptr, 't'},
      {"verbose", no_argument, nullptr, 'v'},
      {nullptr, 0, nullptr, 0},
    };
    while ((c = getopt_long(
                argc, argv, "ab:c:d:e:f:g:hi:l:m:n:o:r:s:t:vx:",
                long_options, nullptr)) != -1) {
      switch (c) {
      case '?':
        PrintUsage(stderr, program_name, EXIT_FAILURE);
        break;
      case 'a':
        config.as_drawn = true;
        break;
      case 'b':
        config.blend_tolerance = atof(optarg);
        break;
      case 'c':
        ParseColorElevationPair(optarg, &config.color_to_elevation);
        break;
      case 'd':
        config.diameter = atof(optarg);
        break;
      case 'e':
        hole_color = optarg;
        break;
      case 'f':
        config.feed_rate = atof(optarg);
        break;
      case 'g':
        config.output_svg_path = optarg;
        break;
      case 'h':
        PrintUsage(stderr, program_name, EXIT_SUCCESS);
        break;
      case 'i':
        config.ignore_colors.insert(ParseColor(optarg));
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
    if (!hole_color.empty()) {
      if (config.as_drawn) {
        const uint32_t color = ParseColor(hole_color.c_str());
        config.color_to_elevation[color] = config.through_elevation;
      } else {
        fprintf(stderr, "--hole-color requires --as-drawn.\n");
        exit(EXIT_FAILURE);
      }
    }
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

  void CutToPolygon(const Paths &cut, double radius, Paths *result) {
    // Offset a path to a polygon representing the removed material.
    ClipperOffset co;
    co.ArcTolerance = kQuantaPerInch / 1000;
    co.AddPaths(cut, jtRound, etClosedLine);
    co.Execute(*result, InchesToQuanta(radius));
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

      cuts->insert(cuts->end(), cut.begin(), cut.end());

      if (!remaining.empty()) {
        Paths mask;
        CutToPolygon(cut, config.diameter * .5, &mask);
        if (!SubtractFrom(mask, &remaining))
          return false;
      }
    }
    return true;
  }

  struct CutLoop;
  struct Segment {
    CutLoop *parent;
    // The open path for this particular cut.
    Path path;
  };

  struct CutLoop {
    double elevation;
    // Closed path. Full loop.
    Path path;
    // Union of all polygons with elevations above this. May be shard with other paths.
    std::shared_ptr<Paths> mask;
    // List of other CutLoops blocked by this one.
    std::vector<CutLoop*> blocks;
    // Sub-segments of this loop that actually need to be cut.
    std::vector<Segment> segments;
  };

  void ToCutLoops(const Paths &in_cuts,
                  double elevation,
                  const Paths &higher_layers,
                  std::vector<CutLoop> *out_cuts) {
    std::shared_ptr<Paths> mask(new Paths(higher_layers));
    for (const Path &c : in_cuts) {
      out_cuts->push_back(CutLoop{elevation, c, mask, {}, {}});
    }
  }

  void ComputeCuts(const Config &config,
                   const std::map<double, Paths> &layers,
                   std::vector<CutLoop> *cuts) {
    Paths higher_layer_union;
    for (auto iter = layers.rbegin(); iter != layers.rend(); ++iter) {
      const double elevation = iter->first;
      const Paths &polygons = iter->second;
      if (elevation < config.material_thickness) {
        Paths mill_area = polygons;
        Paths tmp_cuts;
        if (!higher_layer_union.empty()) {
          // Mill the edges of the layer above.
          OR_DIE(MillEdges(config, higher_layer_union, &tmp_cuts));
          Paths tmp;
          OR_DIE(ComputeOffset(higher_layer_union, config.diameter, &tmp));
          OR_DIE(SubtractFrom(tmp, &mill_area));
        }
        OR_DIE(MillSurface(config, mill_area, &tmp_cuts));
        ToCutLoops(tmp_cuts, elevation, higher_layer_union, cuts); 
      }
      OR_DIE(UnionInto(polygons, &higher_layer_union));
    }
    Paths tmp_cuts;
    OR_DIE(MillEdges(config, higher_layer_union, &tmp_cuts));
    ToCutLoops(tmp_cuts, config.through_elevation, higher_layer_union, cuts); 
  }

  void InsertCutLoop(std::map<CutLoop*, std::vector<CutLoop*> > *roots,
                     CutLoop *cl) {
    std::vector<CutLoop*> tmp;
    std::vector<CutLoop*> to_erase;
    for (auto &iter : *roots) {
      Clipper c;
      Paths intersection;
      OR_DIE(c.AddPath(cl->path, ptSubject, true) &&
             c.AddPath(iter.first->path, ptClip, true) &&
             c.Execute(ctIntersection, intersection));
      if (intersection.empty())
        continue;
      const double cla = fabs(Area(cl->path));
      const double pa = fabs(Area(iter.first->path));
      if (pa >= cla) {
        OR_DIE(tmp.empty());
        iter.second.push_back(cl);
        return;
      }
      tmp.insert(tmp.end(), iter.second.begin(), iter.second.end());
      tmp.push_back(iter.first);
      to_erase.push_back(iter.first);
    }
    for (auto cl : to_erase)
      roots->erase(cl);
    (*roots)[cl] = tmp;
  }

  bool IsDelayed(const Path &path,
                 const double radius,
                 const std::vector<IntPoint> &delay_points) {
    const Paths paths{path};
    Paths approx_orig;
    OR_DIE(ComputeOffset(paths, radius, &approx_orig));
    OR_DIE(approx_orig.size() == 1);
    for (const auto &pt : delay_points) {
      if (PointInPolygon(pt, approx_orig[0])) {
        return true;
      }
    }
    return false;
  }

  std::vector<CutLoop*> ComputeBlocks(
      std::vector<CutLoop*> *cut_loops,
      const double radius,
      const std::vector<IntPoint> &delay_points) {
    // Generate a map of disjoint cuts and a list of cuts inside each.
    std::map<CutLoop*, std::vector<CutLoop*> > roots;
    for (CutLoop *cl : *cut_loops) {
      InsertCutLoop(&roots, cl);
    }
    // Recursively find the disjoint cuts within each child. Update children
    // to block their parents.
    for (auto &iter : roots) {
      for (auto &grand_child
             : ComputeBlocks(&iter.second, radius, delay_points)) {
        grand_child->blocks.push_back(iter.first);
      }
    }
    std::vector<CutLoop*> delayed;
    std::vector<CutLoop*> nondelayed;
    for (auto &iter : roots) {
      (iter.second.empty() && IsDelayed(iter.first->path, radius, delay_points)
       ? delayed : nondelayed).push_back(iter.first);
    }
    for (auto nd : nondelayed) {
      nd->blocks.insert(nd->blocks.end(), delayed.begin(), delayed.end());
    }
    // Return a list of the children.
    std::vector<CutLoop*> children;
    for (auto &iter : roots) {
      children.push_back(iter.first);
    }
    return children;
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

  bool TrimCutToBoundingBox(cInt, cInt, const Path &in, Paths *out) MUST_USE_RESULT;
  bool TrimCutToBoundingBox(cInt width, cInt height, const Path &in, Paths *out) {
    out->clear();
    const Path bbox{{0, 0}, {width, 0}, {width, height}, {0, height}};
    // This test is done to work around a clipper bug where segments are
    // dropped from open paths with the final point the same as the starting
    // point in the *output*. This only happens with paths that aren't clipped
    // at all.
    if (OnOrInside(in, bbox)) {
      out->push_back(in);
      return true;
    }
    Clipper c;
    if (!c.AddPath(in, ptSubject, false))
      return false;
    PolyTree solution;
    if (!c.AddPath(bbox, ptClip, true) ||
        !c.Execute(ctIntersection, solution))
      return false;
    OpenPathsFromPolyTree(solution, *out);
    return true;    
  }

  void ComputeSegments(cInt width, cInt height, std::vector<CutLoop*> *cut_loops) {
    const auto ElevationLess = [](const CutLoop *a, const CutLoop *b) -> bool {
      return a->elevation < b->elevation;
    };
    std::sort(cut_loops->begin(), cut_loops->end(), ElevationLess);

    UniqSeg uniq_seg;
    for (CutLoop *cl : *cut_loops) {
      const Path &path = cl->path;
      Paths segments;
      for (size_t i = 0; i < path.size(); ++i) {
        Paths tmp;
        uniq_seg.RemoveRedundant(path[i], path[(i + 1) % path.size()], &tmp);
        for (const Path &p : tmp) {
          // TODO: encode this in the data type.
          OR_DIE(p.size() == 2);
          if (!segments.empty() && segments.back().back() == p.front()) {
            segments.back().push_back(p.back());
          } else {
            segments.push_back(p);
          }
        }
      }
      if (segments.size() > 1 &&
          segments.front().front() == segments.back().back()) {
        // The end of the last segment connects back to the beginning of the
        // first segment. Merge them.
        segments.back().insert(segments.back().end(),
                               segments.front().begin() + 1,
                               segments.front().end());
        segments.front().swap(segments.back());
        segments.pop_back();
      }
      for (const auto &segment : segments) {
        Paths tmp;
        OR_DIE(TrimCutToBoundingBox(width, height, segment, &tmp));
        for (const Path &p : tmp) {
          cl->segments.push_back(Segment{cl, p});
        }
      }
    }
  }

  void CollapseBlocks(const CutLoop *cl, std::vector<CutLoop*> *new_blocks) {
    for (const auto target : cl->blocks) {
      if (target->segments.empty()) {
        CollapseBlocks(target, new_blocks);
      } else {
        new_blocks->push_back(target);
      }
    }
  }

  void DropEmptyCuts(std::vector<CutLoop*> *cut_loops) {
    std::vector<CutLoop*> tmp;
    for (const auto cl : *cut_loops) {
      if (!cl->segments.empty()) {
        std::vector<CutLoop*> new_blocks;
        CollapseBlocks(cl, &new_blocks);
        cl->blocks.swap(new_blocks);
        tmp.push_back(cl);
      }
    }
    cut_loops->swap(tmp);
  }

  cInt DistanceSquared(const IntPoint &a, const IntPoint &b) {
    const cInt dx = a.X - b.X;
    const cInt dy = a.Y - b.Y;
    return dx*dx + dy*dy;
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

  Segment *RemoveClosestSegment(const IntPoint &last_pos,
                                std::set<Segment*> *segs,
                                size_t *min_index) {
    cInt min_d2 = std::numeric_limits<cInt>::max();
    Segment *min_seg = nullptr;
    for (const auto &seg : *segs) {
      size_t index;
      const cInt d2 = ClosestPointInPath(last_pos, seg->path, &index);
      if (d2 < min_d2) {
        min_d2 = d2;
        *min_index = index;
        min_seg = seg;
      }
    }
    segs->erase(min_seg);
    return min_seg;
  }

  IntPoint LastPosition(std::vector<const Segment*> &history) {
    if (history.empty()) {
      return IntPoint{0, 0};
    } else {
      return history.back()->path.back();
    }
  }

  void ComputeOrdering(const std::vector<CutLoop*> &cut_loops,
                       std::vector<const Segment*> *segments) {
    std::map<const CutLoop*, int> block_count;
    for (const auto cl : cut_loops) {
      for (const auto target : cl->blocks) {
        if (!cl->segments.empty()) {
          block_count[target] += cl->segments.size();
        }
      }
    }
    std::set<Segment*> open;
    for (auto cl : cut_loops) {
      if (!block_count.count(cl)) {
        for (auto &seg : cl->segments) {
          open.insert(&seg);
        }
      }
    }
    while (!open.empty()) {
      size_t index;
      Segment *next =
        RemoveClosestSegment(LastPosition(*segments), &open, &index);
      if (index) {
        Path tmp(next->path.begin() + index, next->path.end());
        tmp.insert(tmp.end(), next->path.begin() + 1, next->path.begin() + index + 1);
        next->path.swap(tmp);
      }
      segments->push_back(next);
      for (CutLoop *target : next->parent->blocks) {
        if (!(--block_count[target])) {
          for (Segment &seg : target->segments) {
            open.insert(&seg);
          }
        }
      }
    }
    for (const auto &iter : block_count) {
      OR_DIE(!iter.second);
    }
  }

  bool CanSlide(const Segment &prev, const Segment &seg, double radius) {
    const auto &start = prev.path.back();
    const auto &end = seg.path.front();
    if (start == end)
      return true;
    const double kMaxSlideDistance = radius * 4;
    if (QuantaToInches(sqrt(DistanceSquared(start, end))) > kMaxSlideDistance)
      return false;
    const Paths &layers_above = prev.parent->elevation > seg.parent->elevation ?
      *prev.parent->mask : *seg.parent->mask;
    if (layers_above.empty())
      return true;
    Paths swept_path;
    {  // Create a swept path for the tool. Ideal would use 'jtRound' with
      // exactly radius. Instead etOpenButt is used. We assume that the start
      // and end are already safe and make the path wider (kEpsilon) to give
      // a bit of extra space. Meatspace is not as precise as math.
      ClipperOffset co;
      const Path line{start, end};
      co.AddPath(line, jtSquare /* doesn't actually matter */, etOpenButt);
      const double kEpsilon = 0.005;
      co.Execute(swept_path, InchesToQuanta(radius + kEpsilon));
    }
    Clipper c;
    Paths intersection;
    OR_DIE(c.AddPaths(layers_above, ptSubject, true) &&
           c.AddPaths(swept_path, ptClip, true) &&
           c.Execute(ctIntersection, intersection));
    return intersection.empty();
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

  double PathLength(const Path &cut) {
    double total = 0;
    for (size_t i = 0; i + 1 < cut.size(); ++i) {
      total += sqrt(DistanceSquared(cut[i], cut[i+1]));
    }
    return QuantaToInches(total);
  }

  struct CutPath {
    double elevation;
    Path path;
    bool slide;
  };

  double QuantaToSvg(cInt x) {
    return QuantaToInches(x) * kSvgUnitsPerInch;
  }

  void WriteSvgPath(FILE *fp, const Path &path, const char *color_str,
                    double diameter, double opacity, bool arrow, cInt height) {
    fprintf(fp,
            "  <path style=\"stroke-width:%.5f;stroke:%s;fill:none;"
            "stroke-linecap:round;stroke-linejoin:round;stroke-opacity:%.3f\" "
            "%sd=\"M", diameter * kSvgUnitsPerInch, color_str, opacity,
            arrow ? "marker-end=\"url(#Triangle)\" " : "");
    for (const auto &pt : path) {
      fprintf(fp, " %.5f,%.5f", QuantaToSvg(pt.X), QuantaToSvg(height - pt.Y));
    }
    fprintf(fp, "\" />\n");
  }

  void WriteCutsToSvg(const std::string &path,
                      const std::vector<CutPath> &ordered_cuts,
                      const cInt width, const cInt height,
                      double diameter) {
    FILE *fp = fopen(path.c_str(), "w");
    if (!fp) {
      perror(path.c_str());
      exit(EXIT_FAILURE);
    }

    const int svg_width = int(QuantaToSvg(width));
    const int svg_height = int(QuantaToSvg(height));
    fprintf(fp,
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
            "<svg width=\"%d\" height=\"%d\" version=\"1.1\">\n"
            "  <rect style=\"fill:white\" width=\"%d\" height=\"%d\" />\n"
            "  <defs>\n"
            "    <marker id=\"Triangle\"\n"
            "      viewBox=\"0 0 6 6\" refX=\"6\" refY=\"3\"\n"
            "      markerUnits=\"strokeWidth\"\n"
            "      markerWidth=\"4\" markerHeight=\"3\"\n"
            "      orient=\"auto\">\n"
            "      <path d=\"M 0,0 6,3 0,6 z\" />\n"
            "    </marker>\n"
            "  </defs>\n",
            svg_width, svg_height, svg_width, svg_height);

    IntPoint last(0, 0);
    double color = 0;
    double dc = 255. / (ordered_cuts.size() - 1);
    for (const auto &cp : ordered_cuts) {
      // Move
      if (last != cp.path.front()) {
        Path move{last, cp.path.front()};
        WriteSvgPath(fp, move, cp.slide ? "yellow" : "red", .01, 1., true, height);
      }
      // Cut
      char hex_color[10];
      sprintf(hex_color, "#00%02X%02X", int(color), 255 - int(color));
      WriteSvgPath(fp, cp.path, hex_color, diameter, .3, false, height);
      WriteSvgPath(fp, cp.path, hex_color, .01, 1., false, height);
      last = cp.path.back();
      color += dc;
    }

    fprintf(fp, "</svg>\n");
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
            "G20 (imperial)\n", config.svg_path.c_str());
    if (config.blend_tolerance > 0.) {
      fprintf(fp, "G64 P%f (blend path mode)\n",
              config.blend_tolerance);
    } else {
      fprintf(fp, "G61 (exact path mode)\n");
    }
    fprintf(fp,
            "F%.3f (feed rate)\n"
            "S%.1f (spindle speed)\n"
            "M3 (start spindle)\n"
            "G04 p1 (wait for 1 second)\n",
            config.feed_rate, config.spindle_speed);

    const double safe_elevation =
      config.material_thickness + config.clearance_space;
    IntPoint last{0, 0};
    double last_elevation = 0;
    for (const auto &cp : ordered_cuts) {
      const auto &first = cp.path.front();
      if (first != last) {
        fprintf(fp, "G%s Z%f\n", cp.slide ? "1" : "0" , cp.slide ?
                std::max(last_elevation, cp.elevation) + .02 : safe_elevation);
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
      last_elevation = cp.elevation;
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

  std::map<double, Paths> layers;
  std::vector<IntPoint> delay_points;
  cInt width, height;
  if (!SvgToPolygons(config.svg_path.c_str(),
                     config.color_to_elevation, config.ignore_colors,
                     config.as_drawn, &layers, &delay_points,
                     &width, &height)) {
    return EXIT_FAILURE;
  }

  std::vector<CutLoop> cut_loops_storage;
  ComputeCuts(config, layers, &cut_loops_storage);

  std::vector<CutLoop*> cut_loops;
  for (CutLoop &cl : cut_loops_storage)
    cut_loops.push_back(&cl);

  ComputeBlocks(&cut_loops, config.diameter / 2, delay_points);
  ComputeSegments(width, height, &cut_loops);
  DropEmptyCuts(&cut_loops);

  std::vector<const Segment*> segments;
  ComputeOrdering(cut_loops, &segments);

  std::vector<CutPath> all_ordered_cuts;
  const Segment *prev = nullptr;
  for (const Segment *seg : segments) {
    CutPath cp;
    cp.path = seg->path;
    cp.slide = prev && CanSlide(*prev, *seg, config.diameter / 2);
    for (const double elevation : PassElevations(config, seg->parent->elevation)) {
      cp.elevation = elevation;
      all_ordered_cuts.push_back(cp);
      cp.slide = true;
    }
    prev = seg;
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

  if (!config.output_svg_path.empty()) {
    WriteCutsToSvg(config.output_svg_path, all_ordered_cuts,
                   width, height, config.diameter);
  }

  if (!config.output_nc_path.empty()) {
    WriteCutsToNc(config, all_ordered_cuts);
  }

  return EXIT_SUCCESS;
}
