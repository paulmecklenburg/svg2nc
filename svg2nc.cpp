// /usr/local/include/polyclipping/clipper.hpp

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

extern "C" {
#include "svgtiny.h"
}

namespace {
  struct Config {
    double diameter = -1.;
    std::map<uint32_t, double> color_to_height;
    // std::set<uint32_t> separators;
    std::string svg_path;
  };

  void PrintUsage(FILE *stream, const char *program_name, int exit_code) {
    fprintf(stream, "Usage:  %s [OPTION]... [FILE]\n", program_name);
    fprintf(stream,
            "  -h  --help             Display this usage information.\n"
            "  -d --diameter=DIAMETER Set the tool DIAMETER.\n"
            "  -c --color-height=<hex color>:<height>\n"
            "  -v  --verbose          Print verbose messages.\n");
    exit(exit_code);
  }

  void ParseColorHeightPair(const char *optarg,
                            std::map<uint32_t, double> *color_to_height) {
    uint32_t color;
    double height;
    if (2 == sscanf(optarg, "%x:%lg", &color, &height)) {
      (*color_to_height)[color] = height;
    } else {
      fprintf(stderr,
              "'%s' does not match the expected format "
              "<hex color>:<decimal height>",
              optarg);
      exit(EXIT_FAILURE);
    }
  }

  Config ParseArgs(int argc, char *argv[]) {
    const char *program_name = argv[0];
    Config config;
    int c;
    static const struct option long_options[] = {
      {"diameter", required_argument, NULL, 'd'},
      {"color-height", required_argument, NULL, 'c'},
      {"help", no_argument, NULL, 'h'},
      {NULL, 0, NULL, 0},
    };
    while ((c = getopt_long(
                argc, argv, "hd:c:", long_options, nullptr)) != -1) {
      switch (c) {
      case 'h':        
        PrintUsage(stderr, program_name, EXIT_SUCCESS);
        break;
      case '?':
        PrintUsage(stderr, program_name, EXIT_FAILURE);
        break;
      case 'd':
        config.diameter = atof(optarg);
        break;
      case 'c':
        ParseColorHeightPair(optarg, &config.color_to_height);
        break;
      }
    }

    if (optind != argc - 1) {
      PrintUsage(stderr, program_name, EXIT_FAILURE);
    }
    config.svg_path = argv[optind];
    return config;
  }

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
}  // namespace

int main(int argc, char *argv[]) {
  const Config config = ParseArgs(argc, argv);
  struct svgtiny_diagram *diagram = LoadSvg(config.svg_path.c_str());
  if (!diagram)
    return EXIT_FAILURE;

  svgtiny_free(diagram);
  return EXIT_SUCCESS;
}
