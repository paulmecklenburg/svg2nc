#ifndef _COMMON_H
#define _COMMON_H

#include <stdio.h>

#include "polyclipping/clipper.hpp"

#define MUST_USE_RESULT __attribute__((warn_unused_result))

const ClipperLib::cInt kQuantaPerInch = 100000;

inline ClipperLib::cInt InchesToQuanta(double x) {
  return x * kQuantaPerInch;
}

inline double QuantaToInches(ClipperLib::cInt x) {
  return double(x) / kQuantaPerInch;
}

#define OR_DIE(exp) \
  if (!(exp)) { \
    fprintf(stderr, "DIE: %s:%d %s\n", __FILE__, __LINE__, #exp); \
    exit(EXIT_FAILURE); \
  }

#endif  // _COMMON_H
