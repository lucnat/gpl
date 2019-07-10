#include <ginac/ginac.h>
using namespace GiNaC;
#include <cln/cln.h>
#include <stdio.h>
#include <iostream>

typedef struct {double r,i;} complex_t;
extern "C"{
complex_t geval_(complex_t * z, int* n);
};
complex_t geval_(complex_t * z, int* n) {
  cln::cl_inhibit_floating_point_underflow = true;
  lst w;
  for(long i=0;i<(*n)-1;i++)
  {
    w.append((z->r)+(z->i)*I);
    z++;
  }
  ex ans = G(w,z->r).evalf();
  return {
    .r = ex_to<numeric>(evalf(real_part(ans))).to_double(),
    .i = ex_to<numeric>(evalf(imag_part(ans))).to_double()
  };

}
