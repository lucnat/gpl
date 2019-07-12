:Evaluate:
  gpl`args2r[a_]:=Re[N[a/.SubPlus|SubMinus->Identity]];
  gpl`args2i[a_]:=Im[N[a/.SubPlus|SubMinus->Identity]];
  gpl`args2e[a_]:=Switch[Head[#], SubPlus, 1, SubMinus, -1, _, 1]& /@ a;
:Begin:
:Function: gpl
:Pattern: gG[a__]/;And @@ (NumberQ /@ ({a} /. SubPlus | SubMinus -> Identity))
:Arguments: { gpl`args2r[{a}], gpl`args2i[{a}], gpl`args2e[{a}] }
:ArgumentTypes: {RealList,RealList,IntegerList}
:ReturnType: Manual
:End:

#include "mathlink.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

typedef struct {double r,i;} complex;
typedef struct {complex c; signed char i0;} inum;

extern complex __gpl_module_MOD_g_superflatn(inum*,long*);

void gpl(double * re, long nr, double * im, long ni, int*ieps, long ne)
{
    assert(nr==ni);
    assert(nr==ne);
    inum input[nr];
    complex ans;
    for(long i=0;i<nr;i++)
    {
        input[i].c.r = *(re+i);
        input[i].c.i = *(im+i);
        input[i].i0  = *(ieps+i);
    }
    ans = __gpl_module_MOD_g_superflatn(&input[0],&nr);
    
    
    if(ans.i == 0)
        MLPutReal(stdlink, ans.r);
    else
    {
        MLPutFunction(stdlink, "Complex", 2);
        MLPutReal(stdlink, ans.r);
        MLPutReal(stdlink, ans.i);
    }
}


int main(int argc, char **argv)
{
    return MLMain(argc, argv);
}


