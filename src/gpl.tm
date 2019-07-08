:Begin:
:Function: gpl
:Pattern: G[a__]/;And @@ (NumberQ /@ ({a} /. SubPlus | SubMinus -> Identity))
:Arguments: {Re@N[{a}], Im@N[{a}]}
:ArgumentTypes: {RealList,RealList}
:ReturnType: Manual
:End:

#include "mathlink.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

typedef struct {double r,i;} complex;

extern complex __gpl_module_MOD_g_superflatn(complex*,long*);

void gpl(double * re, long nr, double * im, long ni)
{
    assert(nr==ni);
    complex input[nr], ans;
    for(long i=0;i<nr;i++)
    {
        input[i].r = *(re+i);
        input[i].i = *(im+i);
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


