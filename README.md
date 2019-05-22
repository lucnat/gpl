
## Evaluate G(...) function

Open `src/eval.f90` and modify the line   

```res = GPL(cmplx([1,2,5]))```

This example computes G(1,2,5). Then write

```$ make eval```

```$ ./eval ```

and the result will be printed in standard output. To see what is going on increase verbosity

```$ ./eval -verb 100```

## Test 

Before running it is advised to run the tests

```$ make test```

```$ ./test ```
