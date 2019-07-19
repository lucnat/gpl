
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


## Verbosity legend

10 tell if we need to make convergent

30 give some details about how it is made convergent

50 print whenever G(...) is called

70 print on what polylogs it is reduced

## License
GPL v3

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

