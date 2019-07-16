GetChenGs[]:=Block[Evaluate[ToExpression["F"<>ToString[#]]&/@Range[40]],
  Import["https://arxiv.org/e-print/1801.01033", {"TAR",{"results.m"}}];
  Complement[
    DeleteDuplicates[Cases[ToExpression["F"<>ToString[#]]&/@Range[40],_G,-1]],{
    G[1,1],
    (* ginac has a problem with these *)
    G[y, y + Sqrt[-1 + y^2], y, 0, 1], 
    G[y - Sqrt[-1 + y^2], -1, y, 0, 1], 
    G[y - Sqrt[-1 + y^2], -1, y, 0, x], 
    G[y - Sqrt[-1 + y^2], 1, y, 0, x], 
    G[y - Sqrt[-1 + y^2], y + Sqrt[-1 + y^2], y, 0, x]

  }]
]

GetMuoneGs[]:=DeleteDuplicates[
 Cases[{
    Import["https://arxiv.org/src/1709.07435v2/anc/arXiv-fam1.m"],
    Import["https://arxiv.org/src/1709.07435v2/anc/arXiv-fam2.m"]
    }, _G, -1] /. G[{a__}, b_] :> G[a, b]
]

GetMuoneNPGs[]:=Block[{MInp, weightsRule},
  Import["https://arxiv.org/src/1806.08241v2/anc/MInp.m"];
  DeleteDuplicates@Cases[MInp/.weightsRule,_G,-1] /. G[{a__}, b_] :> G[a, b]
]

el2f[ 0]:="zero"
el2f[-1]:="mone"
el2f[+1]:="one"
el2f[ I]:="im"
el2f[-I]:="mim"
el2f[(1-I*Sqrt[3])/2 |-(-1)^(2/3)]:="degm60"
el2f[(1+I*Sqrt[3])/2 | (-1)^(1/3)]:="deg60"
el2f[n_]:=ToString[n,FortranForm]

G2Code[G[a__,b_]]:="(/ "<>StringRiffle[el2f/@{a},", "]<>",  "<>el2f[b]<>" /)"

SetBody[gs_]:=Table[
  "  args("<>ToString[i]<>",1:"<>ToString[Length[gs[[i]]]]<>") = "<>G2Code[gs[[i]]],
  {i,Length[gs]}
]
SetRoutine[gs_,name_,vars_]:=StringRiffle[Join[{
    "MODULE gtest"<>name,
    "contains",
    "  FUNCTION ARGS("<>StringRiffle[vars,","]<>")",
    "  use globals, only: prec",
    "  implicit none",
    "  complex(kind=prec), parameter :: zero = ( 0., 0.)",
    "  complex(kind=prec), parameter ::  one = ( 1., 0.)",
    "  complex(kind=prec), parameter :: mone = (-1., 0.)",
    "  complex(kind=prec), parameter ::   im = ( 0., 1.)",
    "  complex(kind=prec), parameter ::  mim = ( 0.,-1.)",
    "  complex(kind=prec), parameter :: deg60 = ("
          <>ToString@N@Re[(-1)^(1/3)]<>","<>ToString@N[Im[(-1)^(1/3)],10]<>")",
    "  complex(kind=prec), parameter :: degm60 = ("
          <>ToString@N@Re[-(-1)^(2/3)]<>","<>ToString@N[Im[-(-1)^(2/3)],10]<>")",
    "  complex(kind=prec) d_qnan",
    "  complex(kind=prec) :: "<>StringRiffle[vars,","],
    "  complex(kind=prec) :: args("<>ToString@Length[gs]<>","<>ToString@Max[Length/@gs]<>")",
    "  d_qnan = 0.",
    "  d_qnan = d_qnan/0",
    "  args = D_QNAN",
    ""
    },SetBody[gs],{
    "  END FUNCTION",
    "END MODULE"
  }],"\n"]


MakeFile[getter_, name_, vars_]:=Block[{gs,filename},
  filename="checks/test-"<>name<>".f90";
  WriteString["stdout", "Getting G functions for "<>name<>"..."];
  gs = AbsoluteTiming[getter[]];
  WriteString["stdout", " done, took "<>ToString[gs[[1]]]<>"s\n"];
  WriteString["stdout", "Writing "<>ToString[Length[gs[[2]]]]<>" G functions to file "<>filename<>"\n\n"];
  Export[filename, SetRoutine[gs[[2]],name,vars], "Text"];
]

Switch[Last[$CommandLine],
  "checks/test-chen.f90",    MakeFile[GetChenGs   , "chen"   , {x, y}],
  "checks/test-muone.f90",   MakeFile[GetMuoneGs  , "muone"  , {x, y}],
  "checks/test-muoneNP.f90", MakeFile[GetMuoneNPGs, "muoneNP", {w, z}],
  _, Print["Availble are chen, muone, muonenp"];
];
