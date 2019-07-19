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

GetChenFFGs[]:=Uncompress["1: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"]/.z->y/.xchen->x

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
  "checks/test-chenff.f90",  MakeFile[GetChenFFGs , "chenff" , {x, y}],
  "checks/test-muone.f90",   MakeFile[GetMuoneGs  , "muone"  , {x, y}],
  "checks/test-muoneNP.f90", MakeFile[GetMuoneNPGs, "muoneNP", {w, z}],
  _, Print["Availble are chen, muone, muonenp"];
];
