(* ::Package:: *)

Install["LoopTools"]; (*Using LoopTools, it can automatically handle PaVe functions and change them into the LoopTool convention*)


(*
Reading the Simplied Couplings Data, and return the whole expression
Or you can just use Get[file] to get the form stored in the file
*)
ReadingSimplifiedCouplings[file_]:=Block[{expori,expFree,expRePart,expPVPart,expReFunc,expReCoeff,expPVFunc,expPVCoeff},
expori=Get[file];
expFree=expori[[1]];
expRePart = expori[[2]];
expPVPart = expori[[3]];
expReFunc=expRePart[[1]];
expReCoeff=expRePart[[2]];
expPVFunc=expPVPart[[1]];
expPVCoeff=expPVPart[[2]];
expFree+expReFunc.expReCoeff+expPVFunc.expPVCoeff
];


(*
Seperate the original expressions into three parts
*)
ExtractPatternCoeff[expin_,pattern_List]:=Block[{elems,findings,freerules,frees,coeff},
elems=Union[Level[expin,Infinity]];
findings=Union@@(Cases[elems,#]&/@pattern);
freerules=(#->0)&/@findings;
frees=expin/.freerules;
coeff=Table[Coefficient[expin,findings[[i]]],{i,1,Length[findings]}];
{frees,findings,coeff}
];
ExtractPatternCoeff[expin_,patternFirst_List,patternSecond_List]:=Block[{first,second},
	first=ExtractPatternCoeff[expin,patternFirst];
  second=ExtractPatternCoeff[first[[1]],patternSecond];
  (*{second[[1]],Join[first[[2]],second[[2]]],Join[first[[3]],second[[3]]]}*)
  {second[[1]],{first[[2]],first[[3]]},{second[[2]],second[[3]]}}
];
ExtractPatternCoeffWithSimplify[expin_,pattern_List]:=Block[{elems,findings,freerules,frees,coeff},
elems=Union[Level[expin,Infinity]];
findings=Union@@(Cases[elems,#]&/@pattern);
freerules=(#->0)&/@findings;
frees=expin/.freerules;
coeff=Table[Simplify[Coefficient[expin,findings[[i]]]],{i,1,Length[findings]}];
{frees,findings,coeff}
];


(*
The original expression has explicit dependence on D,
Use following fuctions to
Limit the dimension to 4
*)
PVReplaceRule={
A0[m_]:>m/EPS+A0[m],
A00[m_]:>m^2/(4EPS)+A00[m],
B0[args___]:>1/EPS+B0[args],
B0i[bb0,args___]:>1/EPS+B0i[bb0,args],
B0i[bb1,args___]:>-1/(2 EPS)+B0i[bb1,args],
B0i[bb00,p_,m1_,m2_]:>((m1+m2)/4-p/12)/EPS+B0i[bb00,p,m1,m2],
B0i[bb11,args___]:>1/(3EPS)+B0i[bb11,args],
B0i[bb001,p_,m1_,m2_]:>(p-2m1-4m2)/(24EPS)+B0i[bb001,p,m1,m2],
B0i[bb111,args___]:>-1/(4 EPS)+B0i[bb111,args],
B0i[dbb00,args___]:>-1/(12 EPS)+B0i[dbb00,args],
C0i[cc00,args___]:>1/(4EPS)+C0i[cc00,args],
C0i[cc001,args___]:>-1/(12EPS)+C0i[cc001,args],
C0i[cc002,args___]:>-1/(12EPS)+C0i[cc002,args],
C0i[cc0000,p1_,p2_,p1p2_,m1_,m2_,m3_]:>-(p1+p2+p1p2-4(m1+m2+m3))/(96EPS)+C0i[cc0000,p1,p2,p1p2,m1,m2,m3],
C0i[cc0011,args___]:>1/(24EPS)+C0i[cc0011,args],
C0i[cc0022,args___]:>1/(24EPS)+C0i[cc0022,args],
C0i[cc0012,args___]:>1/(48EPS)+C0i[cc0012,args],
D0i[dd0000,args___]:>1/(24EPS)+D0i[dd0000,args],
D0i[dd00001,args___]:>-1/(96EPS)+D0i[dd00001,args],
D0i[dd00002,args___]:>-1/(96EPS)+D0i[dd00002,args],
D0i[dd00003,args___]:>-1/(96EPS)+D0i[dd00003,args]
};

DimensionTo4[expin_]:=Block[{expDReplace,expDEXP,expReEXP},
expDReplace=expin/.PVReplaceRule/.{D->4-2EPS};
expReEXP=If[FreeQ[expDReplace,_Re],expDReplace,ComplexExpand[expDReplace,{_A0,_A00,_B0i,_C0i,_D0i,_E0i,_B0,_DB0,_C0,_D0,_PaVe,_PaVeSC}]];
expDEXP=Series[EPS expReEXP,{EPS,0,1}];
If[expDEXP===0,{0,0},(*CoefficientList[expDEXP,EPS]*){Coefficient[expDEXP,EPS,0],Coefficient[expDEXP,EPS,1]}]
];


Finite=1;
PaVeSCtoPaVeRule={PaVeSC[p__,m1_,m2_,m3_]:>PaVe[0,0,{p},{m1,m2,m3}]};
PaVeSCtoLoopToolRule={PaVeSC[args___]:>C0i[cc00,args]};
AlignmentRule={alp->beta-\[Pi]/2};
ConstantSimplifyRule={
(*MZ2 -> MZ^2, MW2 -> MW^2, 
MU2 -> MU^2, MD2 -> MD^2, MC2 -> MC^2, MS2 -> MS^2, MB2 -> MB^2, MT2 -> MT^2, 
ME2 -> ME^2, ML2 -> ML^2, MM2 -> MM^2,Mh02->MH^2,Mh0->MH,MH2->MH^2,
Alfa2->Alfa^2,Alfa->EL^2/(4 Pi),EL2->EL^2*)
};
(*Only use the following ConstantValueRule when trying to get numerical result*)
ConstantValueRule={MZ -> 91.1876, MW -> 80.38, 
MU -> 7.356*10^(-2), MD -> MU, 
MC -> 1.27500,  MS -> 9.5*10^(-2), 
MB -> 4.6600, MT -> 173.21, 
ME -> 5.11*10^(-4), ML -> 1.77684, MM -> 0.10566,
MZ2 -> MZ^2, MW2 -> MW^2, 
MU2 -> MU^2, MD2 -> MD^2, MC2 -> MC^2, MS2 -> MS^2, MB2 -> MB^2, MT2 -> MT^2, 
ME2 -> ME^2, ML2 -> ML^2, MM2 -> MM^2,
CW -> MW/MZ, CW2 -> CW^2, SW2 -> 1 - CW2, SW -> Sqrt[SW2],w->ArcCos[CW], 
Alfa -> 1/137.035999074, EL -> Sqrt[4 \[Pi] Alfa], Alfa2 -> Alfa^2,
Mh->125.0,Mh0 -> 125.0, MH -> 125.0, Mh02 -> Mh0^2, MH2 -> MH^2,(* Alfas -> 0.1032728*)Alfas->0};


(*Get Coupling using Rule to simply*)
CouplingsWithRule[file_,rule_]:=Block[{expin,Rules,
expFreeTerm,expReCoeff,expReFuncs,expPVCoeff,expPVFuncs
},
expin=Get[file];
Rules=Join[rule,ConstantSimplifyRule];
expFreeTerm=expin[[1]]//.Rules;
expReCoeff=Simplify[#//.Rules]&/@expin[[2]][[2]];
expReFuncs=expin[[2]][[1]]/.{Mh02->MH2};
expPVCoeff=Simplify[#//.Rules]&/@expin[[3]][[2]];
expPVFuncs=expin[[3]][[1]]/.{Mh02->MH2}//.PaVeSCtoLoopToolRule;
If[Length[expReFuncs]==0,expReFuncs={0};expReCoeff={0};];
If[Length[expPVFuncs]==0,expPVFuncs={0};expPVCoeff={0};];
expFreeTerm+expReCoeff.expReFuncs+expPVCoeff.expPVFuncs
];
RemoveLightFermion[coupling_List,LightFermion_List:{MU2,MD2,MC2,MS2,ME2,MM2,ML2}]:=Block[{cons,ReFunc,ReCoeff,PVFunc,PVCoeff,
ReNoLightFunc,ReNoLightCoeff,
PVNoLightFunc,PVNoLightCoeff,
lenRe,lenPV,i,FreeCheck
},
cons=coupling[[1]];
ReFunc=coupling[[2]][[1]];
ReCoeff=coupling[[2]][[2]];
PVFunc=coupling[[3]][[1]];
PVCoeff=coupling[[3]][[2]];
ReNoLightFunc={};
ReNoLightCoeff={};
PVNoLightFunc={};
PVNoLightCoeff={};
lenRe=Length[ReFunc];
lenPV=Length[PVFunc];
FreeCheck[exp_]:=And@@(FreeQ[exp,#]&/@LightFermion);
For[i=1,i<=lenRe,i++,
If[FreeCheck[ReFunc[[i]]],ReNoLightFunc={ReNoLightFunc,ReFunc[[i]]};ReNoLightCoeff={ReNoLightCoeff,ReCoeff[[i]]};,ReNoLightFunc=ReNoLightFunc;ReNoLightCoeff=ReNoLightCoeff;];
];
For[i=1,i<=lenPV,i++,
If[FreeCheck[PVFunc[[i]]],PVNoLightFunc={PVNoLightFunc,PVFunc[[i]]};PVNoLightCoeff={PVNoLightCoeff,PVCoeff[[i]]};,PVNoLightFunc=PVNoLightFunc;PVNoLightCoeff=PVNoLightCoeff;];
];
{cons,{Flatten[ReNoLightFunc],Flatten[ReNoLightCoeff]},{Flatten[PVNoLightFunc],Flatten[PVNoLightCoeff]}}
];


(*Get Coupling in Dimension 4*)
(*Using following function if the expression in the file depend on D*)
CouplingsIND4WithRule[file_,rule_]:=Block[{expin,Rules,
expFreeTerm,expReCoeff,expReFuncs,expPVCoeff,expPVFuncs,
expReDimExpan,expPVDimExpan,
UVDivRe,UVDivPV,UVDivAll
},
expin=Get[file];
Rules=Join[rule,ConstantSimplifyRule];
expFreeTerm=expin[[1]]//.Rules;
expReCoeff=Simplify[#//.Rules]&/@expin[[2]][[2]];
expReFuncs=expin[[2]][[1]];
expPVCoeff=Simplify[#//.Rules]&/@expin[[3]][[2]];
expPVFuncs=expin[[3]][[1]]//.PaVeSCtoLoopToolRule;
If[Length[expReFuncs]==0,expReFuncs={0};expReCoeff={0};];
If[Length[expPVFuncs]==0,expPVFuncs={0};expPVCoeff={0};];
expReDimExpan=DimensionTo4/@(expReCoeff expReFuncs);
expPVDimExpan=DimensionTo4/@(expPVCoeff expPVFuncs);
UVDivRe=Plus@@expReDimExpan;
UVDivPV=Plus@@expPVDimExpan;
UVDivAll=Simplify[UVDivRe[[1]]+UVDivPV[[1]]];
If[UVDivAll===0,Print["UV GOOD!"],Print["UV Divergence: ",UVDivAll]];
expFreeTerm+UVDivRe[[2]]+UVDivPV[[2]]//.{Mh02->MH2}
];
CouplingsAlignIND4[file_]:=CouplingsIND4WithRule[file,AlignmentRule];


(*Expand the expin around alignment limit to order order_Integer*)
ExpandAroundAlign[expin_,order_Integer]:=Block[{expexpand,expout},
	expexpand=Series[expin//.{alp->beta-Pi/2+delta},{delta,0,order}];
	expout=Table[Coefficient[expexpand,delta,i]//Simplify,{i,0,order}];
	expout
];


(*Get 2HDM ONLY Expressions*)
SubtractExp[expin2HDM_,expinSM_]:=Block[{semi,sep,Recoeff,PVcoeff,Cons},
semi=expin2HDM-expinSM;
sep=ExtractPatternCoeff[semi,{_Re},{_A0,_A0i,_B0i,_C0i,_D0i}]; (*{a,{Re,Coeff},{PV,Coeff}}*)
Recoeff=Simplify/@sep[[2]][[2]];
PVcoeff=Simplify/@sep[[3]][[2]];
Cons=sep[[1]]//Simplify;
Cons+Recoeff.sep[[2]][[1]]+PVcoeff.sep[[3]][[1]]
];
SubtractExpByParts[expin2HDM_,expinSM_]:=Block[{semi,sep,Recoeff,PVcoeff,Cons},
semi=expin2HDM-expinSM;
sep=ExtractPatternCoeff[semi,{_Re},{_A0,_A0i,_B0i,_C0i,_D0i}]; (*{a,{Re,Coeff},{PV,Coeff}}*)
Recoeff=Simplify/@sep[[2]][[2]];
PVcoeff=Simplify/@sep[[3]][[2]];
Cons=sep[[1]]//Simplify;
{Cons,{sep[[2]][[1]],Recoeff},{sep[[3]][[1]],PVcoeff}}
];


(**)
Couplings2HDM[exp_,mH_,mA_,mHp_,betain_,M2in_,m12in_,m22in_,m32in_]:=Block[{result},
result=exp//.{MHH2->mH^2,MA02->mA^2,MHp2->mHp^2,beta->betain,M2->M2in,m12->m12in,m22->m22in,m32->m32in};
Re[result]
];
CouplingsSM[exp_,m12in_,m22in_,m32in_]:=Block[{result},
result=exp//.{m12->m12in,m22->m22in,m32->m32in};
Re[result]
];
