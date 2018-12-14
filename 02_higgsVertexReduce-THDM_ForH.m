#!/usr/local/bin/MathematicaScript -script

(* ::Package:: *)

(*Needs["FeynCalc`"];*)
TypeName={"Type-I","Type-II","Type-LS","Type-FL"};
type=ToExpression[$ScriptCommandLine[[2]]]; 
Type=TypeName[[type]];
prefixFrom="/home/ycwu/Workings/EWPM/CouplingsData/CouplingsAllDiagramsWithQCD_ForH/";
prefixTo="/home/ycwu/Workings/EWPM/CouplingsData/CouplingsAllDiagramsWithQCDCollection_ForH/";
FermionFileName[ff_]:=Block[{NLOnamelist,LOnamelist,lenNLO,lenLO,i,j,lenAll,nameAllFrom,nameAllTo},
	NLOnamelist={"inTHDMNLOLT_Scalar.dat","inTHDMNLOLT_Pseudo.dat"};
	LOnamelist={"inTHDMTree_Scalar.dat","inTHDMTree_Pseudo.dat"};
	lenNLO=Length[NLOnamelist];
	lenLO=Length[LOnamelist];
	nameAllFrom={NLOnamelist,LOnamelist};
	nameAllTo={NLOnamelist,LOnamelist};
	lenAll={lenNLO,lenLO};
	For[i=1,i<=2,i++,
		For[j=1,j<=lenAll[[i]],j++,
			nameAllFrom[[i]][[j]]=prefixFrom<>Type<>"/"<>ff<>nameAllFrom[[i]][[j]];
			nameAllTo[[i]][[j]]=prefixTo<>Type<>"/FC"<>ff<>nameAllTo[[i]][[j]];
		];
	];
	{nameAllFrom,nameAllTo}
];
VectorFileName[vv_]:=Block[{NLOnamelist,LOnamelist,lenNLO,lenLO,i,j,lenAll,nameAllFrom,nameAllTo},
	NLOnamelist={"inTHDMNLOLT_Metric.dat"};
	LOnamelist={"inTHDMTree_Metric.dat"};
	lenNLO=Length[NLOnamelist];
	lenLO=Length[LOnamelist];
	nameAllFrom={NLOnamelist,LOnamelist};
	nameAllTo={NLOnamelist,LOnamelist};
	lenAll={lenNLO,lenLO};
	For[i=1,i<=2,i++,
		For[j=1,j<=lenAll[[i]],j++,
			nameAllFrom[[i]][[j]]=prefixFrom<>Type<>"/"<>vv<>nameAllFrom[[i]][[j]];
			nameAllTo[[i]][[j]]=prefixTo<>Type<>"/FC"<>vv<>nameAllTo[[i]][[j]];
		];
	];
	{nameAllFrom,nameAllTo}
];
SMRules={MH2->Mh02,MZ->MW/CW};
(*DB00Rule={DB00[k_,m1_,m2_]:>1/36(-3 PaVe[0, {k}, {m1, m1}]-3(k-4m1)DB0[k,m1,m1]-2)/;m1===m2,DB00[k_,m1_,m2_]:>(-3*(m1-m2)^2*PaVe[0, {0}, {m1, m2}]-3*(k^2-(m1-m2)^2)*PaVe[0, {k}, {m1, m2}]+k*(-2*k-3*(k^2+(m1-m2)^2-2*k*(m1+m2))*DB0[k,m1,m2]))/(36*k^2)/;m1=!=m2};
(*SetOptions[PaVe,PaVeAutoReduce\[Rule]True];*)
SetOptions[PaVeReduce,A0ToB0->True];
SetOptions[A0,A0ToB0->True];
ToOldBRules={(*B0i[bb0,args__]->B0[args],*)B0i[bb1,args__]->B1[args],B0i[bb00,args__]->B00[args],B0i[bb11,args__]->B11[args],B0i[bb001,args__]->B001[args],B0i[bb111,args__]->B111[args],B0i[dbb0,args__]->DB0[args],B0i[dbb1,args__]->DB1[args],B0i[dbb00,args__]->DB00[args],B0i[dbb11,args__]->DB11[args]};
ToFeynCalc[expr_,tempfile_]:= Block[{A0i,B0i,C0i,D0i,E0i,F0i,PaVeSC,PaVeS},
C0i[cc0,args___]:=C0[args];
C0i[cc00,args___]:=PaVeSC[args];
D0i[dd0,args___]:=D0[args];
E0i[ee0,args___]:=E0[args];
F0i[ff0,args___]:=F0[args];
B0i[bb0,p__,m1_,m2_]:=PaVeS[0,{p},{m1,m2}];
C0i[i_,p__,m1_,m2_,m3_]:=PaVeS[Sequence@@(ToExpression/@Drop[Characters[ToString[i]],2]),{p},{m1,m2,m3}];
D0i[i_,p__,m1_,m2_,m3_,m4_]:=PaVeS[Sequence@@(ToExpression/@Drop[Characters[ToString[i]],2]),{p},{m1,m2,m3,m4}];
E0i[i_,p__,m1_,m2_,m3_,m4_,m5_]:=PaVeS[Sequence@@(ToExpression/@Drop[Characters[ToString[i]],2]),{p},{m1,m2,m3,m4,m5}];
F0i[i_,p__,m1_,m2_,m3_,m4_,m5_,m6_]:=PaVeS[Sequence@@(ToExpression/@Drop[Characters[ToString[i]],2]),{p},{m1,m2,m3,m4,m5,m6}];
Put[expr/.ToOldBRules,tempfile];
];*)
(*******************
The ExtractPatternCoeff Function:
The pattern_List is a list that contains the pattern for which you want to be isolate out.
Then the expin will be separated into three parts:
expin = part1 + part2.part3;
where part 1 is the part that is free of any pattern in the pattern_List, 
part2 is a list, and it contains all the part in the expin that match the pattern in the pattern_List,
part3 is also a list, same length with part2, and is the 1-order coefficient of the corresponding term in part2
********************)
ExtractPatternCoeff[expin_,pattern_List]:=Block[{elems,findings,freerules,frees,coeff},
elems=Union[Level[expin,Infinity]];
findings=Union@@(Cases[elems,#]&/@pattern);
freerules=(#->0)&/@findings;
frees=expin/.freerules;
coeff=Table[Coefficient[expin,findings[[i]]],{i,1,Length[findings]}];
{frees,findings,coeff}
];
CollectionPattern={_A0,_B0i,_C0i,_D0i};
ReduceProcedure[expin_]:=Block[{tempname,step1,step2,step3,step4,step4ReFree,step4ReFunc,step4ReCoeff,step4ReFuncExpand,step4ReAll,step4ReFreeDecom,step4ReAllDecom,step5,expout,LoopFuncs,elems,LoopCoeff,CoeffConst,FinalConst,FinalReFunc,FinalReCoeff,FinalPVFunc,FinalPVCoeff},
		(*Print["ToFeynCalc Convention......"];
		tempname=prefixTo<>Type<>"/temp.dat";
		ToFeynCalc[expin,tempname];
        Print["Reading ToFeynCalc Temp File......"];
		PaVeS[args___]:=PaVe[args];
		step1=Get[tempname];
		Print["Reduce the Tensor coefficient to scalar one......"];
		step2=PaVeReduce[step1];
		Print["Reduce DB00......"];
		step3=step2//.DB00Rule;*)
		step3=expin;
		Print["Seperate Expr into two parts: 1. Containing Re[], 2. Free of Re[]......"];
		step4=ExtractPatternCoeff[step3,{_Re}];
		step4ReFree=step4[[1]]; (*The Part that is free of Re[]*)
		step4ReFunc=step4[[2]]; (*The Re[] functions in the original expr*)
		Print[step4ReFunc];
		step4ReCoeff=step4[[3]]; (*The Coefficient for above Re[] functions in the original expr*)
        Print["Expanding the Re[] now......"];
		step4ReFuncExpand=ComplexExpand[#,CollectionPattern]& /@ step4ReFunc; (*Only expand the part containing Re[] functions*)
		(*Print["Remove possible Im[] function......"];
		step4ReFuncExpandRe=step4ReFuncExpand//.{Im[F_[args___]]:>-I(F[args]-Re[F[args]]),Finite->1};*)
        Print["Collect the Re[] Part, now inside Re[], it just contain the scalar integral......"];
        step4ReAll=step4ReCoeff.step4ReFuncExpand;
		Print["Seperate ReAll agagin into two Parts: 1. Free of Re[], 2. Containing Re[]......"];
		step4ReAllDecom=ExtractPatternCoeff[step4ReAll,{_Re}];
		Print["Seperate ReFree part into two Parts: 1. Containing scalar integrals, 2. Free of these integrals......"];
		step4ReFreeDecom=ExtractPatternCoeff[step4ReFree,CollectionPattern];
		FinalConst=step4ReAllDecom[[1]]+step4ReFreeDecom[[1]]//.SMRules; (*This is the constant part that doesn't contain any Re[] or Scalar integrals, it could be really complicated, so better not to simplify it without any assumptions*)
		Print["Simplify the coefficient for each Re[] part, and Scalar integral part......"];
		FinalReFunc=step4ReAllDecom[[2]];
		Print[FinalReFunc];
		FinalReCoeff=Simplify[#//.SMRules]&/@step4ReAllDecom[[3]];
		FinalPVFunc=step4ReFreeDecom[[2]];
		Print[FinalPVFunc];
		FinalPVCoeff=Simplify[#//.SMRules]&/@step4ReFreeDecom[[3]];
		Print["Output"];
		Clear[PaVeS];
		{FinalConst,{FinalReFunc,FinalReCoeff},{FinalPVFunc,FinalPVCoeff}}
];
FermionReduce[ff_]:=Block[{tempname,filename,filesfrom,filesto,filesfromNLO,filesfromTree,filestoNLO,filestoTree,lenTree,lenNLO,i,expin,step1,step2,step3,step4,step5,expout},
	filename = FermionFileName[ff];
	filesfrom=filename[[1]];
	filesto=filename[[2]];
	filesfromNLO=filesfrom[[1]];
	filesfromTree=filesfrom[[2]];
	filestoNLO=filesto[[1]];
	filestoTree=filesto[[2]];
	lenTree=Length[filesfromTree];
	For[i=1,i<=lenTree,i++,
		CopyFile[filesfromTree[[i]],filestoTree[[i]]];
	];
	lenNLO=Length[filesfromNLO];
	For[i=1,i<=lenNLO,i++,
		Print["Process the file: ",filesfromNLO[[i]],"......"];
		expin=Get[filesfromNLO[[i]]];
		expout=ReduceProcedure[expin];
		Print["Output"];
		Put[expout,filestoNLO[[i]]];
	];
];
VectorReduce[vv_]:=Block[{filename,filesfrom,filesto,filesfromNLO,filesfromTree,filestoNLO,filestoTree,lenTree,lenNLO,i,expin,step1,step2,step3,step4,step5,expout},
	filename = VectorFileName[vv];
	filesfrom=filename[[1]];
	filesto=filename[[2]];
	filesfromNLO=filesfrom[[1]];
	filesfromTree=filesfrom[[2]];
	filestoNLO=filesto[[1]];
	filestoTree=filesto[[2]];
	lenTree=Length[filesfromTree];
	For[i=1,i<=lenTree,i++,
		CopyFile[filesfromTree[[i]],filestoTree[[i]]];
	];
	lenNLO=Length[filesfromNLO];
	For[i=1,i<=lenNLO,i++,
		Print["Process the file: ",filesfromNLO[[i]],"......"];
		expin=Get[filesfromNLO[[i]]];
		expout=ReduceProcedure[expin];
		Print["Output"];
		Put[expout,filestoNLO[[i]]];
	];
];
FermionReduce["CHu"];
FermionReduce["CHd"];
FermionReduce["CHc"];
FermionReduce["CHs"];
FermionReduce["CHt"];
FermionReduce["CHb"];
FermionReduce["CHele"];
FermionReduce["CHmuon"];
FermionReduce["CHtau"];
VectorReduce["CHGA"];
VectorReduce["CHZ"];
VectorReduce["CHG"];
VectorReduce["CHW"];
VectorReduce["CHGAZ"];
