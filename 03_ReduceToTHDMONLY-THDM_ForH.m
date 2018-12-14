#!/usr/local/bin/MathematicaScript -script
(* ::Package:: *)

TypeName={"Type-I","Type-II","Type-LS","Type-FL"};
type=ToExpression[$ScriptCommandLine[[2]]]; 
Type=TypeName[[type]];
InputDir = "/home/ycwu/Workings/EWPM/CouplingsData/CouplingsAllDiagramsWithQCDCollection_ForH/";
OutDir = InputDir <> "../CouplingsAllDiagramsWithQCDCollectionTHDMONLYNoLightFermion_AllRe_ForH/";

SetDirectory[InputDir<>"../"];
<< CouplingsHandling`;

FermionFileName[ff_] := 
  Block[{THDMnamelist, SMnamelist, lenTHDM, lenSM, i, j, lenAll, 
    nameAllFrom, nameAllTo},
   	THDMnamelist = {"inTHDMNLOLT_Scalar.dat", 
     "inTHDMNLOLT_Pseudo.dat"};
   	lenTHDM = Length[THDMnamelist];
   	nameAllFrom = {THDMnamelist};
   	nameAllTo = {THDMnamelist};
   	lenAll = {lenTHDM};
   	For[i = 1, i <= 1, i++,
    		For[j = 1, j <= lenAll[[i]], j++,
      			nameAllFrom[[i]][[j]] = 
       InputDir <> Type <> "/FC" <> ff <> nameAllFrom[[i]][[j]];
      			nameAllTo[[i]][[j]] = 
       OutDir <> Type <>"/FC" <> ff <> nameAllTo[[i]][[j]];
      		];
    	];
   	{nameAllFrom, nameAllTo}
   ];
VectorFileName[vv_] := 
  Block[{THDMnamelist, SMnamelist, lenTHDM, lenSM, i, j, lenAll, 
    nameAllFrom, nameAllTo},
   	THDMnamelist = {"inTHDMNLOLT_Metric.dat"};
   	lenTHDM = Length[THDMnamelist];
   	nameAllFrom = {THDMnamelist};
   	nameAllTo = {THDMnamelist};
   	lenAll = {lenTHDM};
   	For[i = 1, i <= 1, i++,
    		For[j = 1, j <= lenAll[[i]], j++,
      			nameAllFrom[[i]][[j]] = 
       InputDir <> Type <> "/FC" <> vv <> nameAllFrom[[i]][[j]];
      			nameAllTo[[i]][[j]] = 
       OutDir <> Type <>"/FC" <> vv <> nameAllTo[[i]][[j]];
      		];
    	];
   	{nameAllFrom, nameAllTo}
   ];

ReduceProcedure[fileTHDM_, fileSM_] := Block[{cthdm, csm},
   cthdm = CouplingsAlignIND4[fileTHDM];
   csm = CouplingsAlignIND4[fileSM];
   SubtractExp[cthdm, csm]
   ];

FermionReduce[ff_] := 
  Block[{filename, filesfrom, filesto, filesfromTHDM, filesfromSM, 
    filestoTHDM, filestoSM,
    i, lenSM, lenTHDM,expsm,expsmByParts,expthdm, expthdmonly
    },
   filename = FermionFileName[ff];
   filesfrom = filename[[1]];
   filesto = filename[[2]];
   filesfromTHDM = filesfrom[[1]];
   filestoTHDM = filesto[[1]];
   lenTHDM = Length[filesfromTHDM];
   For[i = 1, i <= lenTHDM, i++,
		expthdm=CouplingsWithRule[filesfromTHDM[[i]],{}];
        expthdmonly = SubtractExpByParts[expthdm,0];
        expthdmonly = RemoveLightFermion[expthdmonly];
		expthdmonly = expthdmonly[[1]]+(ComplexExpand[#]&/@expthdmonly[[2]][[1]]).expthdmonly[[2]][[2]]+expthdmonly[[3]][[1]].expthdmonly[[3]][[2]];
		expthdmonly = ExtractPatternCoeffWithSimplify[expthdmonly,{_A0i,_B0i,_C0i,_D0i}];
        Put[expthdmonly, filestoTHDM[[i]]];
    ];
   ];
VectorReduce[vv_] := 
  Block[{filename, filesfrom, filesto, filesfromTHDM, filesfromSM, 
    filestoTHDM, filestoSM,
    i, lenSM, lenTHDM,expsm, expsmByParts,expthdm, expthdmonly
    },
   filename = VectorFileName[vv];
   filesfrom = filename[[1]];
   filesto = filename[[2]];
   filesfromTHDM = filesfrom[[1]];
   filestoTHDM = filesto[[1]];
   lenTHDM = Length[filesfromTHDM];
   For[i = 1, i <= lenTHDM, i++,
		expthdm=CouplingsWithRule[filesfromTHDM[[i]],{}];
        expthdmonly = SubtractExpByParts[expthdm,0];
		expthdmonly = RemoveLightFermion[expthdmonly];
		expthdmonly = expthdmonly[[1]]+(ComplexExpand[#]&/@expthdmonly[[2]][[1]]).expthdmonly[[2]][[2]]+expthdmonly[[3]][[1]].expthdmonly[[3]][[2]];
		expthdmonly = ExtractPatternCoeffWithSimplify[expthdmonly,{_A0i,_B0i,_C0i,_D0i}];
        Put[expthdmonly, filestoTHDM[[i]]];
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
