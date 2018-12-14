#!/usr/local/bin/MathematicaScript -script
(* ::Package:: *)

Needs["FeynArts`"];
Needs["FormCalc`"];
CKM = IndexDelta;
type=ToExpression[$ScriptCommandLine[[2]]];
TypeName={"Type-I","Type-II","Type-LS","Type-FL"};
YukawaFile={"THDM_Yukawa_TypeI_NLO","THDM_Yukawa_TypeII_NLO","THDM_Yukawa_TypeLS_NLO","THDM_Yukawa_TypeFL_NLO"};
SetOptions[InsertFields,GenericModel -> Lorentznew, InsertionLevel -> {Particles}];
SetOptions[CreateFeynAmp, Truncated -> False];
SetOptions[CalcFeynAmp, FermionChains -> VA];
ifonshell = False;
protop = 1 -> 2;
prefix="/home/ycwu/Workings/EWPM/CouplingsData/CouplingsAllDiagramsWithQCD_ForH/"<>TypeName[[type]]<>"/";
THDMaddon = {"THDM_SSV_SVV_NLO", "THDM_SSVV_NLO", "THDM_SSS_SB_NLO", 
   "THDM_SSSS_SB_NLO", YukawaFile[[type]], "THDM_RC"};
THDMNLO = Join[{"THDM_EW_BasicQCD"}, THDMaddon];
AngleRules = {CA -> Cos[alp], SA -> Sin[alp], CB -> Cos[beta], 
   SB -> Sin[beta], TB -> Tan[beta], C2B -> Cos[2 beta], 
   S2B -> Sin[2 beta], C2A -> Cos[2 alp], S2A -> Sin[2 alp], 
   CAB -> Cos[alp + beta], SAB -> Sin[alp + beta], 
   CBA -> Cos[beta - alp], SBA -> Sin[beta - alp], 
   SBA2 -> Sin[beta - alp]^2, CBA2 -> Cos[beta - alp]^2, 
   SB2 -> Sin[beta]^2, TB2 -> Tan[beta]^2, CA2 -> Cos[alp]^2, 
   SA2 -> Sin[alp]^2, CB2 -> Cos[beta]^2};
VertexCouplingTree[tops_, process_, model_] := 
  Block[{toptree, instree, amptree, topCTs, insCTs, ampCTs, toploop, 
    insloop, amploop, ampTreeAndLoop, tmp, tmp1, ren, uv, step1, 
    step2, uvsimple},
   SetOptions[InsertFields, Model -> model];
   ClearProcess[];
   toptree = CreateTopologies[0, tops];
   instree = InsertFields[toptree, process];
   amptree = CreateFeynAmp[instree];
   tmp = Plus @@ CalcFeynAmp[amptree, OnShell -> ifonshell];
   tmp1 = tmp //. Abbr[] //. Subexpr[] //. Abbr[];
   tmp1
   ];
VertexCoupling[tops_, process_, model_] := 
  Block[{toptree, instree, amptree, topCTs, insCTs, ampCTs, toploop, 
    insloop, amploop, ampTreeAndLoop, tmp, tmp1, ren, uv, step1, 
    step2, uvsimple},
   SetOptions[InsertFields, Model -> model];
   ClearProcess[];
   toptree = CreateTopologies[0, tops];
   instree = InsertFields[toptree, process];
   amptree = CreateFeynAmp[instree];
   topCTs = 
    CreateCTTopologies[1, tops, 
     ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs}];
   insCTs = InsertFields[topCTs, process];
   ampCTs = CreateFeynAmp[insCTs];
   toploop = 
    CreateTopologies[1, tops, 
     ExcludeTopologies -> {Tadpoles, WFCorrections}];
   insloop = InsertFields[toploop, process];
   Paint[insloop];
   amploop = CreateFeynAmp[insloop];
   tmp = CalcFeynAmp[amptree, amploop, ampCTs, OnShell -> ifonshell];
   tmp1 = tmp //. Abbr[] //. Subexpr[] //. Abbr[];
   (*Print[tmp1];*)
   ren = CalcRenConst[tmp1];
   ampTreeAndLoop = tmp1 //. List@@ren[[1]];
   uv = UVDivergentPart[ampTreeAndLoop];
   step1 = uv[[1]] // Expand // Simplify // SMSimplify;
   step2 = 
    step1 //. {CA -> Cos[\[Alpha]], SA -> Sin[\[Alpha]], 
      CB -> Cos[\[Beta]], SB -> Sin[\[Beta]], TB -> Tan[\[Beta]], 
      C2B -> Cos[2 \[Beta]], S2B -> Sin[2 \[Beta]], 
      C2A -> Cos[2 \[Alpha]], S2A -> Sin[2 \[Alpha]], 
      CAB -> Cos[\[Alpha] + \[Beta]], SAB -> Sin[\[Alpha] + \[Beta]], 
      CBA -> Cos[\[Beta] - \[Alpha]], SBA -> Sin[\[Beta] - \[Alpha]], 
      SBA2 -> Sin[\[Beta] - \[Alpha]]^2, 
      CBA2 -> Cos[\[Beta] - \[Alpha]]^2, SB2 -> Sin[\[Beta]]^2, 
      TB2 -> Tan[\[Beta]]^2, CA2 -> Cos[\[Alpha]]^2, 
      SA2 -> Sin[\[Alpha]]^2, CB2 -> Cos[\[Beta]]^2};
   uvsimple = step2 // Simplify // SMSimplify;
   If[uvsimple =!= 0, 
    Print["UV Divergence! Checking your files!"]; {uvsimple, 
     ampTreeAndLoop}, {uvsimple, ampTreeAndLoop}]
   (*{uvsimple,ampTreeAndLoop}*)
   ];
hVVMetricPart[vertex_, \[Mu]1_, \[Mu]2_, \[Mu]3_] := Block[{result},
   result = 
    vertex //. {Pair[_e | _ec, _k] -> 0, Pair[_e, _ec] -> 1, 
      Pair[k[1], k[1]] -> \[Mu]1, Pair[k[2], k[2]] -> \[Mu]2, 
      Pair[k[3], k[3]] -> \[Mu]3, Pair[_ec, _ec] -> 1, 
      Eps[_, _, _, _] -> 0, SUNT[_, _, _, _] -> 1, Mat[1] -> 1};
   result
   ];
hffScalarPart[vertex_, \[Mu]1_, \[Mu]2_, \[Mu]3_] := Block[{result},
   result = 
    vertex //. {DiracChain[_, _, _, _, _] -> 0, 
      DiracChain[_, _, _, _] -> 0, DiracChain[__, 5, __] -> 0, 
      Mat[0] -> 0, DiracChain[__, 1, __] -> 1, Mat[SUNT[__, __]] -> 1,
       Mat[1] -> 1, Pair[k[1], k[1]] -> \[Mu]1, 
      Pair[k[2], k[2]] -> \[Mu]2, Pair[k[3], k[3]] -> \[Mu]3};
   result
   ];
hffPseudoPart[vertex_, \[Mu]1_, \[Mu]2_, \[Mu]3_] := Block[{result},
   result = 
    vertex //. {DiracChain[_, _, _, _, _] -> 0, 
      DiracChain[_, _, _, _] -> 0, DiracChain[__, 1, __] -> 0, 
      Mat[0] -> 0, DiracChain[__, 5, __] -> 1, Mat[SUNT[__, __]] -> 1,
       Mat[1] -> 1, Pair[k[1], k[1]] -> \[Mu]1, 
      Pair[k[2], k[2]] -> \[Mu]2, Pair[k[3], k[3]] -> \[Mu]3};
   result
   ];
hffGeneralVertexAll[name_, f_, model_,modelSM_] := 
  Block[{process, VinTHDM, VinSM, VinTHDMTree ,VinSMTree, CinTHDMS, CinSMS, 
    CinTHDMTreeS,CinSMTreeS, CinTHDMP, CinSMP, CinTHDMTreeP, CinSMTreeP},
   process = S[2] -> {f, -f};
   VinTHDM = VertexCoupling[protop, process, model];
   VinTHDMTree = VertexCouplingTree[protop,process,model];
   CinTHDMS = 
    hffScalarPart[Plus @@ VinTHDM[[2]], m12, m22, m32] //. 
     AngleRules;
   CinTHDMTreeS = hffScalarPart[VinTHDMTree,0,0,0] //. AngleRules;
   CinTHDMP = 
    hffPseudoPart[Plus @@ VinTHDM[[2]], m12, m22, m32] //. 
     AngleRules;
   CinTHDMTreeP = hffPseudoPart[VinTHDMTree,0,0,0] //. AngleRules;
   Put[CinTHDMS, prefix<>name <> "inTHDMNLOLT_Scalar.dat"];
   Put[CinTHDMTreeS, prefix <> name <> "inTHDMTree_Scalar.dat"];
   Put[CinTHDMP, prefix<>name <> "inTHDMNLOLT_Pseudo.dat"];
   Put[CinTHDMTreeP, prefix <> name <> "inTHDMTree_Pseudo.dat"];
   ];
hVVGeneralVertexMetric[name_, V_,model_,modelSM_] := 
  Block[{process, VinTHDM, VinSM, VinTHDMTree, VinSMTree, CinTHDM, CinSM, 
    CinTHDMTree, CinSMTree},
   process = S[2] -> {V, V};
   VinTHDM = VertexCoupling[protop, process, model];
   VinTHDMTree = VertexCouplingTree[protop,process,model];
   CinTHDM = 
    hVVMetricPart[Plus @@ VinTHDM[[2]], m12, m22, m32] //. 
     AngleRules;
   CinTHDMTree = hVVMetricPart[VinTHDMTree,0,0,0] //. AngleRules;
   Put[CinTHDM, prefix<>name <> "inTHDMNLOLT_Metric.dat"];
   Put[CinTHDMTree, prefix<>name<>"inTHDMTree_Metric.dat"];
   ];
hVVGeneralVertexMetric[name_, V_, V1_,model_,modelSM_] := 
  Block[{process, VinTHDM, VinSM, VinTHDMTree, VinSMTree, CinTHDM, CinSM, 
    CinTHDMTree, CinSMTree},
   process = S[2] -> {V, V1};
   VinTHDM = VertexCoupling[protop, process, model];
   VinTHDMTree = VertexCouplingTree[protop,process,model];
   CinTHDM = 
    hVVMetricPart[Plus @@ VinTHDM[[2]], m12, m22, m32] //. 
     AngleRules;
   CinTHDMTree = hVVMetricPart[VinTHDMTree,0,0,0] //. AngleRules;
   Put[CinTHDM, prefix<>name <> "inTHDMNLOLT_Metric.dat"];
   Put[CinTHDMTree, prefix<>name<>"inTHDMTree_Metric.dat"];
   ];
hffGeneralVertexAll["CHu", F[3,{1}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHd", F[4,{1}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHc", F[3,{2}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHs", F[4,{2}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHt", F[3,{3}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHb", F[4,{3}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHele", F[2,{1}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHmuon", F[2,{2}], THDMNLO,SMQCD];
hffGeneralVertexAll["CHtau", F[2,{3}], THDMNLO,SMQCD];
hVVGeneralVertexMetric["CHGA", V[1],THDMNLO,SMQCD];
hVVGeneralVertexMetric["CHZ", V[2],THDMNLO,SMQCD];
hVVGeneralVertexMetric["CHG", V[5],THDMNLO,SMQCD];
hVVGeneralVertexMetric["CHW", V[3],-V[3],THDMNLO,SMQCD];
hVVGeneralVertexMetric["CHGAZ", V[1],V[2],THDMNLO,SMQCD];
Print["Good!"];
