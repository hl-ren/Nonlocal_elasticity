(* ::Package:: *)

Needs["NDSolve`FEM`"]
AbortOnMessage[t_:True] := If[t, messageHandler = If[Last[#1], Abort[]] & ; 
     Internal`AddHandler["Message", messageHandler], 
    Internal`RemoveHandler["Message", messageHandler]]
CheckTextCommand[path_String] := Block[{res}, If[FileExistsQ[path], 
     res = ReadList[path, String, 1][[1]]; If[StringStartsQ[res, "1"], 
       Print["Abort command was found in file: ", path]; Abort[]]]]
CreateFolder[mydir_] := Block[{sn = NotebookDirectory[], dirname}, 
    dirname = FileNameJoin[{sn, mydir}]; Switch[FileType[dirname], None, 
      CreateDirectory[dirname], Directory, Null, File, 
      Print["File with same name already exists!!"]]; 
     StringJoin[dirname, StringTake[sn, {-1, -1}]]]
FindPoints[v1_List, xmin_List, xmax_List, show_:False] := 
   Module[{ilist = {}, y, r, ndim = Length[v1[[1]]], len = Length[v1]}, 
    Do[If[LessThan2[v1[[i]], xmax] && LessThan2[xmin, v1[[i]]], 
       AppendTo[ilist, i]], {i, len}]; 
     If[show, If[Length[v1[[1]]] == 2, Print[ListPlot[v1[[ilist]], 
          DataRange -> Automatic]]]; If[Length[v1[[1]]] == 3, 
        Print[ListPointPlot3D[v1[[ilist]]], DataRange -> Full]]]; ilist]
LessThan2[v1_List, v2_List] := Module[{p = True}, 
    Do[If[v1[[i]] > v2[[i]], p = False; Break[]], {i, Length[v1]}]; Return[p]]
FindPointsC[v1_List, xmin_List, xmax_List] := 
   Flatten[Switch[Length[v1[[1]]], 0, FindPointsC1[v1, xmin, xmax], 2, 
     FindPointsC2[v1, xmin, xmax], 3, FindPointsC3[v1, xmin, xmax]]]
Attributes[ForceListAll] = {HoldFirst}
 
ForceListAll[neivol_, neiList_, neicoord_, neicoord0_, stress_, fmat_, 
    invKList_, trKList_, hgPen_, neiWeiL_, hgtol_] := 
   Module[{res, forceL, ndim = Length[invKList[[1]]], neivol2, 
     num = Length[trKList], hgPbond, volSign}, 
    forceL = ConstantArray[0., {num, ndim}]; 
     res = ForceAccumulateHgC[neicoord, neicoord0, stress, fmat, invKList, 
       trKList, neivol, hgPen, neiWeiL]; hgPbond = res[[All,All,-1]]; 
     {volSign, neivol2} = CheckDamageByBondHourglass[neivol, hgPbond, hgtol]; 
     res = res[[All,All,1 ;; ndim]]; Do[forceL[[neiList[[i]]]] += res[[i]]; , 
      {i, num}]; {hgPbond, volSign, forceL, neivol2}]
CheckDamageByBondHourglass[neivol_, hgPbond_, tol_] := 
   Block[{volSign = Table[0, {Length[neivol]}], neivol2 = neivol}, 
    Do[Do[If[hgPbond[[i,j]] > tol, neivol2[[i,j]] = 0.; volSign[[i]] = 1; ], 
       {j, 2, Length[neivol[[i]]]}], {i, Length[neivol]}]; {volSign, neivol2}]
GridNdim[xmin_List, xmax_List, dx_] := 
   Module[{idx = {"i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", 
       "t"}, id2 = {}, st = "Table[", al, nx = Ceiling[(xmax - xmin)/dx], 
     ndx, dxL, ndim, res}, ndim = Length[xmin]; idx = idx[[1 ;; ndim]]; 
     ndx = (xmax - xmin)/nx; st = StringJoin[st, ToString[idx]]; 
     Do[id2 = StringJoin[id2, ",", ToString[{idx[[i]], 0, nx[[i]]}]], 
      {i, ndim}]; al = Flatten[ToExpression[StringJoin[st, id2, "]"]]]; 
     res = ArrayReshape[al, {Length[al]/ndim, ndim}]; 
     Do[res[[All,i]] *= ndx[[i]], {i, ndim}]; Do[res[[i]] += xmin, 
      {i, Length[res]}]; res]
KmatrixL[neicoord_, neiWeiL_, neivol_] := 
   Block[{res1, trk, num = Length[neicoord], ndim = Length[neicoord[[1,1]]], 
     invk, x1}, res1 = KmatrixC[neicoord, neiWeiL, neivol]; 
     trk = res1[[All,1]]; invk = ConstantArray[0., {num, ndim, ndim}]; 
     Do[x1 = res1[[i,2 ;; All]]; invk[[i]] = ArrayReshape[x1, 
         {ndim, ndim}]; , {i, num}]; {trk, invk}]
Attributes[KmatrixLPartial] = {HoldAll}
 
KmatrixLPartial[trK_, invK_, neicoord_, neiWeiL_, neivol_, volSign0_] := 
   Block[{res1, num, ndim = Length[neicoord[[1,1]]], x1, vp}, 
    vp = Flatten[Position[volSign0, 1]]; num = Length[vp]; 
     If[num > 0, res1 = KmatrixC[neicoord[[vp]], neiWeiL[[vp]], 
         neivol[[vp]]]; trK[[vp]] = res1[[All,1]]; 
       Do[x1 = res1[[i,2 ;; All]]; invK[[vp[[i]]]] = ArrayReshape[x1, 
           {ndim, ndim}]; , {i, num}]]; ]
MyDelaunayMesh[coord_] := Block[{mm, mm2, ndim}, ndim = Length[coord[[1]]]; 
     If[ndim < 2 || ndim > 3, 
      Print["Error, delaunayMesh only applies for 2D and 3D"]; Abort[]; ]; 
     mm = DelaunayMesh[coord]; mm2 = If[ndim == 2, Cases[MeshCells[mm, 2], 
        Polygon[x_] :> x, -1], Cases[MeshCells[mm, 3], Tetrahedron[x_] :> x, 
        -1]]; Abaqus2Mesh[coord, mm2]]
 
MyDelaunayMesh[coord_, uvw_] := Block[{mm, mm2, ndim}, 
    ndim = Length[coord[[1]]]; mm = DelaunayMesh[coord]; 
     mm2 = If[ndim == 2, Cases[MeshCells[mm, 2], Polygon[x_] :> x, -1], 
       Cases[MeshCells[mm, 3], Tetrahedron[x_] :> x, -1]]; 
     Abaqus2Mesh[coord + uvw, mm2]]
Abaqus2Mesh[coord_, mesh_] := Block[{ec, el = {}}, 
    ec = ElementClassify[mesh, Length[coord[[1]]]]; 
     If[Length[ec[[1]]] > 0, AppendTo[el, TriangleElement[mesh[[ec[[1]]]]]]]; 
     If[Length[ec[[2]]] > 0, AppendTo[el, QuadElement[mesh[[ec[[2]]]]]]]; 
     If[Length[ec[[3]]] > 0, AppendTo[el, TetrahedronElement[
        mesh[[ec[[3]]]]]]]; If[Length[ec[[4]]] > 0, 
      AppendTo[el, HexahedronElement[mesh[[ec[[4]]]]]]]; 
     ToElementMesh["Coordinates" -> coord, "MeshElements" -> el]]
ElementClassify[mesh_, ndim_] := Block[{pos = ConstantArray[0, Length[mesh]], 
     toqh = {{}, {}, {}, {}}}, 
    Do[pos[[i]] = ElementClassifyEc[Length[mesh[[i]]], ndim]; , 
      {i, Length[mesh]}]; toqh[[1]] = Flatten[Position[pos, 1]]; 
     toqh[[2]] = Flatten[Position[pos, 2]]; 
     toqh[[3]] = Flatten[Position[pos, 3]]; 
     toqh[[4]] = Flatten[Position[pos, 4]]; toqh]
ElementClassifyEc[len_, n_] := If[n == 2, If[len == 3 || len == 6, 1, 2], 
    If[len == 4 || len == 10, 3, 4]]
NeighborPartitionC[coord_List, Nei_List] := 
   Switch[Length[Dimensions[coord[[1]]]], 0, NeighborPartitionC1[coord, Nei], 
    1, NeighborPartitionC2[coord, Nei], 2, NeighborPartitionC3[coord, Nei]]
NeighborUpdateByPlane[nei_List, coord_List, geoM_List, SegmentListCut2P_, 
    xmin_, xmax_] := Block[{Nei, pin, broNei = ConstantArray[{}, 
       Length[coord]], ni, xi, xiN, broI}, 
    Nei = nei; pin = FindPoints[coord, xmin, xmax]; 
     Do[ni = nei[[i]]; xi = coord[[i]]; xiN = coord[[ni]]; broI = {}; 
       Do[If[SegmentListCut2P[geoM, xi, xiN[[j]]], AppendTo[broI, ni[[j]]]], 
        {j, Length[ni]}]; broNei[[i]] = broI; , {i, pin}]; 
     Do[Nei[[i]] = ComplementOrder[nei[[i]], broNei[[i]]], {i, pin}]; 
     {Nei, broNei}]
 
NeighborUpdateByPlane[nei_List, coord_List, geoM_List, SegmentListCut2P_, 
    checkList_] := Block[{Nei, pin, broNei = ConstantArray[{}, 
       Length[coord]], ni, xi, xiN, broI}, 
    Nei = nei; Do[ni = nei[[i]]; xi = coord[[i]]; xiN = coord[[ni]]; 
       broI = {}; Do[If[SegmentListCut2P[geoM, xi, xiN[[j]]], 
         AppendTo[broI, ni[[j]]]], {j, Length[ni]}]; broNei[[i]] = broI; , 
      {i, checkList}]; Do[Nei[[i]] = ComplementOrder[nei[[i]], broNei[[i]]], 
      {i, checkList}]; {Nei, broNei}]
SegmentListCut2P[segL_List, p1_, p2_] := Block[{fd = False}, 
    Do[If[SegmentCut2P[segL[[i]], p1, p2], fd = True; Break[]; ], 
      {i, Length[segL]}]; fd]
SegmentCut2P[{x0_, y0_, x1_, y1_}, {x2_, y2_}, {x3_, y3_}] := 
   Block[{t1, t2, t3, t4}, t1 = (x2 - x0)*(y1 - y0) - (y2 - y0)*(x1 - x0); 
     t2 = (x3 - x0)*(y1 - y0) - (y3 - y0)*(x1 - x0); 
     t3 = (x0 - x2)*(y3 - y2) - (y0 - y2)*(x3 - x2); 
     t4 = (x1 - x2)*(y3 - y2) - (y1 - y2)*(x3 - x2); 
     If[t1*t2 < 0 && t3*t4 < 0, Return[True], Return[False]]]
ComplementOrder[vs_List, vi_] := Block[{v2 = ConstantArray[True, Length[vs]], 
     vi2}, If[MemberQ[{List}, Head[vi]], vi2 = vi, vi2 = {vi}]; 
     Do[Do[If[vs[[j]] == vi2[[i]], v2[[j]] = False; Break[]], 
       {j, Length[vs]}], {i, Length[vi2]}]; vs[[Flatten[Position[v2, True]]]]]
OvitoOutput[path_String, fname_String, istep_, tag_, x__] := 
   Block[{file, str, xd = {x}, len, num, al, als, strl0, strl}, 
    file = FileNameJoin[{path, StringJoin[fname, ToString[istep], ".txt"]}]; 
     len = Length[Flatten[xd[[1 ;; All,1]]]]; num = Length[xd[[1]]]; 
     strl = {"ITEM: TIMESTEP", ToString[istep], "ITEM: NUMBER OF ATOMS", 
       ToString[num], StringJoin["ITEM: ATOMS id type ", tag]}; 
     strl0 = Table[StringJoin[ToString[i], " 1 ", DoubleListToString[
         Flatten[xd[[1 ;; All,i]]]]], {i, num}]; 
     Export[file, Join[strl, strl0], "Text"]; Print["Output: ", file]; ]
DoubleListToString[a_List, ndigit_:5] := Block[{as = ""}, 
    Do[as = StringJoin[as, " ", Internal`DoubleToString[a[[i]], False, 
          ndigit]]; , {i, Length[a]}]; StringReplace[as, "*^" -> "e"]]
PlaneStress[Es_, mu_, flag_:True] := If[flag, (Es/(1. - mu^2))*{mu, 1. - mu}, 
    (Es/(1. + mu)/(1. - 2*mu))*{mu, 1. - 2*mu}]
Plot2DField[vals_List, mesh_ElementMesh] := 
   Block[{vi, vm, color = "Rainbow", ntick = 10}, {vm, vi} = MinMax[vals]; 
     Legended[Graphics[ElementMeshToGraphicsComplex[mesh, 
        VertexColors -> ColorData[color] /@ ((1./(vi - vm + 10.^(-20)))*
           (vals - vm))]], BarLegend[{color, {vm, vi}}, ntick]]]
 
Plot2DField[vals_List, mesh_ElementMesh, uvw_List] := 
   Plot2DField[vals, DeformMesh[mesh, uvw]]
 
Plot2DField[vals_List, mesh_List] := 
   Block[{}, Plot2DField[vals, MyDelaunayMesh[mesh]]]
 
Plot2DField[vals_List, mesh_List, uvw_List] := 
   Block[{}, Plot2DField[vals, MyDelaunayMesh[mesh, uvw]]]
DeformMesh[mesh_ElementMesh, uvw_List] := 
   ToElementMesh["Coordinates" -> mesh["Coordinates"] + uvw, 
    "MeshElements" -> {mesh["MeshElements"][[1]]}]
StrainEnergyCal[stress_, deform_, vol_] := 
   Block[{res = 0., I3 = IdentityMatrix[Length[deform[[1]]]]}, 
    Do[res += 0.5*vol[[i]]*Total[stress[[i]]*(deform[[i]] - I3), -1], 
      {i, Length[vol]}]; res]
StressList[bulkK_, shearG_, fmat_] := 
   Block[{emat, sl, imat = IdentityMatrix[Length[fmat[[1]]]]}, 
    sl = Table[0., {Length[fmat]}, {Length[imat]}, {Length[imat]}]; 
     Do[emat = 0.5*(fmat[[i]] + Transpose[fmat[[i]]]) - imat; 
       sl[[i]] = shearG*emat + (bulkK*Tr[emat])*imat; , {i, Length[fmat]}]; 
     sl]
FindPointsC1 = CompiledFunction[{11, 12., 5468}, {{_Real, 1}, {_Real, 1}, 
     {_Real, 1}}, {{3, 1, 0}, {3, 1, 1}, {3, 1, 2}, {2, 2, 4}}, 
    {{0, {2, 0, 5}}, {1, {2, 0, 3}}}, {0, 8, 3, 0, 6}, 
    {{33, 0, 6}, {6, 5, 7}, {35, 6, 2, 3}, {6, 5, 2}, {3, 2}, 
     {36, 7, 5, 2, 3}, {4, 2, 6, -1}, {33, 0, 7}, {6, 5, 2}, {3, 2}, 
     {46, Function[{v1, xmin, xmax}, If[v1[[i]] > xmin && v1[[i]] < xmax, 
        vt[[i]] = 1]], {i, 2, 0, 2, Block}, {vt, 2, 1, 3, Block}, 3, 1, 0, 3, 
      1, 1, 3, 1, 2, 6, 0, 17}, {4, 2, 7, -1}, {42, "Position", 2, 1, 3, -2, 
      0, 3, 2, 2, 4}, {1}}, Function[{v1, xmin, xmax}, 
     Block[{y, r, vt = Table[0, {Length[v1]}]}, 
      Do[If[v1[[i]] > xmin && v1[[i]] < xmax, vt[[i]] = 1]; , 
        {i, Length[v1]}]; Position[vt, 1]]], Evaluate]
FindPointsC2 = CompiledFunction[{11, 12., 5468}, {{_Real, 2}, {_Real, 1}, 
     {_Real, 1}}, {{3, 2, 0}, {3, 1, 1}, {3, 1, 2}, {2, 2, 5}}, 
    {{False, {1, 0, 4}}, {0, {2, 0, 5}}, {2, {2, 0, 7}}, {1, {2, 0, 3}}, 
     {7., {3, 0, 2}}}, {6, 8, 3, 0, 6}, {{33, 0, 6}, {6, 5, 2}, 
     {35, 6, 2, 3}, {6, 5, 4}, {3, 2}, {36, 2, 5, 2, 3}, {4, 4, 6, -1}, 
     {33, 0, 2}, {6, 5, 4}, {3, 28}, {38, 0, 0, 4, 0, 3, 0, 0}, 
     {38, 1, 0, 3, 0, 1}, {27, 7, 2, 0, 1, 1}, {2, 1, 20}, 
     {38, 0, 0, 4, 0, 7, 0, 0}, {38, 1, 0, 7, 0, 1}, {27, 7, 2, 0, 1, 0}, 
     {2, 0, 13}, {38, 0, 0, 4, 0, 3, 0, 0}, {38, 2, 0, 3, 0, 1}, 
     {27, 3, 2, 0, 1, 2}, {2, 2, 6}, {38, 0, 0, 4, 0, 7, 0, 0}, 
     {38, 2, 0, 7, 0, 1}, {27, 3, 2, 0, 1, 5}, {5, 5, 3}, {3, 2}, {5, 4, 3}, 
     {5, 3, 2}, {3, 2}, {5, 4, 2}, {5, 2, 0}, {3, 2}, {5, 4, 0}, {2, 0, 3}, 
     {39, 3, 0, 4, 0, 3}, {3, 1}, {4, 4, 2, -27}, {42, "Position", 2, 1, 3, 
      -2, 0, 3, 2, 2, 5}, {1}}, Function[{v1, xmin, xmax}, 
     Block[{y, r, vt = Table[0, {Length[v1]}]}, 
      Do[If[v1[[i,1]] > xmin[[1]] && v1[[i,2]] > xmin[[2]] && 
           v1[[i,1]] < xmax[[1]] && v1[[i,2]] < xmax[[2]], vt[[i]] = 1]; , 
        {i, Length[v1]}]; Position[vt, 1]]], Evaluate]
FindPointsC3 = CompiledFunction[{11, 12., 5468}, {{_Real, 2}, {_Real, 1}, 
     {_Real, 1}}, {{3, 2, 0}, {3, 1, 1}, {3, 1, 2}, {2, 2, 5}}, 
    {{False, {1, 0, 6}}, {0, {2, 0, 5}}, {2, {2, 0, 7}}, {1, {2, 0, 3}}, 
     {7., {3, 0, 2}}, {3, {2, 0, 8}}}, {8, 9, 3, 0, 6}, 
    {{33, 0, 6}, {6, 5, 2}, {35, 6, 2, 3}, {6, 5, 4}, {3, 2}, 
     {36, 2, 5, 2, 3}, {4, 4, 6, -1}, {33, 0, 2}, {6, 5, 4}, {3, 42}, 
     {38, 0, 0, 4, 0, 3, 0, 0}, {38, 1, 0, 3, 0, 1}, {27, 7, 2, 0, 1, 1}, 
     {2, 1, 34}, {38, 0, 0, 4, 0, 7, 0, 0}, {38, 1, 0, 7, 0, 1}, 
     {27, 7, 2, 0, 1, 0}, {2, 0, 27}, {38, 0, 0, 4, 0, 8, 0, 0}, 
     {38, 1, 0, 8, 0, 1}, {27, 7, 2, 0, 1, 2}, {2, 2, 20}, 
     {38, 0, 0, 4, 0, 3, 0, 0}, {38, 2, 0, 3, 0, 1}, {27, 3, 2, 0, 1, 3}, 
     {2, 3, 13}, {38, 0, 0, 4, 0, 7, 0, 0}, {38, 2, 0, 7, 0, 1}, 
     {27, 3, 2, 0, 1, 4}, {2, 4, 6}, {38, 0, 0, 4, 0, 8, 0, 0}, 
     {38, 2, 0, 8, 0, 1}, {27, 3, 2, 0, 1, 7}, {5, 7, 5}, {3, 2}, {5, 6, 5}, 
     {5, 5, 4}, {3, 2}, {5, 6, 4}, {5, 4, 3}, {3, 2}, {5, 6, 3}, {5, 3, 2}, 
     {3, 2}, {5, 6, 2}, {5, 2, 0}, {3, 2}, {5, 6, 0}, {2, 0, 3}, 
     {39, 3, 0, 4, 0, 3}, {3, 1}, {4, 4, 2, -41}, {42, "Position", 2, 1, 3, 
      -2, 0, 3, 2, 2, 5}, {1}}, Function[{v1, xmin, xmax}, 
     Block[{y, r, vt = Table[0, {Length[v1]}]}, 
      Do[If[v1[[i,1]] > xmin[[1]] && v1[[i,2]] > xmin[[2]] && 
           v1[[i,3]] > xmin[[3]] && v1[[i,1]] < xmax[[1]] && 
           v1[[i,2]] < xmax[[2]] && v1[[i,3]] < xmax[[3]], vt[[i]] = 1]; , 
        {i, Length[v1]}]; Position[vt, 1]]], Evaluate]
ForceAccumulateHgC = CompiledFunction[{11, 12., 5788}, 
    {{_Real, 2}, {_Real, 2}, {_Real, 2}, {_Real, 2}, {_Real, 2}, _Real, 
     {_Real, 1}, _Real, {_Real, 1}}, {{3, 2, 0}, {3, 2, 1}, {3, 2, 2}, 
     {3, 2, 3}, {3, 2, 4}, {3, 0, 0}, {3, 1, 5}, {3, 0, 1}, {3, 1, 6}, 
     {3, 2, 8}}, {{0, {2, 0, 10}}, {1., {3, 0, 4}}, {4, {2, 0, 13}}, 
     {2, {2, 0, 12}}, {-1, {2, 0, 14}}, {1, {2, 0, 1}}, {0., {3, 0, 2}}}, 
    {0, 15, 7, 0, 16}, {{33, 1, 2}, {38, 1, 0, 1, 1, 8}, {33, 8, 0}, 
     {6, 2, 6}, {12, 0, 1, 7}, {6, 10, 9}, {35, 6, 7, 3, 8}, {6, 10, 3}, 
     {3, 5}, {6, 10, 11}, {3, 2}, {36, 9, 2, 3, 8}, {4, 11, 7, -1}, 
     {4, 3, 6, -4}, {42, "CopyTensor", 3, 2, 3, 3, 2, 11}, {6, 0, 5}, 
     {6, 10, 8}, {3, 4}, {38, 11, 0, 8, 0, 8, 0, 5}, 
     {41, 258, 3, 0, 5, 3, 0, 4, 3, 0, 3}, {39, 11, 0, 8, 0, 8, 0, 3}, 
     {4, 8, 5, -3}, {6, 2, 5}, {6, 1, 8}, {3, 34}, {38, 1, 0, 8, 1, 9}, 
     {38, 1, 0, 1, 1, 14}, {40, 43, 3, 1, 14, 3, 1, 12}, {44, 9, 12, 9}, 
     {38, 5, 0, 1, 0, 3}, {38, 5, 0, 8, 0, 5}, {38, 6, 0, 8, 0, 6}, 
     {16, 3, 5, 6, 3}, {42, "Dot", 3, 2, 11, 3, 1, 9, 2, 0, 13, 3, 1, 12}, 
     {38, 0, 0, 8, 1, 14}, {38, 0, 0, 1, 1, 10}, {40, 43, 3, 1, 10, 3, 1, 
      13}, {44, 14, 13, 14}, {40, 43, 3, 1, 14, 3, 1, 13}, {44, 12, 13, 12}, 
     {42, "Dot", 3, 2, 4, 3, 1, 9, 2, 0, 13, 3, 1, 13}, 
     {42, "Dot", 3, 2, 2, 3, 1, 13, 2, 0, 13, 3, 1, 14}, 
     {40, 60, 3, 0, 0, 3, 0, 5}, {16, 1, 5, 6}, {41, 259, 3, 0, 6, 3, 1, 12, 
      3, 1, 13}, {40, 43, 3, 1, 13, 3, 1, 10}, {44, 14, 10, 14}, 
     {41, 259, 3, 0, 3, 3, 1, 14, 3, 1, 10}, {42, "DotVV", 3, 1, 12, 3, 1, 
      12, 2, 0, 13, 3, 0, 6}, {19, 6, 5}, {34, 1, 1, 14, 2, 14}, 
     {34, 1, 1, 14, 0, 13}, {42, "Insert", 3, 1, 10, -3, 0, 5, 2, 2, 13, 3, 
      1, 15}, {38, 8, 0, 1, 1, 10}, {44, 10, 15, 10}, {39, 8, 0, 1, 1, 10}, 
     {40, 43, 3, 1, 15, 3, 1, 14}, {39, 8, 0, 8, 1, 14}, {4, 8, 5, -33}, 
     {1}}, Function[{uvwi, coordi, stress, Fmat, invK, trK, voli, hgPen, 
      WeiL}, Block[{fi, fl, num = Length[coordi], ndim = Length[coordi[[1]]], 
       xi, ri, vi, vw, NablaU}, fl = Table[0., {num}, {ndim + 1}]; 
       NablaU = Fmat; Do[NablaU[[i,i]] -= 1., {i, ndim}]; 
       Do[xi = coordi[[j]] - coordi[[1]]; vw = voli[[1]]*voli[[j]]*WeiL[[j]]; 
         vi = NablaU . xi - (uvwi[[j]] - uvwi[[1]]); 
         fi = Append[vw*(stress . (invK . xi) - (hgPen/trK)*vi), -vi . vi]; 
         fl[[1]] += fi; fl[[j]] = -fi; , {j, 2, num}]; fl], Listable], 
    Evaluate]
KmatrixC = CompiledFunction[{11, 12., 5788}, {{_Real, 2}, {_Real, 1}, 
     {_Real, 1}}, {{3, 2, 0}, {3, 1, 1}, {3, 1, 2}, {3, 1, 6}}, 
    {{0, {2, 0, 7}}, {1., {3, 0, 6}}, {{2, -1, 1, 0}, {2, 1, 8}}, 
     {2, {2, 0, 9}}, {1, {2, 0, 1}}, {3, {2, 0, 11}}, {0., {3, 0, 0}}}, 
    {1, 12, 25, 0, 13}, {{33, 0, 2}, {38, 0, 0, 1, 1, 9}, {33, 9, 0}, 
     {7, 0, 3}, {7, 0, 11}, {7, 0, 12}, {7, 0, 8}, {7, 0, 9}, {7, 0, 2}, 
     {7, 0, 4}, {7, 0, 7}, {7, 0, 1}, {7, 0, 5}, {6, 0, 10}, {6, 0, 8}, 
     {6, 7, 6}, {35, 10, 8, 3, 9}, {6, 7, 3}, {3, 5}, {6, 7, 4}, {3, 2}, 
     {36, 6, 0, 3, 9}, {4, 4, 8, -1}, {4, 3, 10, -4}, {15, 0, 0, 6}, 
     {12, 1, 6, 10}, {6, 7, 8}, {35, 10, 3, 6}, {6, 7, 6}, {3, 2}, 
     {36, 8, 0, 3, 6}, {4, 6, 10, -1}, {6, 2, 4}, {6, 1, 3}, {3, 22}, 
     {38, 0, 0, 3, 1, 10}, {38, 0, 0, 1, 1, 7}, {40, 43, 3, 1, 7, 3, 1, 11}, 
     {44, 10, 11, 10}, {38, 2, 0, 3, 0, 13}, {38, 1, 0, 3, 0, 14}, 
     {16, 13, 14, 13}, {6, 0, 8}, {6, 7, 6}, {3, 11}, {6, 0, 10}, {6, 7, 5}, 
     {3, 7}, {38, 9, 0, 6, 0, 5, 0, 14}, {38, 10, 0, 6, 0, 17}, 
     {38, 10, 0, 5, 0, 15}, {16, 13, 17, 15, 19}, {13, 14, 19, 14}, 
     {39, 9, 0, 6, 0, 5, 0, 14}, {4, 5, 10, -6}, {4, 6, 8, -10}, 
     {4, 3, 4, -21}, {38, 9, 0, 1, 0, 1, 0, 14}, {7, 14, 11}, 
     {38, 9, 0, 9, 0, 9, 0, 14}, {7, 14, 9}, {38, 9, 0, 1, 0, 9, 0, 14}, 
     {7, 14, 12}, {38, 9, 0, 9, 0, 1, 0, 14}, {7, 14, 8}, {24, 0, 11, 0}, 
     {2, 0, 12}, {38, 9, 0, 1, 0, 11, 0, 14}, {7, 14, 2}, 
     {38, 9, 0, 11, 0, 1, 0, 14}, {7, 14, 4}, {38, 9, 0, 9, 0, 11, 0, 14}, 
     {7, 14, 7}, {38, 9, 0, 11, 0, 9, 0, 14}, {7, 14, 1}, 
     {38, 9, 0, 11, 0, 11, 0, 14}, {7, 14, 5}, {3, 1}, {24, 0, 9, 0}, 
     {2, 0, 9}, {19, 12, 14}, {16, 14, 8, 14}, {16, 11, 9, 19}, 
     {13, 14, 19, 14}, {40, 60, 3, 0, 14, 3, 0, 19}, {16, 6, 19, 14}, 
     {7, 14, 3}, {3, 14}, {19, 2, 14}, {16, 14, 9, 4, 14}, 
     {16, 12, 7, 4, 19}, {16, 2, 8, 1, 17}, {16, 11, 7, 1, 15}, {19, 15, 18}, 
     {16, 12, 8, 5, 15}, {19, 15, 20}, {16, 11, 9, 5, 15}, 
     {13, 14, 19, 17, 18, 20, 15, 14}, {40, 60, 3, 0, 14, 3, 0, 19}, 
     {16, 6, 19, 14}, {7, 14, 3}, {24, 0, 9, 0}, {2, 0, 9}, {13, 11, 9, 14}, 
     {39, 6, 0, 1, 0, 14}, {19, 12, 19}, {19, 8, 17}, 
     {34, 1, 4, 9, 19, 17, 11, 3, 11}, {41, 259, 3, 0, 3, 3, 1, 11, 3, 1, 7}, 
     {39, 6, 3, 8, 1, 7}, {3, 42}, {13, 11, 9, 5, 14}, {39, 6, 0, 1, 0, 14}, 
     {19, 2, 19}, {7, 19, 10}, {19, 7, 19}, {16, 19, 1, 19}, {16, 9, 5, 17}, 
     {13, 19, 17, 19}, {16, 2, 1, 17}, {16, 12, 5, 18}, {19, 18, 20}, 
     {13, 17, 20, 17}, {16, 10, 9, 20}, {16, 12, 7, 18}, {13, 20, 18, 20}, 
     {16, 7, 4, 18}, {16, 8, 5, 15}, {19, 15, 16}, {13, 18, 16, 18}, 
     {16, 10, 4, 16}, {16, 11, 5, 15}, {13, 16, 15, 16}, {16, 2, 8, 15}, 
     {16, 11, 7, 22}, {19, 22, 21}, {13, 15, 21, 15}, {19, 9, 21}, 
     {16, 21, 4, 21}, {16, 8, 1, 22}, {13, 21, 22, 21}, {16, 12, 4, 22}, 
     {16, 11, 1, 24}, {19, 24, 23}, {13, 22, 23, 22}, {19, 12, 23}, 
     {16, 23, 8, 23}, {16, 11, 9, 24}, {13, 23, 24, 23}, 
     {34, 1, 9, 19, 17, 20, 18, 16, 15, 21, 22, 23, 3, 11}, 
     {41, 259, 3, 0, 3, 3, 1, 11, 3, 1, 12}, {39, 6, 3, 8, 1, 12}, {1}}, 
    Function[{coordi, WeiF, voli}, Block[{Compile`$63}, 
      Block[{num = Length[coordi], ndim = Length[coordi[[1]]], kmat, vw, xi, 
        ri, ik$det = 0., k$11 = 0., k$12 = 0., k$21 = 0., k$22 = 0., 
        k$13 = 0., k$31 = 0., k$23 = 0., k$32 = 0., k$33 = 0., invkv}, 
       kmat = Table[0., {ndim}, {ndim}]; invkv = Table[0., {1 + ndim*ndim}]; 
        Do[xi = coordi[[j]] - coordi[[1]]; vw = voli[[j]]*WeiF[[j]]; 
          Do[kmat[[i,jj]] += vw*xi[[i]]*xi[[jj]], {i, ndim}, {jj, ndim}]; , 
         {j, 2, num}]; k$11 = kmat[[1,1]]; k$22 = kmat[[2,2]]; 
        k$12 = kmat[[1,2]]; k$21 = kmat[[2,1]]; If[ndim == 3, 
         k$13 = kmat[[1,3]]; k$31 = kmat[[3,1]]; k$23 = kmat[[2,3]]; 
          k$32 = kmat[[3,2]]; k$33 = kmat[[3,3]]]; If[ndim == 2, 
         ik$det = 1./((-k$12)*k$21 + k$11*k$22), 
         ik$det = 1./((-k$13)*k$22*k$31 + k$12*k$23*k$31 + k$13*k$21*k$32 - 
            k$11*k$23*k$32 - k$12*k$21*k$33 + k$11*k$22*k$33)]; 
        If[ndim == 2, invkv[[1]] = k$11 + k$22; invkv[[2 ;; -1]] = 
           ik$det*{k$22, -k$12, -k$21, k$11}, 
         invkv[[1]] = k$11 + k$22 + k$33; invkv[[2 ;; -1]] = 
           (Compile`$63 = -k$13; ik$det*{(-k$23)*k$32 + k$22*k$33, 
              k$13*k$32 - k$12*k$33, Compile`$63*k$22 + k$12*k$23, 
              k$23*k$31 - k$21*k$33, Compile`$63*k$31 + k$11*k$33, 
              k$13*k$21 - k$11*k$23, (-k$22)*k$31 + k$21*k$32, 
              k$12*k$31 - k$11*k$32, (-k$12)*k$21 + k$11*k$22})]; invkv]], 
     Listable], Evaluate]
NeighborPartitionC1 = CompiledFunction[{11, 12., 5468}, 
    {{_Real, 1}, {_Integer, 2}}, {{3, 1, 0}, {2, 2, 1}, {3, 2, 2}}, 
    {{0, {2, 0, 6}}, {-1, {2, 0, 5}}, {1, {2, 0, 3}}}, {0, 8, 0, 0, 6}, 
    {{33, 1, 2}, {6, 5, 4}, {35, 2, 4, 3, 2}, {6, 6, 7}, {3, 4}, 
     {38, 1, 0, 7, 1, 4}, {38, 0, 1, 4, 1, 5}, {36, 4, 5, 0, 2}, 
     {4, 7, 2, -3}, {1}}, Function[{coord, Nei}, Table[coord[[Nei[[i]]]], 
      {i, Length[Nei]}]], Evaluate]
NeighborPartitionC2 = CompiledFunction[{11, 12., 5468}, 
    {{_Real, 2}, {_Integer, 2}}, {{3, 2, 0}, {2, 2, 1}, {3, 3, 2}}, 
    {{0, {2, 0, 6}}, {-1, {2, 0, 5}}, {1, {2, 0, 3}}}, {0, 8, 0, 0, 6}, 
    {{33, 1, 2}, {6, 5, 4}, {35, 2, 3, 4, 3, 2}, {6, 6, 7}, {3, 4}, 
     {38, 1, 0, 7, 1, 4}, {38, 0, 1, 4, 2, 5}, {36, 4, 5, 0, 2}, 
     {4, 7, 2, -3}, {1}}, Function[{coord, Nei}, Table[coord[[Nei[[i]]]], 
      {i, Length[Nei]}]], Evaluate]
NeighborPartitionC3 = CompiledFunction[{11, 12., 5468}, 
    {{_Real, 3}, {_Integer, 2}}, {{3, 3, 0}, {2, 2, 1}, {3, 4, 2}}, 
    {{0, {2, 0, 6}}, {-1, {2, 0, 5}}, {1, {2, 0, 3}}}, {0, 8, 0, 0, 6}, 
    {{33, 1, 2}, {6, 5, 4}, {35, 2, 3, 3, 4, 3, 2}, {6, 6, 7}, {3, 4}, 
     {38, 1, 0, 7, 1, 4}, {38, 0, 1, 4, 3, 5}, {36, 4, 5, 0, 2}, 
     {4, 7, 2, -3}, {1}}, Function[{coord, Nei}, Table[coord[[Nei[[i]]]], 
      {i, Length[Nei]}]], Evaluate]
CheckCriticalStretch = CompiledFunction[{11, 12., 5788}, 
    {{_Real, 2}, {_Real, 2}}, {{3, 2, 0}, {3, 2, 1}, {3, 0, 5}}, 
    {{4, {2, 0, 5}}, {-1., {3, 0, 0}}, {2, {2, 0, 0}}, {1, {2, 0, 3}}, 
     {7., {3, 0, 3}}, {0., {3, 0, 2}}}, {1, 6, 7, 0, 6}, 
    {{7, 0, 5}, {33, 0, 2}, {6, 3, 4}, {3, 22}, {38, 0, 0, 4, 1, 4}, 
     {38, 0, 0, 3, 1, 2}, {40, 43, 3, 1, 2, 3, 1, 5}, {44, 4, 5, 4}, 
     {38, 1, 0, 4, 1, 5}, {38, 1, 0, 3, 1, 2}, {40, 43, 3, 1, 2, 3, 1, 3}, 
     {44, 5, 3, 5}, {42, "DotVV", 3, 1, 4, 3, 1, 5, 2, 0, 5, 3, 0, 1}, 
     {27, 7, 3, 1, 2, 0}, {2, 0, 8}, {42, "DotVV", 3, 1, 4, 3, 1, 4, 2, 0, 5, 
      3, 0, 1}, {42, "DotVV", 3, 1, 5, 3, 1, 5, 2, 0, 5, 3, 0, 4}, 
     {40, 60, 3, 0, 4, 3, 0, 6}, {16, 1, 6, 1}, {40, 57, 3, 0, 1, 3, 0, 6}, 
     {7, 6, 1}, {3, 2}, {7, 2, 1}, {42, "MaxR", 3, 0, 5, 3, 0, 1, 3, 0, 6}, 
     {7, 6, 5}, {4, 4, 2, -21}, {1}}, Function[{uvwi, coordi}, 
     Block[{s = -1., uij, rij}, Do[uij = uvwi[[j]] - uvwi[[1]]; 
         rij = coordi[[j]] - coordi[[1]]; 
         s = Max[s, If[uij . rij > 0., Sqrt[uij . uij/rij . rij], 0.]], 
        {j, 2, Length[uvwi]}]; s], Listable], Evaluate]
FmatrixC = CompiledFunction[{11, 12., 5788}, {{_Real, 2}, {_Real, 2}, 
     {_Real, 1}, {_Real, 2}, {_Real, 1}}, {{3, 2, 0}, {3, 2, 1}, {3, 1, 2}, 
     {3, 2, 3}, {3, 1, 4}, {3, 2, 8}}, {{0, {2, 0, 7}}, {4, {2, 0, 11}}, 
     {2, {2, 0, 9}}, {1, {2, 0, 1}}, {0., {3, 0, 0}}}, {0, 12, 6, 0, 11}, 
    {{33, 1, 2}, {38, 1, 0, 1, 1, 8}, {33, 8, 0}, {6, 0, 5}, {6, 0, 3}, 
     {6, 7, 10}, {35, 5, 3, 3, 8}, {6, 7, 4}, {3, 5}, {6, 7, 8}, {3, 2}, 
     {36, 10, 0, 3, 8}, {4, 8, 3, -1}, {4, 4, 5, -4}, {6, 2, 10}, {6, 1, 8}, 
     {3, 26}, {38, 1, 0, 8, 1, 6}, {38, 1, 0, 1, 1, 9}, 
     {40, 43, 3, 1, 9, 3, 1, 7}, {44, 6, 7, 6}, {38, 0, 0, 8, 1, 7}, 
     {38, 0, 0, 1, 1, 9}, {40, 43, 3, 1, 9, 3, 1, 10}, {44, 7, 10, 7}, 
     {38, 2, 0, 8, 0, 1}, {38, 4, 0, 8, 0, 2}, {16, 1, 2, 1}, {6, 0, 4}, 
     {6, 7, 5}, {3, 11}, {6, 0, 3}, {6, 7, 6}, {3, 7}, 
     {38, 8, 0, 5, 0, 6, 0, 2}, {38, 7, 0, 5, 0, 5}, {38, 6, 0, 6, 0, 3}, 
     {16, 1, 5, 3, 4}, {13, 2, 4, 2}, {39, 8, 0, 5, 0, 6, 0, 2}, 
     {4, 6, 3, -6}, {4, 5, 4, -10}, {4, 8, 10, -25}, 
     {42, "Dot", 3, 2, 8, 3, 2, 3, 2, 0, 11, 3, 2, 10}, 
     {42, "CopyTensor", 3, 2, 10, 3, 2, 8}, {6, 0, 10}, {6, 7, 8}, {3, 5}, 
     {38, 8, 0, 8, 0, 8, 0, 2}, {10, 1, 4}, {13, 2, 4, 2}, 
     {39, 8, 0, 8, 0, 8, 0, 2}, {4, 8, 10, -4}, {1}}, 
    Function[{uvwi, coordi, voli, invK, WeiL}, 
     Block[{num = Length[coordi], ndim = Length[coordi[[1]]], fmat, vw, xi, 
       ri, ui}, fmat = Table[0., {ndim}, {ndim}]; 
       Do[xi = coordi[[j]] - coordi[[1]]; ui = uvwi[[j]] - uvwi[[1]]; 
         vw = voli[[j]]*WeiL[[j]]; Do[fmat[[i,j]] += vw*ui[[i]]*xi[[j]], 
          {i, ndim}, {j, ndim}]; , {j, 2, num}]; fmat = fmat . invK; 
       Do[fmat[[i,i]] += 1, {i, ndim}]; fmat], Listable], Evaluate]
WeightByNeighborC2 = CompiledFunction[{11, 12., 5788}, {{_Real, 2}}, 
    {{3, 2, 0}, {3, 1, 3}}, {{0, {2, 0, 5}}, {4, {2, 0, 8}}, {2, {2, 0, 7}}, 
     {1, {2, 0, 0}}, {0., {3, 0, 1}}, {1.*^-10, {3, 0, 0}}}, {0, 9, 4, 0, 8}, 
    {{42, "Dimensions", 3, 2, 0, 2, 1, 1}, {38, 1, 0, 0, 0, 1}, {6, 5, 3}, 
     {35, 1, 3, 3}, {6, 5, 4}, {3, 2}, {36, 3, 0, 3, 3}, {4, 4, 1, -1}, 
     {38, 1, 0, 7, 0, 3}, {6, 5, 1}, {35, 3, 3, 4}, {6, 5, 4}, {3, 2}, 
     {36, 1, 1, 3, 4}, {4, 4, 3, -1}, {33, 0, 6}, {6, 0, 1}, {3, 9}, 
     {38, 0, 0, 1, 1, 6}, {38, 0, 0, 0, 1, 5}, {40, 43, 3, 1, 5, 3, 1, 7}, 
     {44, 6, 7, 6}, {42, "CopyTensor", 3, 1, 6, 3, 1, 4}, 
     {42, "DotVV", 3, 1, 4, 3, 1, 4, 2, 0, 8, 3, 0, 3}, 
     {40, 60, 3, 0, 3, 3, 0, 2}, {39, 3, 0, 1, 0, 2}, {4, 1, 6, -8}, {1}}, 
    Function[{ni}, Block[{res, rii, dm = Dimensions[ni]}, 
      res = Table[1.*^-10, {dm[[1]]}]; rii = Table[0., {dm[[2]]}]; 
       Do[rii = ni[[j]] - ni[[1]]; res[[j]] = 1/rii . rii; (*default weight function WeiF[r_]:=1/r^3 is used. Other 
       othe weight function is applicable as well.*), 
        {j, 2, Length[ni]}]; res], Listable], Evaluate]
        
        
        PlaneListCut2P[p3dL_List,p1_,p2_]:=Block[{fd},Do[If[PlaneCut2P[p3dL[[i]],p1,p2],fd=True;Break[]],{i,Length[p3dL]}];fd];

PlaneCut2P[p3d_List,p4_,p5_]:=Block[{p0,n1,r0,rn,n45,nn,d,p,dis,p1,p2,p3,area},If[Length[p3d]==7,(* ///circle plane in 3d:x0,y0,z0,n1,n2,n3,r0;*)p0=p3d[[1;;3]];n1=p3d[[4;;6]];r0=p3d[[7]];
rn=Sqrt[n1.n1];
n1=1./rn*n1;
n45=p5-p4;
nn=n1.n45;
If[Abs[nn]<10^-5,Return[False]];
d=1/nn*(p0-p4).n1;
If[d>=1||d<0.,Return[False]];
p=d n45+p4;
dis=Norm[p-p0];
If[dis<=r0,Return[True],Return[False]]];
If[Length[p3d]==12,(*//plane3D:x1 y1 z1,x2 y2 z2,x3 y3 z3,x4 y4 z4;*)p0=p3d[[1;;3]];p1=p3d[[4;;6]];p2=p3d[[7;;9]];
p3=p3d[[10;;12]];n1=Cross[(p1-p0),p2-p0];
n1*=1/Sqrt[n1.n1];
n45=p5-p4;
nn=n1.n45;
If[Abs[nn]<10^-5,Return[False]];
d=1/nn*(p0-p4).n1;
If[d>=1.||d<0.,Return[False]];
p=d*n45+p4;
(*area={Area3point[p0,p1,p],Area3point[p1,p2,p],Area3point[p2,p3,p],Area3point[p3,p0,p]};*)(*Cross[p1-p0,p-p0],Cross[p3-p2,p-p2],Cross[p2-p1,p-p1],Cross[p0-p3,p-p3]*)If[Cross[p1-p0,p-p0].Cross[p3-p2,p-p2]>=0.&&Cross[p2-p1,p-p1].Cross[p0-p3,p-p3]>=0.,Return[True],Return[False]]];
Print["Error,the plane is not well defined",p3d];Abort[];];
Area3point[p1_,p2_,p3_]:=0.5 Norm[Cross[p2-p1,p3-p1]];
PointOnTriangle[v1_,p1_,p2_,p3_,eps_]:=(*3D case;-1,-normal direction;0,not on;1,normal direction;*)Block[{x,y,z,d,x1,y1,z1,x2,y2,z2,x3,y3,z3,A1,A2,A,normal},{x2,y2,z2}=p2-p1;{x3,y3,z3}=p3-p1;
normal=Cross[{x2,y2,z2},{x3,y3,z3}];
{x,y,z}=v1-p1;d=normal.{x,y,z};
If[d<=-eps||d>=eps,Return[0]];{x1,y1,z1}={x,y,z}-d*normal;
A1=normal.{y2 z1-y1 z2,z2 x1-z1 x2,x2 y1-x1 y2};
A2=normal.{y1 z3-y3 z1,z1 x3-z3 x1,x1 y3-x3 y1};
A=normal.{y2 z3-y3 z2,z2 x3-z3 x2,x2 y3-x3 y2};
If[Sign[A1]!=Sign[A2],Return[0]];
If[Abs[A1+A2]>Abs[A],Return[0]];
If[d>=0,Return[1],Return[2]]];

SimplexVolume[coord_List]:=Module[{n,v1={}},n=Length[coord];
Do[AppendTo[v1,coord[[i]]-coord[[1]]],{i,2,n}];
Abs[Det[v1//N]]/((n-1)!)];

Area3point[p1_,p2_,p3_]:=0.5 Norm[Cross[p2-p1,p3-p1]];
PointOnTriangle[v1_,p1_,p2_,p3_,eps_]:=(*3D case;-1,-normal direction;0,not on;1,normal direction;*)Block[{x,y,z,d,x1,y1,z1,x2,y2,z2,x3,y3,z3,A1,A2,A,normal},{x2,y2,z2}=p2-p1;{x3,y3,z3}=p3-p1;
normal=Cross[{x2,y2,z2},{x3,y3,z3}];
{x,y,z}=v1-p1;d=normal.{x,y,z};
If[d<=-eps||d>=eps,Return[0]];{x1,y1,z1}={x,y,z}-d*normal;
A1=normal.{y2 z1-y1 z2,z2 x1-z1 x2,x2 y1-x1 y2};
A2=normal.{y1 z3-y3 z1,z1 x3-z3 x1,x1 y3-x3 y1};
A=normal.{y2 z3-y3 z2,z2 x3-z3 x2,x2 y3-x3 y2};
If[Sign[A1]!=Sign[A2],Return[0]];
If[Abs[A1+A2]>Abs[A],Return[0]];
If[d>=0,Return[1],Return[2]]];
Needs["NDSolve`FEM`"]
Plot3DField[vals_List,mesh_ElementMesh,ndim_:3,nnode_:1]:=Block[{len},(*If[ndim\[NotEqual]3,Print["Plot3DField Error, only for 3D mesh"];Return[];];If[Length[vals]\[NotEqual]nnode,Print["Plot3DField Error, the number of nodes is not the same as the length of vals"];Return[]];*)
Legended[ElementMeshSurfacePlot3D[vals,mesh,Boxed->False,Axes->True,AxesLabel->{x,y,z}],BarLegend[{ColorFunction/.Options[ElementMeshSurfacePlot3D]//First,MinMax@vals}]]];
Plot3DField[vals_List,mesh_ElementMesh,title_String]:=Block[{len},(*If[ndim\[NotEqual]3,Print["Plot3DField Error, only for 3D mesh"];Return[];];If[Length[vals]\[NotEqual]nnode,Print["Plot3DField Error, the number of nodes is not the same as the length of vals"];Return[]];*)
Legended[ElementMeshSurfacePlot3D[vals,mesh,Boxed->False,Axes->True,AxesLabel->{x,y,z},PlotLabel->title],BarLegend[{ColorFunction/.Options[ElementMeshSurfacePlot3D]//First,MinMax@vals}]]];
Plot3DField[vals_List,mesh_ElementMesh,uvw_List,color_:"Rainbow",ntick_:10]:=Plot3DField[vals,DeformMesh[mesh,uvw],3,Length[vals]];(*Block[{m2},(*m2=ToElementMesh["Coordinates"\[Rule]mesh["Coordinates"]+uvw,"MeshElements"\[Rule]{mesh["MeshElements"][[1]]}];*)Plot3DField[vals,DeformMesh[mesh,uvw],3,Length[vals]]];*)
Plot3DField[vals_List,mesh_List]:=Block[{},Plot3DField[vals,MyDelaunayMesh[mesh],Length[mesh[[1]]],Length[mesh]]];
Plot3DField[vals_List,mesh_List,uvw_List]:=Block[{},Plot3DField[vals,MyDelaunayMesh[mesh,uvw],Length[mesh[[1]]],Length[mesh]]];
ParseAbaqusFile::usage="{nodes,elements,sets}=ParseAbaqusFile[FilePath];it can only read keywords: *node, *element, *set. Before processing the file, please convert the *elset to *set and *nset to *set in the input file";
ParseAbaqusFile[file_String]:=Module[{i,j,k,node,element,setList,oneSet,starList,starLine,tlineStart,tlineEnd,tline,strings,String2List, i1,i2,len,pair,ei,j1},
If[FileExtension[file]!="inp",Print["Error,the file should be Abaqus keyword file .inp"];Return[];];
strings=ReadList[file,String];If[strings==$Failed,Print["Error,file not found"];Return[]];i=Dimensions[strings][[1]];node={};
element={};
setList={};
starList={};
Do[If[StringStartsQ[strings[[j]],"*"]&&!StringStartsQ[strings[[j]],{"**","$"}],AppendTo[starList,j]],{j,i}];
String2List[str_String]:=Module[{sl,s2},s2=StringReplace[str,{"E","e"}->"*^"];sl=StringSplit[s2,{","," "}];
DeleteCases[ToExpression[sl],Null]];
Do[starLine=ToLowerCase[strings[[starList[[j]]]]];
tlineStart=starList[[j]]+1;tlineEnd=starList[[j+1]]-1;
Which[StringStartsQ[starLine,"*node"],Do[tline=strings[[k]];
If[!StringStartsQ[tline,"*"],AppendTo[node,String2List[tline]]];,{k,tlineStart,tlineEnd}](**),StringStartsQ[starLine,"*element"],Do[tline=strings[[k]];
If[!StringStartsQ[tline,"*"],AppendTo[element,String2List[tline]]];,{k,tlineStart,tlineEnd}](**),StringStartsQ[starLine,"*set"],oneSet={};
Do[tline=strings[[k]];
If[!StringStartsQ[tline,"*"],AppendTo[oneSet,String2List[tline]]];,{k,tlineStart,tlineEnd}];
AppendTo[setList,Flatten[oneSet]](**)];,{j,1,Length[starList]-1}];
len=Length[node];i1=node[[1,1]];i2=node[[len,1]];
If[i1==1 && i2==len,Return[{node[[All,2;;-1]],element[[All,2;;-1]],setList}]];
pair=ConstantArray[0,i2];
Do[pair[[node[[i,1]]]]=i,{i,len}];
len=Length[element];
Do[ei=element[[i]];Do[j1=ei[[j]];ei[[j]]=pair[[j1]];,{j,2,Length[ei]}];element[[i]]=ei;,{i,len}];
{node[[All,2;;-1]],element[[All,2;;-1]],setList}];
ParseAbaqusFile[file_String]:=Module[{j,node,element,setList,oneSet,starList={},starLine,tlineStart,tlineEnd,tline,strings,String2List,i1,i2,len,pair,ei,j1,nodeRange={},eleRange={},setRange={}},If[FileExtension[file]!="inp",Print["Error,the file should be Abaqus keyword file .inp"];
Return[];];
strings=ReadList[file,String];
If[strings==$Failed,Print["Error,file not found"];Return[]];
Do[If[StringStartsQ[strings[[j]],"*"],AppendTo[starList,j]],{j,Length[strings]}];
AppendTo[starList,Length[strings]];
Do[j=starList[[i]];tline=ToLowerCase[strings[[j]]];
If[StringStartsQ[tline,"*node"],AppendTo[nodeRange,Range[j+1,starList[[i+1]]-1]]];
If[StringStartsQ[tline,"*element"],AppendTo[eleRange,Range[j+1,starList[[i+1]]-1]]];
If[StringStartsQ[tline,"*set"],AppendTo[setRange,Range[j+1,starList[[i+1]]-1]]];,{i,Length[starList]-1}];
nodeRange=Flatten[nodeRange];
setRange=Flatten[setRange];
eleRange=Flatten[eleRange];
String2List[str_String]:=Module[{sl,s2},s2=StringReplace[str,{"E","e"}->"*^"];
sl=StringSplit[s2,{","," "}];
DeleteCases[ToExpression[sl],Null]];
SetAttributes[String2List,Listable];
node=If[Length[nodeRange]<100000,String2List[strings[[nodeRange]]],ParallelTable[String2List[strings[[i]]],{i,nodeRange}]];
element=If[Length[eleRange]<100000,String2List[strings[[eleRange]]],ParallelTable[String2List[strings[[i]]],{i,eleRange}]];
setList=String2List[strings[[setRange]]];
len=Length[node];i1=node[[1,1]];i2=node[[len,1]];
If[i1==1&&i2==len,Return[{node[[All,2;;-1]],element[[All,2;;-1]],setList}]];
pair=ConstantArray[0,i2];
Do[pair[[node[[i,1]]]]=i,{i,len}];
len=Length[element];
Do[ei=element[[i]];
Do[j1=ei[[j]];ei[[j]]=pair[[j1]];,{j,2,Length[ei]}];
element[[i]]=ei;,{i,len}];
{node[[All,2;;-1]],element[[All,2;;-1]],setList}];

AreaOf2DElement[nodes_List,elements_List]:=Block[{ndim,nnode,ve},ndim=Length[nodes[[1]]];nnode=Length[elements];
If[ndim==2,ve=Table[Area[Polygon[nodes[[elements[[i]]]]]],{i,nnode}],Print["Error, only 2D elements are allowed"]];ve];
Area2Nodes[nodes_List,elements_List]:=Block[{ve,ne,ei},ve=AreaOf2DElement[nodes,elements];ne=ConstantArray[0,Length[nodes]];Do[ei=elements[[i]];ne[[ei]]+=ve[[i]]/Length[ei],{i,Length[elements]}];ne];
Simplex2Nodes[nodes_List,elements_List]:=Block[{ve,ne,ei},ne=ConstantArray[0,Length[nodes]];Do[ei=elements[[i]];ne[[ei]]+=SimplexVolume[nodes[[ei]]]/Length[ei],{i,Length[elements]}];ne];
ShowMesh[m_,withLable_:False]:=If[withLable,Show[m["Wireframe"["MeshElementIDStyle"->Black]],m["Wireframe"["MeshElement"->"PointElements","MeshElementIDStyle"->Red]]],Show[m["Wireframe"]]];
My3DpointPlot[coord_,f_,size_:10]:=Block[{ri=Min[f],rx=Max[f],data},data={coord[[All,1]],coord[[All,2]],coord[[All,3]],f}\[Transpose];
Graphics3D[{AbsolutePointSize[size],{Blend[{Blue,Red},Rescale[Last[#],{ri,rx}]],Point[Most[#]]}&/@data}]];
