(* ::Package:: *)

(*
<< Combinatorica`*)


(* :Title: GraphPartition *)

(* :Author: Dylan M. Scott *)

(* :Summary:
	This package contains functions pertaining to spectral and
	isoperimetric graph partitioning algorithms.
	The Combinatorica package is required, as this package makes use
	of its representation of graphs.
*)

(* :Mathematica Version: 6.0
   :Combinatorica Version: 2.1 *)

BeginPackage["GraphPartition`"]

CriterionCut::usage =
	"CriterionCut[g,v] gives the component of v which gives a good
	balance between subgraph sizes and edges cut. It attempts to
	find a reasonably balanced cut (in the middle third of the sorted
	vector) that minimizes the objective function:
	Cut(g1,g2) ( (1/Mass(g1)) + (1/Mass(g2))"

DegreeMatrix::usage =
	"DegreeMatrix[g] returns the degree matrix of the graph g."

FiedlerVector::usage =
	"FiedlerVector[g] returns the Fiedler vector of the graph g.
	That is, the eigenvector of the Laplacian matrix of g which
	corresponds to the smallest non-trivial eigenvalue."

GeneralizedCriterionCut::usage =
	"GeneralizedCriterionCut[criterionFn,g,v] finds the component
	of v which produces a cut in g that minimizes the criterion
	function."

IsoperimetricComparison::usage =
	"IsoperimetricComparison[g] runs the isoperimetric partitioner
	twice as in TwiceIsoperimetricVector and prints the values of
	the criterion cut for each vector for comparative purposes."

IsoperimetricCut::usage =
	"IsoperimetricCut[g] partitions the graph into two graphs, or rather two
	disconnected components of a single Mathematica graph."

IsoperimetricVector::usage =
	"IsoperimetricVector[g,subs] returns the vector x which is equal to
	(k/( {1} M y )) y where y is the solution to the linear system
	L . y = M . {1}, and k is the desired mass of G1, in this case
	Mass[g] / 2. The list 'subs' should contain the indices of vertices
	to be grounded to G2 - at least one in each connected component of G
	will ensure that the system has a solution."

JumpCut::usage =
	"JumpCut[v] gives a threshold for cutting v, determined by sorting
	its components, finding the biggest 'jump' between consecutive
	values, and averaging those two values. This implementation examines
	only the middle third of the sorted vector in an attempt to get a
	more balanced cut."

LaplacianMatrix::usage =
	"LaplacianMatrix[g] returns the Laplacian matrix of the graph g."

MedianCut::usage =
	"MedianCut[v] gives the median element of v as the threshold for
	rounding."

PermuteMatrix::usage =
	"PermuteMatrix[m,p] permutes the given matrix using the given
	permutation. Permute matrix works only on square matrices,
	and permutes the rows and columns, consistent with the built-
	in Permute function which operates on lists."

PlotFiedlerVector::usage =
	"PlotFiedlerVector[g] displays a plot of the components of the
	Fiedler vector of g, in both unsorted and sorted order."

PlotIsoperimetricVector::usage =
	"PlotIsoperimetricVector[g] displays a plot of the components of the
	Isoperimetric vector of g, in both unsorted and sorted order."

RowColumnDrop::usage =
	"RowColumnDrop[m,n] returns the matrix m with it's nth row and
	column removed."

ShowIsoperimetricCut::usage =
	"ShowIsoperimetricCut[g] shows the graph g with the edges that would be
	cut by the Isoperimetric graph partitioning algorithm highlighted."

ShowMultipleFiedlerCut::usage =
	"ShowMultipleFielderCut[g] shows the graph g with the edges that would
	be cut by the Multiple Fielder graph partitioning algorithm highlighted."
ShowSpectralCut::usage =
	"ShowSpectralCut[g] shows the graph g with the edges that would be
	cut by the Spectral graph partitioning algorithm highlighted."

SpectralCut::usage =
	"SpectralCut[g] returns the result of partitioning the graph into
	two graphs, or rather two disconnected components of a single
	Mathematica graph."

VectorCut::usage =
	"VectorCut[g] gives a list of the edges in g that would be cut
	using the given vector and cutoff threshold."

MultipleFiedlerCut::usage =
	"MultipleFiedlerCut[g] gives the list of the edges in g that
	would be cut using the given vector and cutoff threshold"

VertexWeights::usage =
	"VertexWeights[g] gives a list of the weights of each vertex in the
	graph g, starting with the first vertex, then the second, etc."

LookupWeight::usage =
	"LookupWeight[g,e] gives the weight of the edge e in the graph g.
	If edge weights are not set in g, the uniform weight 1 is returned."

RoachGraph::usage =
	"RoachGraph[n] returns a graph on 2n vertcies corresponding to the
	special 'Roach graph' that demonstrates the failure of the spectral
	partitioner to find a reasonable cut."

MeshToEdgeList::usage =
	"MeshListToEdgeList[fList,vertices,weightFn] converts the mesh into
	a graph, where each face is a vertex in the graph and edges connect
	adjacent faces (two faces are adjacent if they share two or more
	vertices. Edges are weighted by the argument weight function which
	takes as arguments the coordinates of the vertices shared by the 
	two faces."

ObjToGraph::usage =
	"ObjToGraphics3D[s] converts the OBJ file specified by the
	argument string into a Mathematica graph"

PartitionObj::usage =
	"PartitionObj[s] partitions the surface mesh defined by
	the file specified by the argument string, and returns it is
	a Mathematica Graphics3D object with colored partitions."






(* ---------------------------------------------------------------
		GRAPH OPERATIONS
		Partitioning algorithms, cut selection, and the functions
		that help them.
*)
	(* ----
		Partitioners
			Implementation of the spectral and isoperimetric
			partitioners. *)

	(* SPECTRAL PARTITIONER *)

	SpectralCut[g_Graph] := EdgeDelete[g,SpectralEdges[g]](*ChangeEdges[g,Complement[Edges[g],SpectralEdges[g]]]*)

	Options[SpectralEdges]={IncludeObj->False}
	SpectralEdges[g_Graph, OptionsPattern[]] := With[{v=FiedlerVector[g]},
		With[{cutoff=CriterionCut[g,v]},
			If[OptionValue[IncludeObj],
				{CriterionPartitionFunction[g,v,cutoff], VectorCut[g,v,cutoff]},
			     VectorCut[g,v,cutoff]]]]
	
	(* :FiedlerVector:
		First solves for the eigensystem of the Laplacian matrix of the
		graph, formats to the form {{value1, vector1},{value2, vector2},...}
		then removes eigenvalue/vector pairs for which the eigenvalue is (very
		close to) 0, sorts the pairs by eigenvalue in ascending order, and selects
		the eigenvector corresponding to the smallest remaining eigenvalue. *)
	FiedlerVector[g_Graph] :=
	    Last[
			Sort[
				(*Select[*)
					Thread[Apply[List,Eigensystem[N[LaplacianMatrix[g]]]]],
				(*Abs[First[#]-0]>10^-10&],*)
			First[#1]<First[#2]&][[2]]]




	



	(* ISOPERIMETRIC PARTITIONER *)

	IsoperimetricCut[g_Graph] := EdgeDelete[g,IsoperimetricEdges[g]]

	(* :IsoperimetricEdges
		Note: Grounds a random vertex in each connected component, unless the
		RunTwice option is true, in which case it runs the partitioner twice,
		the second time grounding the vertices from the first run with the highest
		values. *)
	Options[IsoperimetricEdges]={RunTwice->False, IncludeObj->False, Substitutions->{}}
	IsoperimetricEdges[g_Graph,OptionsPattern[]] := With[{v=IsoperimetricVector[g,
				Substitutions->If[OptionValue[Substitutions]=={}, 
					Sort[Map[ #[[Random[Integer,{1,Length[#]}]]]&ConnectedComponents[g]]],
					OptionValue[Substitutions]],
				IsoperimetricMode->If[OptionValue[RunTwice],Twice,Fixed]]},
			With[{cutoff=CriterionCut[g,v]},
				If[OptionValue[IncludeObj],
					{CriterionPartitionFunction[g,v,cutoff], VectorCut[g,v,cutoff]},
					VectorCut[g,v,cutoff]]]]

	(* :IsoperimetricVector:
		To pass a vertex fixing mode to the partitioner, pass the FixingMode option
		in the form FixingMode->{Mode, f, n} where n is the vertex to be fixed on the
		initial run, and f specifies either how many vertices to fix on the second run,
		or how many runs to do.  By default vertex fixing will fix the top and bottom eighths
		based on an initial (random) run of the isoperimetric partitioner. *)
	Options[IsoperimetricVector]={IsoperimetricMode->Default,Substitutions->{},FirstRun->{},FixedMode->{}}
	IsoperimetricVector[g_Graph,OptionsPattern[]] := Module[ {matrix,vector,soln,subs},
		Switch[OptionValue[IsoperimetricMode],
			Default,
			If[OptionValue[Substitutions]=={},
				subs=Sort[Map[ #[[Random[Integer,{1,Length[#]}]]]&, ConnectedComponents[g]]];,
				subs=OptionValue[Substitutions];];
			matrix=RowColumnDrop[LaplacianMatrix[g],Map[List,subs]];
			vector=Delete[VertexWeights[g],Map[List,subs]];
			With[{k=Total[VertexWeights[g]],y=LinearSolve[N[matrix],vector]},
				soln= (k / (vector . y)) y];
			Map[ soln=Insert[soln,0,#]; &,subs];
			soln,

			Twice,
			With[{v1=If[OptionValue[FirstRun]=={},IsoperimetricVector[g],OptionValue[FirstRun]]},
			IsoperimetricVector[g,Substitutions->
			Map[First[First[Sort[#,Abs[Last[#1]]>Abs[Last[#2]]&]]]&,
				ConnectedComponents[g]/.n_Integer:>{n,v1[[n]]}]]],

			Fixed,
			Module[{v,order,fixed,result},
				If[OptionValue[FirstRun]=={},
					v=IsoperimetricVector[g];,
					v=OptionValue[FirstRun];];
				If[OptionValue[FixedMode]=={},
					fixed=Round[V[g]/8];
					order=Map[First,Sort[ MapIndexed[{#2[[1]],#1}&,v], #1[[2]]>#2[[2]]&]];
					result=FixedIsoperimetricVector[g,order[[Round[Length[order]/2]]],order[[1;;fixed]],order[[-fixed;;-1]]];
					result,
					With[{n = Last[OptionValue[FixedMode]],
						  f = First[Rest[OptionValue[FixedMode]]]},
						result = {};
						Switch[First[OptionValue[FixedMode]],
							TwoPass,
							v = FixedIsoperimetricVector[g,
								First[RandomSample[Complement[Range[V[g]],{Last[OptionValue[FixedMode]]}],1]], (* Random vertex that's not fixed *)
								{},
								{Last[OptionValue[FixedMode]]}];

							order=Map[First,Sort[ MapIndexed[{#2[[1]],#1}&,v], #1[[2]]>#2[[2]]&]];
							result={{0.0,FixedIsoperimetricVector[g,order[[Round[Length[order]/2]]],order[[1;;Ceiling[f/2]]],order[[-Floor[f/2];;-1]]];}},

							LogPasses,
							v = FixedIsoperimetricVector[g,
								First[RandomSample[Complement[Range[V[g]],{Last[OptionValue[FixedMode]]}],1]], (* Random vertex that's not fixed *)
								{},
								{Last[OptionValue[FixedMode]]}];
							fixed = 1;
							While[fixed < f,
								PrintIsoperimetricData[g,v];
								result = Join[result, {{N[Last[CriterionCut[g,v,IncludeCut->True]]],v}}];
								IsoperimetricData[g,v];
								order=Map[First,Sort[ MapIndexed[{#2[[1]],#1}&,v], #1[[2]]>#2[[2]]&]];
								v=FixedIsoperimetricVector[g,order[[Round[Length[order]/2]]],order[[1;;fixed]],order[[-fixed;;-1]]];
								fixed = 2 fixed;];,

							LinearPasses,
							v = FixedIsoperimetricVector[g,
								First[RandomSample[Complement[Range[V[g]],{Last[OptionValue[FixedMode]]}],1]], (* Random vertex that's not fixed *)
								{},
								{Last[OptionValue[FixedMode]]}];
							Module[{low,high},
								low=1;
								high=0;
								fixed=1;
								While[fixed <= f,
									PrintIsoperimetricData[g,v];
									result = Join[result, {{N[Last[CriterionCut[g,v,IncludeCut->True]]],v}}];
									fixed = fixed + 1;
									If[high!=low,high=high+1;,low=low+1;];
									order=Map[First,Sort[ MapIndexed[{#2[[1]],#1}&,v], #1[[2]]>#2[[2]]&]];
									v=FixedIsoperimetricVector[g,order[[Round[Length[order]/2]]],order[[1;;low]],order[[-high;;-1]]];];
							];
							Last[First[Sort[result,First[#1]<First[#2]&]]]
						]
					]
				]	
			]
		]
	]

	(* VERTEX FIXING *)

	(* :FixedIsoperimetricVector:
		Note: Chi is the desired mass of G1. *)
	FixedIsoperimetricVector[g_Graph,balance_Integer,ones_List,zeros_List] :=
	Module[{permutation,laplacianM,mass,chi,unfixed,
			Lnn,Mnn,Lqq,Mqq,M11,Lqn,Lq1,L1n,lhs,rhs,xhat,xn},
		(*Print[balance]; Print[ones];Print[zeros];*)
		
		permutation = ConstructPermutation[V[g],balance,ones,zeros];
		laplacianM = PermuteMatrix[LaplacianMatrix[g],permutation];
		mass = Permute[VertexWeights[g],permutation];
		chi = Total[mass] / 2;
		unfixed = V[g]-(1+Length[ones]+Length[zeros]);
		Lnn = Last[Last[laplacianM]];
		Mnn = Last[mass];
		Lqq = Take[Map[Take[#,unfixed]&,laplacianM],unfixed];
		Mqq = DiagonalMatrix[Take[mass,unfixed]];
		M11 = If[Length[ones]>0,DiagonalMatrix[Take[mass,{unfixed+1,unfixed+Length[ones]}]],{{}}];
		Lqn = Take[Map[{Last[#]}&,laplacianM],unfixed];
		Lq1 = Take[Map[Take[#,{unfixed+1,unfixed+Length[ones]}]&,laplacianM],unfixed];
		L1n = Take[Map[{Last[#]}&,laplacianM],{unfixed+1,unfixed+Length[ones]}];
		
		lhs = 2 (Lqq - (1 / Mnn) (Map[PadLeft[#,Length[Lqn],First[#]]&,Lqn]) . Mqq
				- (1 / Mnn) Mqq . Table[First[Transpose[Lqn]],{Length[Lqn]}]
				+ (Lnn / Mnn^2) (Mqq . Table[1,{Length[Mqq]},{Length[Mqq]}] . Mqq));
		rhs = With[{chiexpr=If[Length[ones]>0,(chi - First[First[Transpose[Table[{1},{Length[M11]}]] . M11 . Table[{1},{Length[M11]}]]]),chi]},
				(((2 Lnn chiexpr / Mnn^2) + If[Length[ones]>0,(2 First[First[Transpose[L1n] . Table[{1},{Length[L1n]}]]] / Mnn),0]) Mqq . Table[{1},{Length[Mqq]}])
				- ((2 chiexpr / Mnn) Lqn)
				- If[Length[ones]>0,(2 Lq1 . Table[{1},{Length[First[Lq1]]}]),0] ];

		xhat = Join[Map[First,LinearSolve[N[lhs],rhs]],Table[1,{Length[ones]}],Table[0,{Length[zeros]}]];
		xn = First[(chi - Most[mass] . Transpose[{xhat}]) / Mnn];
		Permute[Append[xhat,xn],
			InversePermutation[permutation]]
	]
		

	PermuteMatrix[m_,p_List] := Permute[Map[Permute[#,p]&,m],p]

	(* :ConstructPermutation:
		Permutes the given matrix in such a way that, going down the
		diagonal, vertices will appear in the order:
		unfixed - fixed to 1 - fixed to 0 - balance constrained vertex *)
	ConstructPermutation[size_Integer,balance_Integer,ones_List,zeros_List] :=
		With[{end=Join[ones, zeros, {balance}]},
			Join[ Complement[Range[size], end], end] ]



	(* MULTIPLE FIEDLER VECTORS *)
	Options[MultipleFiedlerEdges]={IncludeObj->False}
	MultipleFiedlerEdges[g_Graph, OptionsPattern[]] := With[
		{v=MultipleFiedlerVector[g], cutoff=0.5},
		If[OptionValue[IncludeObj],
			{CriterionPartitionFunction[g,v,cutoff], VectorCut[g, v, cutoff]},
			VectorCut[g, v, 0.5]]]
	
	MultipleFiedlerCut[g_Graph] := EdgeDelete[g, MultipleFiedlerEdges[g]]
	MultipleIsoperimetricCut[g_Graph] := EdgeDelete[g, MultipleIsoperimetricEdges[g]]

	EvaluateMultipleFiedlerEdges[g_Graph] := With[{v=MultipleFiedlerVector[g], cutoff=0.5},
		With[{tCutedges = Timing[VectorCut[g, v, cutoff]]},
		{First[tCutedges], CriterionPartitionFunction[g, v, cutoff], tCutedges,v}]]		

	(* :MultipleFiedlerVector:
		Gives an indicator vector generated by using the vector partitioning algorithm
		on the spectral vectors generated for the graph.
		By default cuts the graph into two partitions of equal number of vertices. *)
	Options[MultipleFiedlerVector]={(*Range->{1/5,4/5}, *)NumBasis->5, Basis->{}}
	MultipleFiedlerVector[g_Graph, OptionsPattern[]] := With[
		{basisVectors= If[OptionValue[Basis]=={}, 
		    SpectralVectors[g,OptionValue[NumBasis]], OptionValue[Basis]]},
		VectorPartition[basisVectors, V[g], OptimizationOn->True,
			OptimizationObj-> (CriterionPartitionFunction[g,#,0.5]&)]]


	SpectralVectors[g_Graph, k_] := With[{eigensystem = Eigensystem[N[LaplacianMatrix[g]],-k]},
	Module[{n,eigenvalues,eigenvectors,xi},
		n = V[g];
		eigenvalues = Reverse[First[eigensystem]];
		Print[eigenvalues];
		eigenvectors = Table[Normalize[v],{v,Reverse[Last[eigensystem]]}];
		(* Using Alpert et al.'s recommendation of xi = 2nd eigenval + nth eigenval *)
		xi = eigenvalues[[2]] + eigenvalues[[-1]];
		Table[Table[N[Sqrt[xi - eigenvalues[[i]]]] * eigenvectors[[i]][[j]],{i,k}],{j,n}]
	]]


	(* MULTIPLE ISOPERIMETRIC VECTORS *)
	Options[MultipleIsoperimetricEdges]={IncludeObj->False, NumBasis->5,
	                                     GroundingScheme-> GroundingRandom,ShowGrounds->False}
	MultipleIsoperimetricEdges[g_Graph, OptionsPattern[]] := Module[
		{grounds, v, cutoff=0.5, result},
		 {grounds, v}=MultipleIsoperimetricVector[g, GroundingScheme->OptionValue[GroundingScheme],
		                                              ShowGrounds->True,
		                                              NumBasis->OptionValue[NumBasis]];
		 result = {If[OptionValue[IncludeObj],CriterionPartitionFunction[g,v,cutoff],Nothing],
		           VectorCut[g, v, cutoff],
		           If[OptionValue[ShowGrounds], grounds, Nothing]};
		 If[Length[result]==1, result[[1]], result]]

	Options[MultipleIsoperimetricVector]={NumBasis->5, GroundingScheme-> GroundingRandom, ShowGrounds->False}
	MultipleIsoperimetricVector[g_Graph, OptionsPattern[]] := 
	(*With[{xi=Eigenvalues[N[LaplacianMatrix[g]]][[1]]},*)
	Module[{grounds, basis,v}, 
	       {grounds, basis} = FingerprintVectors[g,OptionValue[NumBasis],OptionValue[GroundingScheme]];
	       v= MultipleFiedlerVector[g, NumBasis->OptionValue[NumBasis], Basis->basis];
		   If[OptionValue[ShowGrounds], {grounds, v}, v]]
	(*Maximum[
		Table[VectorPartition[SpectralVectors[g], mu],
		{mu, Floor[V[g]*OptionValue[Range][[1]]], Floor[V[g]*OptionValue[Range][[2]]]}],
		-CriterionPartitionFunction[g,#,0.5]&
		][[1]]*)

	(* :SpectralVectors:
		Using the k smallest eigenvalues and corresponding vectors,
		Gives the n-element spectral vectors of the graph, where the jth spectral vector:
			s'_j = [ sqrt[xi-lam1]*v_j1, sqrt[xi-lam2]*v_j2, ..., sqrt[xi-lamN]*v_jN ]
		Where lamN is the nth eigenvalue of the graph's Laplacian and v_ji is the jth element
		of the ith eigenvector. *)




	(* :FingerprintVectors:
		Using the k randomly grounded isoperimetric vectors,
		Gives the n-element fingerprint vectors of the graph, where the jth fingerprint vector:
			t'_j = the jth column of J^(1/2)B^TM  is a real valued vector,
			J = B^T(xi*M - L)B (each column of B is a basis vector)
		Where L is the graph's Laplacian, xi is an arbitrary constant, 
		M is the mass matrix of the vertices (For now all vertices have the same weight 1)
		J^(1/2) is the upper triangular square root of J computed by Cholesky decomposition. *)
	(*Next: Random without replacement*)
	Options[FingerprintVectors]={Xi->{}}
	FingerprintVectors[g_Graph, k_Integer, scheme_, OptionsPattern[]] := Module[
		 {grounds, B, M, L, xi, J},
		 {grounds, B} = IsoperimetricBasis[g, k, scheme];
		 M = IdentityMatrix[V[g]];
		 L=N[LaplacianMatrix[g]]; (*epsilonI = OptionValue[Epsilon]*IdentityMatrix[k]},*)
		 xi = If[OptionValue[Xi]=={},Eigenvalues[N[LaplacianMatrix[g]]][[1]], OptionValue[Xi]];
		 J = Transpose[B].(xi*M-L).B;
		 {grounds,Transpose[CholeskyDecomposition[RoundToSym[Re[J]]].Transpose[B].M]}]

	RoundToSym[M_] := (M+Transpose[M])/2
				(*Transpose so that each row of the output is a fingerprint vector*)
	(*Options[IsoperimetricBasis]={Scheme\[Rule]GroundingRandom}*)
	IsoperimetricBasis[g_Graph, k_Integer, (*OptionsPattern[]*)scheme_] := 
	Module[{groundVertices={}, vectors={}, next},
		If[SameQ[scheme,GroundingRandom]|| SameQ[scheme,GroundingPreviousMax]
		||SameQ[scheme,GroundingPreviousMedian],
		    For[i=1,i<=k, i+=1,
            next=scheme[g,groundVertices, If[Length[vectors]>0,vectors[[-1]], {}]];
            groundVertices=Append[groundVertices,next[[1]]];
            (*Print[groundVertices];*)
            vectors=Append[vectors,next[[2]]]];,
		{groundVertices, vectors}=scheme[g,k]];
        (* Print["vectors", vectors];*)
	    {groundVertices,Transpose[Orthogonalize[vectors, Method->"GramSchmidt"]]}]
(*	Module[{groundVertices, vectors, v, prevMax},
	vectors=Switch[OptionValue[GroundingScheme],
			Sub, IsoperimetricVector[g, Substitutions->{#}]&/@OptionValue[Grounds],
			Rand, IsoperimetricVector[g, Substitutions->{#}]&/@RandomSample[Range[V[g]],k],
			PrevMax, prevMax = RandomInteger[V[g]];
			         Table[v=IsoperimetricVector[g, Substitutions->{prevMax};
			               prevMax=Maximum[v,Identity][[2]]; v],k]];
	Print[vectors];*)
	
	Options[GroundingSubset]={Ratio->3}
	GroundingSubset[g_Graph,k_Integer,OptionsPattern[]] := 
		Module[{candidates,groundAndVectors},
			(*Generate Ratio*k candidate isoperimetric vectors*)
			candidates = RandomSample[Range[V[g]],Min[V[g],k*OptionValue[Ratio]]];
			groundAndVectors = Map[{#, IsoperimetricVector[g, Substitutions->{#}]}&, candidates];
			(*Sort by the objective (sparsity) ascendingly*)
			groundAndVectors = SortBy[groundAndVectors, CriterionPartitionFunction[g,#[[1]],cutoff]&];
			(*Choose the best k vectors from Ratio*k candidates*)
			groundAndVectors = Take[groundAndVectors, k];
			(*Return the best k grounds and the corresponding vectors*)
			{Map[#[[1]]&,groundAndVectors], Map[#[[2]]&,groundAndVectors]}]
		
	GroundingMaxDegree[g_Graph,k_Integer] := 
		Module[{degrees,grounds,vectors},
			degrees= MapIndexed[{#1, #2[[1]]}&,DegreeList[g]];
			degrees=Take[SortBy[degrees, -#[[1]]&],k];
			grounds=Map[#[[2]]&,degrees];
			vectors=Map[IsoperimetricVector[g, Substitutions->{#}]&, grounds];
			{grounds, vectors}]
			
	GroundingRandom[g_Graph, prevGrounds_List,prevV_List] := 
		With[{ground=RandomChoice[Complement[Range[V[g]],prevGrounds]]},
		    {ground, IsoperimetricVector[g, Substitutions->{ground}]}]

	GroundingPreviousMax[g_Graph, prevGrounds_List,prevV_List] := 
		Module[{candidates, labeledV, ground},
			ground=If[prevV=={},
		                   RandomInteger[{1,V[g]}], 
		                   (*label the vector with indexes and only keep 
		                   the indexes that didn't appear in previous grounds*)
		                   labeledV = MapIndexed[{#1,#2[[1]]}&,prevV];
		                   labeledV = Select[labeledV, !MemberQ[prevGrounds, #[[2]]]&];
		                   (*[[1,2]] because Maximum here gives {{value, oldindex}, newindex}
		                   we want the old index before Select. *)
		                   Maximum[labeledV,#[[1]]&][[1,2]]];
		    (*Print[prevV];
		    Print[ground];*)
		    {ground, IsoperimetricVector[g, Substitutions->{ground}]}]
	GroundingPreviousMedian[g_Graph, prevGrounds_List,prevV_List] := 
		Module[{candidates, labeledV, ground},
			ground=If[prevV=={},
		                   RandomInteger[{1,V[g]}], 
		                   (*Only keep the indexes that didn't appear in previous grounds*)
		                   labeledV = MapIndexed[{#1,#2[[1]]}&,prevV];
		                   labeledV = Select[labeledV, !MemberQ[prevGrounds, #[[2]]]&];
		                   labeledV = SortBy[labeledV, #[[1]]&];
		                   labeledV[[ Quotient[Length[prevV],2] ]][[2]]
		                   ];
		                   (*medValue = Median[
		                       Table[prevV[[i]], {i, Complement[Range[V[g]], prevGrounds]}]]
		                   medPositions = Flatten[Position[prevV,medValue]]
		                   Complement[medPositions, prevGrounds][[1]]]*)
		    (*Print[prevV];
		    Print[ground];*)
		    {ground, IsoperimetricVector[g, Substitutions->{ground}]}]



			
	Options[VectorMaxOnDirection]={OptimizationOn->False, OptimizationObj->Identity}
	VectorMaxOnDirection[vectors_List, eta_Integer, dir_List, OptionsPattern[]] := Module[
		{partition, sortedVertices, objectives, bestEta, bestObj=Infinity},
	    labelledVectors = MapIndexed[ {#1, #2[[1]]}&, vectors];
		(*Sort the basis vectors descendingly according to the projected lengths of the
		vectors on the given direction*)
		labelledVectors = SortBy[labelledVectors, -Dot[#[[1]],Normalize[dir]]&];
		sortedVertices = Map[#[[2]]&,labelledVectors];
		(*Print[sortedVertices];*)
		bestEta = Min[eta, Length[vectors]];
		If[OptionValue[OptimizationOn],
			(*Try putting at most eta vertices in partition 1,
			choose the one with smallest sparsity*)
		    partition = Table[0,V[g]];
(*		    For[i=1, i\[LessEqual]Length[sortedVertices],i++,
		        partition[[sortedVertices[[i]]]]=1;
				obj=OptionValue[OptimizationObj][partition];
				If[obj<bestObj, bestObj=obj; bestEta=i, Nothing]];*)
			objectives=Table[partition[[sortedVertices[[i]]]]=1;
				OptionValue[OptimizationObj][partition],{i,1,bestEta}];
			(*Print[N[objectives]]*);
			(*use negative objectives to find the minimum *)
		    {bestObj, bestEta} = Maximum[-objectives,Identity];
		    (*Print[bestEta];*)];
		Table[If[MemberQ[Take[sortedVertices, bestEta],i], 1, 0], {i,Length[vectors]}]]

		
	(* :VectorPartition:
		Performs a basic vector partitioning on the given list of vectors.  It attempts to find
		eta vectors from the list whose sum is maximized.  It does this by placing the largest
		vector in one partition, then adds vectors one by one which maximize the sum so far, until
		it has eta vectors. Returns an indicator vector of the same length as the input list, where
		a 1 indicates the vector is part of the new partition, and a 0 if it is not. *)
	Options[VectorPartition]={OptimizationOn->False, OptimizationObj->Identity, Method->Melo}
	VectorPartition[vectors_List,eta_Integer,OptionsPattern[]] := Module[
	    {n,partition,labelledVectors,max,count, meloResult},
		labelledVectors = MapIndexed[ {#1, #2[[1]]}&, vectors];
		n = 0;
		count = 0;
		partition = {};
		bestObj = Infinity;
		bestPartitionV = {};
		While[count < eta,
			count++;
			max = Maximum[labelledVectors, Norm[n + #[[1]]]&];
			(*Print[max]*)
			n = n + max[[1]][[1]];
			partition = Join[partition, {max[[1]][[2]]}];
			(*Print[partition];*)
			If[OptionValue[OptimizationOn],
				With[{partitionV=Table[ If[MemberQ[partition,i], 1, 0], {i, Length[vectors]}]},
					(*Print[partitionV];*)
					With[{obj=OptionValue[OptimizationObj][partitionV]},
						(*Print[obj];*)
						If[obj<=bestObj,
							bestObj=obj;bestPartitionV=partitionV;,
							Nothing]]],
				Nothing];
		    
			labelledVectors = Drop[labelledVectors, {max[[2]]}];
		];
(*		Table[ If[MemberQ[bestPartition,i], 1, 0], {i, Length[vectors]} ] ]*)
		meloResult=If[OptionValue[OptimizationOn],
					bestPartitionV,
					Table[ If[MemberQ[partition,i], 1, 0], {i, Length[vectors]} ] ];
		Switch[OptionValue[Method],
			Melo, meloResult,
			MeloProbe, 
				With[{dir=Normalize[Total[basisVectors*meloResult]]},
					VectorMaxOnDirection[vectors,eta,dir,OptimizationOn->OptionValue[OptimizationOn],
					OptimizationObj->OptionValue[OptimizationObj]]]]		
		]
	
	


	(* VARIOUS UTILITY FUNCTIONS *)

	(* :Maximum:
		Returns the value in the list for which the given function is maximized, as well as the
		index of that element. *)
	Maximum[l_List, p_] := Module[{image,max,maxIndex,n},
		image = p /@ l;
		maxIndex = 1;
		max = image[[1]];
		n = 2;
		While[n <= Length[l],
			If[ image[[n]] > max,
				max = image[[n]];
				maxIndex = n; ];
			n++;
		];	
		{l[[maxIndex]], maxIndex} ]
	

	DegreeMatrix[g_Graph] := System`DiagonalMatrix[DegreeList[g]]
	(* :DegreeList:
		Note: The weird map is because the ToAdjacencyMatrix procedure with
		edge weights has non-existing edges entered as Infinity. This replaces
		the Infinity entries with zero. *)
	DegreeList[g_Graph] := Apply[Plus, 
					Map[If[#==Infinity,0,#]&,WeightedAdjacencyMatrix[g], {2}]]
		
	LaplacianMatrix[g_Graph] := (DegreeMatrix[g] - 
				Map[If[#==Infinity,0,#]&,WeightedAdjacencyMatrix[g], {2}])

	RowColumnDrop[m_, l_] := Map[Delete[#,l]&,Delete[m,l]]

	Options[VectorCut]={IncludeCut->False}
	VectorCut[g_Graph,v_List,cutoff_] := Select[EdgeList[g],
		(v[[First[#]]] >= cutoff && v[[Last[#]]] < cutoff) ||
		(v[[First[#]]] < cutoff && v[[Last[#]]] >= cutoff) &]

	V[g_Graph] := Length[VertexList[g]]
	VertexWeights[g_Graph] := Map[PropertyValue[{g,#},System`VertexWeight]&,VertexList[g]]
(*	EdgeWeights[g_Graph] := Map[PropertyValue[{g,#},System`EdgeWeight]&,EdgeList[g]]*)

(*Map[If[Length[#]==0,1,Last[First[#]]]&,
		Table[GraphOptions[g,i],{i,V[g]}]]*)

	(* :LookupWeight:
		Has to take into account the fact that looking up an edge in a graph
		whose edge weights have not been set returns an empty list. *)
	LookupWeight[g_Graph, e_] := PropertyValue[{g,e},System`EdgeWeight]
(*	First[GetEdgeWeights[g, {e}]]*)

(*With[{edge=GraphOptions[g,e]},
		If[Length[edge]==0,1,Last[First[edge]]]]*)

	(* ---- *)
	(* ----
		Cuts
			Functions pertaining to the selection of a threshold to round
			the vector resulting from the spectral or isoperimetric
			partitioners. *)

	(* CRITERION CUT
		Tries many possible thresholds in order to find one a cut
		over which minimizes some objective function. *)
		(*:Option: IncludeCut is True if the objective function value coressponding to
		the optimal cut is desired*)
	Options[GeneralizedCriterionCut]={IncludeCut->False}
	GeneralizedCriterionCut[criterionFn_,g_Graph,vector_List,cutoffs_List,OptionsPattern[]] :=
		With[{cut=First[Sort[Map[{#,criterionFn[g,vector,#]}&,cutoffs],Last[#1] < Last[#2]&]]},
			If[OptionValue[IncludeCut],
				cut,
				First[cut]]]

	Options[CriterionCut]={IncludeCut->False,Range->{0,1}}
	CriterionCut[g_Graph,vector_List,OptionsPattern[]] :=
		With[{range=Take[Sort[vector],{Max[{1,Round[Length[vector] * First[OptionValue[Range]]]}],Round[Length[vector] * Last[OptionValue[Range]]]}]},
		GeneralizedCriterionCut[CriterionPartitionFunction,g,vector,range,IncludeCut->OptionValue[IncludeCut]]]

	CriterionPartitionFunction[g_Graph,v_List,c_] := With[{partitionWeights=PartitionWeights[GraphPartition`VertexWeights[g],v,c]},
		If[First[partitionWeights]==0 || Last[partitionWeights]==0, Infinity,
		(Total[Map[LookupWeight[g,#]&,VectorCut[g,v,c]]]) * ((1/(First[partitionWeights])) + (1/(Last[partitionWeights])))]]

	PartitionWeights[vertexWeights_List,v_List,cutoff_] := 
		{  Total[ Select[ Table[i,{i,1,Length[v]}], v[[#]]>=cutoff&]/.n_Integer:>vertexWeights[[n]] ],
		 Total[ Select[ Table[i,{i,1,Length[v]}], v[[#]]<cutoff&]/.n_Integer:>vertexWeights[[n]] ]}

	(* JUMP CUT
		Finds the largest 'jump' (difference) between two values in
		the sorted vector, and uses the average of these as the
		threshold. *)
	JumpCut[v_List] := With[{midv=Take[Sort[v],{Round[Length[v]/3],Round[2 Length[v]/3]}]},
		Mean[First[Sort[Partition[Rest[Drop[Sort[Riffle[midv,midv]],-1]],2],
							Abs[First[#1]-Last[#1]]>Abs[First[#2]-Last[#2]]&]]]]

	(* MEDIAN CUT
		Returns the median of the vector. Guarantees an equal number
		of vertices in each partition, but the quality of the cut is
		generally questionable. *)
	MedianCut[v_List] := Median[v]
	(* ---- *)
(* --------------------------------------------------------------- *)



(* ---------------------------------------------------------------
		ANALYSIS AND VISUALIZATION
		Functions to help analyze the partitions produced by the above
		algorithms, as well as visualize partitions in terms of the
		vectors used, the graph being partitioned, or the partition of
		a 3D mesh.
*)
	(* ---- 
		Special graphs
			Graphs which highlight behaviors or shortcomings of the partitioners,
			or simply to test the partitioners or cuts on. *)

	(* ROACH GRAPH
		A graph whose properties highlight a bad partition on the part of
		the spectral partitioner. *)
	RoachGraph[n_Integer] := With[{vList=Apply[Join,Table[Table[{i,j},{i,0,n-1,1}],{j,0,2,2}]],
			eList=Join[Table[{i,i+1},{i,n-1}],Table[{i,i+1},{i,n+1,2n-1}],Table[{i,i+n},{i,Round[n/2]+1,n}]]},
		ChangeEdges[ChangeVertices[EmptyGraph[2n],vList],eList]]
	(* ---- *)

	(* ----
		Cut analysis
			Functions that plot the vectors produced by the partitioners,
			or help compare cuts. *)

	PlotFiedlerVector[g_Graph] := With[{v=FiedlerVector[g]},
		GraphicsColumn[{ListPlot[v], ListPlot[Sort[v,Less]]}]]

	PlotIsoperimetricVector[g_Graph] := With[{v=IsoperimetricVector[g]},
		Print[v];
		GraphicsColumn[{ListPlot[v], ListPlot[Sort[v,Less]]}]]

	IsoperimetricComparison[g_Graph] := Module[{v1,v2},
		v1=IsoperimetricVector[g];
		v2=IsoperimetricVector[g,IsoperimetricMode->Fixed,FirstRun->v1];
		Print["Criterion function value for normal run ", N[Last[CriterionCut[g,v1,IncludeCut->True]]]];
		Print["Criterion function value for fixed run  ", N[Last[CriterionCut[g,v2,IncludeCut->True]]]];
		GraphicsRow[{ListPlot[Sort[v1],PlotLabel->"Normal"],ListPlot[Sort[v2],PlotLabel->"Fixed"]}]
	]
	
	PrintObjectiveFunctionData[g_Graph, v_] := Module[{},
		Print["Criterion function value for run ", N[Last[CriterionCut[g,v,IncludeCut->True]]]];
		Print[ListPlot[Sort[v]]];]

	PrintIsoperimetricData[g_Graph, v_] := Module[{},
		Print["Criterion function value for run ", N[Last[CriterionCut[g,v,IncludeCut->True]]]];
		Print[ListPlot[Sort[v]]];]
		
	EvaluatePartitionAlg[alg_,g_] := Switch[alg,
		SpectralAlg,
		Flatten[Timing[SpectralEdges[g, IncludeObj->True]],1],
		IsoperimetricAlg,
		Flatten[Timing[IsoperimetricEdges[g, IncludeObj->True,RunTwice->False]],1],
		MultipleFiedlerAlg,
		Flatten[Timing[MultipleFiedlerEdges[g, IncludeObj->True]],1],
		MultipleIsoperimetricAlg,
		Flatten[Timing[MultipleIsoperimetricEdges[g, IncludeObj->True]],1]
		]
	(*TODO: need to rename to sth related to multiple to avoid confusion*)
	Options[CompareIsoperimetricAlgs]={ShowAll->False}
	CompareIsoperimetricAlgs[g_,k_,scheme_, OptionsPattern[]] := Module[
	   {timing, multipleResult, grounds, singleMetrics, multipleMetrics, toGrid, grid1, grid2},
	     {timing, multipleResult} = Timing[MultipleIsoperimetricEdges[g,GroundingScheme->scheme,
	      NumBasis->k, IncludeObj->True,ShowGrounds->True]];
	     grounds = multipleResult[[-1]];
	     multipleMetrics = Prepend[multipleResult[[;;-2]], timing];
	     singleMetrics = Map[(Flatten[Timing[
	        IsoperimetricEdges[g,RunTwice->False,Substitutions->{#},IncludeObj->True]],1])&,grounds];
		 (*m1, m2, m3, toGrid, grid1, grid2},
		 m1 = singleEval/@grounds;
		 m2 = multipleEval[grounds];*)
		 toGrid = ({ShowCuttedEdges[g,#[[3]]],N[#[[2]]]}&);
		 grid1 = toGrid/@singleMetrics;
		 grid2 = toGrid@multipleMetrics;
		 (*Print["grounds:",grounds];*)
		 With[{all=Table[Join[Prepend[grid1[[i]],grounds[[i]]], 
			     If[i==1,grid2,{SpanFromAbove, SpanFromAbove}]],
			     {i,1, Length[grid1]}],
			     best={Join[First[Sort[grid1,#1[[2]]<#2[[2]]&]],grid2]}},
			   If[OptionValue[ShowAll],all, best]]]
		 
	

	(* ---- *)







	(* ----
		Graph display
			Functions to display a graph and the partition - in this case
			by highlighting the edges that would be cut by the partitioner. *)

	ShowSpectralCut[g_Graph] := ShowCuttedEdges[g,SpectralEdges[g]]
	(*ShowGraph[Highlight[g,{SpectralEdges[g]},
			HighlightedEdgeColors->{Red},
			HighlightedEdgeStyle->Dashed]]*)

	ShowIsoperimetricCut[g_Graph] := ShowCuttedEdges[g, IsoperimetricEdges[g]]
	(*ShowGraph[Highlight[g,{IsoperimetricEdges[g]},
			HighlightedEdgeColors->{Red},
			HighlightedEdgeStyle->Dashed]]*)

	ShowMultipleFiedlerCut[g_Graph] := ShowCuttedEdges[g, MultipleFiedlerEdges[g]]
	(*ShowGraph[Highlight[g,{MultipleFiedlerEdges[g]},
			HighlightedEdgeColors->{Red},
			HighlightedEdgeStyle->Dashed]]*)

	ShowMultipleIsoperimetricCut[g_Graph] := ShowCuttedEdges[g, MultipleIsoperimetricEdges[g]]
			
	ShowCuttedEdges[g_Graph, e_List] := HighlightGraph[g,e](*ShowGraph[Highlight[g,{e},
			HighlightedEdgeColors->{Red},
			HighlightedEdgeStyle->Dashed],
			VertexLabel->True]*)
	
	ShowVectorCut[g_Graph,v_List,c_] := ShowCuttedEdges[g, VectorCut[g,v,c]]

	EigenShowGraph[g_Graph] := With[{eigensystem=Last/@Sort[Select[Thread[Apply[List,Eigensystem[N[LaplacianMatrix[g]]]]],
										Abs[First[#]-0]>10^-10&], First[#1]<First[#2]&]},
			ShowGraph[ChangeVertices[g,Range[V[g]]/.n_Integer:>{eigensystem[[2]][[n]],eigensystem[[3]][[n]]}]]]
	(* ---- *)

	(* ----
		Mesh visualization
			Converts a surface mesh into a graph where each face of the mesh
			is a vertex which is connected by an edge to the faces it shares
			an edge with.  The partition of this graph (division of faces
			into two partitions) demonstrated by showing the original mesh
			with meshes of the different partitions displayed in different
			colors. *)

	FacesToGraph[fList_List,vertices_List] := With[{elist =MeshToEdgeList[fList,vertices,1&]},
		Graph[Range[Length[fList]],Map[UndirectedEdge[#[[1]][[1]],#[[1]][[2]]]&, elist], 
			(*VertexLabels\[Rule]Automatic,*)VertexWeight->Table[1/Length[fList], Length[fList]],
			EdgeWeight->Automatic](*Map[#[[2]][[2]]&, elist]*)](*AddEdges[EmptyGraph[Length[fList]],
		MeshToEdgeList[fList,vertices,1&]]*)
	FacesToFrame[fList_List, vertices_List] := With[{elist=MeshToFrameEdgeList[fList, vertices]},
	Graph[Range[Length[vertices]],elist, 
			VertexWeight->Table[1/Length[vertices], Length[vertices]],
			EdgeWeight->Automatic]]

     MeshToFrameEdgeList[fList_List,vertices_List] := 
     Flatten[Table[Map[UndirectedEdge[f[[#]],f[[Mod[# , Length[f]]+1]]]&,Range[Length[f]]],{f, fList}]]
	(* :FacesToWeightedGraph:
		Weights the edges of the graph based on the length of the edge
		the faces on the mesh share.  Longer edges have a larger weight. *)
	FacesToWeightedGraph[fList_List,vertices_List] := AddEdges[EmptyGraph[Length[fList]],
		MeshToEdgeList[fList,vertices,EuclideanDistance]]

	(* :FaceListToEdgeList:
		Converts the list of faces and vertices to an edge list, weighting
		edges based on the weighting function, which must accept as arguments
		two length 3 lists, containing the coordinates of the vertices
		shared by the two faces being connected by an edge. *)
	MeshToEdgeList[fList_List,vertices_List,weightFn_Function] := Module[{labeledFList,edgeList},
		labeledFList=MapIndexed[{#1,First[#2]}&,fList];
		edgeList={};
		While[Length[labeledFList]>0,
				With[{face=First[labeledFList]},
					labeledFList=Rest[labeledFList];
					edgeList=Join[edgeList,Apply[Join,Map[ If[ Length[Intersection[First[#],First[face]]]>=2,
						{{{Last[#],Last[face]},
							EdgeWeight->Apply[weightFn,Intersection[First[#],First[face]]/.n_Integer:>vertices[[n]]]}},
								{}]&,labeledFList]]];]];
		edgeList]

	(* :PartitionObj
		NOTE: OBJ parsing code written by Mark McClur, found at:
		http://facstaff.unca.edu/mcmcclur/java/LiveMathematica/real.html *)
	PartitionObj[s_String] := Module[
		{data,vertices,faces,vector,cutoff},

		data=ReadList[s,Word,RecordLists->True];
		If[data=!=$Failed,
			vertices=Map[ToExpression,Drop[#,1]&/@
				Select[data,First[#]=="v"&],{2}];
			faces=Drop[#,1]&/@
				Select[data,First[#]=="f"&];
			If[vertices=!={}&&faces=!={},
				faces=If[DigitQ[faces[[1]][[1]]],
					Map[ToExpression,faces,{2}],
					Map[ToExpression[
						StringTake[#,StringPosition[#,"/"][[1]][[1]]-1]]&,
							faces,{2}]];
				With[{graph=FacesToGraph[faces,vertices]},
					vector=IsoperimetricVector[graph,IsoperimetricMode->Fixed];
					With[{cut=CriterionCut[graph,vector,IncludeCut->True]},
						cutoff=First[cut];
						(*Print["Weight: ", N[Last[cut]]];*)
					];
				];
				Show[Graphics3D[{Blue, Polygon/@(Apply[Join,MapIndexed[If[vector[[First[#2]]]>=cutoff,{#1},{}]&,faces]])/.n_Integer:>vertices[[n]],
								 Red, Polygon/@(Apply[Join,MapIndexed[If[vector[[First[#2]]]<cutoff,{#1},{}]&,faces]])/.n_Integer:>vertices[[n]]}],
					 Boxed->False],
				Print[s<>" does not seem to be an Obj file."]
			],
			Print[s<>" not found."]
		]
	]

	ObjToGraph[s_String] := Module[
		{data,faces,vertices},
		data=ReadList[s,Word,RecordLists->True];
		vertices=Map[ToExpression,Drop[#,1]&/@
				Select[data,First[#]=="v"&],{2}];
		If[data=!=$Failed,
			faces=Drop[#,1]&/@
				Select[data,First[#]=="f"&];
			If[vertices=!={}&&faces=!={},
				faces=If[DigitQ[faces[[1]][[1]]],
					Map[ToExpression,faces,{2}],
					Map[ToExpression[
						StringTake[#,StringPosition[#,"/"][[1]][[1]]-1]]&,
							faces,{2}]];
				FacesToGraph[faces,vertices],
				Print[s<>" does not seem to be an Obj file."]
			],
			Print[s<>" not found."]
		]
	]
	
(* --------------------------------------------------------------- *)
(*End[ ]
*)
EndPackage[ ]



(*<< Combinatorica`*)
(*$Context
ShowGraph[RoachGraph[4]]*)






