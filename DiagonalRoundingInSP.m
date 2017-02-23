(* ::Package:: *)

(* :Title: Diagonal Rounding in Spectral Partitioning *)

(* :Author: Maurya Talisetti *)

(* :Summary:
	This package contains functions pertaining to "Diagonal Rounding in Spectral Partitioning".
	The Combinatorica package is required, as this package makes use of its representation of graphs.

	:Usage:
	In the Mathematica notebook, making sure that this file and the .obj file are in the default folder ("Documents" in most cases),
	<<Combinatorica`
	<<DiagonalRoundingInSP.m
	diagonalRoundingObj["weird_ball.obj"] OR (NO LONGER SUPPORTED)
	diagonalRoundingRandom[10, fiedlerOrder->True] OR
	diagonalRoundingRandom[10, plot->True]

	:Key to variables used:
	Several local variables used follow the notations in Section D.1 - Diagonal Rounding In Spectral Partitioning.
	dhat, mhat - are D, M  with their nth rows and columns removed
	mnn, dnn   - are elements in M[n][n] and D[n][n]
	one        - Identity Matrix
	modl       - as given by equation 33
	modm       - as given by equation 34
	xhat       - x with its nth row and column removed
	xn         - x[n]
	deltal     - unique diagonal matrix
	lhatn      - is the nth column of L with L[n][n] removed

	:Note: cutoff is considered 0 throughout
	:Note: While running this code from the notebook, make sure to use Clear[<graph_variable>] to ensure that previous runs don't affect current runs *)

(* :Mathematica Version: 6.0
   :Combinatorica Version: 2.1 *)

BeginPackage["DiagonalRoundingInSP`", {"GraphPartition`", "Combinatorica`"}]


(*      ----------------------------------------
						DRIVERS
		----------------------------------------     *)

(*  :diagonalRoundingRandom:
	creates a random graph and calls diagonalRoundingMain[] on it.
	:Param: numOfVertices is the number of vertices
	:Option: plot is True if a plot of vectors v1 and v2 is desired.
	:Option: fiedlerOrder is True if fiedler vector ordering of vertices is desired during the display. Vertices are separated into two separate parts of the graph.
	:Option: roundingAlgorithm allows the user to choose an algorithm for rounding the eigen vector that will be computed in the process of partitioning
	:Returns: nothing important
	:Author: Maurya Talisetti *)
Options[diagonalRoundingRandom] = {plot->False, fiedlerOrder->False, roundingAlgorithm->CriterionCut}
diagonalRoundingRandom[numOfVertices_, OptionsPattern[]] := 
	Module[{g},
		g = createRandomGraph[numOfVertices];
		diagonalRoundingMain[g, plot->OptionValue[plot], fiedlerOrder->OptionValue[fiedlerOrder], roundingAlgorithm->OptionValue[roundingAlgorithm]]
	]

(*  :diagonalRoundingMain:
	performs the Diagonal Rounding for Spectral Partitioning algorithm and displays the graph.
	:Param: g is the input graph
	:Option: plot is True if a plot of vectors v1 and v2 is desired.
	:Option: fiedlerOrder is True if fiedler vector ordering of vertices is desired during the display. Vertices are separated into two separate parts of the graph.
	:Option: roundingAlgorithm allows the user to choose an algorithm for rounding the eigen vector that will be computed in the process of partitioning
	:Returns: the resulting graph
	:Author: Maurya Talisetti *)
Options[diagonalRoundingMain] = {plot->False, fiedlerOrder->False, roundingAlgorithm->CriterionCut}
diagonalRoundingMain[g_Graph, OptionsPattern[]] :=
	Module[{result, indicatorVector, unroundedVector, unroundedVector2, graphWithOptions, h}, 
		h = duplicateGraph[g];
(*		h = MakeGraph[Range[V[g]], False &, Type->Undirected];
		h = ChangeVertices[h, Vertices[g]];
		h = AddEdges[h, Edges[g]];*)
		(* if vertex weights are 0, return incorrect input *)

		result = diagonalRounding[g, roundingAlgorithm->OptionValue[roundingAlgorithm]];
		indicatorVector = First[result];
		unroundedVector = Part[result, 2];
		unroundedVector2 = Last[result]; (* vector corresponding to second smallest eigenvalue *)

		outputResults[h, indicatorVector, unroundedVector, unroundedVector2, fiedlerOrder->OptionValue[fiedlerOrder], plot->OptionValue[plot]]
	]

(*  :outputResults:
	computes a graph or plot according to the option specified, sets its edge and vertex options and displays it. Also prints the vectors.
	:Param: g is the input partitioned graph
	:Param: indicatorVector is the rounded vector given by the partitioning algorithm used
	:Param: unroundedVector is the unrounded vector given by the partitioning algorithm used
	:Param: unroundedVector2 is the unrounded vector corresponding to the second smallest eigenvalue given by the partitioning algorithm used
	:Option: plot is True if a plot of vectors v1 and v2 is desired.
	:Option: fiedlerOrder is True if fiedler vector ordering of vertices is desired during the display. Vertices are separated into two separate parts of the graph.
	:Returns: a graph with its vertex and edge options set and ordered according to the option specified
	:Author: Maurya Talisetti *)
Options[outputResults] = {fiedlerOrder->False, plot->False}
outputResults[g_Graph, indicatorVector_List, unroundedVector_List, unroundedVector2_List, OptionsPattern[]] :=
	Module[{},
		If[!OptionValue[fiedlerOrder],
			(Print["\nunrounded vector: ", unroundedVector];
			Print["\nrounded vector: ", indicatorVector, "\n"]),
			(Print["\nsorted unrounded vector: ", Sort[unroundedVector]];
			Print["\nsorted rounded vector: ", Sort[indicatorVector], "\n"])
		];

		If[OptionValue[plot], 
			drawPlot[g, indicatorVector, unroundedVector, unroundedVector2],
			drawGraph[g, indicatorVector, unroundedVector, fiedlerOrder->OptionValue[fiedlerOrder]]
		]
	]

(*  :diagonalRoundingResults:
	performs the Diagonal Rounding for Spectral Partitioning algorithm and returns the vectors and the partitioned graph
	:Param: g is the input graph
	:Option: plot is True if a plot of vectors v1 and v2 is desired.
	:Option: roundingAlgorithm allows the user to choose an algorithm for rounding the eigen vector that will be computed in the process of partitioning
	:Returns: if plot is False, list containing unrounded vector, rounded vector, regular graph, sorted unrounded vector, sorted rounded vector, fiedler order graph
			  else, list containing unrounded vector, rounded vector, plot of vectors v1 and v2
	:Author: Maurya Talisetti *)
Options[diagonalRoundingResults] = {plot->False, roundingAlgorithm->CriterionCut}
diagonalRoundingResults[g_Graph, OptionsPattern[]] :=
	Module[{result, indicatorVector, unroundedVector, unroundedVector2, graphWithOptions, h}, 
		h = duplicateGraph[g];
		(* if vertex weights are 0, return incorrect input *)

		result = diagonalRounding[g, roundingAlgorithm->OptionValue[roundingAlgorithm]];
		indicatorVector = First[result];
		unroundedVector = Part[result, 2];
		unroundedVector2 = Last[result]; (* vector corresponding to second smallest eigenvalue *)

		outputResults2[h, indicatorVector, unroundedVector, unroundedVector2, plot->OptionValue[plot]]
	]

(*  :outputResults2:
	computes a graph or plot according to the option specified, sets its edge and vertex options and returns it along with the vectors
	:Param: g is the input partitioned graph
	:Param: indicatorVector is the rounded vector given by the partitioning algorithm used
	:Param: unroundedVector is the unrounded vector given by the partitioning algorithm used
	:Param: unroundedVector2 is the unrounded vector corresponding to the second smallest eigenvalue given by the partitioning algorithm used
	:Option: plot is True if a plot of vectors v1 and v2 is desired.
	:Option: fiedlerOrder is True if fiedler vector ordering of vertices is desired during the display. Vertices are separated into two separate parts of the graph.
	:Returns: graphs with vertex and edge options set and ordered according to the option specified and the vectors
	:Author: Maurya Talisetti *)
Options[outputResults2] = {plot->False}
outputResults2[g_Graph, indicatorVector_List, unroundedVector_List, unroundedVector2_List, OptionsPattern[]] :=
	Module[{h, regularGraph, fiedlerOrderGraph, vectorPlot},
		If[!OptionValue[plot],
			(h = duplicateGraph[g];
			regularGraph = getPartitionedGraph[g, indicatorVector, unroundedVector];
			fiedlerOrderGraph = getPartitionedGraph[h, indicatorVector, unroundedVector, fiedlerOrder->True];
			Return[{unroundedVector, indicatorVector, regularGraph, Sort[unroundedVector], Sort[indicatorVector], fiedlerOrderGraph}]),
			(* Else *)
			(vectorPlot = plotVectors[g, indicatorVector, unroundedVector, unroundedVector2];
			Return[{unroundedVector, indicatorVector, vectorPlot}])
		]
	]

(*  :duplicateGraph: 
	duplicates the input graph and returns it
	:Param: g is the input graph
	:Returns: a duplicate of the input graph
	:Author: Maurya Talisetti *)
duplicateGraph[g_Graph] :=
	Module[{h},
		h = MakeGraph[Range[V[g]], False &, Type->Undirected];
		h = ChangeVertices[h, Vertices[g]];
		h = AddEdges[h, Edges[g]];
		Return[h]
	]

(*  :partitionGraphDriver[]: 
	Implements the algorithm "Diagonal Rounding in Spectral Partitioning" without attempting to improve the cutsize.
	:Param: g is the graph which is being partitioned
	:Returns: smallest eigenvalue and vectors that correspond to the smallest (represents the optimal cut) and the second smallest eigenvalues
	:Author: Maurya Talisetti *)
partitionGraphDriver[g_Graph] :=
	Module[{modm, modl, mhat, mnn, one, results},
		results = computeMatrices[g];

		modm = results[[1]];
		modl = results[[2]];
		mhat = results[[3]];
		mnn = results[[4]];
		one = results[[5]];

		partitionGraph[g, modm, modl, mhat, mnn, one]
	]


(*      ----------------------------------------
					SAMPLE GRAPHS
		----------------------------------------     *)

(*  :createRandomGraph:
	creates a random graph and returns it.
	:Param: numOfVertices is the number of vertices
	:Returns: the newly created graph
	:Author: Maurya Talisetti *)
createRandomGraph[numOfVertices_] := 
	Module[{g, edgeProb},
		edgeProb = 0.2; (* probability that an edge exists between two vertices *)

		g = MakeGraph[Range[numOfVertices], If[RandomReal[] < edgeProb, True] &, Type->Undirected];

		(* If the randomly generated graph is not connected, then keep generating a new graph until it is connected. *)
		While[!ConnectedQ[g],
			g = MakeGraph[Range[numOfVertices], If[RandomReal[] < edgeProb, True] &, Type->Undirected]];

		g = RemoveSelfLoops[g];
		g = SetVertexWeights[g, Array[1 &, V[g]]];

		g = SetEdgeWeights[g, WeightingFunction->RandomInteger, WeightRange->{1, numOfVertices}];
		Return[g];
	]

		(*      ----------------------------------------
							DOUBLE TREE
				----------------------------------------     *)

(*  :createDoubleTree:
	creates a Double Tree graph with 14 nodes.
	:Returns: the newly computed graph
	:Author: Maurya Talisetti *)
createDoubleTree[] :=
	Module[{g, v, list1, list2, list3, list4, list5, BUNCH, HEIGHT}, 

		v = Table[12i+j, {i, 0, 13}, {j, 12}];
		g = MakeGraph[Range[168], False &, Type->Undirected];
		g = SetVertexWeights[g, Array[1 &, V[g]]];
		BUNCH = 12;
		HEIGHT = 2;

		list4 = xCoordinates[];
		list5 = yCoordinates[HEIGHT, BUNCH];

		list6 = {};
		(*list6 = (#1, #2) &, xCoordinates[], yCoordinates[]);*)
		For[i = 1, i <= Length[list4], i++,
			list6 = Append[list6, {list4[[i]], list5[[i]]}]];
		g = ChangeVertices[g, list6];

		list1 = createTree[Take[v, Length[v]/2]];
		list2 = createTree[Take[v, -Length[v]/2]];
		list3 = createEdges[v[[4]], v[[11]]];
		list4 = internalEdges[v]; (* edges between vertices in a set of vertices *)
		g = AddEdges[g, Join[list1, list2, list3, list4]]
]

(*  :createTree:
	creates a Tree of 7 nodes with v4 being the root. Branches are created at nodes v2, v6, v4.
	:Param: v is a list of list of vertices
	:Returns: the newly computed list of edges
	:Author: Maurya Talisetti *)
createTree[v_List] := 
	Module[{list1, list2, list3},
		list1 = createBranch[v[[1]], v[[2]], v[[3]]];
		list2 = createBranch[v[[5]], v[[6]], v[[7]]];
		list3 = createBranch[v[[2]], v[[4]], v[[6]]];
		Join[list1, list2, list3]]

(*  :createBranch:
	creates a branch with v2 as the node and v1 and v3 having edges to v2. The nth vertex in the first list is connected to the nth vertex in the second list by an 
	edge and the same applies to the pair of lists v3 and v2.
	:Param: v1 is a list of vertices
	:Param: v2 is a list of vertices
	:Param: v2 is a list of vertices
	:Returns: the newly computed list of edges
	:Author: Maurya Talisetti *)
createBranch[v1_List, v2_List, v3_List] :=
	Module[{list1, list2},
		list1 = createEdges[v1, v2];
		list2 = createEdges[v2, v3];
		Join[list1, list2]]

(*  :createEdges:
	creates edges between two lists of vertices. The nth vertex in the first list is connected to the nth vertex in the second list by the edge.
	:Param: v1 is a list of vertices
	:Param: v2 is a list of vertices
	:Returns: the newly computed list of edges
	:Author: Maurya Talisetti *)
createEdges[v1_List, v2_List] :=
	Module[{i, list}, 
		list = {};
		For[i = 1, i <= 12, i++, 
			list = Append[list, {v1[[i]], v2[[i]]}]];
		Return[list]]

(*  :internalEdges:
	given a list of vertex lists, computes the list of edges between vertices in each list of vertices.
	:Param: v is a list of vertex lists
	:Returns: the list of edges
	:Author: Maurya Talisetti *)
internalEdges[v_List] :=
	Module[{i, list, w},
		list = {};

		(* for each vertex set *)
		For[i = 1, i <= Length[v], i++, 
			(w = v[[i]];

				(* for each pair of vertices in a vertex set *)
				For[j = 1, j <= Length[w] - 1, j++, 
					list = Append[list, {w[[j]], w[[j + 1]]}]
				]
			)
		];

		Return[list];
	]

(*  :xCoordinates:
	changes the x coordinates of the vertices of the input graph to reflect a double tree structure.
	:Param: g is the graph whose vertices are being changed
	:Returns: nothing important
	:Author: Maurya Talisetti *)
xCoordinates[] :=
	Module[{xlist, deltaX, prevX, i},
		xlist = Array[0 &, 84];
		deltaX = N[2/84];
		prevX = -1;
		For[i = 1, i <= 84, i++, (xlist[[i]] = prevX + deltaX; prevX = xlist[[i]])];
		Return[Join[xlist, xlist]]
	]

(*  :yCoordinates:
	creates a list of y coordinates for all the vertices in a double tree
	:Param: height is the height of a single tree in the double tree graph
	:Param: bunch is the number of subnodes at each node of the double tree
	:Returns: the list of y coordinates for all the vertices
	:Author: Maurya Talisetti *)
yCoordinates[height_, bunch_] :=
	Module[{incr, coordinatesTree1, coordinatesTree2, coordinates, root1, root2},
		(* how much to increment by depends on the number of levels*)
		incr = 0.9/height;

		root1 = -0.1;
		root2 = 0.1;
		
		coordinatesTree1 = giveCoordinates[height, root1, -incr, bunch];
		coordinatesTree2 = giveCoordinates[height, root2, incr, bunch];

		coordinates = Join[coordinatesTree1, coordinatesTree2];
		Return[coordinates];
	]

(*  :giveCoordinates:
	creates a list of y-coordinates for all the vertices in a subtree of height "height ".
	:Param: height is the height of a the subtree
	:Param: root is the y-coordinate of the root of the subtree
	:Param: incr is the y-coordinate distance between consecutive levels of the nodes
	:Param: bunch is the number of subnodes at each node of the subtree
	:Returns: the list of y-coordinates for all the vertices in the subtree
	:Author: Maurya Talisetti *)
giveCoordinates[height_, root_, incr_, bunch_] := 
	Module[{list1, list2}, 
		list2 = Array[root &, bunch];
		If[height == 0, Return[list2]];
		
		list1 = giveCoordinates[height - 1, root + incr, incr, bunch];
		Return[Join[list1, list2, list1]];
	]


(*      ----------------------------------------
					GRAPH DISPLAY
		----------------------------------------     *)

(*  :drawGraph:
	computes a graph ordered according to the option specified, sets its edge and vertex options and displays it.
	:Param: g is the input partitioned graph
	:Param: indicatorVector is the rounded vector given by the partitioning algorithm used
	:Param: unroundedVector is the unrounded vector given by the partitioning algorithm used
	:Option: plot is True if a plot of vectors v1 and v2 is desired.
	:Option: fiedlerOrder is True if fiedler vector ordering of vertices is desired during the display. Vertices are separated into two separate parts of the graph.
	:Returns: a graph with its vertex and edge options set and ordered according to the option specified
	:Author: Maurya Talisetti *)
Options[drawGraph] = {plot->False, fiedlerOrder->False}
drawGraph[g_Graph,  indicatorVector_List, unroundedVector_List, OptionsPattern[]] :=
	Module[{graphWithOptions},
		graphWithOptions = getPartitionedGraph[g, indicatorVector, unroundedVector, fiedlerOrder->OptionValue[fiedlerOrder]];

		(*	printValues[{"diagonalRoundingMain", result, indicatorVector, unroundedVector, unroundedVector2, g, graphWithOptions}];*)

		ShowGraph[graphWithOptions]
	]

(*  :drawGraph:
	plots the vertices according to the eigenvectors corresponding to the smallest and the second smallest eigenvalues and displays it
	:Param: g is the input partitioned graph
	:Param: indicatorVector is the rounded vector given by the partitioning algorithm used
	:Param: unroundedVector is the unrounded vector given by the partitioning algorithm used
	:Param: unroundedVector2 is the unrounded vector corresponding to the second smallest eigenvalue given by the partitioning algorithm used
	:Returns: the newly created plot
	:Author: Maurya Talisetti *)
drawPlot[g_Graph, indicatorVector_List, unroundedVector_List, unroundedVector2_List] :=
	Module[{graphWithOptions},
		graphWithOptions = plotVectors[g, indicatorVector, unroundedVector, unroundedVector2];
		ShowGraph[graphWithOptions]
	]

(*  :plotVectors:
	plots the vertices according to the eigenvectors corresponding to the smallest and the second smallest eigenvalues.
	:Param: g is the input graph
	:Param: v1 is the eigenvector corresponding to the smallest eigenvalue
	:Param: v2 is the eigenvector corresponding to the second smallest eigenvalue
	:Returns: the newly created plot graph
	:Author: Maurya Talisetti *)
plotVectors[g_Graph, indicatorVector_List, v1_List, v2_List] := 
	Module[{numOfVert, i, vertices, plot},
		numOfVert = V[g];
		vertices = Array[{0, 0} &, numOfVert];
		For[i = 1, i <= numOfVert, i++, vertices[[i]] = {v1[[i]], v2[[i]]}];
		plot = ChangeVertices[g, vertices];

		(* graph options *)
		setGraphOptions[plot, indicatorVector]
	]

(*  :getPartitionedGraph:
	computes a partitioned graph whose vertices are ordered according to the option specified and returns it with its vertex and edge options specified.
	:Param: g is the input graph
	:Param: indicatorVector is the rounded vector given by diagonal rounding
	:Param: unroundedVector is the unrounded vector given by diagonal rounding
	:Option: fiedlerOrder is True if fiedler vector ordering of vertices is desired during the display. Vertices are separated into two separate parts of the graph.
	:Option: permVector is the permutation vector used to reorder vertices according to the fiedler order
	:Returns: a partitioned graph with its vertex and edge options set and ordered according to the option specified
	:Author: Maurya Talisetti *)
Options[getPartitionedGraph] = {fiedlerOrder->False, permVector->{}}
getPartitionedGraph[g_Graph, indicatorVector_List, unroundedVector_List, OptionsPattern[]] := 
	(If[OptionValue[fiedlerOrder],

		(* then, display according to fiedler order *)
		If[Length[OptionValue[permVector]] > 0,
			createFiedlerOrderGraph[g, indicatorVector, permVector->OptionValue[permVector]],
			setFiedlerOrder[g, indicatorVector, unroundedVector]],
	
		(* else, original ordering of vertices *)
		setGraphOptions[g, indicatorVector]])

(*  :setFiedlerOrder:
	computes a new graph with vertices in the order of the fiedler vector and sets vertex and edge options.
	:Param: g is the input graph
	:Param: indicatorVector is the rounded vector given by diagonal rounding
	:Param: unroundedVector is the unrounded vector given by diagonal rounding
	:Returns: a graph with its vertex and edge options set and ordered according to the fiedler vector
	:Author: Maurya Talisetti *)
setFiedlerOrder[g_Graph, indicatorVector_List, unroundedVector_List] := 
	Module[{permutationVector, fiedlerOrderGraph, newEdges}, 
		permutationVector = findPermutationVector[unroundedVector];
		createFiedlerOrderGraph[g, indicatorVector, permutationVector]]

(*  :createFiedlerOrderGraph:
	given the permutation vector, computes a new graph with vertices in the order of the fiedler vector and sets vertex and edge options.
	:Param: g is the input graph
	:Param: indicatorVector is the rounded vector given by diagonal rounding
	:Param: permVector maps the components of the unrounded vector to the sorted unrounded vector
	:Returns: a graph with its vertex and edge options set and ordered according to the fiedler vector
	:Author: Maurya Talisetti *)
createFiedlerOrderGraph[g_Graph, indicatorVector_List, permVector_List] :=
	Module[{fiedlerOrderGraph},
		fiedlerOrderGraph = reorderGraph[g, permVector];
		setGraphOptions[fiedlerOrderGraph, Sort[indicatorVector]]]

(*  :reorderGraph:
	computes a new graph with vertices in the order of the sorted fiedler vector.
	:Param: g is the input graph
	:Param: permutationVector maps the components of the fiedler vector to the sorted fiedler vector
	:Returns: a graph with its vertices ordered according to the sorted fiedler vector
	:Author: Maurya Talisetti *)
reorderGraph[g_Graph, permutationVector_] :=
	Module[{fiedlerOrderGraph, newEdges}, 
		(*construct new graph *)
		fiedlerOrderGraph = MakeGraph[Range[V[g]], False &, Type->Undirected];

		(* fill edges in this graph *)
		newEdges = Map[{permutationVector[[First[#]]], permutationVector[[Last[#]]]} &, Edges[g]];
		fiedlerOrderGraph = AddEdges[fiedlerOrderGraph, newEdges];

		(* set edge weights *)
		fiedlerOrderGraph = SetEdgeWeights[fiedlerOrderGraph, GetEdgeWeights[g]];
		Return[fiedlerOrderGraph];
	]

		(*      ----------------------------------------
							GRAPH OPTIONS
				----------------------------------------     *)

(*  :setGraphOptions:
	sets the color options for vertices according to the partition and style option for the edges which are cut.
	:Param: g is the partitioned graph which is being displayed
	:Param: indicatorVector is the rounded vector given by diagonal rounding
	:Returns: a graph with its options set
	:Author: Maurya Talisetti *)
setGraphOptions[g_Graph, indicatorVector_List] :=
	Module[{i, subset1, subset2, cutoff}, 
		cutoff = 0;
		subset1 = {};
		subset2 = {};

		(* divide vertices based on the partition into two subsets *)
		For[i = 1, i <= V[g], i++, 
			(subset1 = If[indicatorVector[[i]] <= cutoff, Append[subset1, i], subset1]; subset2 = If[indicatorVector[[i]] > cutoff, Append[subset2, i], subset2])];

		(* set vertex options *)
		g = SetGraphOptions[g, {Append[subset1, VertexColor->Blue], Append[subset2, VertexColor->Red]}];

		(* set edge options. this is optional because it takes n^2 running time. *)
		setEdgeOptions[g, subset1, subset2];
		Return[g]]

(*  :setEdgeOptions:
	adds style to the edges crossed by the cut.
	:Param: g is the input graph
	:Param: v1 is the first subset of the partition of the vertices
	:Param: v2 is the second subset of the partition of the vertices
	:Returns: nothing important
	:Author: Maurya Talisetti *)
setEdgeOptions[g_Graph, v1_List, v2_List] := 
	Module[{i, j},
		For[i = 1, i <= Length[v1], i++, 
			For[j = 1, j <= Length[v2], j++, 
				g = SetGraphOptions[g, {{v1[[i]], v2[[j]]}, EdgeStyle->Dashed, EdgeColor->Red}]]]]

(*      ----------------------------------------
					DIAGONAL ROUNDING
		----------------------------------------     *)

(*  :diagonalRounding[]: 
	Implements the algorithm "Diagonal Rounding in Spectral Partitioning"
	:Param: g is the graph which is being partitioned.
	:Option: roundingAlgorithm allows the user to choose an algorithm for rounding the eigen vector that will be computed in the process of partitioning
	:Returns: the vector (both the rounded and the unrounded) that represents the optimum/improved cut and the eigenvector that corresponds to the second smallest
		eigenvalue
	:Author: Maurya Talisetti *)
Options[diagonalRounding] = {roundingAlgorithm->CriterionCut}
diagonalRounding[g_Graph, OptionsPattern[]] := 
	(Module[{mhat, mnn, one, modl, modm, x, x2, lambda, cutDetails, improvedResult, eigenSoln, solution, results},

		results = computeMatrices[g];
		modm = results[[1]];
		modl = results[[2]];
		mhat = results[[3]];
		mnn = results[[4]];
		one = results[[5]];

		(* find x and x2 (vectors corresponding to smallest and second smallest eigenvalues) *)
		solution = partitionGraph[g, modm, modl, mhat, mnn, one];
		lambda = First[solution];
		x = Part[solution, 2];
		x2 = Last[solution];

		(* round x *)
		cutDetails = roundEigenvector[g, x, roundingAlgorithm->OptionValue[roundingAlgorithm]];
		
		(* check if a better cut can by found by repetition *)
		improvedResult = Improvement[g, lambda, modm, {x, x2}, cutDetails, modl];

		Print["Number of iterations performed: ", Last[improvedResult]];

(*		printValues[{"diagonalRounding", g, m, mhat, mnn, one, modl, modm, dhat, x, x2, lambda, cutDetails, numOfVert, improvedResult, eigenSoln, mass}];*)

		Return[Take[improvedResult, 3]]])

(*  :partitionGraph[]: 
	Implements the algorithm "Diagonal Rounding in Spectral Partitioning" without attempting to improve the cutsize.
	:Param: g is the graph which is being partitioned
	:Param: modl is the modified Laplacian matrix as given by equation 33
	:Param: modm is the modified Mass matrix
	:Param: mhat is the Mass Matrix with its last row and column removed.
	:Param: mnn is the element in that falls in the last row and column of the Mass Matrix.
	:Param: one is the Identity vector.
	:Returns: smallest eigenvalue and vectors that correspond to the smallest (represents the optimal cut) and the second smallest eigenvalues
	:Author: Maurya Talisetti *)
partitionGraph[g_Graph, modm_List, modl_List, mhat_List, mnn_, one_List] :=
	Module[{x, x2, lambda, eigenSoln, mass},
		mass = Mass[g];

		(* find x and x2 (vectors corresponding to smallest and second smallest eigenvalues) *)
		eigenSoln = findEigenSoln[modl, mass, {modm, mhat, mnn}, one];
		lambda = First[eigenSoln];
		x = Part[eigenSoln, 2];
		x2 = Last[eigenSoln];

(*		printValues[{"partitionGraph", g, m, mhat, mnn, one, modl, modm, dhat, x, x2, lambda, numOfVert, eigenSoln, mass}];*)

		Return[{lambda, x, x2}]]

(*  :Improvement[]: 
	Repeats the Diagonal Rounding algorithm until it doesn't find a better cut.
	:Param: g is the graph which is being partitioned.
	:Param: prevlambda is the eigenvalue obtained in the first round of execution.
	:Param: modm is the modified Mass Matrix.
	:Param: unroundedPrevx FILL THIS
	:Param: oldCutDetails is a list containing the vector that represents the cut and the cutoff value which were obtained in the first round of execution.
	:Param: prevmodl is the modified Laplacian Matrix, obtained in the first round of execution.
	:Returns: a list that contains the vectors (both rounded and unrounded) that represent the optimum cut and the number of iterations performed.
	:Author: Maurya Talisetti *)
Improvement[g_Graph, prevlambda_, modm_, unroundedPrevx_List, oldCutDetails_, prevmodl_] :=
	Module[{y, list, d, modl, eigenSoln, lambda, mhat, m, one, mnn, x, x2, unroundedOldx, unroundedOldx2, oldx, oldCutSize, oldxhat, oldlambda, oldmodl, 
		numOfIter, numOfVert, i, unroundedx2, unroundedxhat},
		m = MassMatrix[g]; numOfVert = V[g]; mhat = Mhat[m, numOfVert]; one = ConstantArray[1, numOfVert - 1]; mnn = Mnn[m, numOfVert]; oldlambda = prevlambda; 
		oldmodl = prevmodl; numOfIter = 0; unroundedOldx = First[unroundedPrevx]; unroundedOldx2 = Last[unroundedPrevx];

		oldx = First[oldCutDetails];
		oldCutSize = oldCutDetails[[2]];
		oldxhat = Part[oldx, 1;;numOfVert - 1];

		While[True,
			numOfIter++;

			(* Print iteration details *)
			printIterDetails[numOfIter, oldlambda];
	
			(* Find the unique diagonal matrix deltal *)
			d = findUniqDiagMatrix[oldmodl, oldxhat];

			(* Calc new modl *)
			modl = oldmodl + d;

			(* Calc smallest and second smallest eigenvalue and eigenvector *)
			eigenSoln = findEigenSoln[modl, Mass[g], {modm, mhat, mnn}, one];
			lambda = First[eigenSoln];
			x = Part[eigenSoln, 2];
			unroundedxhat = Part[x, 1;;numOfVert - 1];
			unroundedx2 = Last[eigenSoln];

(*			Print["lambda is :", N[lambda]]; (* to be removed *)
			Print["oldlambda is :", N[oldlambda]]; (* to be removed *)
			If[N[lambda, 5] == N[oldlambda, 5], Print["true"], Print["false"]];*)
			
			(* test code - can be safely removed *)
			(*checkEigen[oldxhat, oldlambda, eigenSoln];*)

			If[lambda == 0, 
				(* If the smallest eigenvalue of the updated eigensystem is equal to 0, then optimal cut is found. *)
				(Print["Optimal Cut found!"];
					Print["eigenvalue: ", lambda];
					v1 = oldx;
					v2 = x;
					mm = m;
					Return[{oldx, unroundedOldx, unroundedOldx2, numOfIter}]), 

				(* Else if lambda < 0 *)
				(* first ensure that it is not a false positive *)
				If[vectorsAlmostEqual[oldx, x, m],
					(Print["Optimal Cut found!"];
						Print["eigenvalue: ", lambda];
						Return[{oldx, unroundedOldx, unroundedOldx2, numOfIter}])					
				];

				(* compute the indicator vector found by rounding the eigenvector (associated with lambda) of the updated eigensystem. *)
				Module[{xhat, xn, unroundedx, cutDetails, cutSize},
					(* round x *)
					cutDetails = roundEigenvector[g, x];
					unroundedx = x;
					x = First[cutDetails];
					xhat = Part[x, 1;;numOfVert - 1];

					cutSize = cutDetails[[2]];

					(* to be removed *)
					(*Print["old cut size :", N[oldCutSize]];
					Print["new cut size :", N[cutSize]];*)

					If[cutSize >= oldCutSize, 
						(* If the indicator vector found by rounding the eigenvector x represents a bigger cut than the previously computed cut,
							then, stop looking for a better cut. *)
						(Print["Optimal Cut NOT found!"];
							Print["eigenvalue of the cut found: ", oldlambda, "\n"];
							Return[{oldx, unroundedOldx, unroundedOldx2, numOfIter}]), 

						(* Else, continue iterating until there is no more improvement. *)
						Module[{}, oldlambda = lambda;
							oldmodl = modl;
							oldCutSize = cutSize;
							oldx = x;
							unroundedOldx = unroundedx;
							unroundedOldx2 = unroundedx2;
							oldxhat = xhat]
					];

(*					printValues["Improvement2", xhat, xn, unroundedx, cutDetails, cutSize];*)

				]
			]
		];

(*		printValues["Improvement1 ", g, prevlambda, modm, unroundedPrevx, oldCutDetails, prevmodl, y, list, d, modl, eigenSoln, lambda, mhat, m, one, mnn, x, x2, 
			unroundedOldx, unroundedOldx2, oldx, oldCutSize, oldxhat, oldlambda, oldmodl, numOfIter, numOfVert, i, unroundedx2, unroundedxhat];*)

	]

(*  :findEigenSoln[]: finds the eigen solution and returns the vectors and eigenvalue
	:Param: modl is the modified Laplacian matrix as given by equation 33
	:Param: mass is the mass of the input graph
	:Param: m is a list that contains modm, mhat and mnn where
			modm is the modified Mass Matrix.
			mhat is M  with its nth rows and columns removed
			mnn is M[n][n]
	:Param: one is the Identity Matrix
	:Returns: returns the new vector which includes the nth component
	:Author: Maurya Talisetti *)
findEigenSoln[modl_List, mass_, m_List, one_] := 
	Module[{modm, mhat, mnn, eigenSoln, lambda, xhat, xhat2, x, x2}, 
		modm = First[m];
		mhat = Part[m, 2];
		mnn = Last[m];

		eigenSoln = Eigensystem[{modl, modm}];

		lambda = Last[First[eigenSoln]];
		xhat = Last[Last[eigenSoln]];

		(* find eigenvector2 *)
		xhat2 = Part[Last[eigenSoln], Length[Last[eigenSoln]] - 1];

		(* find x and x2 *)
		x = scaleEigen[xhat, mass, m, one];
		x2 = scaleEigen[xhat2, mass, m, one];

(*		printValues[{"findEigenSoln", modl, mass, m, one, modm, mhat, mnn, eigenSoln, lambda, xhat, xhat2, x, x2}];*)

		Return[{lambda, x, x2}]]

(*  :roundEigenvector: 
	Rounds the supplied eigenvector using criterion cut and returns it and the cut size. Also computes the permutationVector which maps the components of the fiedler
	vector x to the sorted fiedler vector. Will be expanded to include the option of using isoperimetric/ratio cut.
	:Param: g_Graph is the graph which is being partitioned
	:Param: x is the eigenvector which is to be rounded
	:Option: roundingAlgorithm allows the user to choose an algorithm for rounding the input eigen vector
	:Returns: a list containing the rounded vector, the cut size and the permutation vector
	:Author: Maurya Talisetti *)
Options[roundEigenvector] = {roundingAlgorithm->CriterionCut}
roundEigenvector[g_Graph, x_List, OptionsPattern[]] := 
	Module[{cutDetails, i, cutoff, cutSize, a, permutationVector},
		cutDetails = OptionValue[roundingAlgorithm][g, x, IncludeCut->True, Range->{1/Length[x], 1}];
		cutoff = First[cutDetails];
		cutSize = Last[cutDetails];
		a = x;
		For[i = 1, i <= V[g] , i++,
			a[[i]] = If[a[[i]] < cutoff, -1, 1]];

		permutationVector = findPermutationVector[x];

(*		printValues[{"roundEigenvector", g, x, OptionValue[roundingAlgorithm], cutDetails, cutoff, cutSize, a, permutationVector}];*)

		Return[{a, cutSize, permutationVector}];
	]

(*  :reorderVector:
	computes the sorted fiedler vector without actually sorting the fiedler vector.
	:Param: fiedlerVector is the fiedler vector given by any partitioning algorithm
	:Param: permutationVector maps the components of the fiedler vector to the sorted fiedler vector
	:Returns: the sorted fiedler vector
	:Time: O(Length[fiedlerVector])
	:Author: Maurya Talisetti *)
reorderVector[fiedlerVector_, permVector_] :=
	Module[{length, sortedVector}, 
		length = Length[fiedlerVector];
		sortedVector = Array[0 &, length];
		For[i = 1, i <= length, i++,
			sortedVector[[permVector[[i]]]] = fiedlerVector[[i]]
		];

		Return[sortedVector];
	]

(*  :findPermutationVector:
	computes a permutation vector which maps the components of the unrounded vector to the sorted unrounded vector
	:Param: unroundedVector is the unrounded vector given by any partitioning algorithm
	:Returns: permutationVector which maps the components of the unrounded vector to the sorted unrounded vector
	:Author: Maurya Talisetti *)
findPermutationVector[unroundedVector_List] :=
	Module[{fiedlerOrderVector, permutationVector},
		fiedlerOrderVector = Sort[unroundedVector];
		permutationVector = Array[0 &, Length[unroundedVector]];
		For[i = 1, i <= Length[unroundedVector], i++, 
			(permutationVector[[i]] = First[First[Position[fiedlerOrderVector, unroundedVector[[i]]]]]; fiedlerOrderVector[[permutationVector[[i]]]] = $MaxNumber)];
		Return[permutationVector];
	]

(*  :partitionAndRound[]: 
	find the eigenvector that dictates the partition of graph using "Diagonal Rounding In Spectral Partitioning" and rounds it.
	:Param: g is the graph which is being partitioned.
	:Option: roundingAlgorithm allows the user to choose an algorithm for rounding the eigen vector that will be computed in the process of partitioning
	:Returns: a list containing the rounded vector, the cut size and the permutation vector
	:Author: Maurya Talisetti *)
Options[partitionAndRound] = {roundingAlgorithm->CriterionCut}
partitionAndRound[g_Graph] := 
	Module[{results, modm, modl, mhat, mnn, one, eigenSoln}, 
		results = computeMatrices[g];
		modm = results[[1]];
		modl = results[[2]];
		mhat = results[[3]];
		mnn = results[[4]];
		one = results[[5]];

		eigenSoln = partitionGraph[g, modm, modl, mhat, mnn, one];
		Return[roundEigenvector[g, Part[eigenSoln, 2], roundingAlgorithm->OptionValue[roundingAlgorithm]]];
	]

		(*      ----------------------------------------
					DIAGONAL ROUNDING HELPER FUNCTIONS
				----------------------------------------     *)

(*  :scaleEigen[]: scales the eigen vector xhat to satisfy the quadratic constraint (34), computes the nth component and returns the new vector.
	:Param: xhat is the eigen vector being scaled
	:Param: mass is the mass of the input graph
	:Param: m is a list that contains modm, mhat and mnn where
			modm is the modified Mass Matrix.
			mhat is M  with its nth rows and columns removed
			mnn is M[n][n]
	:Param: one is the Identity Matrix
	:Returns: returns the new vector which includes the nth component
	:Author: Maurya Talisetti *)
scaleEigen[xhat_List, mass_, m_List, one_List] := 
	Module[{newxhat, xn, x, modm, mhat, mnn}, 
		modm = First[m];
		mhat = Part[m, 2];
		mnn = Last[m];

		(* scale xhat *)
		newxhat = Flatten[xhat]; (* not necessary? *)
		newxhat = scaleVector[newxhat, mass, modm];
		xn = - (one.mhat.newxhat)/mnn;
		x = Join[newxhat, {xn}];

(*		printValues["scaleEigen", xhat, mass, m, one, newxhat, xn, x, modm, mhat, mnn];*)

		Return[x];
]

(*  :computeMatrices: 
	calculates modm, modl, mhat, mnn and one and returns them
	:Param: g is the graph which is being partitioned
	:Param: modm is the modified Mass matrix
	:Param: modl is the modified Laplacian matrix as given by equation 33
	:Param: mhat is the Mass Matrix with its last row and column removed.
	:Param: mnn is the element in that falls in the last row and column of the Mass Matrix.
	:Param: one is the Identity vector.
	:Returns: modm, modl, mhat, mnn, one as a list
	:Author: Maurya Talisetti *)
computeMatrices[g_Graph] :=
	Module[{m, numOfVert, mhat, mnn, one, modm, modl, keyResults},
		m = MassMatrix[g]; numOfVert = V[g]; mhat = Mhat[m, numOfVert]; mnn = Mnn[m, numOfVert]; one = ConstantArray[1, numOfVert - 1];

		(* calculate dhat, modl, modm *)
		keyResults = setUp[g, m, mhat, mnn, one];
		modl = keyResults[[1]];
		modm = keyResults[[2]];

		Return[{modm, modl, mhat, mnn, one}];
	]

(*  :setUp: 
	calculates dhat, modl and modm and returns them
	:Param: g is the graph which is being partitioned
	:Param: m is the mass matrix
	:Param: mhat is the Mass Matrix with its last row and column removed.
	:Param: mnn is the element in that falls in the last row and column of the Mass Matrix.
	:Param: one is the Identity vector.
	:Returns: modl and modm
	:Author: Maurya Talisetti *)
setUp[g_Graph, m_List, mhat_List, mnn_, one_List] :=
	Module[{dhat, modl, modm},
		(* dhat is zero initially *)
		dhat = initDhat[V[g]];

		(* calculate modl, modm*)
		modl = N[ModifiedLaplacianMatrix[g, m, mhat, mnn, one, dhat]];
		modm = N[ModifiedMassMatrix[mhat, mnn, one]];

		Return[{modl, modm}];
	]

(*  :vectorsAlmostEqual: 
	checks if the input vectors are almost equal
	:Param: x1 is a vector
	:Param: x2 is a vector
	:Param: modm is the modified mass matrix of a graph
	:Returns: returns True if x1 is almost equal to x2 and returns False otherwise
	:Author: Maurya Talisetti *)
vectorsAlmostEqual[x1_List, x2_List, m_List] :=
	Module[{k},
		k = (x2.m.x1)/(Sqrt[(x2.m.x2)(x1.m.x1)]);
		If[k < 0.5, Return[False], Return[True]];
	]

(*  :findUniqDiagMatrix: 
	finds the unique diagonal matrix d that satisfies (modl + d) xhat = 0
	:Param: oldmodl is the previously calculated modl
	:Param: oldxhat is the previously calculated eigenvector
	:Returns: the computed diagonal matrix
	:Author: Maurya Talisetti *)
findUniqDiagMatrix[oldmodl_List, oldxhat_List] := 
	Module[{y, list, d, i},
		y = - oldmodl.oldxhat;
		list = Table[y[[i]]/oldxhat[[i]], {i, Length[oldxhat]}];
		d = DiagonalMatrix[list];

		Return[d];
	]

(*  :scaleVector: 
	scales the eigenvector to satisfy the quadratic constraint (34).
	:Param: xhat is the eigenvector being scaled.
	:Param: mass is the mass of the graph.
	:Param: modm is the modified Mass Matrix.
	:Author: Maurya Talisetti *)
scaleVector[xhat_List, mass_, modm_List] :=
	Module[{theta}, 
		theta = Sqrt[mass/(xhat.modm.xhat)];
		Return[theta xhat]]

(*  :ModifiedLaplacianMatrix: 
	returns the modified Laplacian Matrix.
	:Param: g is the graph which is being partitioned.
	:Param: mass is the mass of the graph.
	:Param: mhat is the Mass Matrix with its last row and column removed.
	:Param: mnn is the element in that falls in the last row and column of the Mass Matrix.
	:Param: one is the Identity vector.
	:Param: dhat is the Diagonal Matrix D with its last row and column removed.
	:Author: Maurya Talisetti *)
ModifiedLaplacianMatrix[g_Graph, m_List, mhat_List, mnn_, one_List, dhat_List] :=
	Module[{n, l, lhatn, Dnn}, n = V[g]; l = LaplacianMatrix[g]; lhatn = Lhatn[l, n]; Dnn = 0;
		Lhat[l, n] + dhat - 1/mnn (lhatn.Transpose[Transpose[{one}]].mhat + mhat.Transpose[{one}].Transpose[lhatn])
			+ ((Lnn[l, n] + Dnn)/mnn^2) mhat.Transpose[{one}].Transpose[Transpose[{one}]].mhat]

(*  :ModifiedMassMatrix: 
	returns the modified Mass Matrix.
	:Param: mhat is the Mass Matrix with its last row and column removed.
	:Param: mnn is the element in that falls in the last row and column of the Mass Matrix.
	:Param: one is the Identity vector.
	:Author: Maurya Talisetti *)
ModifiedMassMatrix[mhat_List, mnn_, one_List] := mhat + (mhat.Transpose[{one}].Transpose[Transpose[{one}]].mhat)/mnn 
(* if possible, look for an elegant way to transpose *)

(*  :chooseLambda: 
	chooses lambda to minimize the expression |lambda modm.xhat - modl.xhat|^2
	:Param: modm is the modified Mass Matrix.
	:Param: modl is the modified Laplacian Matrix.
	:Param: xhat is the eigenvector of the modified Laplacian Matrix with respect to the modified Mass Matrix.
	:Returns: lambda that minimizes the above expression.
	:Author: Maurya Talisetti *)
chooseLambda[modm_List, modl_List, xhat_List] := ((modm.xhat).modl.xhat)/((modm.xhat).modm.xhat)

Lhat[l_, n_] := RowColumnDrop[l, n]

Lhatn[l_List, n_] := Take[l, {1, n - 1}, {n, n}]

Lnn[l_List, n_] := l[[n, n]]

initDhat[n_] := Table[0, {i, n - 1}, {j, n - 1}]

Mnn[m_, numOfVert_] := m[[numOfVert, numOfVert]]
	
Mhat[m_List, numOfVert_] := RowColumnDrop[m, numOfVert]

MassMatrix[g_Graph] := DiagonalMatrix[GetVertexWeights[g]]

Mass[g_Graph] := Total[GetVertexWeights[g]]

		(*      ----------------------------------------
						PRINTING AND DEBUGGING
				----------------------------------------     *)

(*  :printIterDetails: 
	prints computed values during this iteration.
	:Param: numOfIter is the iteration number
	:Param: lambda is the eigenvalue found during this iteration
	:Returns: nothing important
	:Author: Maurya Talisetti *)
printIterDetails[numOfIter_, lambda_] := (Print["Iteration: ", numOfIter, "\n eigenvalue computed: ", lambda, "\n "])

(*  :checkEigen: (Test code)
	checks if the rounded eigen vector/value pair is an eigen vector/value pair of the new eigensystem and prints the appropriate message.
	:Param: oldxhat is the rounded eigenvector in the previous iteration.
	:Param: oldlambda is the eigenvalue in the previous iteration.
	:Param: newEigenSoln is the solution of the new eigen system in this iteration.
	:Returns: nothing important
	:Author: Maurya Talisetti *)
checkEigen[oldxhat_List, oldlambda_, newEigenSoln_List] := 
	(If[MemberQ[Last[newEigenSoln], oldxhat] && MemberQ[First[newEigenSoln], oldlambda], Print["checkEigen passed!"], Print["checkEigen failed!"]])

(*  :printValues:
	prints each of the elements of the list listOfArguments. This function is used primarily for debugging. The variables in each function form the listOfArguments.
	By convention, the first element of the list is the name of the function.
	:Param: listOfArguments is a list of elements
	:Returns: nothing important
	:Author: Maurya Talisetti *)
printValues[listOfArguments_List] :=
	Module[{i},
		For[i = 1, i <= Length[listOfArguments], i++, Print[listOfArguments[[i]]]];
	]


(*  ----------------------------------------------------------------------
							MESH STUFF
	---------------------------------------------------------------------- *)

FaceListToEdgeList[fList_] := Module[{labeledFList,edgeList},
	labeledFList=MapIndexed[{#1,First[#2]}&,fList];
	edgeList={};
	While[Length[labeledFList]>0,
			With[{face=First[labeledFList]},
				labeledFList=Rest[labeledFList];
				edgeList=Join[edgeList,Apply[Join,Map[ If[ Length[Intersection[First[#],First[face]]]>=2,
															{{Last[#],Last[face]}},
															{}]&,labeledFList]]];]];
	edgeList]
	
(* :FacesToGraph:
	Converts a list of faces into a graph where each face is a vertex,
	and edges connect faces that share an edge. *)
FacesToGraph[fList_] := AddEdges[EmptyGraph[Length[fList]],
	FaceListToEdgeList[fList]]

(* :DiagonalRoundingObj:
	NOTE: OBJ parsing code written by Mark McClur, found at:
	http://facstaff.unca.edu/mcmcclur/java/LiveMathematica/real.html *)
DiagonalRoundingObj[s_String] := Module[
	{data,vertices,faces,vector,cutoff},

	data=ReadList[s,Word,RecordLists->True];
	If[data=!=$Failed,
		vertices=Map[ToExpression,Drop[#,1]&/@
			Select[data,First[#]=="v "&],{2}];
		faces=Drop[#,1]&/@
			Select[data,First[#]=="f "&];
		If[vertices=!={}&&faces=!={},
			faces=If[DigitQ[faces[[1]][[1]]],
				Map[ToExpression,faces,{2}],
				Map[ToExpression[
					StringTake[#,StringPosition[#,"/"][[1]][[1]]-1]]&,
						faces,{2}]];
				With[{graph=FacesToGraph[faces]},
				vector=First[diagonalRounding[graph]];
				cutoff=0; (* temporary *)
				(*cutoff=Criterion[graph,vector];*)
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
	{data,faces},

	data=ReadList[s,Word,RecordLists->True];
	If[data=!=$Failed,
		faces=Drop[#,1]&/@
			Select[data,First[#]=="f "&];
		If[vertices=!={}&&faces=!={},
			faces=If[DigitQ[faces[[1]][[1]]],
				Map[ToExpression,faces,{2}],
				Map[ToExpression[
					StringTake[#,StringPosition[#,"/"][[1]][[1]]-1]]&,
						faces,{2}]];
			FacesToGraph[faces],
			Print[s<>" does not seem to be an Obj file."]
		],
		Print[s<>" not found."]
	]
]

EndPackage[ ]
