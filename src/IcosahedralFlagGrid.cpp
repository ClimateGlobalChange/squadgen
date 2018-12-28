///////////////////////////////////////////////////////////////////////////////
///
///	\file    IcosahedralFlagGrid.cpp
///	\author  Paul Ullrich
///	\version February 21, 2012
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "IcosahedralFlagGrid.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

void GenerateQuadrilaterals(
	int nRefineLevel,
	const Edge & edge0,
	const Edge & edge1,
	const Edge & edge2,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	int i;
	int j;
	int k;

	int ixEndNode;

	int ixInt;

	// Edges
	Edge edgeBot;
	Edge edgeTop;

	// Initial bottom edge
	edgeBot.push_back(edge0[0]);

	// Loop over all refined faces
	for (j = 0; j <= nRefineLevel; j++) {

		// Generate top level vertices
		if (j == nRefineLevel) {
			edgeTop = edge2;
		} else {
			GenerateEdgeVertices(
				j+1, edge0[j+1], edge1[j+1], vecNodes, edgeTop);
		}

		// Generate quads
		printf("%i %i\n", edgeBot.size(), edgeTop.size());
		for (i = 2; i < 2*j; i += 2) {
			vecFaces.push_back(Face(
				edgeTop[i], edgeTop[i+1], edgeBot[i/2+1], edgeBot[i/2]));
		}

		// New bottom edge
		edgeBot = edgeTop;
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateFacesFromTriangle(
	int nRefineLevel,
	const Edge & edge0,
	const Edge & edge1,
	const Edge & edge2,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	int i;
	int j;
	int k;

	int ixEndNode;

	int ixInt;

	// Edges
	Edge edgeBot;
	Edge edgeMid;
	Edge edgeTop;

	// Initial bottom edge
	edgeBot.push_back(edge0[0]);

	// Loop over all refined faces
	for (j = 0; j < nRefineLevel; j++) {

		// Generate mid level vertices
		GenerateEdgeVertices(
			2*j+1, edge0[2*j+1], edge1[2*j+1], vecNodes, edgeMid);

		// Generate top level vertices
		if (j == nRefineLevel-1) {
			edgeTop = edge2;
		} else {
			GenerateEdgeVertices(
				2*j+2, edge0[2*j+2], edge1[2*j+2], vecNodes, edgeTop);
		}

		// Generate faces
		for (i = 0; i < 2*j+1; i++) {
			// Downward pointing faces
			if (i % 2 == 0) {
				ixInt = InsertTriFaceCentroidNode(
					edgeMid[i], edgeMid[i+1], edgeTop[i+1], vecNodes);

				vecFaces.push_back(Face(
					edgeBot[i], edgeMid[i], ixInt, edgeMid[i+1]));
				vecFaces.push_back(Face(
					edgeTop[i], edgeTop[i+1], ixInt, edgeMid[i]));
				vecFaces.push_back(Face(
					edgeTop[i+2], edgeMid[i+1], ixInt, edgeTop[i+1]));

			// Upward pointing faces
			} else {

				ixInt = InsertTriFaceCentroidNode(
					edgeMid[i], edgeMid[i+1], edgeBot[i], vecNodes);

				vecFaces.push_back(Face(
					edgeBot[i-1], edgeMid[i], ixInt, edgeBot[i]));
				vecFaces.push_back(Face(
					edgeTop[i+1], edgeMid[i+1], ixInt, edgeMid[i]));
				vecFaces.push_back(Face(
					edgeBot[i+1], edgeBot[i], ixInt, edgeMid[i+1]));
			}
		}

		// New bottom edge
		edgeBot = edgeTop;
	}
}

///////////////////////////////////////////////////////////////////////////////

void ConvertFromLonLatToCartesian(
	const LonLatNodeVector & vecLonLatNodes,
	NodeVector & vecNodes
) {
	vecNodes.resize(vecLonLatNodes.size());

	// Loop over all nodes
	int i;
	for (i = 0; i < vecLonLatNodes.size(); i++) {
		vecNodes[i].x =
			sin(vecLonLatNodes[i].lon) * cos(vecLonLatNodes[i].lat);
		vecNodes[i].y =
			cos(vecLonLatNodes[i].lon) * cos(vecLonLatNodes[i].lat);
		vecNodes[i].z =
			sin(vecLonLatNodes[i].lat);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateIcosahedralQuadGrid(
	int nRefineLevel,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	// Latitude of nodes (Northern Hemisphere)
	const double NodeLat = atan(0.5);

	// Store all icosahedral nodes
	LonLatNodeVector vecLonLatNodes;

	vecLonLatNodes.push_back(LonLatNode(0.0,          -0.5*M_PI));
	vecLonLatNodes.push_back(LonLatNode(0.0,          -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.2, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.4, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.6, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.8, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.1, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.3, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.5, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.7, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.9, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(0.0,          +0.5*M_PI));

	// Convert icosahedral nodes to Cartesian geometry
	ConvertFromLonLatToCartesian(vecLonLatNodes, vecNodes);

	// Vector of edges
	EdgeVector vecEdges;
	vecEdges.resize(30);

	// Generate vertices along edges
	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, 0, i+1, vecNodes, vecEdges[i]);
	}

	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, i+1, ((i+1)%5)+1, vecNodes, vecEdges[i+5]);
	}

	GenerateEdgeVertices(2*nRefineLevel, 1, 6, vecNodes, vecEdges[10]);
	GenerateEdgeVertices(2*nRefineLevel, 6, 2, vecNodes, vecEdges[11]);
	GenerateEdgeVertices(2*nRefineLevel, 2, 7, vecNodes, vecEdges[12]);
	GenerateEdgeVertices(2*nRefineLevel, 7, 3, vecNodes, vecEdges[13]);
	GenerateEdgeVertices(2*nRefineLevel, 3, 8, vecNodes, vecEdges[14]);
	GenerateEdgeVertices(2*nRefineLevel, 8, 4, vecNodes, vecEdges[15]);
	GenerateEdgeVertices(2*nRefineLevel, 4, 9, vecNodes, vecEdges[16]);
	GenerateEdgeVertices(2*nRefineLevel, 9, 5, vecNodes, vecEdges[17]);
	GenerateEdgeVertices(2*nRefineLevel, 5, 10, vecNodes, vecEdges[18]);
	GenerateEdgeVertices(2*nRefineLevel, 10, 1, vecNodes, vecEdges[19]);

	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, i+6, ((i+1)%5)+6, vecNodes, vecEdges[i+20]);
	}

	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, i+6, 11, vecNodes, vecEdges[i+25]);
	}

	// Generate south polar faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i],
			vecEdges[(i+1)%5],
			vecEdges[i+5],
			vecNodes,
			vecFaces
		);
	}

	// Generate south equatorial faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[2*i+10],
			vecEdges[i+5],
			vecEdges[2*i+11],
			vecNodes,
			vecFaces
		);
	}

	// Generate north equatorial faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i+20],
			vecEdges[2*i+11],
			vecEdges[2*((i+1)%5)+10].Flip(),
			vecNodes,
			vecFaces
		);
	}

	// Generate north polar faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i+25],
			vecEdges[i+20],
			vecEdges[((i+1)%5)+25].Flip(),
			vecNodes,
			vecFaces
		);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateOctahedralQuadGrid1(
	int nRefineLevel,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	// Store all octahedral nodes
	LonLatNodeVector vecLonLatNodes;

	vecLonLatNodes.push_back(LonLatNode(0.0,      -0.5*M_PI));
	vecLonLatNodes.push_back(LonLatNode(0.0,       0.0));
	vecLonLatNodes.push_back(LonLatNode(0.5*M_PI,  0.0));
	vecLonLatNodes.push_back(LonLatNode(1.0*M_PI,  0.0));
	vecLonLatNodes.push_back(LonLatNode(1.5*M_PI,  0.0));
	vecLonLatNodes.push_back(LonLatNode(0.0,      +0.5*M_PI));

	// Convert octahedral nodes to Cartesian geometry
	ConvertFromLonLatToCartesian(vecLonLatNodes, vecNodes);

	// Vector of edges
	EdgeVector vecEdges;
	vecEdges.resize(12);

	// Generate vertices along edges
	for (int i = 0; i < 4; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, 0, i+1, vecNodes, vecEdges[i]);
	}

	for (int i = 0; i < 4; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, i+1, ((i+1)%4)+1, vecNodes, vecEdges[i+4]);
	}

	for (int i = 0; i < 4; i++) {
		GenerateEdgeVertices(
			2*nRefineLevel, i+1, 5, vecNodes, vecEdges[i+8]);
	}

	// Generate south polar faces
	for (int i = 0; i < 4; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i],
			vecEdges[(i+1)%4],
			vecEdges[i+4],
			vecNodes,
			vecFaces
		);
	}

	// Generate north polar faces
	for (int i = 0; i < 4; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i+8],
			vecEdges[i+4],
			vecEdges[((i+1)%4)+8].Flip(),
			vecNodes,
			vecFaces
		);
	}

}

///////////////////////////////////////////////////////////////////////////////

void GenerateOctahedralQuadGrid2(
	int nRefineLevel,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	// Store refine level as Nr
	const int nR = nRefineLevel;

	// Store all icosahedral nodes
	LonLatNodeVector vecLonLatNodes;

	double dDeltaLat = 0.5 * M_PI / static_cast<double>(nR);

	// South polar nodes
	for (int k = 0; k < nR; k++) {
		int nLevelNodes = 4 * (k+1);

		double dLonOffset = M_PI/static_cast<double>(nLevelNodes)*k;

		for (int i = 0; i < nLevelNodes; i++) {
			vecLonLatNodes.push_back(
				LonLatNode(
					2.0*M_PI/nLevelNodes*i + dLonOffset,
					-0.5*M_PI + dDeltaLat*(k+1)));
		}
	}

	// North polar nodes
	for (int k = nR-2; k >= 0; k--) {
		int nLevelNodes = 4 * (k+1);

		double dLonOffset = M_PI/static_cast<double>(nLevelNodes)*k;

		for (int i = 0; i < nLevelNodes; i++) {
			vecLonLatNodes.push_back(
				LonLatNode(
					2.0*M_PI/nLevelNodes*i + dLonOffset,
					0.5*M_PI - dDeltaLat*(k+1)));
		}
	}

	// Convert icosahedral nodes to Cartesian geometry
	ConvertFromLonLatToCartesian(vecLonLatNodes, vecNodes);

	// South polar face
	vecFaces.push_back(Face(0, 1, 2, 3));

	// Southern hemisphere faces
	for (int k = 0; k < nR-1; k++) {
		for (int l = 0; l < 4; l++) {

			// Square elements
			{
				int i0 = 2*k*(k+1) + l*(k+1);
				int i1 = 2*(k+1)*(k+2) + l*(k+2);
				int i2 = 2*(k+1)*(k+2)+1 + l*(k+2);
				int i3 = 2*k*(k+1)+1 + l*(k+1);

				if ((k == 0) && (l == 3)) {
					i3 = 0;
				}

				vecFaces.push_back(Face(i0, i1, i2, i3));
			}

			// Diamond elements
			if (k != nRefineLevel-2) {
				for (int m = 0; m <= k; m++) {
					int i0 = 2*k*(k+1)+1 + l*(k+1) + m;
					int i1 = 2*(k+1)*(k+2)+1 + l*(k+2) + m;
					int i2 = 2*(k+2)*(k+3)+1 + l*(k+3) + m + 1;
					int i3 = 2*(k+1)*(k+2)+1 + l*(k+2) + m + 1;

					if ((l == 3) && (m == k)) {
						i0 -= 4*(k+1);
						i3 -= 4*(k+2);
					}

					vecFaces.push_back(Face(i0, i1, i2, i3));
				}
			}
		}
	}

	// Northern hemisphere faces
	int iNodeEq = 2 * nR * (nR+1) - 4 * nR;

	for (int k = 0; k < nR-1; k++) {
		for (int l = 0; l < 4; l++) {

			// Square elements
			{
				int i0 = iNodeEq + l*(nR-k);
				int i1 = iNodeEq + 4*(nR-k) + l*(nR-(k+1));
				int i2 = iNodeEq + 4*(nR-k) + l*(nR-(k+1)) + 1;
				int i3 = iNodeEq + l*(nR-k) + 1;

				if ((k == nR-2) && (l == 3)) {
					i2 -= 4*(nR-(k+1));
				}

				vecFaces.push_back(Face(i0, i1, i2, i3));
			}


			// Diamond elements
			if (k != nR-2) {
				for (int m = 1; m < nR-k-1; m++) {
					int i0 = iNodeEq + l*(nR-k) + m + 1;
					int i1 = iNodeEq + 4*(nR-k) + l*(nR-(k+1)) + m;
					int i2 = iNodeEq + 4*(nR-k) + 4*(nR-(k+1)) + l*(nR-(k+2)) + m;
					int i3 = iNodeEq + 4*(nR-k) + l*(nR-(k+1)) + m + 1;

					if ((m == nR-k-2) && (l == 3)) {
						i2 -= 4*(nR-(k+2));
						i3 -= 4*(nR-(k+1));
					}

					vecFaces.push_back(Face(i0, i1, i2, i3));
				}
			}
		}
		iNodeEq += 4*(nRefineLevel-k);
	}

	// South polar face
	int iNodeCount = static_cast<int>(vecNodes.size());
	vecFaces.push_back(
		Face(iNodeCount-4, iNodeCount-3, iNodeCount-2, iNodeCount-1));
}

///////////////////////////////////////////////////////////////////////////////

