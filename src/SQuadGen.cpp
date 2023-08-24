///////////////////////////////////////////////////////////////////////////////
///
///	\file    SQuadGen.cpp
///	\author  Paul Ullrich
///	\version February 13, 2012
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

#include <png.h>
#include "netcdfcpp.h"

#include <string>
#include <cmath>
#include <cstdlib>

#include "GridElements.h"
#include "CubedSphereGrid.h"
#include "IcosahedralFlagGrid.h"
#include "CSRefinementMap.h"
#include "RefineGrid.h"
#include "RefinementTemplateCUBIT.h"
#include "RefinementTemplateLOWCONN.h"
#include "RefinementTemplateLOWCONNOLD.h"
#include "SpringDynamics.h"
#include "Tessellate.h"
#include "MeshUtilities.h"

#include "Exception.h"
#include "CommandLine.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	try {

	// Grid type
	std::string strGridType;

	// Type of refinement
	std::string strRefineType;

	// Initial resolution
	int nResolution;

	// Number of levels of refinement
	int nRefinementLevel;

	// Refinement file
	std::string strRefineFile;

	// Output file
	std::string strOutputFile;

	// Refinement map
	bool fLoadCSRefinementMap;

	// Type of smoothing
	std::string strSmoothType;

	// Transition region smoothing distance
	int nTransSmoothDist;

	// Number of iterations for smoothing
	int nSmoothIterations;

	// Reverse orientation
	bool fReverseOrientation = false;

	// Invert the image
	bool fInvertImage;

	// Image longitude shift
	double dImageLonBase;

	// Image latitude shift
	double dImageLatBase;

	// Grid rotation about the X axis
	double dGridXRotate;

	// Grid rotation about the Y axis
	double dGridYRotate;

	// Reference longitude
	double dReferenceLonDeg;

	// Reference latitude
	double dReferenceLatDeg;

	// Reference orientation
	double dReferenceOrientDeg;

	// Number of tesselations
	int nTessellations;

	// Sub-cell resolution
	int nSubCellResolution;

	// Rectangular refinement commands
	std::string strRefineRect;

	// Block refine
	bool fBlockRefine;

	// Parse the command line
	BeginCommandLine()
		CommandLineStringD(strGridType, "grid_type", "CS",
			"(Options: CS | ICO | OCT1 | OCT2)");
		CommandLineStringD(strRefineType, "refine_type", "LOWCONN",
			"(Options: LOWCONN | CUBIT | LOWCONNOLD)");
		CommandLineInt(nRefinementLevel, "refine_level", 2);
		CommandLineInt(nResolution, "resolution", 10);
		CommandLineString(strRefineFile, "refine_file", "");
		CommandLineString(strRefineRect, "refine_rect", "");
		CommandLineString(strOutputFile, "output", "");
		CommandLineBool(fLoadCSRefinementMap, "loadcsrefinementmap");
		CommandLineStringD(strSmoothType, "smooth_type", "NONE",
			"(Options: NONE | SPRING | PRESSURE)");
		CommandLineIntD(nTransSmoothDist, "smooth_dist", 1,
			"(Smooth distance, -1 = smooth entire mesh)");
		CommandLineInt(nSmoothIterations, "smooth_iter", 10);
		CommandLineDouble(dImageLonBase, "lon_base", -180.0);
		CommandLineDouble(dImageLatBase, "lat_base", 0.0);
		CommandLineDouble(dGridXRotate, "x_rotate", 0.0);
		CommandLineDouble(dGridYRotate, "y_rotate", 0.0);
		CommandLineDouble(dReferenceLonDeg, "lon_ref", 0.0);
		CommandLineDouble(dReferenceLatDeg, "lat_ref", 0.0);
		CommandLineDouble(dReferenceOrientDeg, "orient_ref", 0.0);
		CommandLineInt(nTessellations, "tessellate", 0);
		CommandLineInt(nSubCellResolution, "subcellres", 0);
		CommandLineBool(fInvertImage, "invert");
		CommandLineBool(fBlockRefine, "block_refine");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	printf("----------------------------------------\n");

	// Check for output file
	if (strOutputFile == "") {
		std::cout << argv[0] << ": No output file specified" << std::endl;
		return (-1);
	}

	// Check for odd resolution
	if ((strRefineFile != "") || (strRefineRect != "")) {
		if ((strRefineType == "LOWCONN") || (strRefineType == "LOWCONNOLD")) {
			if ((nResolution % 2) == 1) {
				_EXCEPTIONT("\nERROR: "
					"Only even resolutions are supported by LOWCONN template.");
			}
			nResolution /= 2;
		}
	}

	// Uppercase the grid type
	for (int i = 0; i < strGridType.length(); i++) {
		strGridType[i] = toupper(strGridType[i]);
	}

	// Uppercase the refine type
	for (int i = 0; i < strRefineType.length(); i++) {
		strRefineType[i] = toupper(strRefineType[i]);
	}

	// Nodes of the grid
	NodeVector vecNodes;
	FaceVector vecFaces;

	// Generate grid
	printf("Generating mesh\n");
	if (strGridType == "CS") {
		GenerateCubedSphere(nResolution, vecNodes, vecFaces);
	} else if (strGridType == "ICO") {
		GenerateIcosahedralQuadGrid(nResolution, vecNodes, vecFaces);
	} else if (strGridType == "OCT1") {
		GenerateOctahedralQuadGrid1(nResolution, vecNodes, vecFaces);
	} else if (strGridType == "OCT2") {
		GenerateOctahedralQuadGrid2(nResolution, vecNodes, vecFaces);
	} else {
		_EXCEPTION1("Unknown grid type: %s\n", strGridType.c_str());
	}
	printf("..Base mesh complete\n");

	// Perform refinement (if refinement file specified)
	if ((nRefinementLevel > 0) &&
		((strRefineRect != "") || (strRefineFile != "") || (fLoadCSRefinementMap))
	) {
		int nCommandLineCount = 
			  ((strRefineRect != "")?(1):(0))
			+ ((strRefineFile != "")?(1):(0))
			+ (fLoadCSRefinementMap?(1):(0));

		if (nCommandLineCount > 1) {
			_EXCEPTIONT("Only one of --refine_map, --refine_rect or --loadcsrefinementmap allowed");
		}

		CSRefinementMap refmap(nResolution, nRefinementLevel);

		if (fLoadCSRefinementMap) {
			printf("..Loading mesh from refinement file\n");
			refmap.FromFile("refine_map.dat");

		} else if (strRefineRect != "") {
			refmap.InitializeFromRefineRect(
				strRefineRect,
				dReferenceLonDeg,
				dReferenceLatDeg,
				dReferenceOrientDeg);

			refmap.Normalize();
			refmap.ToFile("refine_map.dat");

		} else {
			refmap.InitializeFromPNG(
				strRefineFile,
				dImageLonBase,
				dImageLatBase,
				fInvertImage,
				fBlockRefine,
				dReferenceLonDeg,
				dReferenceLatDeg,
				dReferenceOrientDeg);

			refmap.Normalize();
			refmap.ToFile("refine_map.dat");
		}

		// CUBIT refinement template
		if (strRefineType == "CUBIT") {
			RefinementTemplateCUBIT reftempCUBIT;
			RefineGrid(vecNodes, vecFaces, reftempCUBIT, refmap);

		// LOWCONN refinement template
		} else if (strRefineType == "LOWCONN") {
			RefinementTemplateLOWCONN reftempLOWCONN;
			RefineGrid(vecNodes, vecFaces, reftempLOWCONN, refmap);

		// LOWCONNOLD refinement template
		} else if (strRefineType == "LOWCONNOLD") {
			RefinementTemplateLOWCONNOLD reftempLOWCONNOLD;
			RefineGrid(vecNodes, vecFaces, reftempLOWCONNOLD, refmap);

		} else {
			_EXCEPTIONT("Invalid refinement type");
		}
	}

	// Tessellate
	if ((nTessellations < 0) || (nTessellations > 100)) {
		_EXCEPTIONT("Invalid number of tesselations; expected [0,100]");
	} else if (nTessellations == 0) {
		// Do nothing
	} else {
		for (int n = 0; n < nTessellations; n++) {
			Tessellate(vecNodes, vecFaces);
		}
	}

	// Add sub-cell resolution
	if (nSubCellResolution < 0) {
		_EXCEPTIONT("--subcellres must be >= 0");
	} else if (nSubCellResolution == 0) {
		// Do nothing
	} else {
		RefineEverything(vecNodes, vecFaces, nSubCellResolution+1);
	}

	// Reverse orientation
	if ((fReverseOrientation) && (strGridType != "CS")) {
		ReverseFaceOrientation(vecFaces);
	} else if ((!fReverseOrientation) && (strGridType == "CS")) {
		ReverseFaceOrientation(vecFaces);
	}

	// No smoothing
	if (strSmoothType == "NONE") {

	// Perform smoothing
	} else if (strSmoothType == "SPRING" ) {
		SpringDynamics(
			vecNodes,
			vecFaces,
			nTransSmoothDist,
			nSmoothIterations);

	// Pressure dynamics
	} else if (strSmoothType == "PRESSURE") {
		PressureDynamics(
			vecNodes,
			vecFaces,
			nTransSmoothDist,
			nSmoothIterations);

	} else {
		_EXCEPTIONT("Invalid smoothing type");
	}

	printf("Mesh refinement complete!\n");


	// Rotate around the X axis
	if (dGridXRotate != 0.0) {
	        printf("Rotating grid around X axis.\n");
		double dCosTheta = cos(dGridXRotate * M_PI / 180.0);
		double dSinTheta = sin(dGridXRotate * M_PI / 180.0);

		for (int i = 0; i < vecNodes.size(); i++) {
			double dTempY = vecNodes[i].y;
			double dTempZ = vecNodes[i].z;

			vecNodes[i].y = dCosTheta * dTempY - dSinTheta * dTempZ;
			vecNodes[i].z = dSinTheta * dTempY + dCosTheta * dTempZ;
		}
	}

	// Rotate around the Y axis
	if (dGridYRotate != 0.0) {
	        printf("Rotating grid around Y axis.\n");
		double dCosTheta = cos(dGridYRotate * M_PI / 180.0);
		double dSinTheta = sin(dGridYRotate * M_PI / 180.0);

		for (int i = 0; i < vecNodes.size(); i++) {
			double dTempX = vecNodes[i].x;
			double dTempZ = vecNodes[i].z;

			vecNodes[i].x = dCosTheta * dTempX - dSinTheta * dTempZ;
			vecNodes[i].z = dSinTheta * dTempX + dCosTheta * dTempZ;
		}
	}

	// Rotate nodes in 3D to match reference point
	RotateMeshByOrientLatLon(vecNodes, dReferenceLonDeg, dReferenceLatDeg, dReferenceOrientDeg);

	// Output number of nodes and faces
	printf("----------------------------------------\n");
	printf("Node Count:  %lu\n", vecNodes.size());
	printf("Face Count:  %lu\n", vecFaces.size());
	printf("----------------------------------------\n");

	// Compute mesh quality
	ComputeMeshQuality(vecNodes, vecFaces);
/*
	// Output the nodes and connectivity
	FILE * fpNodes = fopen("nodes.dat", "w");
	for (int i = 0; i < vecNodes.size(); i++) {
		fprintf(fpNodes, "%1.5e %1.5e %1.5e\n",
			vecNodes[i].x, vecNodes[i].y, vecNodes[i].z);
	}
	fclose(fpNodes);

	FILE * fpFaces = fopen("faces.dat", "w");
	for (int i = 0; i < vecFaces.size(); i++) {
		fprintf(fpFaces, "%i %i %i %i %i\n",
			vecFaces[i][0],
			vecFaces[i][1],
			vecFaces[i][2],
			vecFaces[i][3],
			vecFaces[i].nRefineLevel);
	}
	fclose(fpFaces);
*/
	// Write Exodus file
	OutputNetCDFFile(strOutputFile.c_str(), vecNodes, vecFaces);

	// Success
	std::cout << "SQuadGen completed successfully." << std::endl;

	} catch (Exception e) {
		std::cout << e.ToString() << std::endl;
		std::cout << "SQuadGen failed" << std::endl;
	}

	return (0);
}

///////////////////////////////////////////////////////////////////////////////

