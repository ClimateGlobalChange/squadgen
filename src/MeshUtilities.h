///////////////////////////////////////////////////////////////////////////////
///
///	\file    MeshUtilities.h
///	\author  Paul Ullrich
///	\version February 8, 2022
///
///	<remarks>
///		Copyright 2022 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _MESHUTILITIES_H_
#define _MESHUTILITIES_H_

#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality(
	const NodeVector & vecNodes,
	const FaceVector & vecFaces
);

///////////////////////////////////////////////////////////////////////////////

void UnrotateCoordByOrientLatLon(
	double dRotLonDeg,
	double dRotLatDeg,
	double dRotOrientDeg,
	double & dX,
	double & dY,
	double & dZ
);

///////////////////////////////////////////////////////////////////////////////

void RotateMeshByOrientLatLon(
	NodeVector & vecNodes,
	double dRotLonDeg,
	double dRotLatDeg,
	double dRotOrientDeg
);

///////////////////////////////////////////////////////////////////////////////

void OutputNetCDFFile(
	std::string strFilename,
	const NodeVector & vecNodes,
	const FaceVector & vecFaces
);

///////////////////////////////////////////////////////////////////////////////

#endif // _MESHUTILITIES_H_
