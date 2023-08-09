///////////////////////////////////////////////////////////////////////////////
///
///	\file    CSRefinementMap.h
///	\author  Paul Ullrich
///	\version March 12, 2012
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

#include "CSRefinementMap.h"

#include "Exception.h"
#include "MathHelper.h"
#include "CubedSphereTrans.h"
#include "CoordTransforms.h"
#include "MeshUtilities.h"

#include "lodepng.h"

///////////////////////////////////////////////////////////////////////////////

template <typename T>
T clamp(const T& n, const T& lower, const T& upper) {
  return std::max(lower, std::min(n, upper));
}

///////////////////////////////////////////////////////////////////////////////

CSRefinementMap::CSRefinementMap(
	int nBaseResolution,
	int nMaxRefineLevel
) :
	m_nBaseResolution(nBaseResolution),
	m_nMaxRefineLevel(nMaxRefineLevel)
{
	int nMaxResolution = nBaseResolution * IntPow(2, nMaxRefineLevel-1);
	m_nMap.Initialize(6, nMaxResolution, nMaxResolution);
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::SetRefineLevel(
	int iPanel,
	int iA,
	int iB,
	int nRefineLevel
) {
	if (nRefineLevel > m_nMaxRefineLevel) {
		_EXCEPTIONT("RefineLevel exceed maximum refinement level.");
	}

	m_nMap[iPanel][iA][iB] = nRefineLevel;
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::SetMinimumRefineLevel(
	int iPanel,
	int iA,
	int iB,
	int iRefineLevel,
	int nActiveFineRatio
) {
	int iABegin = iA * nActiveFineRatio;
	int iBBegin = iB * nActiveFineRatio;
	int iAEnd   = (iA+1) * nActiveFineRatio;
	int iBEnd   = (iB+1) * nActiveFineRatio;

	for (int i = iABegin; i < iAEnd; i++) {
	for (int j = iBBegin; j < iBEnd; j++) {
		if (m_nMap[iPanel][i][j] < iRefineLevel) {
			m_nMap[iPanel][i][j] = iRefineLevel;
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::InitializeFromRefineRect(
	const std::string & strRefineRect,
	double dReferenceLonDeg,
	double dReferenceLatDeg,
	double dReferenceOrientDeg
) {
	int nActivePanel = (-1);
	size_t pos;
	std::string strRefineRectCopy(strRefineRect);
	while (strRefineRectCopy.length() != 0) {

		// Break up by semicolon
		pos = strRefineRectCopy.find(";");
		std::string strActiveRect(strRefineRectCopy.substr(0, pos));
		if (pos == std::string::npos) {
			strRefineRectCopy = "";
		} else {
			strRefineRectCopy = strRefineRectCopy.substr(pos+1);
		}

		// Break up by commas
		int iStage = 0;

		double dLonDeg0;
		double dLonDeg1;
		double dLatDeg0;
		double dLatDeg1;
		int nRefineLevel;

		while ((pos = strActiveRect.find(",")) != std::string::npos) {
			if (iStage == 0) {
				dLonDeg0 = std::stod(strActiveRect.substr(0,pos));
			} else if (iStage == 1) {
				dLatDeg0 = std::stod(strActiveRect.substr(0,pos));
			} else if (iStage == 2) {
				dLonDeg1 = std::stod(strActiveRect.substr(0,pos));
			} else if (iStage == 3) {
				dLatDeg1 = std::stod(strActiveRect.substr(0,pos));
			}
			iStage++;
			strActiveRect.erase(0, pos+1);
		}

		if (iStage != 4) {
			_EXCEPTION1("Malformed --refine_rect string: Insufficient entries in \"%s\"", strRefineRect.c_str());
		}
		nRefineLevel = std::stoi(strActiveRect.substr(0,pos));

		if (nRefineLevel > m_nMaxRefineLevel) {
			_EXCEPTION1("Malformed --refine_rect string: Refinement level cannot exceed --refine_level %i\n", m_nMaxRefineLevel);
		}
		if (nRefineLevel < 1) {
			_EXCEPTION1("Malformed --refine_rect string: Refinement level must be greater than 0\n", m_nMaxRefineLevel);
		}

		//printf("%1.5f %1.5f %1.5f %1.5f %i\n", dLonDeg0, dLonDeg1, dLatDeg0, dLatDeg1, nRefineLevel);

		// Convert coordinates to ABP
		double dX0, dY0, dZ0;
		double dX1, dY1, dZ1;
		RLLtoXYZ_Deg(dLonDeg0, dLatDeg0, dX0, dY0, dZ0);
		RLLtoXYZ_Deg(dLonDeg1, dLatDeg1, dX1, dY1, dZ1);

		UnrotateCoordByOrientLatLon(
			-dReferenceLonDeg,
			-dReferenceLatDeg,
			-dReferenceOrientDeg,
			dX0, dY0, dZ0);

		UnrotateCoordByOrientLatLon(
			-dReferenceLonDeg,
			-dReferenceLatDeg,
			-dReferenceOrientDeg,
			dX1, dY1, dZ1);

		double dUnrotLonRad0, dUnrotLatRad0;
		double dUnrotLonRad1, dUnrotLatRad1;

		XYZtoRLL_Rad(dX0, dY0, dZ0, dUnrotLonRad0, dUnrotLatRad0);
		XYZtoRLL_Rad(dX1, dY1, dZ1, dUnrotLonRad1, dUnrotLatRad1);

		//printf("%1.5f %1.5f %1.5f %1.5f\n",
		//	RadToDeg(dUnrotLonRad0),
		//	RadToDeg(dUnrotLonRad1),
		//	RadToDeg(dUnrotLatRad0),
		//	RadToDeg(dUnrotLatRad1));

		double dA0, dB0;
		double dA1, dB1;
		int nP0, nP1;

		CubedSphereTrans::ABPFromRLL(dUnrotLonRad0, dUnrotLatRad0, dA0, dB0, nP0);
		CubedSphereTrans::ABPFromRLL(dUnrotLonRad1, dUnrotLatRad1, dA1, dB1, nP1);

		if (nActivePanel == (-1)) {
			nActivePanel = nP0;
		}
		if ((nP0 != nP1) || (nP0 != nActivePanel)) {
			_EXCEPTIONT("At present all --refine_rect coordinates must be on the same cubed sphere panel");
		}

		// Convert equiangular ABP coordinate to [0,1]x[0,1]
		dA0 = (dA0 + 0.25 * M_PI) / (0.5 * M_PI);
		dA1 = (dA1 + 0.25 * M_PI) / (0.5 * M_PI);
		dB0 = (dB0 + 0.25 * M_PI) / (0.5 * M_PI);
		dB1 = (dB1 + 0.25 * M_PI) / (0.5 * M_PI);

		dA0 = clamp(dA0, 0.0, 1.0);
		dA1 = clamp(dA1, 0.0, 1.0);
		dB0 = clamp(dB0, 0.0, 1.0);
		dB1 = clamp(dB1, 0.0, 1.0);

		if (dA1 < dA0) {
			std::swap(dA0, dA1);
		}
		if (dB1 < dB0) {
			std::swap(dB0, dB1);
		}

		// Refinement level
		int nPatchSize = IntPow(2, m_nMaxRefineLevel-nRefineLevel);
		int nActiveResolution = m_nBaseResolution * IntPow(2, nRefineLevel-1);

		// Identify integer refinement region
		int iA0 = static_cast<int>(dA0 * static_cast<double>(nActiveResolution));
		int iA1 = static_cast<int>(dA1 * static_cast<double>(nActiveResolution));
		if (iA0 > nActiveResolution-1) {
			iA0 = nActiveResolution-1;
		}
		if (iA1 > nActiveResolution-1) {
			iA1 = nActiveResolution-1;
		}
		if (iA1 < nActiveResolution-1) {
			iA1++;
		}

		int iB0 = static_cast<int>(dB0 * static_cast<double>(nActiveResolution));
		int iB1 = static_cast<int>(dB1 * static_cast<double>(nActiveResolution));
		if (iB0 > nActiveResolution-1) {
			iB0 = nActiveResolution-1;
		}
		if (iB1 > nActiveResolution-1) {
			iB1 = nActiveResolution-1;
		}
		if (iB1 < nActiveResolution-1) {
			iB1++;
		}

		//printf("%1.5f %1.5f %i : %1.5f %1.5f %i\n", dA0, dB0, nP0, dA1, dB1, nP1);
		//printf("%i %i %i : %i %i %i : %i\n", iA0, iB0, nP0, iA1, iB1, nP1, nPatchSize);
		//_EXCEPTION();

		int nMaxResolution = m_nBaseResolution * IntPow(2, m_nMaxRefineLevel-1);
		_ASSERT(iA0*nPatchSize < nMaxResolution);
		_ASSERT(iA1*nPatchSize < nMaxResolution);
		_ASSERT(iB0*nPatchSize < nMaxResolution);
		_ASSERT(iB1*nPatchSize < nMaxResolution);

		// Refine
		for (int i = iA0*nPatchSize; i <= iA1*nPatchSize; i++) {
		for (int j = iB0*nPatchSize; j <= iB1*nPatchSize; j++) {
			m_nMap[nActivePanel][i][j] = nRefineLevel;
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::InitializeFromPNG(
	const std::string & strPNGFile,
	double dLonBase,
	double dLatBase,
	bool fInvert,
	bool fBlockRefine,
	double dReferenceLonDeg,
	double dReferenceLatDeg,
	double dReferenceOrientDeg
) {
	std::vector<unsigned char> image;
	unsigned int nWidth;
	unsigned int nHeight;

	unsigned int iError =
		lodepng::decode(image, nWidth, nHeight, strPNGFile.c_str());

	if (iError) {
		_EXCEPTION2("\n  Decoder error %i: %s",
			iError, lodepng_error_text(iError));
	}

	printf("----------------------------------------\n");
	printf("PNG loaded successfully\n");
	printf("... Dimensions: %i x %i\n", nWidth, nHeight);

	// Convert arguments to radians
	dLonBase *= M_PI / 180.0;
	dLatBase *= M_PI / 180.0;

	double dReferenceLonRad = dReferenceLonDeg * M_PI / 180.0;
	double dReferenceLatRad = dReferenceLatDeg * M_PI / 180.0;
	double dReferenceOrientRad = dReferenceOrientDeg * M_PI / 180.0;

	// Calculate image greyscale levels
	printf("... Levels: ");
	int iLastLevel = (-1);
	for (int iC = 0; iC < 255; iC++) { 

		double dC = static_cast<double>(iC) / 256.0;
		int iCurrentLevel = floor(dC * (m_nMaxRefineLevel + 1));

		if (iCurrentLevel != iLastLevel) {
			printf("%i ", iC);
			iLastLevel = iCurrentLevel;
		}
	}
	printf("\n");

	// Loop throught all cubed-sphere faces
	for (int iP = 0; iP < 6; iP++) {
	for (int iA = 0; iA < m_nMap.GetColumns(); iA++) {
	for (int iB = 0; iB < m_nMap.GetColumns(); iB++) {
		int iMaxRefineLevel = 0;

		for (int iS = 0; iS < 2; iS++) {
		for (int iT = 0; iT < 2; iT++) {
			double dA = -0.25 * M_PI + 0.5 * M_PI
				* ((static_cast<double>(iA + iS)) / m_nMap.GetColumns());
			double dB = -0.25 * M_PI + 0.5 * M_PI
				* ((static_cast<double>(iB + iT)) / m_nMap.GetColumns());

			double dLonRad;
			double dLatRad;

			CubedSphereTrans::RLLFromABP(dA, dB, iP, dLonRad, dLatRad);

			double dX0;
			double dY0;
			double dZ0;

			RLLtoXYZ_Rad(dLonRad, dLatRad, dX0, dY0, dZ0);

			// Rotate around X axis by +dReferenceOrientRad
			double dX1 = dX0;
			double dY1 = + cos(dReferenceOrientRad) * dY0 - sin(dReferenceOrientRad) * dZ0;
			double dZ1 = + sin(dReferenceOrientRad) * dY0 + cos(dReferenceOrientRad) * dZ0;

			// Rotate around Y axis by -dReferenceLatRad
			double dX2 = + cos(dReferenceLatRad) * dX1 - sin(dReferenceLatRad) * dZ1;
			double dY2 = dY1;
			double dZ2  = + sin(dReferenceLatRad) * dX1 + cos(dReferenceLatRad) * dZ1;

			// Rotate around Z axis by +dReferenceLonRad
			double dX3 = + cos(dReferenceLonRad) * dX2 - sin(dReferenceLonRad) * dY2;
			double dY3 = + sin(dReferenceLonRad) * dX2 + cos(dReferenceLonRad) * dY2;
			double dZ3 = dZ2;

			XYZtoRLL_Rad(dX3, dY3, dZ3, dLonRad, dLatRad);

			// Sample image
			dLatRad += dLatBase;
			dLonRad += dLonBase;
			dLonRad = dLonRad - 2.0 * M_PI * floor(dLonRad / (2.0 * M_PI));
			if ((dLonRad < 0.0) || (dLonRad > 2.0 * M_PI)) {
				_EXCEPTIONT("Invalid longitude");
			}

			int iLon = floor((dLonRad / (2.0 * M_PI)) * nWidth);
			int iLat = floor((dLatRad + 0.5 * M_PI) / M_PI * nHeight);

			if (iLon < 0) {
				iLon = 0;
			}
			if (iLat < 0) {
				iLat = 0;
			}
			if (iLon >= nWidth) {
				iLon = nWidth - 1;
			}
			if (iLat >= nHeight) {
				iLat = nHeight - 1;
			}

			int iImageIx = ((nHeight - 1 - iLat) * nWidth + iLon) * 4;

			// Calculate luminance
			double dY =
				+ 0.2126 * static_cast<double>(image[iImageIx]) / 256.0
				+ 0.7152 * static_cast<double>(image[iImageIx+1]) / 256.0
				+ 0.0722 * static_cast<double>(image[iImageIx+2]) / 256.0;

			// Find active fine ratio
			int iRefineLevel =
				floor(dY * static_cast<double>(m_nMaxRefineLevel + 1));

			// Invert the image
			if (fInvert) {
				iRefineLevel = m_nMaxRefineLevel - iRefineLevel;
			}

			// Store the maximum refinement level among all nodes
			if (iRefineLevel > iMaxRefineLevel) {
				iMaxRefineLevel = iRefineLevel;
			}
		}
		}

		if (iMaxRefineLevel != 0) {
			int nActiveFineRatio;
			if (fBlockRefine) {
				nActiveFineRatio =
					IntPow(2, m_nMaxRefineLevel - 1);
			} else {
				nActiveFineRatio =
					IntPow(2, m_nMaxRefineLevel - iMaxRefineLevel);
			}

			SetMinimumRefineLevel(
				iP, iA / nActiveFineRatio, iB / nActiveFineRatio,
				iMaxRefineLevel, nActiveFineRatio);
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::Normalize() {
	int iActiveLevel = m_nMaxRefineLevel;
	int nActiveResolution = m_nMap.GetColumns();

	// Loop over all refinement levels
	while (iActiveLevel > 0) {
		int nActiveFineRatio = m_nMap.GetColumns() / nActiveResolution;

		// Smooth the refinement map at this level
		bool fChanges = true;
		while (fChanges) {
			fChanges = false;

			for (int p = 0; p < 6; p++) {
			for (int i = 0; i < nActiveResolution; i++) {
			for (int j = 0; j < nActiveResolution; j++) {

				int ix = i * nActiveFineRatio;
				int jx = j * nActiveFineRatio;

				if (m_nMap[p][ix][jx] >= iActiveLevel) {
					continue;
				}

				// Find neighbors in the eight cardinal coordinate directions
				int iE_src = ix + nActiveFineRatio;
				int iW_src = ix - nActiveFineRatio;
				int jN_src = jx + nActiveFineRatio;
				int jS_src = jx - nActiveFineRatio;

				int pE, pW, pN, pS;
				int iE, iW, iN, iS;
				int jE, jW, jN, jS;

				bool fNE, fNW, fSE, fSW;
				int pNE, pNW, pSE, pSW;
				int iNE, iNW, iSE, iSW;
				int jNE, jNW, jSE, jSW;

				CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, iE_src, jx, pE, iE, jE);
				CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, iW_src, jx, pW, iW, jW);
				CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, ix, jN_src, pN, iN, jN);
				CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, ix, jS_src, pS, iS, jS);

				fNE = CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, iE_src, jN_src, pNE, iNE, jNE);
				fNW = CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, iW_src, jN_src, pNW, iNW, jNW);
				fSE = CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, iE_src, jS_src, pSE, iSE, jSE);
				fSW = CubedSphereTrans::RelativeCoord(
					m_nMap.GetColumns(), p, iW_src, jS_src, pSW, iSW, jSW);

				// Check for cross-over corners (although these are not
				// necessarily problematic).  A better refinement algorithm
				// could negate the necessity of refining these points.
				if ((m_nMap[pE][iE][jE] >= iActiveLevel) &&
					(m_nMap[pN][iN][jN] >= iActiveLevel) &&
					(fNE) && (m_nMap[pNE][iNE][jNE] < iActiveLevel)
				) {
					SetMinimumRefineLevel(
						p, i, j, iActiveLevel, nActiveFineRatio);
					fChanges = true;
					continue;
				}

				if ((m_nMap[pE][iE][jE] >= iActiveLevel) &&
					(m_nMap[pS][iS][jS] >= iActiveLevel) &&
					(fSE) && (m_nMap[pSE][iSE][jSE] < iActiveLevel)
				) {
					SetMinimumRefineLevel(
						p, i, j, iActiveLevel, nActiveFineRatio);
					fChanges = true;
					continue;
				}

				// Check for exterior corner intersections
				if (fNE && (m_nMap[pNE][iNE][jNE] >= iActiveLevel)) {
					if (((m_nMap[pE][iE][jE] < iActiveLevel) &&
							(m_nMap[pS][iS][jS] >= iActiveLevel)) ||
						((m_nMap[pN][iN][jN] < iActiveLevel) &&
							(m_nMap[pW][iW][jW] >= iActiveLevel))
					) {
						SetMinimumRefineLevel(
							p, i, j, iActiveLevel, nActiveFineRatio);
						fChanges = true;
						continue;
					}
				}

				if (fNW && (m_nMap[pNW][iNW][jNW] >= iActiveLevel)) {
					if (((m_nMap[pW][iW][jW] < iActiveLevel) &&
							(m_nMap[pS][iS][jS] >= iActiveLevel)) ||
						((m_nMap[pN][iN][jN] < iActiveLevel) &&
							(m_nMap[pE][iE][jE] >= iActiveLevel))
					) {
						SetMinimumRefineLevel(
							p, i, j, iActiveLevel, nActiveFineRatio);
						fChanges = true;
						continue;
					}
				}

				if (fSE && (m_nMap[pSE][iSE][jSE] >= iActiveLevel)) {
					if (((m_nMap[pE][iE][jE] < iActiveLevel) &&
							(m_nMap[pN][iN][jN] >= iActiveLevel)) ||
						((m_nMap[pS][iS][jS] < iActiveLevel) &&
							(m_nMap[pW][iW][jW] >= iActiveLevel))
					) {
						SetMinimumRefineLevel(
							p, i, j, iActiveLevel, nActiveFineRatio);
						fChanges = true;
						continue;
					}
				}

				if (fSW && (m_nMap[pSW][iSW][jSW] >= iActiveLevel)) {
					if (((m_nMap[pW][iW][jW] < iActiveLevel) &&
							(m_nMap[pN][iN][jN] >= iActiveLevel)) ||
						((m_nMap[pS][iS][jS] < iActiveLevel) &&
							(m_nMap[pE][iE][jE] >= iActiveLevel))
					) {
						SetMinimumRefineLevel(
							p, i, j, iActiveLevel, nActiveFineRatio);
						fChanges = true;
						continue;
					}
				}

				// Check for faces bounded on opposite sides
				if ((m_nMap[pE][iE][jE] >= iActiveLevel) &&
					(m_nMap[pW][iW][jW] >= iActiveLevel)
				) {
					SetMinimumRefineLevel(
						p, i, j, iActiveLevel, nActiveFineRatio);
					fChanges = true;
					continue;
				}

				if ((m_nMap[pN][iN][jN] >= iActiveLevel) &&
					(m_nMap[pS][iS][jS] >= iActiveLevel)
				) {
					SetMinimumRefineLevel(
						p, i, j, iActiveLevel, nActiveFineRatio);
					fChanges = true;
					continue;
				}

				if (fNE && fNW && fSE && fSW) {
					if ((m_nMap[pNE][iNE][jNE] >= iActiveLevel) &&
						(m_nMap[pSW][iSW][jSW] >= iActiveLevel) &&
						(m_nMap[pSE][iSE][jSE] < iActiveLevel) &&
						(m_nMap[pNW][iNW][jNW] < iActiveLevel)
					) {
						SetMinimumRefineLevel(
							p, i, j, iActiveLevel, nActiveFineRatio);
						fChanges = true;
						continue;
					}
					if ((m_nMap[pNE][iNE][jNE] < iActiveLevel) &&
						(m_nMap[pSW][iSW][jSW] < iActiveLevel) &&
						(m_nMap[pSE][iSE][jSE] >= iActiveLevel) &&
						(m_nMap[pNW][iNW][jNW] >= iActiveLevel)
					) {
						SetMinimumRefineLevel(
							p, i, j, iActiveLevel, nActiveFineRatio);
						fChanges = true;
						continue;
					}
				}
			}
			}
			}
		}

		if (iActiveLevel == 1) { 
			break;
		}

		// Loop over all elements at this refinement level
		for (int p = 0; p < 6; p++) {
		for (int i = 0; i < nActiveResolution; i++) {
		for (int j = 0; j < nActiveResolution; j++) {
			int iref = i * nActiveFineRatio;
			int jref = j * nActiveFineRatio;

			if (m_nMap[p][iref][jref] < iActiveLevel) {
				continue;
			}

			// Loop over all neighbors
			for (int ix = -1; ix <= 1; ix++) {
			for (int jx = -1; jx <= 1; jx++) {
				if ((ix == 0) && (jx == 0)) {
					continue;
				}

				int i_src = (i+ix+2)/2-1;
				int j_src = (j+jx+2)/2-1;

				int p_dest;
				int i_dest;
				int j_dest;

				bool fNeighborExists =
					CubedSphereTrans::RelativeCoord(
						nActiveResolution/2, p, i_src, j_src,
						p_dest, i_dest, j_dest);

				if (!fNeighborExists) {
					continue;
				}

				SetMinimumRefineLevel(
					p_dest, i_dest, j_dest,
					iActiveLevel-1, nActiveFineRatio*2);
			}
			}
		}
		}
		}

		iActiveLevel--;
		nActiveResolution /= 2;
	}
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::ToFile(const char * szFile) const {

	FILE * fp = fopen(szFile, "w");

	for (int p = 0; p < 6; p++) {
	for (int j = m_nMap.GetColumns()-1; j >= 0; j--) {
		for (int i = 0; i < m_nMap.GetColumns(); i++) {
			fprintf(fp, "%i ", m_nMap[p][i][j]);
		}
		fprintf(fp, "\n");
	}
	}

	fclose(fp);
}

///////////////////////////////////////////////////////////////////////////////

void CSRefinementMap::FromFile(const char * szFile) {

	FILE * fp = fopen(szFile, "r");
	if (fp == NULL) {
		_EXCEPTION1("CSRefinement file \"%s\" not found", szFile);
	}

	for (int p = 0; p < 6; p++) {
	for (int j = m_nMap.GetColumns()-1; j >= 0; j--) {
		for (int i = 0; i < m_nMap.GetColumns(); i++) {
			fscanf(fp, "%i ", &(m_nMap[p][i][j]));
		}
	}
	}

	fclose(fp);

	// Normalize the input refinement map
	Normalize();
}

///////////////////////////////////////////////////////////////////////////////

