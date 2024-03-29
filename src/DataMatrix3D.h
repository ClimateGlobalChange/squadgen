///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataMatrix3D.h
///	\author  Paul Ullrich
///	\version July 26, 2010
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

#ifndef _DATAMATRIX3D_H_
#define _DATAMATRIX3D_H_

///////////////////////////////////////////////////////////////////////////////

#include "DataVector.h"
#include "Exception.h"

#include <iostream>
#include <cstdlib>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A 3D data matrix is a datatype that stores data in a 3D structure.
///		Arithmatic operations are not supported for this datatype.
///	</summary>

template <typename DataType>
class DataMatrix3D {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		DataMatrix3D() :
			m_data(NULL)
		{
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
		}

		///	<summary>
		///		Constructor.
		///	</summary>
		DataMatrix3D(
			size_t sRows,
			size_t sColumns,
			size_t sSubColumns
		) :
			m_data(NULL)
		{
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;

			Initialize(sRows, sColumns, sSubColumns);
		}

		///	<summary>
		///		Copy constructor.
		///	</summary>
		DataMatrix3D(const DataMatrix3D<DataType> & dm)
			: m_data(NULL)
		{
			Assign(dm);
		}

		///	<summary>
		///		Destructor.
		///	</summary>
		virtual ~DataMatrix3D() {
			if (m_data != NULL) {
				free(reinterpret_cast<void*>(m_data));
			}
		}

	public:
		///	<summary>
		///		Determine if this DataMatrix3D is initialized.
		///	</summary>
		bool IsInitialized() const {
			if (m_data == NULL) {
				return false;
			} else {
				return true;
			}
		}

	public:
		///	<summary>
		///		Deallocate data for this object.
		///	</summary>
		void Deinitialize() {
			if (m_data != NULL) {
				free(reinterpret_cast<void*>(m_data));
			}

			m_data = NULL;
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
		}

		///	<summary>
		///		Allocate data for this object.
		///	</summary>
		void Initialize(
			size_t sRows,
			size_t sColumns,
			size_t sSubColumns,
			bool fAutoZero = true
		) {
			size_t sI;
			size_t sJ;

			// Check for zero size
			if ((sRows == 0) || (sColumns == 0) || (sSubColumns == 0)) {
				Deinitialize();
				return;
			}

			// No need to reallocate memory if this matrix already has
			// the correct dimensions.
			if ((m_sSize[0] == sRows) &&
				(m_sSize[1] == sColumns) &&
				(m_sSize[2] == sSubColumns)
			) {
				// Auto zero
				if (fAutoZero) {
					Zero();
				}

				return;
			}

			// Deinitialize existing content
			Deinitialize();

			// Calculate the per-row footprint
			size_t sRowPtrFootprint =
				sRows * sizeof(DataType **);
			size_t sColumnPtrFootprint =
				sColumns * sizeof(DataType *);
			size_t sColumnFootprint =
				sSubColumns * sizeof(DataType);

			// Allocate memory
			char *rawdata = reinterpret_cast<char*>(
				malloc(
					sRowPtrFootprint +
					sRows * sColumnPtrFootprint +
					sRows * sColumns * sColumnFootprint
				));

			if (rawdata == NULL) {
				_EXCEPTIONT("Out of memory.");
			}

			// Assign memory pointers
			char *pColumnPtrs = rawdata + sRowPtrFootprint;
			char *pDataStart = pColumnPtrs + sRows * sColumnPtrFootprint;

			m_data = reinterpret_cast<DataType***>(rawdata);

			for (sI = 0; sI < sRows; sI++) {
				m_data[sI] = reinterpret_cast<DataType**>(
					pColumnPtrs + sI * sColumnPtrFootprint
				);

				for (sJ = 0; sJ < sColumns; sJ++) {
					m_data[sI][sJ] = reinterpret_cast<DataType*>(
						pDataStart + (sI * sColumns + sJ) * sColumnFootprint
					);
				}
			}

			// Assign dimensions
			m_sSize[0] = sRows;
			m_sSize[1] = sColumns;
			m_sSize[2] = sSubColumns;

			// Auto zero
			if (fAutoZero) {
				Zero();
			}
		}

	public:
		///	<summary>
		///		Assignment operator.
		///	</summary>
		void Assign(const DataMatrix3D<DataType> & dm) {

			// Check initialization status
			if (!dm.IsInitialized()) {
				Deinitialize();
				return;
			}

			// Allocate memory
			Initialize(dm.m_sSize[0], dm.m_sSize[1], dm.m_sSize[2]);

			// Copy data
			size_t sRowPtrFootprint = m_sSize[0] * sizeof(DataType **);
			size_t sColPtrFootprint = m_sSize[1] * sizeof(DataType *);

			memcpy(
				reinterpret_cast<char*>(m_data) +
					sRowPtrFootprint + m_sSize[0] * sColPtrFootprint,
				reinterpret_cast<char*>(dm.m_data) +
					sRowPtrFootprint + m_sSize[0] * sColPtrFootprint,
				m_sSize[0] * m_sSize[1] * m_sSize[2] * sizeof(DataType)
			);
		}

		///	<summary>
		///		Assignment operator.
		///	</summary>
		DataMatrix3D & operator= (const DataMatrix3D<DataType> & dm) {
			Assign(dm);
			return (*this);
		}

		///	<summary>
		///		Zero the data content of this object.
		///	</summary>
		void Zero() {

			// Check initialization status
			if (!IsInitialized()) {
				_EXCEPTIONT(
					"Attempted operation on uninitialized DataMatrix3D.");
			}

			// Set content to zero
			memset(
				&(m_data[0][0][0]),
				0,
				m_sSize[0] * m_sSize[1] * m_sSize[2] * sizeof(DataType)
			);
		}

	public:
		///	<summary>
		///		Copy the data of this object into a DataVector.
		///	</summary>
		void Vectorize(DataVector<DataType> & vec) const {

			// Check initialization status
			if (!IsInitialized()) {
				_EXCEPTIONT(
					"Attempted operation on uninitialized DataMatrix3D.");
			}

			// Allocate memory
			vec.Initialize(m_sSize[0] * m_sSize[1] * m_sSize[2]);

			// Copy data
			memcpy(
				reinterpret_cast<char*>(&(vec[0])),
				reinterpret_cast<char*>(&(m_data[0][0][0])),
				m_sSize[0] * m_sSize[1] * m_sSize[2] * sizeof(DataType));
		}

	public:
		///	<summary>
		///		Get the number of rows in this matrix.
		///	</summary>
		inline size_t GetRows() const {
			return m_sSize[0];
		}

		///	<summary>
		///		Get the number of columns in this matrix.
		///	</summary>
		inline size_t GetColumns() const {
			return m_sSize[1];
		}

		///	<summary>
		///		Get the number of columns in this matrix.
		///	</summary>
		inline size_t GetSubColumns() const {
			return m_sSize[2];
		}

		///	<summary>
		///		Get the number of rows in this matrix.
		///	</summary>
		inline size_t GetSize(int dim) const {
			return m_sSize[dim];
		}

		///	<summary>
		///		Get the total number of elements in this matrix.
		///	</summary>
		inline size_t GetTotalElements() const {
			return m_sSize[0] * m_sSize[1] * m_sSize[2];
		}

	public:
		///	<summary>
		///		Cast to an array.
		///	</summary>
		inline operator DataType***() const {
			return m_data;
		}

	private:
		///	<summary>
		///		The number of elements in each dimension of this matrix.
		///	</summary>
		size_t m_sSize[3];

		///	<summary>
		///		A pointer to the data associated with this matrix.
		///	</summary>
		DataType*** m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

