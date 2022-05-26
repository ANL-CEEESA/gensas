//
// Copyright (C) 2022, UChicago Argonne, LLC. All rights reserved.
//
// OPEN SOURCE LICENSE
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 
// 
// ******************************************************************************************************
// DISCLAIMER
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ***************************************************************************************************
// 
#ifndef _Che_Comp_Util_H_
#define _Che_Comp_Util_H_

#include "CheDataFormat.h"
#include "CheState.h"
#include "CheYMatrix.h"
#include "CheEvents.h"
#include <list>
#include <vector>
#include <assert.h>
#include <sstream>

using namespace che::io;
using namespace std;

namespace che {
	namespace util {
		class CheCompUtil {
		public:
			static CheStateIdx getCheStateIdx(const chedata::PsatDataSet &cheData);
			static CheYMatrix getCheYMatrix(const chedata::PsatDataSet &cheData, const list<Fault>& faultList = list<Fault>());
			static list<chedata::PsatDataSet> splitIslands(const chedata::PsatDataSet &cheData, const uvec& islands);
			static uvec searchIslands(const chedata::PsatDataSet &cheData);
			static umat spgetseq(int n, int d);
			
			/**
		* @brief sp_submatrix Function to extract columns and rows from a sparse matrix
		* @param A Sparse matrix pointer
		* @param rows vector of the rown indices to keep (must be sorted)
		* @param cols vector of the clumn indices to keep (must be sorted)
		* @return Sparse matrix of the indicated indices
		*/
			//template<typename T> static arma::SpMat<T> sp_submatrix(const arma::SpMat<T> &A, const arma::uvec &rows, const arma::uvec &cols) {

			//	// declare reduced matrix
			//	std::size_t n_rows = rows.size();
			//	std::size_t n_cols = cols.size();

			//	bool found = false;
			//	std::size_t n = 0;
			//	std::size_t p = 0;
			//	std::size_t found_idx = 0;

			//	arma::Col<T> new_val(A.n_nonzero);
			//	arma::uvec new_row_ind(A.n_nonzero);
			//	arma::uvec new_col_ptr(n_cols + 1);

			//	new_col_ptr(p) = 0;

			//	for (auto const& j : cols) { // for every column in the cols vector

			//		for (std::size_t k = A.col_ptrs[j]; k < A.col_ptrs[j + 1]; k++) {  // k is the index of the "values" and "row_indices" that corresponds to the column j

			//			// search row_ind[k] in rows
			//			found = false;
			//			found_idx = 0;
			//			while (!found && found_idx < n_rows) {
			//				if (A.row_indices[k] == rows.at(found_idx))
			//					found = true;
			//				found_idx++;
			//			}

			//			// store the values if the row was found in rows
			//			if (found) { // if the row index is in the designated rows...
			//				new_val(n) = A.values[k]; // store the value
			//				new_row_ind(n) = found_idx - 1;  // store the index where the original index was found inside "rows"
			//				n++;
			//			}
			//		}

			//		p++;
			//		new_col_ptr(p) = n;
			//	}
			//	new_col_ptr(p) = n;

			//	// reshape the vectors to the actual number of elements
			//	new_val.reshape(n, 1);
			//	new_row_ind.reshape(n, 1);

			//	// return new object
			//	return arma::SpMat<T>(new_row_ind, new_col_ptr, new_val, n_rows, n_cols);
			//}

			template<typename T> static arma::SpMat<T> sp_submatrix(const arma::SpMat<T> &A, const arma::uvec &rows, const arma::uvec &cols) {

				// declare reduced matrix
				std::size_t n_rows = rows.size();
				std::size_t n_cols = cols.size();

				bool found = false;
				std::size_t n = 0;
				std::size_t p = 0;
				std::size_t found_idx = 0;

				arma::Col<T> new_val(A.n_nonzero);
				arma::uvec new_row_ind(A.n_nonzero);
				arma::uvec new_col_ptr(n_cols + 1);

				arma::uvec rowOldToNew(A.n_rows);
				arma::uvec colOldToNew(A.n_cols);
				rowOldToNew.fill(-1);
				colOldToNew.fill(-1);
				rowOldToNew(rows)= regspace<uvec>(0, n_rows - 1);
				colOldToNew(cols) = regspace<uvec>(0, n_cols - 1);

				new_col_ptr(p) = 0;

				for (auto const& j : cols) { // for every column in the cols vector

					for (std::size_t k = A.col_ptrs[j]; k < A.col_ptrs[j + 1]; k++) {  // k is the index of the "values" and "row_indices" that corresponds to the column j

						// search row_ind[k] in rows
						//found = false;
						//found_idx = 0;
						//while (!found && found_idx < n_rows) {
						//	if (A.row_indices[k] == rows.at(found_idx))
						//		found = true;
						//	found_idx++;
						//}

						//// store the values if the row was found in rows
						//if (found) { // if the row index is in the designated rows...
						//	new_val(n) = A.values[k]; // store the value
						//	new_row_ind(n) = found_idx - 1;  // store the index where the original index was found inside "rows"
						//	n++;
						//}

						int newRowIdx = rowOldToNew(A.row_indices[k]);
						if (newRowIdx != -1) {
							new_val(n) = A.values[k]; // store the value
							new_row_ind(n) = newRowIdx;  // store the index where the original index was found inside "rows"
							n++;
						}
					}

					p++;
					new_col_ptr(p) = n;
				}
				new_col_ptr(p) = n;

				// reshape the vectors to the actual number of elements
				new_val.reshape(n, 1);
				new_row_ind.reshape(n, 1);

				// return new object
				return arma::SpMat<T>(new_row_ind, new_col_ptr, new_val, n_rows, n_cols);
			}
		};
			   
		template<typename T> string armaToString(const Mat<T>& mat) {
			stringstream ss;
			ss << mat;
			return ss.str();
		}

		template<typename T> string armaToString(const Col<T>& mat) {
			stringstream ss;
			ss << mat;
			return ss.str();
		}

		template<typename T> string armaToString(const Row<T>& mat) {
			stringstream ss;
			ss << mat;
			return ss.str();
		}

		template<typename T> string armaToString(const SpMat<T>& mat) {
			stringstream ss;
			ss << mat;
			return ss.str();
		}
	}
}

#endif