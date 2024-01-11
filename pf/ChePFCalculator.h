//
// Copyright (C) 2022, UChicago Argonne, LLC. All rights reserved.
//
// Software Name: Generic Semi-Analytical Simulation Tool (GenSAS)
// By: Argonne National Laboratory
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
#ifndef _Che_ChePFCalculator_H_
#define _Che_ChePFCalculator_H_

#include "util/AbstractCheCalculator.h"

using namespace che::util;

namespace che {
	namespace core {
		class ChePfEmbedSystem :public CheSingleEmbedSystem {
		public:
			ChePfEmbedSystem(const chedata::PsatDataSet &sys);

			ChePfEmbedSystem(const chedata::PsatDataSet &sys, const CheState &st, double alpha=0);

			virtual vec calcEqBalance(CheSolution* sol, double alpha);

			virtual CheSingleEmbedSystem* getNewEmbeddedSystem(const CheState &st, double alpha);
		};

		class ChePfCalculator :public AbstractCheCalculator {
		public:
			vec paraEf;
			vec paraPm;
			vec pShare;
			uvec islands;
			int* perm_c;
			int* perm_r;
			int* etree;
			// int* perm_ci;
			// int* perm_ri;
			// int* etreei;

			ChePfCalculator(const chedata::PsatDataSet &sys,
				const CheCompOptions& compOpt,
				const vec &ef = vec(1).fill(1.2), const vec &pm = vec(1).fill(0.0));

			ChePfCalculator(const chedata::PsatDataSet &sys,
				const CheCompOptions& compOpt, const uvec& islands,
				const vec &ef = vec(1).fill(1.2), const vec &pm = vec(1).fill(0.0));

			virtual CheSingleEmbedSystem* getInitSystem(const chedata::PsatDataSet &sys);

			virtual chedata::PsatDataSet regulateIsland(const chedata::PsatDataSet &sys);

			virtual int calc();

			virtual CheState exportResult();

			virtual ~ChePfCalculator();

			virtual void writeMatFile(const char *);

			virtual void writeMatFile(const char *,double);

			virtual CheSingleEmbedSystem* getNewStage();

			virtual CheSolution* getCheSolution();
		};
	}
}

#if defined(ARMA_USE_SUPERLU)

extern "C"
{
	extern void arma_wrapper(dgstrs)(superlu::trans_t, superlu::SuperMatrix *, superlu::SuperMatrix *,int *, int *, superlu::SuperMatrix *,	superlu::SuperLUStat_t *, int *);
	extern void arma_wrapper(dgsisx)(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, int*, char*, double*, double*, superlu::SuperMatrix*, superlu::SuperMatrix*, void*, int, superlu::SuperMatrix*, superlu::SuperMatrix*, double*, double*, superlu::GlobalLU_t*, superlu::mem_usage_t*, superlu::SuperLUStat_t*, int*);
	extern void arma_wrapper(dgsitrf)(superlu::superlu_options_t*, superlu::SuperMatrix*,int,int, int*, void* , int, int*, int*, superlu::SuperMatrix*, superlu::SuperMatrix*,superlu::GlobalLU_t*,superlu::SuperLUStat_t*, int*);
	extern void arma_wrapper(ilu_set_default_options)(superlu::superlu_options_t*);
}

#endif

#endif