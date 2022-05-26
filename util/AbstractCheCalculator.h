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
#ifndef _Che_AbstractCheCalculator_H_
#define _Che_AbstractCheCalculator_H_

#include "CheState.h"
#include "CheDataFormat.h"
#include <list>
#include "CheYMatrix.h"

using namespace che::util;
using namespace che::io;
using namespace arma;
using namespace std;

namespace che {
	namespace core {		
		enum CheSolutionType{ CHESOL_NONE, CHESOL_PS, CHESOL_PADE};
		enum PadeSolverType { PADE_LU, PADE_LEVINSON };

		class CheCompOptions {
		public:
			int nLvl;
			double maxAlpha;
			double alphaTol;
			double segLen;
			double diffTol;
			double diffTolMax;

			CheCompOptions(int nLvl, double maxAlpha, double alphaTol, double segLen, double diffTol, double diffTolMax);
		
			CheCompOptions(const CheCompOptions&) = default;

			CheCompOptions& operator=(const CheCompOptions&) = default;
		};

		class CheSolution {
		public:
			int nState;
			int nLvl;
			mat solution;
			int type;

			CheSolution(int nState, int nLvl);

			virtual vec getSolValue(double alpha)=0;

			virtual ~CheSolution();

			CheSolution(const CheSolution&) = delete;

			CheSolution& operator=(const CheSolution&) = delete;
		};

		class CheSolutionPowerSeries: public CheSolution {
		public:
			CheSolutionPowerSeries(int nState, int nLvl);

			virtual vec getSolValue(double alpha);

			virtual ~CheSolutionPowerSeries();
		};

		class CheSolutionPade : public CheSolution{
		public:
			int num;
			int den;
			mat numerator;
			mat denomenator;
			bool ready;
			PadeSolverType solver;

			CheSolutionPade(int nState, int num, int den, PadeSolverType sol = PADE_LEVINSON);

			CheSolutionPade(int nState, int nLvl, PadeSolverType sol = PADE_LEVINSON);

			virtual vec getSolValue(double alpha);

			virtual ~CheSolutionPade();

		private:
			bool genPadeCoeffLU();

			bool genPadeCoeffLevinson();
		};

		class CheSolutionFactory {
		public:
			static CheSolution* makeInitCheSol(int type, int nState, int nLvl);
			static CheSolution* makeCopyCheSol(CheSolution* sol);
			static CheSolution* makeDerivativeSol(CheSolution* sol);
		};

		class CheSingleEmbedSystem {
		public:
			CheState initState;
			chedata::PsatDataSet baseSys;
			double startAlpha;
			CheYMatrix yMatrix;

			CheSingleEmbedSystem(const CheState&, const chedata::PsatDataSet&, double);

			virtual vec calcEqBalance(CheSolution* sol,double alpha) = 0;

			virtual CheSingleEmbedSystem* getNewEmbeddedSystem(const CheState &st, double alpha) = 0;
			
			CheSingleEmbedSystem(const CheSingleEmbedSystem&)= delete;

			CheSingleEmbedSystem& operator=(const CheSingleEmbedSystem&) = delete;
		};

		class AbstractCheCalculator {
		public:
			CheCompOptions compOpt;
			bool reachesMaxAlpha;
			chedata::PsatDataSet baseSys;
			list<CheSingleEmbedSystem*> cheList;
			list<CheSolution*> solList;

			AbstractCheCalculator(const chedata::PsatDataSet &psat, const CheCompOptions& compOpt);

			virtual CheSingleEmbedSystem* getInitSystem(const chedata::PsatDataSet &sys) = 0;

			virtual chedata::PsatDataSet regulateIsland(const chedata::PsatDataSet &sys);

			virtual int calc();

			virtual ~AbstractCheCalculator();

			virtual CheState exportResult() = 0;

			virtual void writeMatFile(const char *) =0;
			
			virtual void writeMatFile(const char *, double) =0;

		//protected:
			virtual CheSingleEmbedSystem* getNewStage() = 0;

			virtual CheSolution* getCheSolution() = 0;
		};
	}
}

#endif