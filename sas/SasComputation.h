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
#ifndef SAS_COMPUTATION_H
#define SAS_COMPUTATION_H

#include "SasConfig.h"
#include "SasInput.h"
#include "SasLexico.h"
#include <list>
#include <vector>
#include "AbstractCheCalculator.h"

namespace che {
	namespace core {

		class IdReplaceTable {
		public:
			AstNode* ori = NULL;
			AstNode* rep = NULL;
			IdReplaceTable(AstNode* ori, AstNode* rep);
		};

		class SasComputationModel;

		class SasComputationOptions {
		public:
			double endTime = 10.0;
			double segment = 1.0;
			double alphaTol = 1e-4;
			double errorTol = 1e-6;
			int nLvl = 15;
		};

		class SasSolution {
		public:
			double absStart;
			double absEnd;
			CheSolution* solution = NULL;
			SasComputationModel* linkCompModel=NULL;

			SasSolution(SasComputationModel* linkCompModel, int type, int nState, int nLvl);

			virtual ~SasSolution();
		};

		class SasSolutionSet {
		public:
			SasSolutionSet();

			std::list<std::shared_ptr<SasSolution>> solutionLink;

			std::list<std::shared_ptr<SasSolution>> findSasSolution(double time);

			std::shared_ptr<SasSolution> findSasSolution(std::shared_ptr<SasComputationModel>&, double time);

			void writeMatFile(const char* fileName,double interval);

			void writeMatFile(const char* fileName);

			void writeJSONFile(const char* fileName,double interval);

			virtual ~SasSolutionSet();
		};

		class SasComputationModel {
		public:
			//Model storage
			std::list<std::shared_ptr<StrPearl>> strTable;
			std::list<std::shared_ptr<SasModel>> topModels;
			//DAE
			std::vector<shared_ptr<IdPearl>> idTable;
			std::vector<shared_ptr<AstTree>> eqnTable;
			int nDE = 0;
			int nAE = 0;
			int nX = 0;
			int nY = 0;

			SasComputationModel();

			void generateCompModel(GlobalPool&);

			void preprocessModels();

			void generateDAEs();

			SasSolutionSet* solve(const SasComputationOptions& options);

			virtual ~SasComputationModel();
		private:
			int tmpVarCnt = 1;

			AstNode* copyAstNodeItself(std::shared_ptr<SasModel>& pModel, AstNode* ori);

			AstNode* copyAstTree(std::shared_ptr<SasModel>& pModel, AstNode* ori);

			std::shared_ptr<StrPearl> generateNewTempVar(std::shared_ptr<SasModel>& pModel);

			bool compareTwoSubTrees(AstNode* aTree, AstNode* bTree);

			void mergeConst(AstNode* aTree);
			
			void modelFlatten();

			void doReplacement(AstNode* p, bool waiveNext, int depth, std::shared_ptr<SasModel>& currentModel,
				std::vector<std::list<std::shared_ptr<IdReplaceTable>>>& repTable);

			void processInitialValues();

			vec getStartStateIVP();

			vec calcDiff(CheSolution* solution, double alpha);

			shared_ptr<SasSolution> solveSegment(const vec& init, double curAlpha, double seg, const SasComputationOptions& options);
		};

		enum IdLocatorType{LOCATOR_INIT,LOCATOR_STATE};

		class IdLocator {
		public:
			int type= LOCATOR_INIT;
			std::list<std::shared_ptr<SasModel>> searchDomain;
			arma::vec state;
			arma::vec der;

			IdLocator(const std::list<std::shared_ptr<SasModel>>&);

			void setStates(arma::vec&, arma::vec&);

			void unsetStates();

			static std::shared_ptr<IdPearl> findId(std::list<std::shared_ptr<SasModel>>& searchDomain, const char* iden, int len);

			double findVarValue(const char* iden, int len);

			double findDerValue(const char* iden, int len);
		};
	}
}

#endif