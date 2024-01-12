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
#include "util/AbstractCheCalculator.h"
#include "util/CheCompUtil.h"

using namespace che::util;
using namespace che::io;

namespace che {
	namespace core {
		static vec getPowerSeriesValue(const mat &coeff, double alpha) {
			vec sol = coeff.col(coeff.n_cols - 1);
			for (int i = coeff.n_cols - 2; i >= 0; i--) {
				sol = alpha * sol + coeff.col(i);
			}
			return sol;
		}

		static mat solveToepLU(const mat &ct, const mat &y) {
			int d = ct.n_rows;
			int n = y.n_cols;

			mat x(d, n, fill::zeros);

			for (int i = 0; i < d; i++) {
				mat toep(n, n);
				for (int j = -n + 1; j < n; j++) {
					toep.diag(j).fill(ct(i, n - 1 - j));
				}
				//toep.print("toep:");
				try {
					x.row(i) = solve(toep, y.row(i).t(), solve_opts::no_approx).t();
				}
				catch (const std::runtime_error& e) {
					x.row(i).fill(datum::nan);
				}
			}
			return x;
		}

		static mat solveToepLevinson(const mat &ct, const mat &y) {
			int d = ct.n_rows;
			int n = y.n_cols;

			mat f(d, n, fill::zeros);
			mat b(d, n, fill::zeros);
			mat temp(d, n, fill::zeros);
			mat x(d, n, fill::zeros);
			vec epsf(d);
			vec epsb(d);
			vec epsx(d);
			vec alphaf(d);
			vec betaf(d);
			vec alphab(d);
			vec betab(d);

			f.col(0) = ones<vec>(d) / ct.col(n - 1);
			b.col(n - 1) = f.col(0);
			x.col(0) = y.col(0) % f.col(0);

			for (int i = 1; i < n; i++) {
				epsf = sum(f.cols(0, i - 1) % fliplr(ct.cols(n, n + i - 1)), 1);
				epsb = sum(b.cols(n - i, n - 1) % fliplr(ct.cols(n - i - 1, n - 2)), 1);
				epsx = sum(x.cols(0, i - 1) % fliplr(ct.cols(n, n + i - 1)), 1);

				alphaf = ones<vec>(d) / (ones<vec>(d) - epsf % epsb);
				betaf = -epsf % alphaf;
				alphab = -epsb % alphaf;
				betab = alphaf;

				temp.cols(0, i - 1) = f.cols(0, i - 1);
				f.cols(0, i) = temp.cols(0, i).each_col() % alphaf + b.cols(n - i - 1, n - 1).each_col() % betaf;
				b.cols(n - i - 1, n - 1) = temp.cols(0, i).each_col() % alphab + b.cols(n - i - 1, n - 1).each_col() % betab;

				x.cols(0, i) += b.cols(n - i - 1, n - 1).each_col() % (y.col(i) - epsx);
			}

			return x;
		}

		CheCompOptions::CheCompOptions(int nLvl, double maxAlpha, double alphaTol, double segLen, double diffTol, double diffTolMax) {
			this->nLvl = nLvl;
			this->maxAlpha = maxAlpha;
			this->alphaTol = alphaTol;
			this->segLen = segLen;
			this->diffTol = diffTol;
			this->diffTolMax = diffTolMax;
		}

		CheSolution::CheSolution(int nState, int nLvl) :solution(nState, nLvl, fill::zeros) {
			this->type = CHESOL_NONE;
			this->nState = nState;
			this->nLvl = nLvl;
		}

		CheSolution::~CheSolution(){}

		CheSolutionPowerSeries::CheSolutionPowerSeries(int nState, int nLvl) :CheSolution(nState, nLvl) {
			type = CHESOL_PS;
		}

		vec CheSolutionPowerSeries::getSolValue(double alpha) {
			return getPowerSeriesValue(solution, alpha);
		}

		CheSolutionPowerSeries::~CheSolutionPowerSeries() {}

		CheSolutionPade::CheSolutionPade(int nState, int num, int den, PadeSolverType sol) :CheSolution(nState, num + den) {
			this->num = num;
			this->den = den;
			ready = false;
			solver = sol;
			this->type = CHESOL_PADE;
		}

		CheSolutionPade::CheSolutionPade(int nState, int nLvl, PadeSolverType sol) 
			:CheSolutionPade(nState, nLvl - nLvl / 2, nLvl / 2,  sol) {}

		vec CheSolutionPade::getSolValue(double alpha) {
			if (!ready) {
				if (solver == PADE_LEVINSON) {
					ready = genPadeCoeffLevinson();
				}
				else {
					ready = genPadeCoeffLU();
				}
			}
			if (!ready) {
				return getPowerSeriesValue(solution, alpha);
			}
			else {
				vec sol = getPowerSeriesValue(numerator, alpha) / getPowerSeriesValue(join_rows(ones<mat>(denomenator.n_rows, 1), denomenator), alpha);
				uvec nfrows = find_nonfinite(sol);
				if (!nfrows.is_empty()) {
					sol.rows(nfrows) = getPowerSeriesValue(solution.rows(nfrows), alpha);
					numerator.rows(nfrows)=solution.submat(nfrows,regspace<uvec>(0,this->num-1));
					denomenator.rows(nfrows).fill(0.0);
				}
				return sol;
			}
		}

		bool CheSolutionPade::genPadeCoeffLU() {
			mat augTc = join_rows(mat(nState, den, fill::zeros), solution);
			int nT = 2 * den - 1;
			mat ttxc = augTc.cols(augTc.n_cols - nT - 1, augTc.n_cols - 2);
			mat ytxc = augTc.tail_cols(den);

			this->denomenator = solveToepLU(ttxc, ytxc);
			this->numerator = augTc.cols(den, num + den - 1);

			for (int i = 0; i < den; i++) {
				this->numerator += augTc.cols(den - i - 1, num + den - i - 2).each_col() % this->denomenator.col(i);
			}
			return true;
		}

		bool CheSolutionPade::genPadeCoeffLevinson() {
			mat augTc = join_rows(mat(nState, den, fill::zeros), solution);
			int nT = 2 * den - 1;
			mat ttxc = augTc.cols(augTc.n_cols - nT - 1, augTc.n_cols - 2);
			mat ytxc = augTc.tail_cols(den);

			this->denomenator = -solveToepLevinson(ttxc, ytxc);
			this->numerator = augTc.cols(den, num + den - 1);

			/*cout.precision(16);
			cout.setf(ios::scientific);

			this->denomenator.raw_print(cout, "den");
			this->numerator.raw_print(cout, "num");*/

			for (int i = 0; i < den; i++) {
				this->numerator += augTc.cols(den - i - 1, num + den - i - 2).each_col() % this->denomenator.col(i);
				//this->numerator.raw_print(cout, "num");
			}
			return true;
		}

		CheSolutionPade::~CheSolutionPade() {}

		CheSingleEmbedSystem::CheSingleEmbedSystem(const CheState& init, const chedata::PsatDataSet& baseSys, double startAlpha = 0)
			:initState(init), baseSys(baseSys) {
			this->startAlpha = startAlpha;
			yMatrix = CheCompUtil::getCheYMatrix(baseSys);
		}

		CheSolution* CheSolutionFactory::makeInitCheSol(int type, int nState, int nLvl) {
			CheSolution* pSol = NULL;
			if (type == CHESOL_PS) {
				pSol = new CheSolutionPowerSeries(nState, nLvl);
			}
			else if (type == CHESOL_PADE) {
				pSol = new CheSolutionPade(nState, nLvl, PADE_LEVINSON);
			}
			return pSol;
		}

		CheSolution* CheSolutionFactory::makeCopyCheSol(CheSolution* sol) {
			if (sol == NULL)
				return NULL;
			if (sol->type == CHESOL_PS) {
				CheSolutionPowerSeries* pSol = new CheSolutionPowerSeries(sol->nState, sol->nLvl);
				CheSolutionPowerSeries* pOri = dynamic_cast<CheSolutionPowerSeries*>(sol);
				if (pOri != nullptr) {
					pSol->solution = pOri->solution;
				}
				return pSol;
			}
			else if (sol->type == CHESOL_PADE) {
				CheSolutionPade* pSol = new CheSolutionPade(sol->nState, sol->nLvl, PADE_LEVINSON);
				CheSolutionPade* pOri = dynamic_cast<CheSolutionPade*>(sol);
				if (pOri != nullptr) {
					pSol->num = pOri->num;
					pSol->den = pOri->num;
					pSol->numerator = pOri->numerator;
					pSol->denomenator = pOri->denomenator;
					pSol->ready = pOri->ready;
				}
				return pSol;
			}
			return NULL;
		}

		mat getDer(const mat& c) {
			mat d(c.n_rows, c.n_cols, fill::zeros);
			for (int i = 1; i < c.n_cols; i++) {
				d.col(i - 1) = i * c.col(i);
			}
			return d;
		}

		CheSolution* CheSolutionFactory::makeDerivativeSol(CheSolution* sol) {
			if (sol == NULL)
				return NULL;
			if (sol->type == CHESOL_PS) {
				CheSolutionPowerSeries* pSol = new CheSolutionPowerSeries(sol->nState, sol->nLvl);
				CheSolutionPowerSeries* pOri = dynamic_cast<CheSolutionPowerSeries*>(sol);
				if (pOri != nullptr) {					
					pSol->solution = getDer(pOri->solution);
				}
				return pSol;
			}
			else if (sol->type == CHESOL_PADE) {
				CheSolutionPade* pSol = new CheSolutionPade(sol->nState, sol->nLvl, PADE_LEVINSON);
				CheSolutionPade* pOri = dynamic_cast<CheSolutionPade*>(sol);
				if (pOri != nullptr) {
					pSol->solution = getDer(pOri->solution);
					pSol->num = pOri->num;
					pSol->den = pOri->num;
					pSol->ready = false;
				}
				return pSol;
			}
			return NULL;
		}

		AbstractCheCalculator::AbstractCheCalculator(
			const chedata::PsatDataSet &sys, const CheCompOptions& compOpt)
			: baseSys(regulateIsland(sys)),compOpt(compOpt) {
			this->cheList = list<CheSingleEmbedSystem*>();
			this->solList = list<CheSolution*>();
			this->reachesMaxAlpha = false;
		}

		chedata::PsatDataSet AbstractCheCalculator::regulateIsland(const chedata::PsatDataSet &sys) {
			if (!sys.isFormatted) {
				chedata::PsatDataSet newSys(sys);
				newSys.renumberBuses();
				return newSys;
			}
			else {
				return sys;
			}
		}

		int AbstractCheCalculator::calc() {
			double alpha = 0;
			CheSingleEmbedSystem *pinit = getInitSystem(baseSys);
			cheList.push_back(pinit);
			CheSingleEmbedSystem* currSys;
			while ((currSys = cheList.back()) != NULL) {
				CheSolution* cheSol = getCheSolution();
				solList.push_back(cheSol);
				if (cheSol == NULL) {
					break;
				}
				cheList.push_back(getNewStage());
				if (reachesMaxAlpha) {
					break;
				}
			}
			return 0;
		}

		CheState AbstractCheCalculator::exportResult(){
			return cheList.back()->initState;
		}

		AbstractCheCalculator::~AbstractCheCalculator() {
			for (auto&&che : cheList) {
				if (che != NULL) delete che;
			}
			cheList.clear();
			for (auto&&sol : solList) {
				if (sol != NULL) delete sol;
			}
			solList.clear();
		}
	}
}