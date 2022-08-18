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
#include "CheCompUtil.h"

using namespace arma;
using namespace std;

namespace che {
	namespace util {
		static int nchoosek(int n, int k) {
			if (n >= 0 && k >= 0 && k <= n) {
				if (n == 0) {
					return 1;
				}
				if (k > n / 2) {
					return nchoosek(n, n - k);
				}
				else {
					long numer = n;
					long denom = 1;
					for (int i = 1; i < k; i++) {
						numer *= n - i;
						denom *= i + 1;
					}
					return numer / denom;
				}
			}
			else {
				return 0;
			}
		}

		static double nchoosek(double n, double k) {
			if (n >= 0 && k >= 0 && k <= n) {
				if (n == 0) {
					return 1;
				}
				if (k > n / 2) {
					return nchoosek(n, n - k);
				}
				else {
					double res = n;
					for (int i = 1; i < k; i++) {
						res *= n - i;
						res /= i + 1;
					}
					return res;
				}
			}
			else {
				return 0;
			}
		}

		int round(double x) {
			return (int)(x + 0.5 - (x < 0));
		}

		umat CheCompUtil::spgetseq(int n, int d) {
			unsigned int nrow = round(nchoosek((double)(n + d - 1), (double)(d - 1)));
			umat seq(nrow, d, fill::zeros);
			seq(0, 0) = n;
			unsigned int max = n;
			unsigned int sum = 0;
			unsigned int temp = 0;
			for (int k = 1; k < nrow; k++) {
				if (seq(k - 1, 0) > 0) {
					seq(k, 0) = seq(k - 1, 0) - 1;
					for (int l = 1; l < d; l++) {
						if (seq(k - 1, l) < max) {
							seq(k, l) = seq(k - 1, l) + 1;
							for (int m = l + 1; m < d; m++) {
								seq(k, m) = seq(k - 1, m);
							}
							break;
						}
					}
				}
				else {
					sum = 0;
					for (int l = 1; l < d; l++) {
						if (seq(k - 1, l) < max) {
							seq(k, l) = seq(k - 1, l) + 1;
							sum += seq(k, l);
							for (int m = l + 1; m < d; m++) {
								seq(k, m) = seq(k - 1, m);
								sum += seq(k, m);
							}
							break;
						}
						else {
							temp = 0;
							for (int m = l + 2; m < d; m++) {
								temp += seq(k - 1, m);
							}
							max = n - temp;
							seq(k, l) = 0;
						}
					}
					seq(k, 0) = n - sum;
					max = n - sum;
				}
			}
			return seq;
		}

		static uvec generateIdx(int& sCount, int num) {
			if (num > 0) {
				int start = sCount;
				int end = sCount + num - 1;
				sCount += num;
				return regspace<uvec>(start, end);
			}
			else {
				return uvec();
			}
		}

		CheStateIdx CheCompUtil::getCheStateIdx(const chedata::PsatDataSet &cheData) {
			CheStateIdx stateIdx = CheStateIdx();
			int sCount = 0;
			stateIdx.vrIdx = generateIdx(sCount, cheData.nBus);
			stateIdx.viIdx = generateIdx(sCount, cheData.nBus);
			stateIdx.qIdx = generateIdx(sCount, cheData.nBus);
			stateIdx.pIdx = generateIdx(sCount, cheData.nBus);
			
			stateIdx.sIdx = generateIdx(sCount, cheData.nInd);
			stateIdx.indEr1Idx = generateIdx(sCount, cheData.nInd);
			stateIdx.indEm1Idx = generateIdx(sCount, cheData.nInd);
			stateIdx.indEr2Idx = generateIdx(sCount, cheData.nInd);
			stateIdx.indEm2Idx = generateIdx(sCount, cheData.nInd);
			
			stateIdx.mDeltaIdx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mOmegaIdx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mEq1Idx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mEq2Idx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mEd1Idx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mEd2Idx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mPsidIdx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mPsiqIdx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mPgIdx = generateIdx(sCount, cheData.nSyn);
			stateIdx.mEfIdx = generateIdx(sCount, cheData.nSyn);

			int nExc1 = 0;
			int nExc2 = 0;
			int nExc3 = 0;
			if (cheData.nExc > 0) {
				uvec excType = cheData.get_excs_excType_vec();
				uvec exc1Idx = find(excType == 1);
				uvec exc2Idx = find(excType == 2);
				uvec exc3Idx = find(excType == 3);
				nExc1 = exc1Idx.n_rows;
				nExc2 = exc2Idx.n_rows;
				nExc3 = exc3Idx.n_rows;
			}
			stateIdx.avr1mIdx = generateIdx(sCount, nExc1);
			stateIdx.avr1r1Idx = generateIdx(sCount, nExc1);
			stateIdx.avr1r2Idx = generateIdx(sCount, nExc1);
			stateIdx.avr1rIdx = generateIdx(sCount, nExc1);
			stateIdx.avr1fIdx = generateIdx(sCount, nExc1);
			
			stateIdx.avr2mIdx = generateIdx(sCount, nExc2);
			stateIdx.avr2r1Idx = generateIdx(sCount, nExc2);
			stateIdx.avr2rIdx = generateIdx(sCount, nExc2);
			stateIdx.avr2r2Idx = generateIdx(sCount, nExc2);
			stateIdx.avr2fIdx = generateIdx(sCount, nExc2);

			stateIdx.avrmIdx = generateIdx(sCount, nExc3);
			stateIdx.avrrIdx = generateIdx(sCount, nExc3);
			stateIdx.avrfIdx = generateIdx(sCount, nExc3);
			stateIdx.avrrefIdx = generateIdx(sCount, nExc3);

			int nTg1 = 0;
			int nTg2 = 0;
			if (cheData.nTg > 0) {
				uvec tgType = cheData.get_tgs_tgType_vec();
				uvec tg1Idx = find(tgType == 1);
				uvec tg2Idx = find(tgType == 2);
				nTg1 = tg1Idx.n_rows;
				nTg2 = tg2Idx.n_rows;
			}
			stateIdx.tg1inIndx = generateIdx(sCount, nTg1);
			stateIdx.tg11Idx = generateIdx(sCount, nTg1);
			stateIdx.tg12Idx = generateIdx(sCount, nTg1);
			stateIdx.tg13Idx = generateIdx(sCount, nTg1);

			stateIdx.tg2gIdx = generateIdx(sCount, nTg2);
			stateIdx.tg2mIdx = generateIdx(sCount, nTg2);

			stateIdx.tmechIdx = generateIdx(sCount, cheData.nTg);

			stateIdx.fIdx= generateIdx(sCount, cheData.nBus);
			stateIdx.qpltIdx = generateIdx(sCount, cheData.nBus);
			stateIdx.vgIdx = generateIdx(sCount, cheData.nBus);

			stateIdx.nState = sCount;

			return stateIdx;
		}

		CheYMatrix CheCompUtil::getCheYMatrix(const chedata::PsatDataSet &cheData, const list<Fault>& faultList) {
			CheYMatrix yMatrix = CheYMatrix();

			int nBus = cheData.nBus;

			vec r = cheData.get_lines_r_vec();
			vec x = cheData.get_lines_x_vec();
			cx_vec z(r, x);
			vec b = cheData.get_lines_b_vec();
			vec status = conv_to<vec>::from(cheData.get_lines_status_vec());
			cx_vec chrg1(0.0*b, 0.5*status%b);
			cx_vec chrg2 = chrg1;
			cx_vec y = status / z;

			if (!faultList.empty()) {
				vector<int> idxVec = vector<int>();
				vector<double> posVec = vector<double>();
				vector<complex<double>> zVec = vector<complex<double>>();
				for (auto&&fault : faultList) {
					if (fault.fType == Fault::FAULT_3P || fault.fType == Fault::FAULT_3PG) {
						idxVec.push_back(fault.lineIdx-1);
						posVec.push_back(fault.pos);
						cx_double zz = 1e-6;
						if (!fault.faultY.empty()) zz = ((double)(fault.faultY.n_rows)) / trace(fault.faultY);
						zVec.push_back(zz);
					}
				}
				uvec lineIdx = conv_to<uvec>::from(idxVec);
				vec pos = conv_to<vec>::from(posVec);
				cx_vec zf = conv_to<cx_vec>::from(zVec);

				cx_vec zf1 = z(lineIdx) % pos;
				cx_vec zf2 = z(lineIdx) % (ones<vec>(pos.n_rows) - pos);
				cx_vec zden = zf1 % zf2 + zf2 % zf + zf % zf1;
				cx_vec yffr = zf2 / zden;
				cx_vec yfto = zf1 / zden;
				cx_vec yftr = zf / zden;

				cx_vec rzf2 = ones<vec>(zf2.n_rows) / zf2;
				cx_vec rzf1 = ones<vec>(zf1.n_rows) / zf1;
				uvec izf1g_zero = find(((zf1 == 0.) % (zf == 0.)) == 1);
				uvec izf2g_zero= find(((zf2 == 0.) % (zf == 0.)) == 1);

				yfto(izf1g_zero) = rzf2(izf1g_zero);
				yftr(izf1g_zero) = rzf2(izf1g_zero);
				yffr(izf2g_zero) = rzf1(izf2g_zero);
				yftr(izf2g_zero) = rzf1(izf2g_zero);

				y(lineIdx) = yftr;
				chrg1(lineIdx) += yffr;
				chrg2(lineIdx) += yfto;
			}

			vec k = cheData.get_lines_k_vec();
			k(find(k == 0)).fill(1.0);
			vec angInArc = datum::pi / 180.0*cheData.get_lines_ang_vec();
			cx_vec ts = k % cx_vec(cos(angInArc), sin(angInArc));
			vec ts2 = abs(ts)%abs(ts);

			yMatrix.ytrfr = y / conj(ts);
			yMatrix.ytrto = y / ts;
			yMatrix.yshfr = (y + chrg1) / ts2 - yMatrix.ytrfr;
			yMatrix.yshto = y + chrg2 - yMatrix.ytrto;

			uvec ifr = cheData.get_lines_fromBus_vec()-1;
			uvec ito = cheData.get_lines_toBus_vec()-1;

			umat loc = join_cols(join_rows(ifr, ito), join_rows(ito, ifr), join_rows(ifr, ifr), join_rows(ito, ito)).t();
			cx_vec val = join_cols(-y / conj(ts), -y / ts, (y + chrg1) / ts2, y + chrg2);

			yMatrix.Y = sp_cx_mat(true,loc, val,nBus,nBus);
			yMatrix.Ysh = sum(yMatrix.Y, 1);
			yMatrix.Ytr = yMatrix.Y;
			yMatrix.Ytr.diag() -= yMatrix.Ysh;

			return yMatrix;
		}

		static chedata::PsatDataSet getCheSubSet(const chedata::PsatDataSet &cheData, const uvec &busIdx) {
			assert(cheData.isFormatted);
			uvec busTag(cheData.nBus, fill::zeros);
			busTag(busIdx).fill(1);
			uvec isw = find(busTag(cheData.get_sws_busNumber_vec() - 1) == 1);
			uvec ipv = find(busTag(cheData.get_pvs_busNumber_vec() - 1) == 1);
			uvec ipq = find(busTag(cheData.get_pqs_busNumber_vec() - 1) == 1);
			uvec ishunt = find(busTag(cheData.get_shunts_busNumber_vec() - 1) == 1);
			uvec iline= find((busTag(cheData.get_lines_fromBus_vec() - 1) == 1)% (busTag(cheData.get_lines_toBus_vec() - 1) == 1));
			uvec ipl = find(busTag(cheData.get_pls_busNumber_vec() - 1) == 1);
			uvec iind = find(busTag(cheData.get_inds_busNumber_vec() - 1) == 1);
			uvec isyn= find(busTag(cheData.get_syns_busNumber_vec() - 1) == 1);
			uvec synTag(cheData.nSyn, fill::zeros);
			synTag(isyn).fill(1);
			uvec itg = find(synTag(cheData.get_tgs_synNumber_vec() - 1) == 1);
			uvec iexc = find(synTag(cheData.get_excs_synNumber_vec() - 1) == 1);

			chedata::PsatDataSet newCheData;
			newCheData.nBus = busTag.n_rows;
			if (newCheData.nBus > 0) {
				newCheData.buses = new chedata::Bus[newCheData.nBus];
				for (int i = 0; i < newCheData.nBus; i++) {
					newCheData.buses[i] = chedata::Bus(cheData.buses[busIdx(i)]);
				}
			}
			newCheData.nSw = isw.n_rows;
			if (newCheData.nSw > 0) {
				newCheData.sws = new chedata::SW[newCheData.nSw];
				for (int i = 0; i < newCheData.nSw; i++) {
					newCheData.sws[i] = chedata::SW(cheData.sws[isw(i)]);
				}
			}
			newCheData.nPv = ipv.n_rows;
			if (newCheData.nPv > 0) {
				newCheData.pvs = new chedata::PV[newCheData.nPv];
				for (int i = 0; i < newCheData.nPv; i++) {
					newCheData.pvs[i] = chedata::PV(cheData.pvs[ipv(i)]);
				}
			}
			newCheData.nPq = ipq.n_rows;
			if (newCheData.nPq > 0) {
				newCheData.pqs = new chedata::PQ[newCheData.nPq];
				for (int i = 0; i < newCheData.nPq; i++) {
					newCheData.pqs[i] = chedata::PQ(cheData.pqs[ipq(i)]);
				}
			}
			newCheData.nShunt = ishunt.n_rows;
			if (newCheData.nShunt > 0) {
				newCheData.shunts = new chedata::Shunt[newCheData.nShunt];
				for (int i = 0; i < newCheData.nShunt; i++) {
					newCheData.shunts[i] = chedata::Shunt(cheData.shunts[ishunt(i)]);
				}
			}
			newCheData.nLine = iline.n_rows;
			if (newCheData.nLine > 0) {
				newCheData.lines = new chedata::Line[newCheData.nLine];
				for (int i = 0; i < newCheData.nLine; i++) {
					newCheData.lines[i] = chedata::Line(cheData.lines[iline(i)]);
				}
			}
			newCheData.nPl = ipl.n_rows;
			if (newCheData.nPl > 0) {
				newCheData.pls = new chedata::Pl[newCheData.nPl];
				for (int i = 0; i < newCheData.nPl; i++) {
					newCheData.pls[i] = chedata::Pl(cheData.pls[ipl(i)]);
				}
			}
			newCheData.nSyn = isyn.n_rows;
			if (newCheData.nSyn > 0) {
				newCheData.syns = new chedata::Syn[newCheData.nSyn];
				for (int i = 0; i < newCheData.nSyn; i++) {
					newCheData.syns[i] = chedata::Syn(cheData.syns[isyn(i)]);
				}
			}
			newCheData.nInd = iind.n_rows;
			if (newCheData.nInd > 0) {
				newCheData.inds = new chedata::Ind[newCheData.nInd];
				for (int i = 0; i < newCheData.nInd; i++) {
					newCheData.inds[i] = chedata::Ind(cheData.inds[iind(i)]);
				}
			}
			newCheData.nTg = itg.n_rows;
			if (newCheData.nTg > 0) {
				newCheData.tgs = new chedata::Tg[newCheData.nTg];
				for (int i = 0; i < newCheData.nTg; i++) {
					newCheData.tgs[i] = chedata::Tg(cheData.tgs[itg(i)]);
				}
			}
			newCheData.nExc = iexc.n_rows;
			if (newCheData.nExc > 0) {
				newCheData.excs = new chedata::Exc[newCheData.nExc];
				for (int i = 0; i < newCheData.nExc; i++) {
					newCheData.excs[i] = chedata::Exc(cheData.excs[iexc(i)]);
				}
			}
			newCheData.renumberBuses();
			return newCheData;
		}

		void dfsSearch(std::list<int>& subgraph,int v, bool visited[],std::vector<std::list<int>>& adjacent){
			visited[v]=true;

			subgraph.push_back(v);
			for (std::list<int>::iterator iter=adjacent[v].begin();iter!=adjacent[v].end();iter++){
				if (!visited[*iter]){
					dfsSearch(subgraph, *iter, visited, adjacent);
				}
			}
		}

		std::vector<std::list<int>> connectedComponents(std::vector<std::list<int>>& adjacent){
			std::vector<std::list<int>> res;
			size_t nv=adjacent.size();
			bool* visited = new bool[nv];
			for (size_t v = 0; v < nv; v++){
				visited[v] = false;
			}

			for (size_t v=0;v<nv;v++){
				if (!visited[v]){
					std::list<int> conncomp;
					dfsSearch(conncomp,v,visited,adjacent);
					res.push_back(conncomp);
				}				
			}
			delete [] visited;

			return res;
		}

		uvec CheCompUtil::searchIslands(const chedata::PsatDataSet &cheData) {
			assert(cheData.isFormatted);
			// uvec islandCode = cheData.get_buses_busNumber_vec() - 1;
			// uvec iOnLine = find(cheData.get_lines_status_vec() == 1);
			// uvec ifr = cheData.get_lines_fromBus_vec()(iOnLine) - 1;
			// uvec ito = cheData.get_lines_toBus_vec()(iOnLine) - 1;
			// int nLine = iOnLine.n_rows;
			// for (int i = 0; i < nLine; i++) {
			// 	int frIsland = islandCode(ifr(i));
			// 	int toIsland = islandCode(ito(i));

			// 	if (frIsland < toIsland) {
			// 		uvec iToIsland = find(islandCode == toIsland);
			// 		islandCode(iToIsland).fill(frIsland);
			// 	}
			// 	else if (frIsland > toIsland) {
			// 		uvec iFrIsland = find(islandCode == frIsland);
			// 		islandCode(iFrIsland).fill(toIsland);
			// 	}
			// }
			// uvec cntIsland = islandCode - (cheData.get_buses_busNumber_vec() - 1);
			// uvec islandBusNumber = find(cntIsland == 0);
			// unsigned int nIslands = islandBusNumber.n_rows;
			// uvec islands = islandCode;
			// for (int i = 0; i < nIslands; i++) {
			// 	uvec busesInIsland = find(islandCode == islandBusNumber(i));
			// 	islands(busesInIsland).fill(i);
			// }

			// return islands;

			uvec iOnLine = find(cheData.get_lines_status_vec() == 1);
			uvec ifr = cheData.get_lines_fromBus_vec()(iOnLine) - 1;
			uvec ito = cheData.get_lines_toBus_vec()(iOnLine) - 1;
			size_t nLine = iOnLine.n_rows;

			std::vector<std::list<int>> adjacent(cheData.nBus);
			for (size_t i = 0; i < nLine; i++) {
				adjacent[ifr(i)].push_back(ito(i));
				adjacent[ito(i)].push_back(ifr(i));
			}

			std::vector<std::list<int>> connComps=connectedComponents(adjacent);

			uvec islands(cheData.nBus);
			islands.fill(-1);

			for (size_t i=0; i<connComps.size();i++){
				for (std::list<int>::iterator iter=connComps[i].begin();iter!=connComps[i].end();iter++){
					islands(*iter)=i;
				}
			}

			return islands;
		}

		list<chedata::PsatDataSet> CheCompUtil::splitIslands(const chedata::PsatDataSet &cheData, const uvec& islands) {
			assert(cheData.isFormatted);
			int nIslands = islands.max() + 1;

			list< chedata::PsatDataSet> islandDataList;
			for (int i = 0; i < nIslands; i++) {
				uvec busesInIsland = find(islands == i);
				islandDataList.push_back(getCheSubSet(cheData, busesInIsland));
			} 

			return islandDataList;
		}

	}
}