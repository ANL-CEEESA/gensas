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
#include "ChePFCalculator.h"
#include "CheCompUtil.h"
#include "matio.h"
//#include "slu_ddefs.h"

//#define DEBUG_VERBOSE
#ifdef DEBUG_VERBOSE
#define QN(x) #x
#define QUOTE(x) QN(x)
#define DEBUG_PRINT_MAT(m) (m).print(QUOTE(m));
#else
#define DEBUG_PRINT_MAT(m)
#endif


namespace che {
	namespace core {
		ChePfEmbedSystem::ChePfEmbedSystem(const chedata::PsatDataSet &sys)
			:CheSingleEmbedSystem(CheState(sys), sys, 0.0) {
			initState.state(initState.stateIdx.vrIdx).fill(1.0);
			initState.state(initState.stateIdx.mEfIdx).fill(1.0);
		}

		ChePfEmbedSystem::ChePfEmbedSystem(const chedata::PsatDataSet &sys, const CheState &st, double alpha)
			: CheSingleEmbedSystem(st, sys, alpha) {}

		vec ChePfEmbedSystem::calcEqBalance(CheSolution* sol, double alpha) {
			double absA = startAlpha + alpha;

			vec solVal = sol->getSolValue(alpha);
			CheState curState(baseSys, solVal);

			uvec busType(baseSys.nBus, fill::zeros);
			busType(C_IDX(baseSys.get_pvs_busNumber_vec())).fill(1);
			busType(C_IDX(baseSys.get_sws_busNumber_vec())).fill(2);
			vec pVec(baseSys.nBus, fill::zeros);
			vec qVec(baseSys.nBus, fill::zeros);

			pVec(C_IDX(baseSys.get_pvs_busNumber_vec())) += absA * baseSys.get_pvs_P_vec();
			pVec(C_IDX(baseSys.get_pqs_busNumber_vec())) -= absA * baseSys.get_pqs_P_vec();
			qVec(C_IDX(baseSys.get_pqs_busNumber_vec())) -= absA * baseSys.get_pqs_Q_vec();
			if (baseSys.nPl > 0) {
				pVec(C_IDX(baseSys.get_pls_busNumber_vec())) -= absA * baseSys.get_pls_P_vec() %
					conv_to<vec>::from(baseSys.get_pls_status_vec());
				qVec(C_IDX(baseSys.get_pls_busNumber_vec())) -= absA * baseSys.get_pls_Q_vec() %
					conv_to<vec>::from(baseSys.get_pls_status_vec());
			}
			cx_vec V = cx_vec(curState.getSubVec(curState.stateIdx.vrIdx), curState.getSubVec(curState.stateIdx.viIdx));
			sp_cx_mat Y = yMatrix.Ytr;
			cx_vec Ysh = yMatrix.Ysh;
			cx_vec YshShunt(baseSys.get_shunts_g_vec(), baseSys.get_shunts_b_vec());
			
			Ysh(C_IDX(baseSys.get_shunts_busNumber_vec())) += YshShunt;
			if (baseSys.nPl > 0) {
				Ysh(C_IDX(baseSys.get_pls_busNumber_vec())) += cx_vec(baseSys.get_pls_g_vec(), baseSys.get_pls_b_vec());
			}
			Y.diag() += absA * Ysh;
			cx_vec IInj = Y * V;
			IInj(C_IDX(baseSys.get_pls_busNumber_vec())) += absA *
				cx_vec(baseSys.get_pls_Ip_vec(), -baseSys.get_pls_Iq_vec()) %
				conv_to<vec>::from(baseSys.get_pls_status_vec()) %
				V(C_IDX(baseSys.get_pls_busNumber_vec())) / abs(V(C_IDX(baseSys.get_pls_busNumber_vec())));
			cx_vec SInjRHS = V % conj(IInj);
			if (baseSys.nSyn > 0) {
				uvec synIdx = C_IDX(baseSys.get_syns_busNumber_vec());
				vec Efq = curState.getSubVec(curState.stateIdx.mEfIdx);
				vec Efd = 0 * Efq;
				vec d = curState.getSubVec(curState.stateIdx.mDeltaIdx);
				vec cosd = cos(d);
				vec sind = sin(d);
				vec Cm = real(V(synIdx));
				vec Dm = imag(V(synIdx));
				vec Vd = sind % Cm - cosd % Dm;
				vec Vq = cosd % Cm + sind % Dm;
				vec Rs = baseSys.get_syns_ra_vec();
				vec Xq = baseSys.get_syns_xq_vec();
				vec Xd = baseSys.get_syns_xd_vec();
				vec Id = (Rs % (Efd - Vd) + Xq % (Efq - Vq)) / (Rs%Rs + Xq % Xd);
				vec Iq = (-Xd % (Efd - Vd) + Rs % (Efq - Vq)) / (Rs%Rs + Xq % Xd);
				cx_vec Ig(sind % Id + cosd % Iq, -cosd % Id + sind % Iq);
				SInjRHS(synIdx) -= V(synIdx) % conj(Ig);
			}
			vec diffTind;
			if (baseSys.nInd > 0) {
				// single-cage
				// TODO: deal with double-cage
				vec s = curState.getSubVec(curState.stateIdx.sIdx);
				uvec indIdx = C_IDX(baseSys.get_inds_busNumber_vec());
				vec R1 = baseSys.get_inds_rs_vec();
				vec X1 = baseSys.get_inds_xs_vec();
				vec Xm = baseSys.get_inds_xm_vec();
				cx_vec Ym(0 / Xm, -absA / Xm);
				vec R2 = baseSys.get_inds_rr1_vec();
				vec X2 = baseSys.get_inds_xr1_vec();
				vec T0 = baseSys.get_inds_Ta_vec() + baseSys.get_inds_Tb_vec() + baseSys.get_inds_Tc_vec();
				vec T1 = -baseSys.get_inds_Tb_vec() - 2 * baseSys.get_inds_Tc_vec();
				vec T2 = baseSys.get_inds_Tc_vec();
				cx_vec Y2 = s / (cx_vec(R2, s%X2));
				cx_vec Ytotal = (Ym + Y2) / (cx_vec(R1, X1) % (Ym + Y2) + 1.0);
				cx_vec ILind = V(indIdx) % Ytotal;
				SInjRHS(indIdx) += V(indIdx) % conj(ILind);

				cx_vec Veind = V(indIdx) - ILind % cx_vec(R1, X1);
				cx_vec IRs = Veind % Y2;
				diffTind = real(IRs%conj(IRs) % R2 - absA * (T0 + s % (T1 + s % T2)) % s);
			}
			else {
				diffTind = vec();
			}

			//qVec(find(busType > 0)).fill(0.);
			vec Q = curState.getSubVec(curState.stateIdx.qIdx);
			SInjRHS -= cx_vec(pVec, qVec + Q);
			SInjRHS(find(busType == 2)).fill(0.);

			vec diffV;
			if (baseSys.nPv > 0) {
				uvec pvIdx = C_IDX(baseSys.get_pvs_busNumber_vec());
				diffV = absA * (baseSys.get_pvs_vMag_vec() % baseSys.get_pvs_vMag_vec() - 1) + 1 - real(V(pvIdx) % conj(V(pvIdx)));
			}
			else {
				diffV = vec();
			}

			return join_cols(real(SInjRHS), imag(SInjRHS), diffV, diffTind);
		}

		CheSingleEmbedSystem* ChePfEmbedSystem::getNewEmbeddedSystem(const CheState &st, double alpha) {
			CheSingleEmbedSystem* embSys = new ChePfEmbedSystem(baseSys, st, startAlpha + alpha);
			return embSys;
		}

		ChePfCalculator::ChePfCalculator(const chedata::PsatDataSet &sys,
			const CheCompOptions& compOpt,
			const vec &ef, const vec &pm) :
			ChePfCalculator(sys,compOpt,CheCompUtil::searchIslands(sys),ef,pm){			
		}

		ChePfCalculator::ChePfCalculator(const chedata::PsatDataSet &sys,
			const CheCompOptions& compOpt, const uvec& islands,
			const vec &ef, const vec &pm) :
			AbstractCheCalculator(sys, compOpt) {
			if (ef.n_rows == sys.nSyn) {
				this->paraEf = ef;
			}
			else if (!ef.empty()) {
				this->paraEf = vec(sys.nSyn > 0 ? sys.nSyn : 0).fill(ef(0));
			}
			else {
				this->paraEf = vec(sys.nSyn > 0 ? sys.nSyn : 0).fill(1.2);
			}
			DEBUG_PRINT_MAT(pm)
			if (pm.n_rows == sys.nSyn) {
				this->paraPm = pm;
			}
			else if (!pm.empty()) {
				this->paraPm = vec(sys.nSyn > 0 ? sys.nSyn : 0).fill(pm(0));
			}
			else {
				this->paraPm = vec(sys.nSyn > 0 ? sys.nSyn : 0).fill(0.0);
			}
			// Temporary
			if (!this->paraPm.empty()) {
				this->pShare = vec(paraPm.n_rows).fill(1.0 / paraPm.n_rows);
			}
			else {
				this->pShare = vec();
			}

			if (islands.n_rows!=this->baseSys.nBus){
				this->islands=CheCompUtil::searchIslands(this->baseSys);
			}else{
				this->islands=islands;
			}

			// this->baseSys = regulateIsland(this->baseSys);

			perm_c=NULL;
			perm_r=NULL;
			etree=NULL;
		}

		CheSingleEmbedSystem* ChePfCalculator::getInitSystem(const chedata::PsatDataSet &sys) {
			CheSingleEmbedSystem* embSys = new ChePfEmbedSystem(sys);
			return embSys;
		}

		chedata::PsatDataSet ChePfCalculator::regulateIsland(const chedata::PsatDataSet &sys) {
			chedata::PsatDataSet newSys(sys);
			if (!newSys.isFormatted) {
				newSys.renumberBuses();
			}
			int nbus = newSys.nBus;
			int nSyn = newSys.nSyn > 0 ? newSys.nSyn : 0;
			if (nSyn > 0) {
				uvec synTag(nbus, fill::zeros);
				uvec synIdx = C_IDX(newSys.get_syns_busNumber_vec());
				synTag(synIdx).fill(1);
				uvec swTag = synTag(C_IDX(newSys.get_sws_busNumber_vec()));
				uvec pvTag = synTag(C_IDX(newSys.get_pvs_busNumber_vec()));
				int nSw= newSys.nSw > 0 ? newSys.nSw : 0;
				int nPv = newSys.nPv > 0 ? newSys.nPv : 0;
				int nSwNew = nSw-sum(swTag);
				int nPvNew = nPv-sum(pvTag);
				uvec swIdxNew = find(swTag == 0);
				uvec pvIdxNew = find(pvTag == 0);
				
				vec extraP(nbus, fill::zeros);
				extraP(C_IDX(newSys.get_pvs_busNumber_vec())) += newSys.get_pvs_P_vec();
				uvec nm(nbus, fill::zeros);
				nm(synIdx) += 1;
				vec extraPg = extraP(synIdx) / nm(synIdx);
				this->paraPm += extraPg;
				DEBUG_PRINT_MAT(this->paraPm)

				if (newSys.nSw > 0) {
					if (nSwNew <= 0) {
						newSys.nSw = -1;
						delete[] newSys.sws;
					}
					else if (nSwNew < newSys.nSw) {
						chedata::SW *psw = new chedata::SW[nSwNew];
						for (int i = 0; i < nSwNew; i++) {
							psw[i] = chedata::SW(newSys.sws[swIdxNew(i)]);
						}
						newSys.nSw = nSwNew;
						delete[] newSys.sws;
						newSys.sws = psw;
					}
				}
				if (newSys.nPv > 0) {
					if (nPvNew <= 0) {
						newSys.nPv = -1;
						delete[] newSys.pvs;
					}
					else if (nPvNew < newSys.nPv) {
						chedata::PV *ppv = new chedata::PV[nPvNew];
						for (int i = 0; i < nPvNew; i++) {
							ppv[i] = chedata::PV(newSys.pvs[pvIdxNew(i)]);
						}
						newSys.nPv = nPvNew;
						delete[] newSys.pvs;
						newSys.pvs = ppv;
					}
				}
			}
			return newSys;
		}

		static mat getTaylorPolynomialsSin(const vec &d0, int n) {
			if (n > 4) n = 4;
			if (n < 0) n = 0;
			vec sind0 = sin(d0);
			vec cosd0 = cos(d0);
			mat sinp(d0.n_rows, n + 1,fill::zeros);
			if (n >= 0) sinp.col(0) = sind0;
			if (n >= 1) {
				sinp.col(0) += -cosd0%d0;
				sinp.col(1) += cosd0;
			}
			if (n >= 2) {
				sinp.col(0) += -sind0%d0%d0/2.0;
				sinp.col(1) += sind0%d0;
				sinp.col(2) += -sind0/2.0;
			}
			if (n >= 3) {
				sinp.col(0) += cosd0 % d0%d0%d0/6.0;
				sinp.col(1) += -cosd0 % d0%d0 / 2.0;
				sinp.col(2) += cosd0 % d0 / 2.0;
				sinp.col(3) += -cosd0 / 6.0;
			}
			if (n >= 4) {
				sinp.col(0) += sind0 % d0%d0%d0%d0 / 24.0;
				sinp.col(1) += -sind0 % d0%d0%d0 / 6.0;
				sinp.col(2) += sind0 % d0%d0 / 4.0;
				sinp.col(3) += -sind0 % d0 / 6.0;
				sinp.col(4) += sind0 / 24.0;
			}
			return sinp;
		}

		static mat getTaylorPolynomialsCos(const vec &d0, int n) {
			if (n > 4) n = 4;
			if (n < 0) n = 0;
			vec sind0 = sin(d0);
			vec cosd0 = cos(d0);
			mat cosp(d0.n_rows, n + 1, fill::zeros);
			if (n >= 0) cosp.col(0) = cosd0;
			if (n >= 1) {
				cosp.col(0) += sind0 % d0;
				cosp.col(1) += -sind0;
			}
			if (n >= 2) {
				cosp.col(0) += -cosd0 % d0%d0 / 2.0;
				cosp.col(1) += cosd0 % d0;
				cosp.col(2) += -cosd0 / 2.0;
			}
			if (n >= 3) {
				cosp.col(0) += -sind0 % d0%d0%d0 / 6.0;
				cosp.col(1) += sind0 % d0%d0 / 2.0;
				cosp.col(2) += -sind0 % d0 / 2.0;
				cosp.col(3) += sind0 / 6.0;
			}
			if (n >= 4) {
				cosp.col(0) += cosd0 % d0%d0%d0%d0 / 24.0;
				cosp.col(1) += -cosd0 % d0%d0%d0 / 6.0;
				cosp.col(2) += cosd0 % d0%d0 / 4.0;
				cosp.col(3) += -cosd0 % d0 / 6.0;
				cosp.col(4) += cosd0 / 24.0;
			}
			return cosp;
		}

		static void bicgstab(const sp_mat& A, vec& x, vec& b, superlu::trans_t trans, 
			superlu::SuperMatrix *L, superlu::SuperMatrix *U,int * perm_c, int *perm_r, 
			superlu::SuperLUStat_t *stat, int *info, int max_it, double tol, double* err, int* iter, int* flag){
			
			*iter =0;
			*flag=0;
			double bnrm2=norm(b);
			// A.print("A");
			vec kk=A*x;
			// b.print("b");
			// kk.print("kk");
			// x.print("x");
			// b.print("b");
			vec r(b);
			r-=kk;
			// r.print("r");
			*err=norm(r)/bnrm2;
			if (*err<tol){
				return;
			}

			double omega=1.0;
			vec r_tld=r;
			double beta=0.0;
			double rho_1=1.0;
			double rho=1.0;
			double alpha=0.0;
			vec p,v,s,t,p_hat,s_hat;
			
			superlu_opts superlu_opts_default;
			superlu::superlu_options_t  options;
			sp_auxlib::set_superlu_opts(options, superlu_opts_default);
			options.IterRefine=superlu::NOREFINE;
			options.RefineInitialized=superlu::NO;
			
			for(;*iter<max_it;*iter+=1){
				// iteration of the algorithm
				rho = dot(r_tld,r);
				if (rho==0){
					break;
				}
				if((*iter)>0){
					beta=(rho/rho_1)*(alpha/omega);
					p=r+beta*(p-omega*v);
				}else{
					p=r;
				}
				// p.print("p");
				
				vec xx=p;
				superlu::SuperMatrix superX1;  arrayops::inplace_set(reinterpret_cast<char*>(&superX1), char(0), sizeof(superlu::SuperMatrix));

				const bool status_x1 = sp_auxlib::wrap_to_supermatrix(superX1, xx);
				arma_wrapper(dgstrs)(options.Trans, L, U, perm_c, perm_r, &superX1, stat, info);
				p_hat=xx;
				// p_hat.print("p_hat");
				sp_auxlib::destroy_supermatrix(superX1);

				v=A*p_hat;
				// v.print("v");
				alpha=rho/dot(r_tld,v);
				x=x+alpha*p_hat;
				// x.print("x");
				s=r-alpha*v;
				// s.print("s");
				*err=norm(s)/bnrm2;
				if (*err<tol){					
					break;
				}

				xx=s;
				superlu::SuperMatrix superX2;  arrayops::inplace_set(reinterpret_cast<char*>(&superX2), char(0), sizeof(superlu::SuperMatrix));

				const bool status_x2 = sp_auxlib::wrap_to_supermatrix(superX2, xx);
				arma_wrapper(dgstrs)(options.Trans, L, U, perm_c, perm_r, &superX2, stat, info);
				s_hat=xx;
				sp_auxlib::destroy_supermatrix(superX2);

				t=A*s_hat;
				// t.print("t");
				omega=dot(t,s)/dot(t,t);
				x=x+omega*s_hat;
				// x.print("x");
				r=s-omega*t;
				// r.print("r");
				*err=norm(r)/bnrm2;

				if((*err)<=tol){
					break;
				}
				if (omega==0.0){
					break;
				}
				rho_1=rho;
			}
			// x.print("x_internal");

			vec chk=b-A*x;
			// chk.print("chk");

			if((*err)<=tol||all(s<=tol)){
				if (all(s<=tol)){
					*err=norm(s)/bnrm2;
				}
				*flag=0;
			}else if(omega==0.0){
				*flag=-2;
			}else if (rho=0.0){
				*flag=-1;
			}else{
				*flag=1;
			}
			
		}


		// static void bicgstab2(const sp_mat& A, vec& x, vec& b, superlu::trans_t trans, 
		// 	superlu::SuperMatrix *L, superlu::SuperMatrix *U,int * perm_c, int *perm_r, 
		// 	superlu::SuperLUStat_t *stat, int *info, int max_it, double tol, double* err, int* iter, int* flag){
			
		// 	*iter =0;
		// 	*flag=0;
		// 	double bnrm2=norm(b);
		// 	vec kk=A*x;
		// 	// b.print("b");
		// 	// kk.print("kk");
		// 	x.print("x");
		// 	vec r(b);
		// 	r-=kk;
		// 	r.print("r");
		// 	*err=norm(r)/bnrm2;
		// 	if (*err<tol){
		// 		return;
		// 	}

		// 	double omega=1.0;
		// 	vec r_tld=r;
		// 	double beta=0.0;
		// 	double rho_1=1.0;
		// 	double rho=1.0;
		// 	double alpha=0.0;
		// 	vec p,v,s,t;
			
		// 	superlu_opts superlu_opts_default;
		// 	superlu::superlu_options_t  options;
		// 	sp_auxlib::set_superlu_opts(options, superlu_opts_default);
		// 	options.IterRefine=superlu::NOREFINE;
		// 	options.RefineInitialized=superlu::NO;
			
		// 	for(;*iter<max_it;*iter+=1){
		// 		// iteration of the algorithm
		// 		rho = sum(r_tld%r);
		// 		if (rho==0){
		// 			break;
		// 		}
		// 		if((*iter)>0){
		// 			beta=(rho/rho_1)*(alpha/omega);
		// 			p=r+beta*(p-omega*v);
		// 		}else{
		// 			p=r;
		// 		}
				
		// 		vec xx=p;
		// 		superlu::SuperMatrix superX1;  arrayops::inplace_set(reinterpret_cast<char*>(&superX1), char(0), sizeof(superlu::SuperMatrix));

		// 		const bool status_x1 = sp_auxlib::wrap_to_supermatrix(superX1, xx);
		// 		arma_wrapper(dgstrs)(options.Trans, L, U, perm_c, perm_r, &superX1, stat, info);
		// 		vec p_hat=p;
		// 		sp_auxlib::destroy_supermatrix(superX1);

		// 		v=A*p_hat;
		// 		alpha=rho/sum(r_tld%v);
		// 		s=r-alpha*v;
		// 		if (norm(s)<tol){
		// 			x+=alpha*p_hat;
		// 			*err=norm(s)/bnrm2;
		// 			break;
		// 		}

		// 		xx=s;
		// 		superlu::SuperMatrix superX2;  arrayops::inplace_set(reinterpret_cast<char*>(&superX2), char(0), sizeof(superlu::SuperMatrix));

		// 		const bool status_x2 = sp_auxlib::wrap_to_supermatrix(superX2, xx);
		// 		arma_wrapper(dgstrs)(options.Trans, L, U, perm_c, perm_r, &superX2, stat, info);
		// 		vec s_hat=s;
		// 		sp_auxlib::destroy_supermatrix(superX2);

		// 		t=A*s_hat;
		// 		omega=sum(t%s)/sum(t%t);
		// 		x+=alpha*p_hat+omega*s_hat;
		// 		r=s-omega*t;
		// 		*err=norm(r)/bnrm2;

		// 		if((*err)<=tol){
		// 			break;
		// 		}
		// 		if (omega==0.0){
		// 			break;
		// 		}
		// 		rho_1=rho;
		// 	}
		// 	x.print("x_internal");

		// 	vec chk=b-A*x;
		// 	chk.print("chk");

		// 	if((*err)<=tol||all(s<=tol)){
		// 		if (all(s<=tol)){
		// 			*err=norm(s)/bnrm2;
		// 		}
		// 		*flag=0;
		// 	}else if(omega==0.0){
		// 		*flag=-2;
		// 	}else if (rho=0.0){
		// 		*flag=-1;
		// 	}else{
		// 		*flag=1;
		// 	}
			
		// }

		CheSolution* ChePfCalculator::getCheSolution() {
			int nlvl = compOpt.nLvl;

			CheSingleEmbedSystem* embSys = cheList.back();
			chedata::PsatDataSet baseSys(embSys->baseSys);
			double sAlpha = embSys->startAlpha;
			CheStateIdx stateIdx = embSys->initState.stateIdx;
			int nbus = baseSys.nBus;

			// uvec islands = CheCompUtil::searchIslands(baseSys);
			int nIslands = this->islands.max() + 1;

			uvec busType(nbus, fill::zeros);
			busType(C_IDX(baseSys.get_pvs_busNumber_vec())).fill(1);
			busType(C_IDX(baseSys.get_sws_busNumber_vec())).fill(2);
			uvec isw = find(busType == 2);
			uvec ipv = find(busType == 1);
			uvec ipq = find(busType == 0);
			int npq = ipq.n_rows;
			int npv = ipv.n_rows;
			DEBUG_PRINT_MAT(busType)

			sp_cx_mat Y = embSys->yMatrix.Ytr;
			cx_vec Ysh = embSys->yMatrix.Ysh;
			cx_vec YshShunt(baseSys.get_shunts_g_vec(), baseSys.get_shunts_b_vec());
			DEBUG_PRINT_MAT(YshShunt)

			Ysh(C_IDX(baseSys.get_shunts_busNumber_vec())) += YshShunt;
			if (baseSys.nPl > 0) {
				Ysh(C_IDX(baseSys.get_pls_busNumber_vec())) += cx_vec(baseSys.get_pls_g_vec(), baseSys.get_pls_b_vec());
			}
			Y.diag() += sAlpha * Ysh;

			vec pVec(baseSys.nBus, fill::zeros);
			vec qVec(baseSys.nBus, fill::zeros);
			vec vMagVec(nbus, fill::zeros);

			pVec(C_IDX(baseSys.get_pvs_busNumber_vec())) += baseSys.get_pvs_P_vec();
			pVec(C_IDX(baseSys.get_pqs_busNumber_vec())) -= baseSys.get_pqs_P_vec();
			qVec(C_IDX(baseSys.get_pqs_busNumber_vec())) -= baseSys.get_pqs_Q_vec();
			vMagVec(C_IDX(baseSys.get_pvs_busNumber_vec())) = baseSys.get_pvs_vMag_vec() % baseSys.get_pvs_vMag_vec();
			vMagVec(C_IDX(baseSys.get_sws_busNumber_vec())) = baseSys.get_sws_vMag_vec() % baseSys.get_sws_vMag_vec();

			if (baseSys.nPl > 0) {
				pVec(C_IDX(baseSys.get_pls_busNumber_vec())) -= baseSys.get_pls_P_vec() %
					conv_to<vec>::from(baseSys.get_pls_status_vec());
				qVec(C_IDX(baseSys.get_pls_busNumber_vec())) -= baseSys.get_pls_Q_vec() %
					conv_to<vec>::from(baseSys.get_pls_status_vec());
			}
			vec pVec0 = sAlpha * pVec;
			vec qVec0 = sAlpha * qVec;
			//qVec0(ipv) = embSys->initState.getSubVec(stateIdx.qIdx)(ipv);
			vec VspSq2(nbus, fill::zeros);
			VspSq2(ipv) = vMagVec(ipv) - 1;
			VspSq2(isw) = vMagVec(isw) - 1;

			cx_mat V(nbus, nlvl+1, fill::zeros);
			cx_vec V0(embSys->initState.getSubVec(stateIdx.vrIdx), embSys->initState.getSubVec(stateIdx.viIdx));
			V.col(0) = V0;
			V(isw, uvec(1).fill(1)) = cx_vec(baseSys.get_sws_vMag_vec() % cos(datum::pi / 180.0*baseSys.get_sws_vAng_vec()),
				baseSys.get_sws_vMag_vec() % sin(datum::pi / 180.0*baseSys.get_sws_vAng_vec())) - 1;
			cx_mat W(nbus, nlvl+1, fill::zeros);
			W.col(0) = 1 / V0;
			mat P(nbus, nlvl+1, fill::zeros);
			P.col(0) = pVec0;
			mat Q(nbus, nlvl+1, fill::zeros);
			mat Qxtra(nbus, nlvl + 1, fill::zeros);
			Q.col(0) = embSys->initState.getSubVec(stateIdx.qIdx);
			Qxtra.col(0) = qVec0;
			P.col(1) = pVec;
			Qxtra.col(1) = qVec;
			//Q(find(busType != 0), regspace<uvec>(1, nlvl)).fill(0.);

			vec C0 = real(V.col(0));
			vec D0 = imag(V.col(0));
			vec E0 = real(W.col(0));
			vec F0 = imag(W.col(0));
			vec Pspec = P.col(0);
			vec Qspec = Q.col(0) + Qxtra.col(0);
			vec V0sq = C0 % C0 + D0 % D0;

			sp_mat C0M(true, join_rows(regspace<uvec>(0, nbus - 1), regspace<uvec>(0, nbus - 1)).t(), C0, nbus, nbus);
			sp_mat D0M(true, join_rows(regspace<uvec>(0, nbus - 1), regspace<uvec>(0, nbus - 1)).t(), D0, nbus, nbus);
			sp_mat E0M(true, join_rows(regspace<uvec>(0, nbus - 1), regspace<uvec>(0, nbus - 1)).t(), E0, nbus, nbus);
			sp_mat F0M(true, join_rows(regspace<uvec>(0, nbus - 1), regspace<uvec>(0, nbus - 1)).t(), F0, nbus, nbus);
			sp_mat P0M(true, join_rows(regspace<uvec>(0, nbus - 1), regspace<uvec>(0, nbus - 1)).t(), Pspec, nbus, nbus);
			sp_mat Q0M(true, join_rows(regspace<uvec>(0, nbus - 1), regspace<uvec>(0, nbus - 1)).t(), Qspec, nbus, nbus);

			sp_mat G = real(Y);
			sp_mat B = imag(Y);

			DEBUG_PRINT_MAT(VspSq2)
			DEBUG_PRINT_MAT(G)
			DEBUG_PRINT_MAT(B)
			DEBUG_PRINT_MAT(P0M)
			DEBUG_PRINT_MAT(Q0M)

			//Ind
			int nInd = baseSys.nInd > 0 ? baseSys.nInd : 0;
			uvec indIdx = C_IDX(baseSys.get_inds_busNumber_vec());
			vec s0 = embSys->initState.getSubVec(stateIdx.sIdx);
			mat s(nInd, nlvl + 1, fill::zeros);
			s.col(0) = s0;
			cx_mat IL(nInd, nlvl + 1, fill::zeros);
			cx_mat IR(nInd, nlvl + 1, fill::zeros);
			cx_mat Vm(nInd, nlvl + 1, fill::zeros);

			vec R1 = baseSys.get_inds_rs_vec();
			vec X1 = baseSys.get_inds_xs_vec();
			vec Xm = baseSys.get_inds_xm_vec();
			vec R2 = baseSys.get_inds_rr1_vec();
			vec X2 = baseSys.get_inds_xr1_vec();
			vec T0 = baseSys.get_inds_Ta_vec() + baseSys.get_inds_Tb_vec() + baseSys.get_inds_Tc_vec();
			vec T1 = -baseSys.get_inds_Tb_vec() - 2 * baseSys.get_inds_Tc_vec();
			vec T2 = baseSys.get_inds_Tc_vec();

			cx_vec Z1(R1, X1);
			cx_vec Ym(0 / Xm, -sAlpha / Xm);
			cx_vec cIndTemp = R2 % Ym + s0 % (Ym%cx_vec(0.*X2, X2) + 1);
			IL.col(0) = V0(indIdx) % cIndTemp / (cx_vec(R2, s0%X2) + Z1 % cIndTemp);
			Vm.col(0) = V0(indIdx) - IL.col(0) % Z1;
			IR.col(0) = Vm.col(0) % s0 / cx_vec(R2, s0%X2);

			vec J0 = real(IR.col(0));
			vec K0 = imag(IR.col(0));
			vec JL0 = real(IL.col(0));
			vec KL0 = imag(IL.col(0));

			cx_vec Yeind0(0. / Xm, -sAlpha / Xm);
			cx_vec Yeind1(0. / Xm, -1. / Xm);
			cx_vec Ye1ind0 = Yeind0 % Z1;
			cx_vec Ye1ind1 = Yeind1 % Z1;
			vec Ge = real(Yeind0);
			vec Be = imag(Yeind0);
			vec kg1e = real(Ye1ind0);
			vec kb1e = imag(Ye1ind0);
			DEBUG_PRINT_MAT(Ge)
			DEBUG_PRINT_MAT(Be)
			DEBUG_PRINT_MAT(kg1e)
			DEBUG_PRINT_MAT(kb1e)
			
			vector<mat> LHS_MatInd_Full;
			vector<mat> RHS_C_Shr;
			mat LHS_MatInd_Shr(nInd, 4, fill::zeros);
			mat LHS_MatInd_Bus(nbus, 4, fill::zeros);

			for (int i = 0; i < nInd; i++) {
				mat LHS_MatInd(5, 7);
				LHS_MatInd << R2(i) << -X2(i)*s0(i) << R1(i)*s0(i) << -X1(i)*s0(i) << -K0(i)*X2(i) - C0(indIdx(i)) + JL0(i)*R1(i) - KL0(i)*X1(i) << -s0(i) << 0.0 << endr
					<< X2(i)*s0(i) << R2(i) << X1(i)*s0(i) << R1(i)*s0(i) << J0(i)*X2(i) - D0(indIdx(i)) + JL0(i)*X1(i) + KL0(i)*R1(i) << 0.0 << -s0(i) << endr
					<< C0(indIdx(i)) - JL0(i)*R1(i) + KL0(i)*X1(i) << D0(indIdx(i)) - KL0(i)*R1(i) - JL0(i)*X1(i) << -J0(i)*R1(i) - K0(i)*X1(i) << -K0(i)*R1(i) + J0(i)*X1(i) << -sAlpha * (T1(i) + 2. * T2(i)*s0(i)) << J0(i) << K0(i) << endr
					<< -1. << 0. << 1. + kg1e(i) << -kb1e(i) << 0. << -Ge(i) << Be(i) << endr
					<< 0. << -1. << kb1e(i) << 1. + kg1e(i) << 0. << -Be(i) << -Ge(i) << endr;
				mat MatIndA = LHS_MatInd(span::all, span(0, 4));
				mat MatIndB = LHS_MatInd(span::all, span(5, 6));
				mat MatInvA = inv(MatIndA);
				mat MadCDtoRest = -MatInvA * MatIndB;

				DEBUG_PRINT_MAT(LHS_MatInd)
				DEBUG_PRINT_MAT(MatIndA)
				DEBUG_PRINT_MAT(MatIndB)
				DEBUG_PRINT_MAT(MatInvA)
				DEBUG_PRINT_MAT(MadCDtoRest)

				RHS_C_Shr.push_back(MatInvA);
				LHS_MatInd_Full.push_back(MadCDtoRest);
				LHS_MatInd_Shr.row(i) = MadCDtoRest.rows(span(2, 3)).as_row();

			}
			LHS_MatInd_Bus.rows(indIdx) += LHS_MatInd_Shr;

			DEBUG_PRINT_MAT(J0)
			DEBUG_PRINT_MAT(K0)
			DEBUG_PRINT_MAT(JL0)
			DEBUG_PRINT_MAT(KL0)
			DEBUG_PRINT_MAT(LHS_MatInd_Shr)
			DEBUG_PRINT_MAT(LHS_MatInd_Bus)

			//ZIP
			int nZip = baseSys.nPl > 0 ? baseSys.nPl : 0;
			uvec zipIdx = C_IDX(baseSys.get_pls_busNumber_vec());
			cx_mat IiL(nZip, nlvl + 1, fill::zeros);
			mat BiL(nZip, nlvl + 1, fill::zeros);
			vec Bi0 = abs(V0(zipIdx));
			vec JI = baseSys.get_pls_Ip_vec();
			vec KI = -baseSys.get_pls_Iq_vec();
			cx_vec Ii0L = sAlpha * cx_vec(JI, KI) % V0(zipIdx) / Bi0;
			vec Ji0L = real(Ii0L);
			vec Ki0L = imag(Ii0L);

			IiL.col(0) = Ii0L;
			BiL.col(0) = Bi0;
			vec Ci0 = real(V0(zipIdx));
			vec Di0 = imag(V0(zipIdx));
			mat LHS_MatZip = join_rows(sAlpha*JI / Bi0 - Ci0 % Ji0L / Bi0 / Bi0,
				-sAlpha * KI / Bi0 - Di0 % Ji0L / Bi0 / Bi0,
				sAlpha*KI / Bi0 - Ci0 % Ki0L / Bi0 / Bi0,
				sAlpha*JI / Bi0 - Di0 % Ki0L / Bi0 / Bi0);
			mat Mat_BZip = join_rows(Ci0 / Bi0, Di0 / Bi0);

			DEBUG_PRINT_MAT(Bi0)
			DEBUG_PRINT_MAT(Ji0L)
			DEBUG_PRINT_MAT(Ki0L)
			DEBUG_PRINT_MAT(LHS_MatZip)
			DEBUG_PRINT_MAT(Mat_BZip)


			//Syn
			int nSyn= baseSys.nSyn > 0 ? baseSys.nSyn : 0;
			uvec synIdx = C_IDX(baseSys.get_syns_busNumber_vec());
			vec Rs = baseSys.get_syns_ra_vec();
			vec Xd = baseSys.get_syns_xd_vec();
			vec Xq = baseSys.get_syns_xq_vec();
			vec Ms = baseSys.get_syns_M_vec();

			mat d(nSyn, nlvl + 1, fill::zeros);
			mat JG(nSyn, nlvl + 1, fill::zeros);
			mat KG(nSyn, nlvl + 1, fill::zeros);
			mat Cd(nSyn, nlvl + 1, fill::zeros);
			mat Sd(nSyn, nlvl + 1, fill::zeros);
			mat Ef(nSyn, nlvl + 1, fill::zeros);
			vec d0 = embSys->initState.getSubVec(stateIdx.mDeltaIdx);
			vec Ef0 = embSys->initState.getSubVec(stateIdx.mEfIdx);
			vec Efd = 0.0*Ef0;
			vec Efq = Ef0;
			vec cosd = cos(d0);
			vec sind = sin(d0);
			vec Cg = C0(synIdx);
			vec Dg = D0(synIdx);
			vec Vd = sind % Cg - cosd % Dg;
			vec Vq = cosd % Cg + sind % Dg;
			vec Id = (Rs % (Efd - Vd) + Xq % (Efq - Vq)) / (Rs%Rs + Xq % Xd);
			vec Iq = (-Xd % (Efd - Vd) + Rs % (Efq - Vq)) / (Rs%Rs + Xq % Xd);
			cx_vec IG0(sind%Id + cosd % Iq, -cosd % Id + sind % Iq);
			d.col(0) = d0;
			JG.col(0) = real(IG0);
			KG.col(0) = imag(IG0);
			Cd.col(0) = cosd;
			Sd.col(0) = sind;
			Ef.col(0) = Ef0;
			Ef.col(1) = this->paraEf - 1;
			vec CG0 = Cg;
			vec DG0 = Dg;
			vec JG0 = JG.col(0);
			vec KG0 = KG.col(0);
			vec Cd0 = Cd.col(0);
			vec Sd0 = Sd.col(0);
			mat Pm(nSyn, nlvl + 1, fill::zeros);
			Pm.col(0) = Efq % Iq + Efd % Id;
			Pm.col(1) = paraPm;
			DEBUG_PRINT_MAT(paraPm)
			DEBUG_PRINT_MAT(this->paraEf)
			DEBUG_PRINT_MAT(Rs)
			DEBUG_PRINT_MAT(Xd)
			DEBUG_PRINT_MAT(Xq)
			DEBUG_PRINT_MAT(Ms)
			DEBUG_PRINT_MAT(d0)
			DEBUG_PRINT_MAT(Ef0)
			DEBUG_PRINT_MAT(Efd)
			DEBUG_PRINT_MAT(Efq)
			DEBUG_PRINT_MAT(cosd)
			DEBUG_PRINT_MAT(sind)
			DEBUG_PRINT_MAT(Cg)
			DEBUG_PRINT_MAT(Dg)
			DEBUG_PRINT_MAT(Vd)
			DEBUG_PRINT_MAT(Vq)
			DEBUG_PRINT_MAT(Id)
			DEBUG_PRINT_MAT(Iq)
			DEBUG_PRINT_MAT(JG)
			DEBUG_PRINT_MAT(KG)
			DEBUG_PRINT_MAT(Cd)
			DEBUG_PRINT_MAT(Sd)
			DEBUG_PRINT_MAT(Ef)
			DEBUG_PRINT_MAT(Pm)

			int nTaylor = 4; //Will be moved to a simulation setting structure
			mat cosp = getTaylorPolynomialsCos(d0, nTaylor);
			mat sinp = getTaylorPolynomialsSin(d0, nTaylor);

			vec A1n(nSyn, fill::zeros);
			vec B1n(nSyn, fill::zeros);
			for (int i = 0; i < nTaylor; i++) {
				A1n = A1n % d0 + (nTaylor - i)*cosp.col(nTaylor - i);
				B1n = B1n % d0 + (nTaylor - i)*sinp.col(nTaylor - i);
			}
			DEBUG_PRINT_MAT(cosp)
			DEBUG_PRINT_MAT(sinp)
			DEBUG_PRINT_MAT(A1n)
			DEBUG_PRINT_MAT(B1n)

			uvec synByIsland = this->islands(synIdx);
			vec sumShare(nIslands, fill::zeros);
			sumShare(synByIsland) += this->pShare;
			vec pSharex = this->pShare / sumShare(synByIsland);
			uvec synRegIdx = regspace<uvec>(0, nSyn - 1);
			sp_mat mAuxIsland(true, join_rows(synRegIdx, synByIsland).t(), pSharex, nSyn, nIslands);
			mat mAuxIslandFull(abs(mAuxIsland));
			uvec idxBal(0,fill::zeros);
			if (!mAuxIslandFull.empty()) {
				idxBal = index_max(mAuxIslandFull, 0).t();
			}
			uvec idxBalSyn = idxBal(synByIsland);

			sp_mat MatG1C(true, join_rows(synRegIdx, synIdx).t(), Cd0, nSyn, nbus);
			sp_mat MatG1D(true, join_rows(synRegIdx, synIdx).t(), Sd0, nSyn, nbus);
			sp_mat MatG2C(true, join_rows(synRegIdx, synIdx).t(), Sd0, nSyn, nbus);
			sp_mat MatG2D(true, join_rows(synRegIdx, synIdx).t(), -Cd0, nSyn, nbus);
			sp_mat MatG5AC(true, join_rows(join_cols(synRegIdx, synRegIdx), join_cols(synIdx, synIdx(idxBalSyn))).t(),
				join_cols(-pSharex(idxBalSyn) % JG0, JG0(idxBalSyn) % pSharex), nSyn, nbus);
			sp_mat MatG5AD(true, join_rows(join_cols(synRegIdx, synRegIdx), join_cols(synIdx, synIdx(idxBalSyn))).t(),
				join_cols(-pSharex(idxBalSyn) % KG0, KG0(idxBalSyn) % pSharex), nSyn, nbus);
			sp_mat MatG5BC(nSyn, nbus);
			sp_mat MatG5BD(nSyn, nbus);
			for (int i = 0; i < idxBal.n_rows; i++) {
				MatG5AC.row(idxBal(i)) = MatG5BC.row(i); // sparse matrix does not allow non-contiguous indexing
				MatG5AD.row(idxBal(i)) = MatG5BD.row(i); // sparse matrix does not allow non-contiguous indexing
			}

			sp_mat MatG1J(true, join_rows(synRegIdx, synRegIdx).t(), Rs%Cd0 + Xd % Sd0, nSyn, nSyn);
			sp_mat MatG1K(true, join_rows(synRegIdx, synRegIdx).t(), Rs%Sd0 - Xd % Cd0, nSyn, nSyn);
			sp_mat MatG1d(true, join_rows(synRegIdx, synRegIdx).t(),
				(CG0 + Rs % JG0 - Xd % KG0) % A1n + (DG0 + Rs % KG0 + Xd % JG0) % B1n, nSyn, nSyn);
			sp_mat MatG2J(true, join_rows(synRegIdx, synRegIdx).t(), Rs%Sd0 - Xq % Cd0, nSyn, nSyn);
			sp_mat MatG2K(true, join_rows(synRegIdx, synRegIdx).t(), -Rs%Cd0 - Xq % Sd0, nSyn, nSyn);
			sp_mat MatG2d(true, join_rows(synRegIdx, synRegIdx).t(),
				(-DG0 - Rs % KG0 - Xq % JG0) % A1n + (CG0 + Rs % JG0 - Xq % KG0) % B1n, nSyn, nSyn);
			sp_mat MatG5AJ(true, join_rows(join_cols(synRegIdx, synRegIdx), join_cols(synRegIdx, idxBalSyn)).t(),
				join_cols(-pSharex(idxBalSyn) % (CG0 + 2.0*JG0%Rs), (CG0(idxBalSyn) + 2.0*JG0(idxBalSyn) % Rs(idxBalSyn)) % pSharex), nSyn, nSyn);
			sp_mat MatG5AK(true, join_rows(join_cols(synRegIdx, synRegIdx), join_cols(synRegIdx, idxBalSyn)).t(),
				join_cols(-pSharex(idxBalSyn) % (DG0 + 2.0*KG0%Rs), (DG0(idxBalSyn) + 2.0*KG0(idxBalSyn) % Rs(idxBalSyn)) % pSharex), nSyn, nSyn);
			sp_mat MatG5Ad(nSyn, nSyn);
			sp_mat MatG5BJ(nIslands, nSyn);
			sp_mat MatG5BK(nIslands, nSyn);
			sp_mat MatG5Bd(true, join_rows(synByIsland, synRegIdx).t(), Ms, nIslands, nSyn);
						
			if (nSyn > 1) {
				for (int i = 0; i < nIslands; i++) {
					MatG5AJ.row(idxBal(i)) = MatG5BJ.row(i);
					MatG5AK.row(idxBal(i)) = MatG5BK.row(i);
					MatG5Ad.row(idxBal(i)) = MatG5Bd.row(i);
				}
			}				

			sp_mat MatGA = join_cols(join_rows(MatG1C, MatG1D),
				join_rows(MatG2C, MatG2D), join_rows(MatG5AC, MatG5AD));
			sp_mat MatGB = join_cols(join_rows(MatG1J, MatG1K, MatG1d),
				join_rows(MatG2J, MatG2K, MatG2d), join_rows(MatG5AJ, MatG5AK, MatG5Ad));
			sp_mat MatGBiA(-spsolve(MatGB, mat(MatGA)));
			sp_mat GTrMat(true, join_rows(synIdx, synRegIdx).t(), vec(nSyn, fill::ones), nbus, nSyn);
			sp_mat MatGTrans = join_cols(join_rows(GTrMat, sp_mat(nbus, nSyn), sp_mat(nbus, nSyn)),
				join_rows(sp_mat(nbus, nSyn), GTrMat, sp_mat(nbus, nSyn)));
			sp_mat MatGCD = MatGTrans * MatGBiA;

			DEBUG_PRINT_MAT(pSharex)
			DEBUG_PRINT_MAT(MatG1C)
			DEBUG_PRINT_MAT(MatG1D)
			DEBUG_PRINT_MAT(MatG2C)
			DEBUG_PRINT_MAT(MatG2D)
			DEBUG_PRINT_MAT(MatG5AC)
			DEBUG_PRINT_MAT(MatG5AD)
			DEBUG_PRINT_MAT(MatG5BC)
			DEBUG_PRINT_MAT(MatG5BD)
			DEBUG_PRINT_MAT(MatG1J)
			DEBUG_PRINT_MAT(MatG1K)
			DEBUG_PRINT_MAT(MatG1d)
			DEBUG_PRINT_MAT(MatG2J)
			DEBUG_PRINT_MAT(MatG2K)
			DEBUG_PRINT_MAT(MatG2d)
			DEBUG_PRINT_MAT(MatG5AJ)
			DEBUG_PRINT_MAT(MatG5AK)
			DEBUG_PRINT_MAT(MatG5Ad)
			DEBUG_PRINT_MAT(MatG5BJ)
			DEBUG_PRINT_MAT(MatG5BK)
			DEBUG_PRINT_MAT(MatG5Bd)
			DEBUG_PRINT_MAT(MatGA)
			DEBUG_PRINT_MAT(MatGB)
			DEBUG_PRINT_MAT(GTrMat)
			DEBUG_PRINT_MAT(MatGTrans)
			DEBUG_PRINT_MAT(MatGCD)

			sp_mat Y11 = -G;
			sp_mat Y12 = B;
			sp_mat Y21 = -B;
			sp_mat Y22 = -G;

			sp_mat YEF11 = P0M;
			sp_mat YEF12 = -Q0M;
			sp_mat YEF21 = -Q0M;
			sp_mat YEF22 = -P0M;

			vec auxDiag0(nbus, fill::zeros);
			auxDiag0(zipIdx) += LHS_MatZip.col(0);
			Y11.diag() -= auxDiag0 + LHS_MatInd_Bus.col(0);
			auxDiag0.fill(0);
			auxDiag0(zipIdx) += LHS_MatZip.col(1);
			Y12.diag() -= auxDiag0 + LHS_MatInd_Bus.col(1);
			auxDiag0.fill(0);
			auxDiag0(zipIdx) += LHS_MatZip.col(2);
			Y21.diag() -= auxDiag0 + LHS_MatInd_Bus.col(2);
			auxDiag0.fill(0);
			auxDiag0(zipIdx) += LHS_MatZip.col(3);
			Y22.diag() -= auxDiag0 + LHS_MatInd_Bus.col(3);

			vec C0i = C0 / V0sq;
			vec D0i = D0 / V0sq;
			vec PCQD = Pspec % C0i + Qspec % D0i;
			vec PDQC = Pspec % D0i - Qspec % C0i;
			Y11.diag() -= PCQD % E0 + PDQC % F0;
			Y12.diag() -= -PCQD % F0 + PDQC % E0;
			Y21.diag() -= PDQC % E0 - PCQD % F0;
			Y22.diag() -= -PDQC % F0 - PCQD % E0;

			sp_mat YLHS = join_cols(join_rows(Y11, Y12), join_rows(Y21, Y22));

			YLHS += MatGCD;

			uvec idxNonSw = find(busType != 2);
			uvec idxStackMat = join_cols(idxNonSw, idxNonSw + nbus);

			sp_mat YLHSx = CheCompUtil::sp_submatrix<double>(YLHS, idxStackMat, idxStackMat);

			/*sp_mat YEFx = join_cols(
				join_rows(
					CheCompUtil::sp_submatrix<double>(YEF11, idxNonSw, idxNonSw),
					CheCompUtil::sp_submatrix<double>(YEF12, idxNonSw, idxNonSw),
					-CheCompUtil::sp_submatrix<double>(F0M, idxNonSw, ipv)),
				join_rows(
					CheCompUtil::sp_submatrix<double>(YEF21, idxNonSw, idxNonSw),
					CheCompUtil::sp_submatrix<double>(YEF22, idxNonSw, idxNonSw),
					-CheCompUtil::sp_submatrix<double>(E0M, idxNonSw, ipv)));
			sp_mat YCDx = join_rows(
				CheCompUtil::sp_submatrix<double>(C0M, ipv, idxNonSw),
				CheCompUtil::sp_submatrix<double>(D0M, ipv, idxNonSw),
				sp_mat(npv, 2 * npq + 3 * npv));
			sp_mat CM = CheCompUtil::sp_submatrix<double>(C0M, idxNonSw, idxNonSw);
			sp_mat DM = CheCompUtil::sp_submatrix<double>(D0M, idxNonSw, idxNonSw);
			sp_mat EM = CheCompUtil::sp_submatrix<double>(E0M, idxNonSw, idxNonSw);
			sp_mat FM = CheCompUtil::sp_submatrix<double>(F0M, idxNonSw, idxNonSw);
			sp_mat EFMx = join_cols(join_rows(EM, -FM), join_rows(FM, EM));
			sp_mat CDMx = join_cols(join_rows(CM, -DM), join_rows(DM, CM));*/

			sp_mat YEFx = join_cols(-CheCompUtil::sp_submatrix<double>(F0M, idxNonSw, ipv), -CheCompUtil::sp_submatrix<double>(E0M, idxNonSw, ipv));
			sp_mat YCDx= join_rows(
				CheCompUtil::sp_submatrix<double>(C0M, ipv, idxNonSw),
				CheCompUtil::sp_submatrix<double>(D0M, ipv, idxNonSw),
				sp_mat(npv, npv));

			sp_mat LHS_mat = join_cols(
				join_rows(YLHSx,YEFx),
				YCDx);

			// superlu settings
			superlu_opts superlu_opts_default;
			superlu::superlu_options_t  options;
			sp_auxlib::set_superlu_opts(options, superlu_opts_default);
			options.IterRefine=superlu::NOREFINE;
			options.RefineInitialized=superlu::NO;
					
			const sp_mat& A =   LHS_mat;
			
			superlu::SuperMatrix superX;  arrayops::inplace_set(reinterpret_cast<char*>(&superX), char(0), sizeof(superlu::SuperMatrix));
			superlu::SuperMatrix superA;  arrayops::inplace_set(reinterpret_cast<char*>(&superA), char(0), sizeof(superlu::SuperMatrix));
			superlu::SuperMatrix superB;  arrayops::inplace_set(reinterpret_cast<char*>(&superB), char(0), sizeof(superlu::SuperMatrix));
			const bool status_a = sp_auxlib::copy_to_supermatrix(superA, A);
			superlu::SuperMatrix superL;  arrayops::inplace_set(reinterpret_cast<char*>(&superL), char(0), sizeof(superlu::SuperMatrix));
			superlu::SuperMatrix superU;  arrayops::inplace_set(reinterpret_cast<char*>(&superU), char(0), sizeof(superlu::SuperMatrix));

			int* perm_c = (int*)superlu::malloc((A.n_cols + 1) * sizeof(int));  
			int* perm_r = (int*)superlu::malloc((A.n_rows + 1) * sizeof(int));
    		int* etree  = (int*) superlu::malloc( (A.n_cols+1) * sizeof(int) );
			double* R    = (double*) superlu::malloc( (A.n_rows+1) * sizeof(double) );
			double* C    = (double*) superlu::malloc( (A.n_cols+1) * sizeof(double) );
			double* ferr = (double*) superlu::malloc( (1+1) * sizeof(double) );
			double* berr = (double*) superlu::malloc( (1+1) * sizeof(double) );
			arrayops::inplace_set(perm_c, 0, A.n_cols + 1);
			arrayops::inplace_set(perm_r, 0, A.n_rows + 1);
			arrayops::inplace_set(etree, 0, A.n_cols + 1);

			arrayops::inplace_set(R,    double(0), A.n_rows+1);
			arrayops::inplace_set(C,    double(0), A.n_cols+1);
			arrayops::inplace_set(ferr, double(0), 1+1);
			arrayops::inplace_set(berr, double(0), 1+1);
    
			superlu::GlobalLU_t glu;
			arrayops::inplace_set(reinterpret_cast<char*>(&glu), char(0), sizeof(superlu::GlobalLU_t));
			
			superlu::mem_usage_t  mu;
			arrayops::inplace_set(reinterpret_cast<char*>(&mu), char(0), sizeof(superlu::mem_usage_t));
    
			superlu::SuperLUStat_t stat;
			superlu::init_stat(&stat);

			char equed[8];       // extra characters for paranoia
			double rpg   = double(0);
			double rcond = double(0);
			int superInfo = 0; // Return code.
			
			char  work[8];
			int  lwork = int(0);  // 0 means superlu will allocate memory
    
			DEBUG_PRINT_MAT(LHS_mat)
			// LOOP Body
			for (int lvl = 0; lvl < nlvl; lvl++) {
				umat seq2 = CheCompUtil::spgetseq(lvl + 1, 2);
				umat seq2p = CheCompUtil::spgetseq(lvl + 2, 2);
				umat seq2m = CheCompUtil::spgetseq(lvl, 2);
				umat seq3 = CheCompUtil::spgetseq(lvl + 1, 3);
				uvec idxSeq2 = any(seq2 == lvl + 1, 1);
				uvec idxSeq3 = any(seq3 == lvl + 1, 1);
				umat seq2R = seq2.rows(find(idxSeq2 == 0));
				umat seq3R = seq3.rows(find(idxSeq3 == 0));

				vec RHSILr(nbus, fill::zeros);
				vec RHSILi(nbus, fill::zeros);

				//LOOP-Ind
				mat rhsBus(nInd, 5, fill::zeros);
				cx_vec rhsM = sum(Vm.cols(seq2R.col(0)) % s.cols(seq2R.col(1)), 1) +
					cx_vec(0.0*X2, -X2) % sum(IR.cols(seq2R.col(0)) % s.cols(seq2R.col(1)), 1);
				vec rhsImod = T1 % s.col(lvl) +
					T2 % sum(s.cols(seq2m.col(0)) % s.cols(seq2m.col(1)), 1) +
					sAlpha * T2%sum(s.cols(seq2R.col(0)) % s.cols(seq2R.col(1)), 1) -
					real(sum(V(indIdx, seq2R.col(0)) % conj(IR.cols(seq2R.col(1))), 1)) +
					real(sum(IL.cols(seq2R.col(0)) % conj(IR.cols(seq2R.col(1))), 1) % Z1);
				if (lvl == 0) rhsImod += T0;
				cx_vec rhsIL = V(indIdx, uvec(1).fill(lvl)).as_col()%Yeind1-
					IL.col(lvl) % Ye1ind1;
				for (int i = 0; i < nInd; i++) {
					vec tempRhsInd(5);
					tempRhsInd << real(rhsM(i)) << imag(rhsM(i)) << rhsImod(i) << real(rhsIL(i)) << imag(rhsIL(i));
					rhsBus.row(i) = (RHS_C_Shr[i] * tempRhsInd).t();
				}
				RHSILr(indIdx) += rhsBus.col(2);
				RHSILi(indIdx) += rhsBus.col(3);
				//DEBUG_PRINT_MAT(RHSILr)
				DEBUG_PRINT_MAT(rhsM)  //TODO: DEBUG this
				DEBUG_PRINT_MAT(rhsImod)  //TODO: DEBUG this
				DEBUG_PRINT_MAT(rhsIL)  //TODO: DEBUG this
				DEBUG_PRINT_MAT(RHSILi)  //TODO: DEBUG this

				//LOOP-Zip
				vec RHSIiLr(nbus, fill::zeros);
				vec RHSIiLi(nbus, fill::zeros);

				vec RHS_BZip = (real(sum(V(zipIdx, seq2R.col(0)) % conj(V(zipIdx, seq2R.col(1))), 1)) - 
					sum(BiL.cols(seq2R.col(0)) % BiL.cols(seq2R.col(1)), 1)) / Bi0 / 2.0;
				cx_vec RHZ_BIConv = sum(IiL.cols(seq2R.col(0)) % BiL.cols(seq2R.col(1)), 1);
				vec RHSIiLr_full = (JI%real(V(zipIdx, uvec(1).fill(lvl)).as_col()) - KI % imag(V(zipIdx, uvec(1).fill(lvl)).as_col()))/Bi0 -
					real(RHZ_BIConv) / Bi0 - Ji0L % RHS_BZip / Bi0;
				vec RHSIiLi_full = (KI%real(V(zipIdx, uvec(1).fill(lvl)).as_col()) + JI % imag(V(zipIdx, uvec(1).fill(lvl)).as_col()))/Bi0 -
					imag(RHZ_BIConv) / Bi0 - Ki0L % RHS_BZip / Bi0;
				RHSIiLr(zipIdx) += RHSIiLr_full;
				RHSIiLi(zipIdx) += RHSIiLi_full;
				DEBUG_PRINT_MAT(RHS_BZip)
				DEBUG_PRINT_MAT(RHSIiLr)
				DEBUG_PRINT_MAT(RHSIiLi)
				DEBUG_PRINT_MAT(RHSIiLr_full)
				DEBUG_PRINT_MAT(RHSIiLi_full)

				//LOOP-Syn
				vec RHSIGr(nbus, fill::zeros);
				vec RHSIGi(nbus, fill::zeros);

				vec AG0(nSyn, fill::zeros);
				vec BG0(nSyn, fill::zeros);
				vec tempCD(nSyn, fill::zeros);
				if (nTaylor >= 2) {
					tempCD = sum(d.cols(seq2R.col(0)) % d.cols(seq2R.col(1)), 1);
					AG0 += cosp.col(2) % tempCD;
					BG0 += sinp.col(2) % tempCD;
				}
				if (nTaylor >= 3) {
					tempCD = sum(d.cols(seq3R.col(0)) % d.cols(seq3R.col(1)) % d.cols(seq3R.col(2)), 1);
					AG0 += cosp.col(3) % tempCD;
					BG0 += sinp.col(3) % tempCD;
				}
				if (nTaylor >= 4) {
					umat seq4 = CheCompUtil::spgetseq(lvl + 1, 4);
					uvec idxSeq4 = any(seq4 == lvl + 1, 1);
					umat seq4R = seq4.rows(find(idxSeq4 == 0));

					tempCD = sum(d.cols(seq4R.col(0)) % d.cols(seq4R.col(1)) % d.cols(seq4R.col(2)) % d.cols(seq4R.col(3)), 1);
					AG0 += cosp.col(4) % tempCD;
					BG0 += sinp.col(4) % tempCD;
				}

				vec CCr = sum(real(V(synIdx, seq2R.col(0))) % Cd.cols(seq2R.col(1)), 1);
				vec DCr = sum(imag(V(synIdx, seq2R.col(0))) % Cd.cols(seq2R.col(1)), 1);
				vec CSr = sum(real(V(synIdx, seq2R.col(0))) % Sd.cols(seq2R.col(1)), 1);
				vec DSr = sum(imag(V(synIdx, seq2R.col(0))) % Sd.cols(seq2R.col(1)), 1);
				vec JCr = sum(JG.cols(seq2R.col(0)) % Cd.cols(seq2R.col(1)), 1);
				vec KCr = sum(KG.cols(seq2R.col(0)) % Cd.cols(seq2R.col(1)), 1);
				vec JSr = sum(JG.cols(seq2R.col(0)) % Sd.cols(seq2R.col(1)), 1);
				vec KSr = sum(KG.cols(seq2R.col(0)) % Sd.cols(seq2R.col(1)), 1);

				vec RHSIG1 = Ef.col(lvl + 1) - (CCr + DSr + Rs % (JCr + KSr) + Xd % (JSr - KCr)) -
					(CG0 + Rs % JG0 - Xd % KG0) % AG0 - (DG0 + Rs % KG0 + Xd % JG0) % BG0;
				vec RHSIG2 = - (CSr - DCr + Rs % (JSr - KCr) - Xq % (JCr + KSr)) -
					(-DG0 - Rs % KG0 - Xq % JG0) % AG0 - (CG0 + Rs % JG0 - Xq % KG0) % BG0;
				vec RHSIG3temp = -Pm.col(lvl + 1) +
					sum(real(V(synIdx, seq2R.col(0)) % cx_mat(JG.cols(seq2R.col(1)), -KG.cols(seq2R.col(1)))), 1) +
					(sum(JG.cols(seq2R.col(0)) % JG.cols(seq2R.col(1)), 1) +
						sum(KG.cols(seq2R.col(0)) % KG.cols(seq2R.col(1)), 1)) % Rs;
				vec RHSIG3 = pShare(idxBalSyn) % RHSIG3temp - RHSIG3temp(idxBalSyn) % pShare;
				RHSIG3(idxBal).fill(0.);
				vec RHSIG = spsolve(MatGB, join_cols(RHSIG1, RHSIG2, RHSIG3));
				vec RHSIGJK = MatGTrans * RHSIG;
				RHSIGr = RHSIGJK(span(0, nbus - 1));
				RHSIGi = RHSIGJK(span(nbus, 2 * nbus - 1));

				DEBUG_PRINT_MAT(P)
				DEBUG_PRINT_MAT(Q)
				DEBUG_PRINT_MAT(W)
				DEBUG_PRINT_MAT(Ysh)
				cx_vec RHS1 = sum(cx_mat(-P.cols(seq2.col(0)), Q.cols(seq2.col(0))+ Qxtra.cols(seq2.col(0))) % conj(W.cols(seq2.col(1))), 1) +
					Ysh %V.col(lvl);
				vec RHS2 = -0.5*real(sum(V.cols(seq2.col(0)) % conj(V.cols(seq2.col(1))), 1));
				cx_vec RHS3 = sum(-W.cols(seq2.col(0)) % V.cols(seq2.col(1)), 1);
				/*DEBUG_PRINT_MAT(AG0)
				DEBUG_PRINT_MAT(BG0)
				DEBUG_PRINT_MAT(CCr)
				DEBUG_PRINT_MAT(DCr)
				DEBUG_PRINT_MAT(CSr)
				DEBUG_PRINT_MAT(DSr)
				DEBUG_PRINT_MAT(JCr)
				DEBUG_PRINT_MAT(KCr)
				DEBUG_PRINT_MAT(JSr)
				DEBUG_PRINT_MAT(KSr)
				DEBUG_PRINT_MAT(RHSIGr)
				DEBUG_PRINT_MAT(RHSIGi)*/
				vec RHS3r = real(RHS3);
				vec RHS3i = imag(RHS3);
				vec RHS3xr = PCQD % RHS3r + PDQC % RHS3i;
				vec RHS3xi = PDQC % RHS3r - PCQD % RHS3i;

				if (lvl == 0) RHS2 += 0.5*VspSq2;

				cx_vec compactRHS1 = RHS1(idxNonSw);
				compactRHS1 += CheCompUtil::sp_submatrix<cx_double>(Y, idxNonSw, isw)*V(isw, uvec(1).fill(lvl + 1));
				/*vec RHS = join_cols(
					join_cols(real(compactRHS1) + RHSILr(idxNonSw) + RHSIiLr(idxNonSw) - RHSIGr(idxNonSw),
						imag(compactRHS1) + RHSILi(idxNonSw) + RHSIiLi(idxNonSw) - RHSIGi(idxNonSw)),
					RHS2(ipv),
					join_cols(real(RHS3(idxNonSw)), imag(RHS3(idxNonSw))));*/
				vec RHS= join_cols(
					join_cols(real(compactRHS1) + RHSILr(idxNonSw) + RHSIiLr(idxNonSw) - RHSIGr(idxNonSw)-RHS3xr(idxNonSw),
						imag(compactRHS1) + RHSILi(idxNonSw) + RHSIiLi(idxNonSw) - RHSIGi(idxNonSw) - RHS3xi(idxNonSw)),
					RHS2(ipv));

				DEBUG_PRINT_MAT(RHS1)
				DEBUG_PRINT_MAT(RHS2)
				DEBUG_PRINT_MAT(RHS3)
				DEBUG_PRINT_MAT(RHS)

				superlu_opts opts;
				opts.allow_ugly = true;
				opts.equilibrate = true;
				opts.refine = superlu_opts::REF_NONE;
				//vec x = spsolve(LHS_mat, RHS, "superlu",opts);
				
				vec cb = RHS;
				const mat& B=cb;
				vec x=RHS;
				// x.zeros(A.n_cols, B.n_cols); 
				// vec xx=spsolve(LHS_mat,RHS);
				// xx.print("xx");

				const bool status_x = sp_auxlib::wrap_to_supermatrix(superX, x);
				const bool status_b = sp_auxlib::wrap_to_supermatrix(superB, B);
				bool use_iter_solver=0;
				if (lvl == 0) {
					if(this->perm_c==NULL){
						options.ColPerm=superlu::COLAMD;
					}else{
						options.ColPerm=superlu::MY_PERMC;
						arrayops::copy(perm_c,this->perm_c,A.n_cols + 1);
					}
					options.Fact=superlu::DOFACT;
					if(this->perm_r==NULL||this->etree==NULL||true){
						options.Fact=superlu::DOFACT;
					}else{
						options.Fact=superlu::SamePattern;
						arrayops::copy(perm_r,this->perm_r,A.n_rows + 1);
						arrayops::copy(etree,this->etree,A.n_cols + 1);
					}
					if (use_iter_solver){
						// arma_wrapper(dgssvx)(&options, &superA, perm_c, perm_r, etree, equed, R, C, &superL, &superU, &work[0], lwork, &superB, &superX, &rpg, &rcond, ferr, berr, &glu, &mu, &stat, &superInfo);
						ilu_set_default_options(&options);
						options.Equil=superlu::NO;
						lwork=0;
						char equedx[1]={'N'};
						arma_wrapper(dgsisx)(&options, &superA, perm_c, perm_r, etree, equedx, R, C, &superL, &superU, NULL, 0, &superB, &superX, &rpg, &rcond, &glu, &mu, &stat, &superInfo);
						// arma_wrapper(dgsitrf)(&options, &superA, 2, 10, etree, &work[0], lwork, perm_c, perm_r, &superL, &superU, &glu, &stat, &superInfo);
						// arma_wrapper(dgstrs)(options.Trans, &superL, &superU, perm_c, perm_r, &superX, &stat, &superInfo);
						
						vec xx=RHS;
						superlu::SuperMatrix superX2;  arrayops::inplace_set(reinterpret_cast<char*>(&superX2), char(0), sizeof(superlu::SuperMatrix));

						const bool status_x2 = sp_auxlib::wrap_to_supermatrix(superX2, xx);
						arma_wrapper(dgstrs)(options.Trans, &superL, &superU, perm_c, perm_r, &superX2, &stat, &superInfo);
						x=xx;
						sp_auxlib::destroy_supermatrix(superX2);

						int max_it=5000;
						double tol=1.0e-10;
						double err=0.0;
						int iter=0;
						int flag=0;
						// x.fill(0.0);
						// LHS_mat.print("A");
						vec cbx=RHS;
						bicgstab(LHS_mat, x, cbx, options.Trans, &superL, &superU, perm_c, perm_r, &stat, &superInfo, max_it, tol, &err, &iter, &flag);
					}else{						
						arma_wrapper(dgssvx)(&options, &superA, perm_c, perm_r, etree, equed, R, C, &superL, &superU, &work[0], lwork, &superB, &superX, &rpg, &rcond, ferr, berr, &glu, &mu, &stat, &superInfo);
					}
					// arma_wrapper(dgssvx)(&options, &superA, perm_c, perm_r, etree, equed, R, C, &superL, &superU, &work[0], lwork, &superB, &superX, &rpg, &rcond, ferr, berr, &glu, &mu, &stat, &superInfo);
					
					if(this->perm_c==NULL){
						this->perm_c=new int[A.n_cols + 1];
					}
					if(this->perm_r==NULL){
						this->perm_r=new int[A.n_rows + 1];
					}
					if(this->etree==NULL){
						this->etree=new int[A.n_cols + 1];
					}
					arrayops::copy(this->perm_c,perm_c,A.n_cols + 1);
					arrayops::copy(this->perm_r,perm_r,A.n_rows + 1);
					arrayops::copy(this->etree,etree,A.n_cols + 1);
				}
				else {
					if (use_iter_solver){
						// x.print("xx");
						// arma_wrapper(dgstrs)(options.Trans, &superL, &superU, perm_c, perm_r, &superX, &stat, &superInfo);
						// x.print("xx");
						ilu_set_default_options(&options);
						options.Equil=superlu::NO;
						// lwork=0;
						// char equedx[1]={'B'};
						// arma_wrapper(dgsisx)(&options, &superA, perm_c, perm_r, etree, equedx, R, C, &superL, &superU, NULL, 0, &superB, &superX, &rpg, &rcond, &glu, &mu, &stat, &superInfo);
						
						int max_it=5000;
						double tol=1.0e-10;
						double err=0.0;
						int iter=0;
						int flag=0;
						
						vec xx=RHS;
						superlu::SuperMatrix superX2;  arrayops::inplace_set(reinterpret_cast<char*>(&superX2), char(0), sizeof(superlu::SuperMatrix));

						const bool status_x2 = sp_auxlib::wrap_to_supermatrix(superX2, xx);
						arma_wrapper(dgstrs)(options.Trans, &superL, &superU, perm_c, perm_r, &superX2, &stat, &superInfo);
						x=xx;
						sp_auxlib::destroy_supermatrix(superX2);

						// x.print("x");
						vec cbx=RHS;
						bicgstab(LHS_mat, x, cbx, options.Trans, &superL, &superU, perm_c, perm_r, &stat, &superInfo, max_it, tol, &err, &iter, &flag);
						// x.print("x");
					}else{						
						arma_wrapper(dgstrs)(options.Trans, &superL, &superU, perm_c, perm_r, &superX, &stat, &superInfo);
					}
					
				}

				// x.print("x");	

				V(idxNonSw, uvec(1).fill(lvl + 1)) = 
					cx_vec(x(span(0, npq + npv - 1)), x(span(npq + npv, 2*(npq + npv) - 1)));
				Q(ipv, uvec(1).fill(lvl + 1)) = x.tail(npv);

				vec Cx = real(V.col(lvl+1));
				vec Dx = imag(V.col(lvl + 1));
				vec RHS3xxr = RHS3r - E0 % Cx + F0 % Dx;
				vec RHS3xxi = RHS3i - F0 % Cx - E0 % Dx;
				vec Ex = C0i % RHS3xxr + D0i % RHS3xxi;
				vec Fx = -D0i % RHS3xxr + C0i % RHS3xxi;
				W(idxNonSw, uvec(1).fill(lvl + 1)) =
					cx_vec(Ex(idxNonSw), Fx(idxNonSw));
				DEBUG_PRINT_MAT(V)
				DEBUG_PRINT_MAT(W)
				DEBUG_PRINT_MAT(Q)

				//Aux Ind
				for (int i = 0; i < nInd; i++) {
					vec tempvi(2);
					tempvi << real(V(indIdx(i), lvl + 1)) << imag(V(indIdx(i), lvl + 1));
					vec tempx = LHS_MatInd_Full[i] * tempvi + rhsBus.row(i).t();
					IL(i, lvl + 1) = cx_double(tempx(2), tempx(3));
					IR(i, lvl + 1) = cx_double(tempx(0), tempx(1));
					s(i, lvl + 1) = tempx(4);
					Vm(i, lvl + 1) = V(indIdx(i), lvl + 1) - IL(i, lvl + 1)*Z1(i);
				}

				//Aux Zip
				IiL.col(lvl + 1) = cx_vec(LHS_MatZip.col(0), LHS_MatZip.col(2)) % real(V(zipIdx, uvec(1).fill(lvl + 1))) +
					cx_vec(LHS_MatZip.col(1), LHS_MatZip.col(3)) % imag(V(zipIdx, uvec(1).fill(lvl + 1))) +
					cx_vec(RHSIiLr_full, RHSIiLi_full);
				BiL.col(lvl + 1) = Mat_BZip.col(0) % real(V(zipIdx, uvec(1).fill(lvl + 1))) +
					Mat_BZip.col(1) % imag(V(zipIdx, uvec(1).fill(lvl + 1))) +
					RHS_BZip;

				//Aux Syn
				vec IGJKd = MatGBiA * join_cols(real(V.col(lvl + 1)), imag(V.col(lvl + 1))) + RHSIG;
				if (nSyn > 0) {
					JG.col(lvl + 1) = IGJKd(span(0, nSyn-1));
					KG.col(lvl + 1) = IGJKd(span(nSyn, 2*nSyn - 1));
					d.col(lvl + 1) = IGJKd(span(2 * nSyn, 3 * nSyn - 1));
					Cd.col(lvl + 1) = A1n % d.col(lvl + 1) + AG0;
					Sd.col(lvl + 1) = B1n % d.col(lvl + 1) + BG0;
				}
				DEBUG_PRINT_MAT(IL)
				DEBUG_PRINT_MAT(IR)
				DEBUG_PRINT_MAT(s)
				DEBUG_PRINT_MAT(Vm)
				DEBUG_PRINT_MAT(IiL)
				DEBUG_PRINT_MAT(BiL)
				DEBUG_PRINT_MAT(JG)
				DEBUG_PRINT_MAT(KG)
				DEBUG_PRINT_MAT(d)
				DEBUG_PRINT_MAT(Cd)
				DEBUG_PRINT_MAT(Sd)
				;
			}

			CheSolution* psol = new CheSolutionPade(stateIdx.nState, nlvl + 1);
			psol->solution.rows(stateIdx.vrIdx) = real(V);
			psol->solution.rows(stateIdx.viIdx) = imag(V);
			psol->solution.rows(stateIdx.qIdx) = Q;
			psol->solution.rows(stateIdx.sIdx) = s;
			psol->solution.rows(stateIdx.mDeltaIdx) = d;
			psol->solution.rows(stateIdx.mPgIdx) = Pm;
			psol->solution.rows(stateIdx.mEfIdx) = Ef;			
			
			superlu::free_stat(&stat);

			superlu::free(berr);
			superlu::free(ferr);
			superlu::free(C);
			superlu::free(R);
			superlu::free(etree);
			superlu::free(perm_c);
			superlu::free(perm_r);

			sp_auxlib::destroy_supermatrix(superU);
			sp_auxlib::destroy_supermatrix(superL);
			sp_auxlib::destroy_supermatrix(superA);
			sp_auxlib::destroy_supermatrix(superB);
			sp_auxlib::destroy_supermatrix(superX);

			return psol;
		}

		CheSingleEmbedSystem* ChePfCalculator::getNewStage() {
			return NULL;
		}

		int ChePfCalculator::calc() {
			double alpha = 0.0;
			double alphax = 0.0;
			double alphaConfirm = 0.0;

			double alphaPreMult = 0.95;
			double diffMul = 1.5;
			double alphaTol = this->compOpt.alphaTol;
			double diffTol = this->compOpt.diffTol;
			double diffTolMax = this->compOpt.diffTolMax;
			double segment = this->compOpt.segLen;
			int maxCount = 10;
			int maxNoMove = 5;
			int noMove = 0;
			CheSolution* pSol = NULL;

			CheSingleEmbedSystem* pInitSystem = this->getInitSystem(baseSys);
			this->cheList.push_back(pInitSystem);
			CheSingleEmbedSystem* lastEmbeddedSys = NULL;

			while (alphaConfirm < 1 - alphaTol / 1000.0) {
				alpha = 1 - alphaConfirm;
				if (alpha>segment){
					alpha=segment;
				}
				alphax = alpha / alphaPreMult;
				CheSingleEmbedSystem* pCurrEmbeddedSys = this->cheList.back();
				if (lastEmbeddedSys != pCurrEmbeddedSys || pSol==NULL) {
					pSol = this->getCheSolution();
				}

				/*CheSolutionPade* pSolX = (CheSolutionPade*)pSol;
				vec solVal = pSol->getSolValue(0);

				cout.precision(16);
				cout.setf(ios::scientific);

				pSolX->solution.raw_print(cout,"Sol");
				pSolX->numerator.raw_print(cout, "Num");
				pSolX->denomenator.raw_print(cout, "Den");*/

				double alphaLeft = 0;
				double alphaRight = alphax;
				double absDiff = 0.0;
				while (true) {
					vec diff = pCurrEmbeddedSys->calcEqBalance(pSol, alphax);
					absDiff = max(abs(diff));
					if (absDiff < diffTol && abs(alphax - alpha / alphaPreMult) < alphaTol / 1000.0) {
						alphaLeft = alphaRight;
					}
					else {
						if (absDiff < diffTol) {
							alphaLeft = alphax;
							alphax = (alphaLeft + alphaRight) / 2.0;
						}
						else {
							alphaRight = alphax;
							alphax = (alphaLeft + alphaRight) / 2.0;
						}
					}
					double alphaTolTemp=0.03*alphaRight;
					if (alphaTolTemp<alphaTol){
						alphaTolTemp=alphaTol;
					}

					if (abs(alphaRight - alphaLeft) < alphaTolTemp){
						break;
					}
				}

				double diffParadigm = absDiff + 1e-12;
				int countCheckAlpha = 1;
				alphax = alphaLeft * alphaPreMult;
				if (alpha - alphax < alphaTol / 10) {
					alphax = alpha;
				}

				while (countCheckAlpha < maxCount){
					vec diff = pCurrEmbeddedSys->calcEqBalance(pSol, alphax);
					absDiff = max(abs(diff));
					if (absDiff < 1.05*diffParadigm) {
						break;
					}
					alphax *= alphaPreMult;
					countCheckAlpha++;
				}
				if (countCheckAlpha >= maxCount) {
					alphax *= alphaPreMult;
				}
				alpha = alphax;
				alphaConfirm += alpha;
				cout << "Step=" <<alphaConfirm<<", added="<<alpha<<", (maxDiff<"<<diffTol<<")."<< endl;

				if (alpha == 0.0) {
					cout << "Step did not move!" << endl;
					noMove++;
					if (noMove >= maxNoMove) {
						cout << "Reached consecutive max no move, exit!" << endl;
						if(pSol!=NULL){
							delete pSol;
							pSol=NULL;
						}
						break;
					}
					if (diffTol >= diffTolMax) {
						cout << "Max DiffTol reached and not move, exit!" << endl;
						if(pSol!=NULL){
							delete pSol;
							pSol=NULL;
						}
						break;
					}
					diffTol *= diffMul;
					if (absDiff > diffTol) {
						diffTol = absDiff;
					}
					cout << "Enlarge tol! (Tol=" <<diffTol<<")."<< endl;
					if(pSol!=NULL){
						delete pSol;
						pSol=NULL;
					}
				}
				else {
					noMove = 0;
					CheState curState(pCurrEmbeddedSys->baseSys, pSol->getSolValue(alpha));
					CheSingleEmbedSystem* nextSystem = pCurrEmbeddedSys->getNewEmbeddedSystem(curState, alpha);
					this->cheList.push_back(nextSystem);
					this->solList.push_back(pSol);
				}

			}

			if (alphaConfirm >= 1 - alphaTol / 1000.0) {
				this->reachesMaxAlpha = true;
				return 0;
			}
			else {
				this->reachesMaxAlpha = false;
				return -1;
			}
		}

		CheState ChePfCalculator::exportResult(){
			return cheList.back()->initState;
		}

		
		void ChePfCalculator::writeMatFile(const char *fileName,double interval){
			this->writeMatFile(fileName);
		}

		void ChePfCalculator::writeMatFile(const char *fileName){
			vec result=cheList.back()->initState.state;

			mat solutionMat(result.n_rows, 1, fill::zeros);
			solutionMat.col(0) = result;

			mat_t *matfp;
			matvar_t *matvar;
			size_t dims[2] = {solutionMat.n_rows, solutionMat.n_cols};
			matfp = Mat_CreateVer(fileName, NULL, MAT_FT_DEFAULT);
			if (NULL == matfp)
			{
				cerr << "Error creating MAT file \"" << fileName << "\"." << endl;
				return;
			}

			matvar = Mat_VarCreate("s", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, solutionMat.memptr(), 0);
			if (NULL == matvar)
			{
				cerr << "Error creating variable for ’s’." << endl;
			}
			else
			{
				Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
				Mat_VarFree(matvar);
			}

			Mat_Close(matfp);
		}

		ChePfCalculator::~ChePfCalculator(){
			if(perm_c!=NULL){
				delete [] perm_c;
				perm_c=NULL;
			}
			if(perm_r!=NULL){
				delete [] perm_r;
				perm_r=NULL;
			}
			if(etree!=NULL){
				delete [] etree;
				etree=NULL;
			}
			
			// if(perm_ci!=NULL){
			// 	delete [] perm_ci;
			// 	perm_ci=NULL;
			// }
			// if(perm_ri!=NULL){
			// 	delete [] perm_ri;
			// 	perm_ri=NULL;
			// }
			// if(etreei!=NULL){
			// 	delete [] etreei;
			// 	etreei=NULL;
			// }
		}
	}
}