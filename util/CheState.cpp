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
#include "CheState.h"
#include "CheCompUtil.h"

namespace che {
	namespace util {

		CheStateIdx::CheStateIdx() {
			nState = 0;
		}

		CheStateIdx::CheStateIdx(const CheStateIdx& s) {
			this->nState = s.nState;
			this->vrIdx = s.vrIdx;
			this->viIdx = s.viIdx;
			this->qIdx = s.qIdx;
			this->pIdx = s.pIdx;
			this->sIdx = s.sIdx;
			this->indEr1Idx = s.indEr1Idx;
			this->indEm1Idx = s.indEm1Idx;
			this->indEr2Idx = s.indEr2Idx;
			this->indEm2Idx = s.indEm2Idx;
			this->mDeltaIdx = s.mDeltaIdx;
			this->mOmegaIdx = s.mOmegaIdx;
			this->mEq1Idx = s.mEq1Idx;
			this->mEq2Idx = s.mEq2Idx;
			this->mEd1Idx = s.mEd1Idx;
			this->mEd2Idx = s.mEd2Idx;
			this->mPsidIdx = s.mPsidIdx;
			this->mPsiqIdx = s.mPsiqIdx;
			this->mPgIdx = s.mPgIdx;
			this->mEfIdx = s.mEfIdx;
			this->avr1mIdx = s.avr1mIdx;
			this->avr1r1Idx = s.avr1r1Idx;
			this->avr1r2Idx = s.avr1r2Idx;
			this->avr1rIdx = s.avr1rIdx;
			this->avr1fIdx = s.avr1fIdx;
			this->avr2mIdx = s.avr2mIdx;
			this->avr2r1Idx = s.avr2r1Idx;
			this->avr2rIdx = s.avr2rIdx;
			this->avr2r2Idx = s.avr2r2Idx;
			this->avr2fIdx = s.avr2fIdx;
			this->avrmIdx = s.avrmIdx;
			this->avrrIdx = s.avrrIdx;
			this->avrfIdx = s.avrfIdx;
			this->avrrefIdx = s.avrrefIdx;
			this->tg1inIndx = s.tg1inIndx;
			this->tg11Idx = s.tg11Idx;
			this->tg12Idx = s.tg12Idx;
			this->tg13Idx = s.tg13Idx;
			this->tg2gIdx = s.tg2gIdx;
			this->tg2mIdx = s.tg2mIdx;
			this->tmechIdx = s.tmechIdx;
			this->fIdx = s.fIdx;
			this->qpltIdx = s.qpltIdx;
			this->vgIdx = s.vgIdx;
		}

		CheState::CheState() {
			state = vec();
		}

		CheState::CheState(int nState) {
			state = vec(nState, fill::zeros);
		}

		CheState::CheState(const CheState & st) {
			this->stateIdx = CheStateIdx(st.stateIdx);
			this->state = st.state;
		}

		CheState::CheState(const chedata::PsatDataSet &psat) {
			stateIdx = CheCompUtil::getCheStateIdx(psat);
			state = vec(stateIdx.nState, fill::zeros);
		}

		CheState::CheState(const chedata::PsatDataSet &psat, const vec & state) {
			stateIdx = CheCompUtil::getCheStateIdx(psat);
			this->state = state;
		}

		vec CheState::getSubVec(uvec idx) const {
			return state(idx);
		}

		void CheState::setSubVec(uvec idx, vec value) {
			state(idx) = value;
		}
	}
}