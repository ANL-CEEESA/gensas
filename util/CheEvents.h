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
#ifndef _Che_Events_H_
#define _Che_Events_H_

#include "SafeArmadillo.h"

using namespace arma;

namespace che {
	namespace util {
		class Fault {
		public:
			enum FaultTypeBit {FAULT_A=0x1,FAULT_B=0x2,FAULT_C=0x4,FAULT_N=0x8,FAULT_G=0x16};
			const static int NO_FAULT = 0;
			const static int FAULT_3P = FAULT_A | FAULT_B | FAULT_C ;
			const static int FAULT_3PG= FAULT_A | FAULT_B | FAULT_C | FAULT_G;

			int lineIdx;
			int pos;
			int fType;
			cx_mat faultY;
			double startT;
			double endT;

			Fault(int lineIdx=0,double pos=0.0,int fType=NO_FAULT,const cx_mat& faultY=cx_mat(),double startT=0.0,double endT=0.0) {
				this->lineIdx = lineIdx;
				this->pos = pos;
				this->fType = fType;
				this->faultY = faultY;
				this->startT = startT;
				this->endT = endT;
			}

			Fault(const Fault& fault) {
				this->lineIdx = fault.lineIdx;
				this->pos = fault.pos;
				this->fType = fault.fType;
				this->faultY = fault.faultY;
				this->startT = fault.startT;
				this->endT = fault.endT;
			}
		};
	}
}

#endif