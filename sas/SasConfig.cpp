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
#include "SasConfig.h"

namespace che {
	namespace core {
		unsigned int SasCompUtil::ELFHashStr(const char *str, int len) {
			unsigned int h = 0;
			unsigned int x = 0;
			int i;

			for (i = 0; i < len; ++i)
			{
				h = (h << 4) + *str++;
				if ((x = h & 0xF0000000) != 0)
				{
					h ^= x >> 24;
					h &= ~x;
				}
			}
			return h;
		}

		unsigned int SasCompUtil::ELFHashInt(int val) {
			unsigned int h1 = val%PRIME1;
			unsigned int h2 = (PRIME2 - val % PRIME2)*PRIME3;
			return h1 + h2;
		}

		unsigned int SasCompUtil::ELFHashDouble(double val) {
			int segs = sizeof(double) / sizeof(int);
			int *pVal = (int*)&val;
			int tempVal = pVal[0];
			unsigned int h = ELFHashInt(tempVal);
			for (int i = 0; i < segs; i++) {
				h = ELFHashOfTwo(ELFHashInt(pVal[1]), h);
			}
			return h;
		}

		unsigned int SasCompUtil::ELFHashOfTwo(unsigned int h1, unsigned int h2) {
			return ((PRIME1 - h1 % PRIME1)*PRIME3 + h2 % PRIME2)%LARGE_PRIME;
		}
	}
}