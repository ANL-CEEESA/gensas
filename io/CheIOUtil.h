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
#ifndef _Che_CHEIOUTIL_H_
#define _Che_CHEIOUTIL_H_

#include <list>
#include <vector>
#include <map>
#include <string>
#include "SafeArmadillo.h"
#include "che_io_defs.h"

using namespace std;
using namespace arma;

namespace che {
	namespace io {
		class CheSimpleMDataReader {
		public:
			map<string, mat> matList;
			list<string> msgs;

			CheSimpleMDataReader() {
				reset();
			}

			void reset() {
				matList.clear();
				msgs.clear();
			}

			int parseMatFile(const char* filePath);

		private:
			int trimStrLeft(string &str) {
				int len = str.length();
				size_t left = str.find_first_not_of(" \n\r\t");
				str.erase(0, left);
				if (left == string::npos) {
					return len;
				}
				else {
					return left;
				}
			}

			void trimStrRight(string &str) {
				size_t right = str.find_last_not_of(" \n\r\t");
				if (right != string::npos) {
					str.erase(right + 1);
				}
			}

			int trimStr(string& str) {
				int offset = trimStrLeft(str);
				trimStrRight(str);
				return offset;
			}

			void removeComment(string& str) {
				size_t found = str.find("%");
				if (found != string::npos) {
					str.erase(found);
				}				
				found = str.find("...");
				if (found != string::npos) {
					str.erase(str.find("...")+3);
				}
			}

			void addMsg(const string& s) {
				msgs.push_back(s);
			}

			void addMsg(const string& s, int line) {
				string head = "[Line";
				msgs.push_back(head+to_string(line)+"]"+s);
			}
		};
	}
}

#endif