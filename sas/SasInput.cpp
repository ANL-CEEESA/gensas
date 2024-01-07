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
#include "sas/SasInput.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

namespace che {
	namespace io {

		BufferedFileReader::BufferedFileReader() {
			buffer = NULL;
			cursur = NULL;
			lineHead = NULL;
			line = 0;
			size = 0;
			filename = NULL;
		}


		bool BufferedFileReader::readFile(const char* filename) {
			resetReader();
			ifstream ifs(filename);
			int nameLength = strlen(filename);
			this->filename = new char[nameLength +1];
			strcpy(this->filename, filename);
			if (ifs) {
				ifs.seekg(0, ifs.end);
				int length = ifs.tellg();
				ifs.seekg(0, ifs.beg);
				buffer = new uchar[length+1];
				ifs.read((char*)buffer, length);
				if (ifs){
					cout << "All " << ifs.gcount() << " characters read successfully."<<endl;
				}else{
					cout << "Only " << ifs.gcount() << " characters could be read"<<endl;
					length = ifs.gcount();
				}
				buffer[length] = char_traits<uchar>::eof();
				size = length;
				cursur = buffer;
				lineHead = buffer;
				return true;
			}
			else {
				return false;
			}
		}

		void BufferedFileReader::resetReader() {
			cursur = NULL;
			lineHead = NULL;
			line = 0;
			size = 0;
			if (buffer != NULL) {
				delete[] buffer;
				buffer = NULL;
			}
			if (filename != NULL) {
				delete[] filename;
				filename = NULL;
			}
		}

		BufferedFileReader::~BufferedFileReader() {
			resetReader();
		}
	}
}