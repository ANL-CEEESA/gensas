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
#include "CheIOUtil.h"


namespace che {
	namespace io {
		int CheSimpleMDataReader::parseMatFile(const char* filePath) {
			ifstream mfile(filePath, ios::in);
			size_t found;
			if (!mfile) { //if file is null, return;
				cout << "can't open file: " << filePath << endl;
				return CHE_IO_FAIL;
			}

			string templine;
			int lineCount = 0;
			list<vector<double>> s;
			
			while (true) {
				do {
					getline(mfile, templine);
					lineCount++;
					trimStr(templine);
					removeComment(templine);
				} while (templine.empty()&&mfile.good());
				if (!mfile.good()) {
					if (mfile.eof()) {
						addMsg("Read to the end of file.", lineCount);
						mfile.close();
						return CHE_IO_SUCCESS;
					}
					else {
						addMsg("Unknown failure", lineCount);
						mfile.close();
						return CHE_IO_FAIL;
					}
				}

				string tempAppend;
				while ((found=templine.find("...")) != string::npos) {
					getline(mfile, tempAppend);
					lineCount++;
					if (!mfile.good()) {
						if (mfile.eof()) {
							addMsg("Read to the end of file.", lineCount);
							mfile.close();
						}
						else {
							addMsg("Unknown failure", lineCount);
							mfile.close();
							return CHE_IO_FAIL;
						}
					}
					trimStr(tempAppend);
					removeComment(tempAppend);
					if (!tempAppend.empty()) {
						templine.erase(found);
						templine += tempAppend;
					}
				}

				found = templine.find("=");
				if (found == string::npos) {
					addMsg("Expect '='. Parsing interrupted.", lineCount);
					mfile.close();
					return CHE_IO_SUCCESS;					
				}
				string varName = templine.substr(0, found);
				trimStr(varName);
				string rest = templine.substr(found + 1);
				trimStr(rest);
				if (rest[0] == '[') {// Matrix.
					rest.erase(0,1);
					vector<mat> valMat;
					vector<double> valLine;
					string sx;
					bool endOfMatrix = false;

					while (true) {
						trimStr(rest);
						found = rest.find_first_of("]");
						if (found != string::npos) {
							rest.erase(found+1);
							trimStr(rest);
						}
						if (rest.find_first_not_of(" ,;\t]") == string::npos&&rest.find("]") != string::npos) {
							endOfMatrix = true;
							break;
						}
						while (true) {
							found = rest.find_first_of(" ,;\t]");
							try {
								sx = rest.substr(0, found);
								trimStr(sx);
								double x = stod(sx);
								valLine.push_back(x);
							}
							catch (const logic_error& e) {
								addMsg("Problem when parsing '" + sx + "' to double. Interrupt.", lineCount);
								mfile.close();
								return CHE_IO_SUCCESS;
							}
							if (found == string::npos || rest[found] == ';' || rest[found] == ']') {
								if (!valLine.empty()) {
									mat row = conv_to<mat>::from(valLine).t();
									valMat.push_back(row);
									valLine.clear();
								}
								if (found != string::npos&&rest[found] == ']') {
									endOfMatrix = true;
									break;
								}
							}
							if (found != string::npos) {
								rest.erase(0, found + 1);
							}
							else {
								rest.erase(0, found);
							}
							if (rest.find_first_not_of(" ,;\t]") == string::npos&&rest.find("]") != string::npos) {
								if (!valLine.empty()) {
									mat row = conv_to<mat>::from(valLine).t();
									valMat.push_back(row);
									valLine.clear();
								}
								endOfMatrix = true;
								break;
							}
							rest.erase(0, rest.find_first_not_of(" ,;\n\r\t"));
							trimStr(rest);
							if (rest.empty()) {
								if (!valLine.empty()) {
									mat row = conv_to<mat>::from(valLine).t();
									valMat.push_back(row);
									valLine.clear();
								}
								break;
							}
							if (!mfile.good()) {
								endOfMatrix = true;
								break;
							}
						}
						if (endOfMatrix) {
							break;
						}
						//read new line
						do {
							getline(mfile, rest);
							lineCount++;
							trimStr(rest);
							removeComment(rest);
						} while (rest.empty() && mfile.good());
						if (!mfile.good()) {
							if (mfile.eof()) {
								addMsg("Read to the end of file.", lineCount);
								mfile.close();
							}
							else {
								addMsg("Unknown failure", lineCount);
								mfile.close();
								return CHE_IO_FAIL;
							}
						}
						while ((found = rest.find("...")) != string::npos) {
							getline(mfile, tempAppend);
							lineCount++;
							if (!mfile.good()) {
								if (mfile.eof()) {
									addMsg("Read to the end of file.", lineCount);
									mfile.close();
								}
								else {
									addMsg("Unknown failure", lineCount);
									mfile.close();
									return CHE_IO_FAIL;
								}
							}
							trimStr(tempAppend);
							removeComment(tempAppend);
							if (!tempAppend.empty()) {
								rest.erase(found);
								rest += tempAppend;
							}
						}
					}
					try {
						if (valMat.empty()) {
							matList.insert(pair<string, mat>(varName, mat()));
							valMat.clear();
						}
						else {
							mat mx = valMat[0];
							for (int i = 1; i < valMat.size(); i++) {
								mx = join_cols(mx, valMat[i]);
							}
							matList.insert(pair<string, mat>(varName, mx));
							valMat.clear();
						}
					}
					catch (const logic_error& e) {
						addMsg("Problem when joining matrices. Interrupt.", lineCount);
						mfile.close();
						return CHE_IO_SUCCESS;
					}
					if (mfile.eof()) {
						addMsg("Read to the end of file.", lineCount);
						mfile.close();
						return CHE_IO_SUCCESS;
					}
				}
				else {
					found = rest.find_first_of(" ,;\t");
					string sx;
					try {
						sx = rest.substr(0, found);
						double x = stod(sx);
						mat mx(1,1);
						mx << x << endr;
						matList.insert(pair<string, mat>(varName, mx));
					}
					catch(const logic_error& e){
						addMsg("Problem when parsing '"+ sx+"' to double. Interrupt.", lineCount);
						mfile.close();
						return CHE_IO_SUCCESS;
					}
					rest.erase(0, found+1);
					found = rest.find_first_not_of(" ,;\n\r\t");
					if (found != string::npos) {
						addMsg("Ignoring the content after the first value.", lineCount);
					}
				}
			}
			mfile.close();
			return CHE_IO_SUCCESS;
		}
	}
}