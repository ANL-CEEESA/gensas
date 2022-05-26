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
#pragma once

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include "MDataFormat.h"

using namespace std;

//#define DEFAULT_READER_BUFFER_SIZE 50
#define MFILE_READER_SUCCESS 0
#define MFILE_READER_FAIL -1

//#define PARSE_START 10
//#define PARSE_FUNC 11
//#define PARSE_EXPR 12
//#define PARSE_NUM_LINE 13
//#define PARSE_OBJ_LINE 14
//#define PARSE_END 15
//#define PARSE_EOF 16
//#define PARSE_SYNTAX_ERR 99

namespace che {
	namespace io {
		namespace mdata {

			const char* tokenStrings[] = {
			#define TOKEN(k,s) s,
			#include "MatlabToken.h"
			#undef TOKEN
			};

			enum token {
			#define TOKEN(k,s) k,
			#include "MatlabToken.h"
			#undef TOKEN
			TK_SYMBOL,
			TK_NUM_INT,
			TK_NUM_FLOAT,
			TK_STR,
			TK_UNDETERMINED,
			};

			inline int trimStrLeft(string &str) {
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

			inline void trimStrRight(string &str) {
				str.erase(str.find_last_not_of(" \n\r\t") + 1);
			}

			inline int trimStr(string& str) {
				int offset=trimStrLeft(str);
				trimStrRight(str);
				return offset;
			}

			inline void removeComment(string& str) {
				str.erase(str.find("%"));
				str.erase(str.find("...")+3);
			}

			int checkNumeric(string s);

			int checkSymbol(string s);

			int checkStr(string s);

			int getTokenType(string str) {
				trimStr(str);
				if (str.empty()) return TK_EMPTY;
#define TOKEN(k,x) if(str.length()==strlen(x)&&str.compare(x)==0) return k;
#include "MatlabToken.h"
#undef TOKEN
				int type = TK_UNDETERMINED;
				type = checkNumeric(str);
				if (type != TK_UNDETERMINED) return type;
				type = checkSymbol(str);
				if (type != TK_UNDETERMINED) return type;
				type = checkStr(str);
				if (type != TK_UNDETERMINED) return type;
				return TK_UNDETERMINED;
			}

			class SyntaxNode {
			public:
				string token;
				int tokenType;
				list<SyntaxNode*> pChildren;
				SyntaxNode* pParent;
				// Form parsing tree
				enum state {
					s_start, s_func, s_func_output, s_func_decl, s_func_name, s_func_input, s_func_end, 
					s_output_list, s_input_list, s_exp_end,
					s_exp_assn, s_exp_assn_name, s_exp_assn_mat, s_exp_assn_cell, s_exp_assn_line
				};
				state st;

				SyntaxNode() {
					token = "";
					tokenType = TK_UNDETERMINED;
					pChildren = list<SyntaxNode*>();
					pParent = NULL;
					st = s_start;
				}

				SyntaxNode(string& token, int tokenType, SyntaxNode *pParent,state st) {
					this->token = token;
					this->tokenType = tokenType;
					this->pParent = pParent;
					this->st = st;
					this->pChildren = list<SyntaxNode*>();
				}
			};

			class BufferedMFileReader {
			public:
				string buffer;
				string currToken;
				int line;
				int cursur;
				SyntaxNode sTree;
				SyntaxNode *pNode;
				MDataSet mData;
				ifstream ifs;
				char* filePath;
				
				BufferedMFileReader(char* filePath) {
					resetReader(filePath);
				}

				void resetReader(char* filePath) {
					if (ifs&&ifs.is_open()) {
						ifs.close();
					}
					ifs = ifstream(filePath, ios::in);
					this->filePath = filePath;
					buffer = "";
					currToken = "";
					line = 0;
					cursur = 0;
					sTree = SyntaxNode();
					pNode = &sTree;
					mData = MDataSet();
				}

				int parseFile();

				string nextToken(vector<string> &delims) {
					cursur+=trimStrLeft(buffer);
					if (buffer.empty()) {
						getNewLine();
						if (!ifs.good()) {
							return "";
						}
						else {
							return "\n";
						}
					}
					size_t pos = string::npos;
					for (int i = 0;i< delims.size(); i++) {
						pos = buffer.find(delims[i].c_str(), 0, pos);
						if (pos == 0) pos = delims[i].length(); //This may match a delimiter
					}
					string token = buffer.substr(0, pos);
					cursur += token.length();
					buffer.erase(0, pos);
					return token;
				}

				inline int getNewLine() {
					string s = "";
					int rdState=0;
					do {
						rdState= getline(ifs, s).rdstate();
						line++;
						cursur = 0;
						trimStr(s);
						removeComment(s);
					} while (s.empty()&&ifs.good());
					buffer.append(s);
					return rdState;
				}
			};
		}
	}
}