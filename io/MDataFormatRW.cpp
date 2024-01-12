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
#include "io/MDataFormatRW.h"

using namespace std;

namespace che {
	namespace io {
		namespace mdata {

			int checkNumeric(string s) {
				trimStr(s);
				int idx = 0;
				enum state { s_start, s_num_sign, s_num_int, s_num_float, s_num_dot, s_exp, s_exp_sign, s_exp_body };
				state st = s_start;
				char c = 0;
				while (idx < s.length()) {
					c = s[idx];
					if (st == s_start) {
						if (c == '+' || c == '-') {
							st = s_num_sign;
						}
						else if (c >= '0'&&c <= '9') {
							st = s_num_int;
						}
						else if (c == '.') {
							st = s_num_dot;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_num_sign) {
						if (c == '+' || c == '-') {

						}
						else if (c >= '0'&&c <= '9') {
							st = s_num_int;
						}
						else if (c == '.') {
							st = s_num_dot;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_num_int) {
						if (c == '+' || c == '-') {
							return TK_UNDETERMINED;
						}
						else if (c >= '0'&&c <= '9') {

						}
						else if (c == '.') {
							st = s_num_float;
						}
						else if (c == 'e' || c == 'E') {
							st = s_exp;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_num_dot) {
						if (c == '+' || c == '-') {
							return TK_UNDETERMINED;
						}
						else if (c >= '0'&&c <= '9') {
							st = s_num_float;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_num_float) {
						if (c == '+' || c == '-') {
							return TK_UNDETERMINED;
						}
						else if (c >= '0'&&c <= '9') {

						}
						else if (c == '.') {
							return TK_UNDETERMINED;
						}
						else if (c == 'e' || c == 'E') {
							return s_exp;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_exp) {
						if (c == '+' || c == '-') {
							st = s_exp_sign;
						}
						else if (c >= '0'&&c <= '9') {
							st = s_exp_body;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_exp_sign) {
						if (c >= '0'&&c <= '9') {
							st = s_exp_body;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_exp_body) {
						if (c >= '0'&&c <= '9') {

						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else {
						return TK_UNDETERMINED;
					}
					idx++;
				}
				if (st == s_num_int) {
					return TK_NUM_INT;
				}
				else if (st == s_num_float || st == s_exp_body) {
					return TK_NUM_FLOAT;
				}
				else {
					return TK_UNDETERMINED;
				}
			}

			int checkSymbol(string s) {
				trimStr(s);
				int idx = 0;
				char c = 0;
				enum state { s_start, s_symbol, s_dot };
				state st = s_start;
				while (idx < s.length()) {
					c = s[idx];
					if (st == s_start) {
						if ((c >= 'A'&&c <= 'Z') || (c >= 'a'&&c <= 'z')) {
							st = s_symbol;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_symbol) {
						if ((c >= 'A'&&c <= 'Z') || (c >= 'a'&&c <= 'z') || (c >= '0'&&c <= '9') || c == '_') {
							st = s_symbol;
						}
						else if (c == '.') {
							st = s_dot;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else if (st == s_dot) {
						if ((c >= 'A'&&c <= 'Z') || (c >= 'a'&&c <= 'z')) {
							st = s_symbol;
						}
						else {
							return TK_UNDETERMINED;
						}
					}
					else {
						return TK_UNDETERMINED;
					}
					idx++;
				}
				if (st == s_symbol) {
					return TK_SYMBOL;
				}
				else {
					return TK_UNDETERMINED;
				}
			}

			int checkStr(string s) {
				trimStr(s);
				if (s[0] != '\'') return TK_UNDETERMINED;
				if (s.find('\'', 1) != s.length() - 1) return TK_UNDETERMINED;
				return TK_STR;
			}

			void message(int line, int cursur, const char* msg, const char* level) {
				cout << "[" << level << "]" << msg << " At line " << line << ", column " << cursur << endl;
			}

			void errMsg(int line, int cursur, const char* msg) {
				message(line, cursur, msg, "ERROR");
			}

#define MACRO_READ_SKIP_TRIDOTS currToken=nextToken(delims);while((tokenType=getTokenType(currToken))==TK_TRIDOT){getNewLine();currToken=nextToken(delims);}

			int BufferedMFileReader::parseFile() {
				if (!ifs) { //if file is null, return;
					cout << "can't open file: " << filePath << endl;
					return MFILE_READER_FAIL;
				}

				sTree = SyntaxNode();
				pNode = &sTree;
				mData = MDataSet();

				while (ifs) {
					int cursurRecord = cursur;
					if (pNode->st == SyntaxNode::s_start) {
						if (pNode->tokenType == TK_EMPTY) {
							vector<string> delims = { " ","=","..." };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_KW_function) {
								string emptyStr = string("");
								SyntaxNode sNode = SyntaxNode(emptyStr, TK_EMPTY, pNode, SyntaxNode::s_func);
								pNode->pChildren.push_back(&sNode);
								pNode = &sNode;

								SyntaxNode sNodef = SyntaxNode(currToken, tokenType, pNode, SyntaxNode::s_func_decl);
								pNode->pChildren.push_back(&sNodef);
								pNode = &sNodef;
							}
							else if (tokenType == TK_SYMBOL) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode, SyntaxNode::s_exp_assn);
								pNode->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect <function> and variables.");
								return MFILE_READER_FAIL;
							}
						}
					}
					else if (pNode->st == SyntaxNode::s_func_decl) {
						if (pNode->tokenType == TK_KW_function) {
							vector<string> delims = { " ","...","[","=","(",";" };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL || tokenType == TK_SQBRKT_L) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode, SyntaxNode::s_func_output);
								pNode->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables or list of variables surrounded by [].");
								return MFILE_READER_FAIL;
							}
						}
					}
					else if (pNode->st == SyntaxNode::s_func_output) {
						if (pNode->tokenType == TK_SYMBOL) {
							vector<string> delims = { " ","...","=","(",";",",", };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_EQUAL) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_name);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else if (tokenType == TK_SEMCOLON || tokenType == TK_COMMA || tokenType == TK_NEWLINE) {
								pNode->st = SyntaxNode::s_func_name;//Encounter end of line. The last symbol is function name
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_exp_end);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode; // END of function declaration.
							}
							else if (tokenType == TK_BRKT_L) {
								pNode->st = SyntaxNode::s_func_name;//Encounter (. The last symbol is function name
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_input);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else if (tokenType == TK_SYMBOL || tokenType == TK_KW_RETURN) {
								pNode->st = SyntaxNode::s_func_name;//Encounter following symbol/return. The last symbol is function name
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_exp_assn);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode; // END of function declaration.
							}
							else if (tokenType == TK_KW_end) {
								pNode->st = SyntaxNode::s_func_name;//Encounter following end. The last symbol is function name
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_func_end);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode; // END of function, return to the function level.
							}
							else {
								errMsg(line, cursurRecord, "Only expect {=;,(}.");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_SQBRKT_L) {
							vector<string> delims = { " ","...","]","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL || tokenType == TK_SQBRKT_R) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode, SyntaxNode::s_output_list);
								pNode->pChildren.push_back(&sNode); // Open new children for output list
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables or ].");
								return MFILE_READER_FAIL;
							}
						}
					}
					else if (pNode->st == SyntaxNode::s_output_list) {
						if (pNode->tokenType == TK_SYMBOL) {
							vector<string> delims = { " ","...","]","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL || tokenType == TK_COMMA || tokenType == TK_SQBRKT_R) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_output_list);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables, comma or ].");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_COMMA) {
							vector<string> delims = { " ","...",")","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_output_list);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables.");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_SQBRKT_R) {
							vector<string> delims = { " ","...","=" };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_EQUAL) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_func_name);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect =.");
								return MFILE_READER_FAIL;
							}
						}
					}
					else if (pNode->st == SyntaxNode::s_func_name) {
						if (pNode->tokenType == TK_EQUAL) {
							vector<string> delims = { " ","...","(",";","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_name);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect function name.");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_SYMBOL) {
							vector<string> delims = { " ","...","(",";","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SEMCOLON || tokenType == TK_COMMA || tokenType == TK_NEWLINE) {
								pNode->st = SyntaxNode::s_func_name;//Encounter end of line. The last symbol is function name
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_exp_end);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode; // END of function declaration, return to the upper level.
							}
							else if (tokenType == TK_BRKT_L) {
								pNode->st = SyntaxNode::s_func_name;//Encounter (. The last symbol is function name
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_input);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else if (tokenType == TK_SYMBOL || tokenType == TK_KW_RETURN) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_exp_assn);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode; // END of function declaration, return to the upper level.
							}
							else if (tokenType == TK_KW_end) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_func_end);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = &sNode; // END of function, return to the function level.
							}
							else {
								errMsg(line, cursurRecord, "Only expect {;,(}.");
								return MFILE_READER_FAIL;
							}
						}
					}
					else if (pNode->st == SyntaxNode::s_func_input) {
						if (pNode->tokenType == TK_BRKT_L)
						{
							vector<string> delims = { " ","...",")","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL) {//After ( expect variable
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode, SyntaxNode::s_func_input);
								pNode->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables.");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_SYMBOL) {
							vector<string> delims = { " ","...",")","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL || tokenType == TK_COMMA || tokenType == TK_BRKT_R) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_input);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables, comma or ).");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_COMMA) {
							vector<string> delims = { " ","...",")","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SYMBOL) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_input);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables.");
								return MFILE_READER_FAIL;
							}
						}
						else if (pNode->tokenType == TK_BRKT_R) {
							vector<string> delims = { " ","...",")","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SEMCOLON || tokenType == TK_COMMA || tokenType == TK_NEWLINE) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent->pParent, SyntaxNode::s_exp_end);
								pNode->pParent->pParent->pChildren.push_back(&sNode);
								pNode = pNode->pParent->pParent;
							}
							else {
								errMsg(line, cursurRecord, "Only expect variables.");
								return MFILE_READER_FAIL;
							}
						}
					}
					else if (pNode->st == SyntaxNode::s_func_end) {
						if (pNode->tokenType == TK_KW_end) {
							vector<string> delims = { " ","=","...",";","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_NEWLINE) { // getting TK_NEWLINE means the end of the line.
								//SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_exp_end);
								//pNode->pParent->pChildren.push_back(&sNode);
								//pNode = &sNode;
							}
							else if (tokenType == TK_KW_function) {

							}
							else if (tokenType == TK_KW_end) {

							}
							else if (tokenType == TK_KW_RETURN) {

							}
							else if (tokenType == TK_SYMBOL) {

							}
						}
					}
					else if (pNode->st == SyntaxNode::s_exp_end) {
						if (pNode->tokenType == TK_SEMCOLON ||
							pNode->tokenType == TK_COMMA ||
							pNode->tokenType == TK_NEWLINE) {

							vector<string> delims = { " ","=","...",";","," };
							int tokenType = TK_UNDETERMINED;
							MACRO_READ_SKIP_TRIDOTS;

							if (tokenType == TK_SEMCOLON || tokenType == TK_COMMA || tokenType == TK_NEWLINE) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_exp_end);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else if (tokenType == TK_SYMBOL || tokenType == TK_KW_RETURN) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_exp_assn);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else if (tokenType == TK_KW_end) {
								SyntaxNode sNode = SyntaxNode(currToken, tokenType, pNode->pParent, SyntaxNode::s_func_end);
								pNode->pParent->pChildren.push_back(&sNode);
								pNode = &sNode;
							}
							else {
								errMsg(line, cursurRecord, "Only expect {;,(}.");
								return MFILE_READER_FAIL;
							}
							// TODO: study the effect of TK_NEWLINE and add branches.
						}
					}
				}

				// Translate to MData

				return MFILE_READER_SUCCESS;
			}
		}
	}
}