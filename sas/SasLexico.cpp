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
#include "SasLexico.h"
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdarg>

using namespace std;

namespace che {
	namespace core {

		static const char *opNames[]{
#define OPINFO(op, prec, name, func, opcode) name,
#include "opinfo.txt"
			NULL,
#undef OPINFO
		};

		int writeSingleTokenWidth(AstNode *p, char* b) {
			int width = 0;
			if (p->op == OP_CONST) {
				if (p->ty == TY_INT) {
					width += sprintf(b, "i(%d)", p->value.i[0]);
				}
				else {
					width += sprintf(b, "d(%5.4g)", p->value.d);
				}
			}
			else if (p->op == OP_STR) {
				width += sprintf(b, "(%s)", (char*)p->value.p);
			}
			else if (p->op == OP_ID) {
				width += sprintf(b, "(%s)", (char*)p->value.p);
			}
			else if (p->op != OP_NONE) {
				width += sprintf(b, "(%s)", opNames[p->op]);
			}
			else {
				width += sprintf(b, "(%s)", "N/A");
			}
			return width;
		}

		void getWidth(AstNode *p, int& width) {
			char b[20];
			if (p->subs[0] != NULL) {
				getWidth(p->subs[0], width);
			}
			if (p->subs[1] != NULL) {
				getWidth(p->subs[1], width);
			}
			width += writeSingleTokenWidth(p, b);
		}

		int _printTree(AstNode *p, int isLeft, int offset, int depth, int parentWidth, char **s) {
			char b[20];

			if (!p) return 0;

			int width = writeSingleTokenWidth(p, b);

			int left = _printTree(p->subs[0], 1, offset, depth + 1, width, s);
			int right = _printTree(p->subs[1], 0, offset + left + width, depth + 1, width, s);

			//COMPACT
			/*for (int i = 0; i < width; i++)
				s[depth][offset + left + i] = b[i];

			if (depth && isLeft) {

				for (int i = 0; i < width + right; i++)
					s[depth - 1][offset + left + width / 2 + i] = '-';

				s[depth - 1][offset + left + width / 2] = '.';

			}
			else if (depth && !isLeft) {

				for (int i = 0; i < left + width; i++)
					s[depth - 1][offset - width / 2 + i] = '-';

				s[depth - 1][offset + left + width / 2] = '.';
			}*/

			//NON-COMPACT
			for (int i = 0; i < width; i++)
				s[2 * depth][offset + left + i] = b[i];

			if (depth && isLeft) {

				for (int i = 0; i < width - width / 2 + right + parentWidth / 2; i++)
					s[2 * depth - 1][offset + left + width / 2 + i] = '-';

				s[2 * depth - 1][offset + left + width / 2] = '+';
				s[2 * depth - 1][offset + left + width + right + parentWidth / 2] = '+';

			}
			else if (depth && !isLeft) {

				for (int i = 0; i < left + width/2 + parentWidth - parentWidth / 2; i++)
					s[2 * depth - 1][offset + parentWidth / 2 - parentWidth + i] = '-';

				s[2 * depth - 1][offset + parentWidth / 2 - parentWidth] = '+';
				s[2 * depth - 1][offset + left + width / 2] = '+';
			}

			return left + width + right;
		}

		int exploreDepth(AstNode* p) {
			if (p == NULL) {
				return 0;
			}
			int ld = 0;
			int rd = 0;
			if (p->subs[0] != NULL) {
				ld = exploreDepth(p->subs[0]);
			}
			if (p->subs[1] != NULL) {
				rd = exploreDepth(p->subs[1]);
			}
			if (rd > ld) {
				ld = rd;
			}
			return ld + 1;
		}

		void AstTree::refreshDepth() {
			depth = exploreDepth(this->pHead);
		}

		void AstTree::printTree(ostream& ost) {
			int width = 0;
			refreshDepth();
			if (pHead != NULL) {
				getWidth(pHead, width);
			}
			width += 3;
			if (depth > 0 && width>0) {
				std::streambuf * old = std::cout.rdbuf(ost.rdbuf());

				char **s = new char*[2 * depth];
				for (int i = 0; i < 2 * depth; i++) {
					s[i] = new char[width + 1];
					memset(s[i], ' ', (width) * sizeof(char));
					s[i][width] = '\0';
				}

				_printTree(pHead, 0, 0, 0, 0, s);

				for (int i = 0; i < 2 * depth; i++) {
					s[i][width] = '\0';
				}
				cout << "**************************" << endl;
				cout << Lexer::getInstance().eqnNames[this->eqnType]<<" ("<<this->eqnIdx<<")"<< endl;
				for (int i = 0; i < 2 * depth; i++) {
					cout << s[i] << endl;
				}
				cout << "**************************" << endl;

				for (int i = 0; i < 2 * depth; i++) {
					delete [] s[i];
				}
				delete[] s;

				std::cout.rdbuf(old);
			}
		}

		void freeAstTreeRecursive(AstNode *p) {
			if (p == NULL) {
				return;
			}
			if (p->subs[0] != NULL) {
				freeAstTreeRecursive(p->subs[0]);
				p->subs[0] = NULL;
			}
			if (p->subs[1] != NULL) {
				freeAstTreeRecursive(p->subs[1]);
				p->subs[1] = NULL;
			}
			if (p->next != NULL) {
				freeAstTreeRecursive(p->next);
				p->next = NULL;
			}
			delete p;
		}

		AstTree::~AstTree() {
			if (pHead != NULL) {
				freeAstTreeRecursive(pHead);
				pHead = NULL;
			}
			depth = 0;
		}

		/*void SasModelParser::error(const char* filename, int line, const char* errorMsg) {
			cerr << "(" << filename << "," << line << ") error: " << errorMsg << endl;
			errCount++;
		}

		void SasModelParser::warn(const char* filename, int line, const char* errorMsg) {
			cerr << "(" << filename << "," << line << ") warning: " << errorMsg << endl;
			warnCount++;
		}*/

		SasModelParser::SasModelParser() :reader() {}

		void SasModelParser::error(const char* formatMsg, ...) {
			va_list ap;
			char buffer[MAX_PRINT_LEN + 1];
			va_start(ap, formatMsg);
			vsprintf(buffer, formatMsg, ap);
			cerr << "(" << coord.filename << "," << coord.line << ") error: " << buffer << endl;
			errCount++;
			va_end(ap);
		}

		void SasModelParser::warn(const char* formatMsg, ...) {
			va_list ap;
			char buffer[MAX_PRINT_LEN + 1];
			va_start(ap, formatMsg);
			vsprintf(buffer, formatMsg, ap);
			cerr << "(" << coord.filename << "," << coord.line << ") warn: " << buffer << endl;
			warnCount++;
			va_end(ap);
		}

		void SasModelParser::parseString(const char* str) {
			reader.resetReader();
			int len = strlen(str);
			reader.buffer = new uchar[len + 1];
			strncpy((char*)reader.buffer, str, len);
			reader.buffer[len] = char_traits<uchar>::eof();
			reader.size = len;
			reader.cursur = reader.buffer;
			reader.lineHead = reader.buffer;

			const char* fn = "String";
			int fnlen = strlen(fn);
			reader.filename = new char[fnlen +1];
			strncpy(reader.filename, fn, fnlen);
			reader.filename[fnlen] = 0;
			coord.filename = reader.filename;
			prevCoord = coord;
			peekCoord = coord;
		}

		void SasModelParser::parseFromFile(const char* filename) {
			reader.readFile(filename);
			coord.filename = reader.filename;
			prevCoord = coord;
			peekCoord = coord;
		}


		int SasModelParser::bufferMargin() {
			return reader.size - (reader.cursur - reader.buffer);
		}

		bool SasModelParser::hasReachedBufferSize() {
			return (reader.size - (reader.cursur - reader.buffer)) <= 0;
		}

		void SasModelParser::beginPeekToken() {
			peekPoint = reader.cursur;
			peekValue = value;
			peekCoord = coord;
		}

		void SasModelParser::endPeekToken() {
			reader.cursur = peekPoint;
			value = peekValue;
			coord = peekCoord;
		}

		SasModelParser::~SasModelParser() {
			if (currentModel != nullptr) {
				currentModel.reset();
				currentModel = nullptr;
			}
			if (globalPool != nullptr) {
				globalPool.reset();
				globalPool = nullptr;
			}
			if (currentId != nullptr) {
				currentId.reset();
				currentId = nullptr;
			}
		}

		void setupLexer(scanner* scanners);
		Lexer::Lexer() {
			setupLexer(scanners);
		}
	}
}