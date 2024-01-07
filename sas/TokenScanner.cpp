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
#include "sas/SasLexico.h"
#include <iostream>
#include <cstring>
#include <limits.h>
using namespace std;

namespace che {
	namespace core {
		static struct keyword keywords_[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsA[] =
		{
			{"algorithm", 9, TK_K_ALGORITHM},
			{"and", 3, TK_K_AND},
			{"annotation", 10, TK_K_ANNOTATION},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsB[] =
		{
			{"block", 5, TK_K_BLOCK},
			{"break", 5, TK_K_BREAK},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsC[] =
		{
			{"class",5,TK_K_CLASS},
			{"connect",7,TK_K_CONNECT},
			{"connector",9,TK_K_CONNECTOR},
			{"constant",8,TK_K_CONSTANT},
			{"constrainedby",13,TK_K_CONSTRAINEDBY},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsD[] =
		{
			{"discrete",8,TK_K_DISCRETE},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsE[] =
		{
			{"each",4,TK_K_EACH},
			{"else",4,TK_K_ELSE},
			{"elseif",6,TK_K_ELSEIF},
			{"elsewhen",8,TK_K_ELSEWHEN},
			{"encapsulated",12,TK_K_ENCAPSULATED},
			{"end",3,TK_K_END},
			{"enumeration",11,TK_K_ENUMERATION},
			{"expandable",10,TK_K_EXPANDABLE},
			{"extends",7,TK_K_EXTENDS},
			{"external",8,TK_K_EXTERNAL},			
			{"equation", 8, TK_K_EQN},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsF[] =
		{
			{"false",5,TK_K_FALSE},
			{"final",5,TK_K_FINAL},
			{"flow",4,TK_K_FLOW},
			{"for",3,TK_K_FOR},
			{"function",8,TK_K_FUNCTION},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsG[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsH[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsI[] =
		{
			{"if",2,TK_K_IF},
			{"import",6,TK_K_IMPORT},
			{"impure",6,TK_K_IMPURE},
			{"in",2,TK_K_IN},
			{"inner",5,TK_K_INNER},
			{"input",5,TK_K_INPUT},
			{"loop",4,TK_K_LOOP},
			{"initial", 7, TK_K_INIT},
			{"ivar", 4, TK_K_IVAR},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsJ[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsK[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsL[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsM[] =
		{
			{"model", 5, TK_K_MDL},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsN[] =
		{
			{"not", 3, TK_K_NOT},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsO[] =
		{
			{"outer",5,TK_K_OUTER},
			{"output",6,TK_K_OUTPUT},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsP[] =
		{
			{"package",7,TK_K_PACKAGE},
			{"partial",7,TK_K_PARTIAL},
			{"protected",9,TK_K_PROTECTED},
			{"public",6,TK_K_PUBLIC},
			{"pure",4,TK_K_PURE},
			{"parameter", 9, TK_K_PARA},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsQ[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsR[] =
		{
			{"record",6,TK_K_RECORD},
			{"redeclare",9,TK_K_REDECLARE},
			{"replaceable",11,TK_K_REPLACEABLE},
			{"return",6,TK_K_RETURN},
			{"Real", 4, TK_K_REAL},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsS[] =
		{
			{"stream",6,TK_K_STREAM},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsT[] =
		{
			{"then",4,TK_K_THEN},
			{"true",4,TK_K_TRUE},
			{"type",4,TK_K_TYPE},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsU[] =
		{
			{"unit", 4, TK_K_UNIT},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsV[] =
		{
			{"var", 3, TK_K_VAR},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsW[] =
		{
			{"when",4,TK_K_WHEN},
			{"while",5,TK_K_WHILE},
			{"within",6,TK_K_WITHIN},
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsX[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsY[] =
		{
			{NULL, 0, TK_ID}
		};

		static struct keyword keywordsZ[] =
		{
			{NULL, 0, TK_ID}
		};
		/**
			classify keywords by their first letter,
			to speed up comparing.
			see function FindKeyword()
		 */
		static struct keyword *keywords[] =
		{
			keywords_, keywordsA, keywordsB, keywordsC,
			keywordsD, keywordsE, keywordsF, keywordsG,
			keywordsH, keywordsI, keywordsJ, keywordsK,
			keywordsL, keywordsM, keywordsN, keywordsO,
			keywordsP, keywordsQ, keywordsR, keywordsS,
			keywordsT, keywordsU, keywordsV, keywordsW,
			keywordsX, keywordsY, keywordsZ
		};

		static func defaultFunctions[] =
		{
			{"der", 3, TK_F_DER},
			{"sin", 3, TK_F_SIN},
			{"cos", 3, TK_F_COS},
			{"exp", 3, TK_F_EXP},
			{"sqrt", 4, TK_F_SQRT},
			{"inv", 3, TK_F_INV},
			{NULL, 0, TK_ID}
		};


		inline bool isDigit(char c) { return c >= '0'&&c <= '9'; }
		inline bool isOctDigit(char c){ return c >= '0'&&c <= '7'; }
		inline bool isHexDigit(char c) { return (c >= '0'&&c <= '9') || (c >= 'A' && c <= 'F') || (c >= 'a' && c <= 'f'); }
		inline bool isLetter(char c) { return (c >= 'a' && c <= 'z') || (c == '_') || (c >= 'A' && c <= 'Z'); }
		inline char toUpper(char c) { return c & ~0x20; }
		inline bool isLetterOrDigit(char c){ return (c >= '0'&&c <= '9')||(c >= 'a' && c <= 'z') || (c == '_') || (c >= 'A' && c <= 'Z'); }
		inline int high4Bit(int v) { return ((v) >> (8 * sizeof(int) - 4) & 0x0f); }
		inline int high3Bit(int v) { return ((v) >> (8 * sizeof(int) - 3) & 0x07); }
		inline int high2Bit(int v) { return ((v) >> (8 * sizeof(int) - 2) & 0x03); }
		inline int high1Bit(int v) { return ((v) >> (8 * sizeof(int) - 1) & 0x01); }

		const char* TokenStrings[] =
		{
		#define TOKEN(k, s) s,
		#include "sas/token.txt"
		#undef  TOKEN
		};

		static void skipWhiteSpace(SasModelParser *parser) {
			int ch = *(parser->reader.cursur);
			while (ch == '\t' || ch == '\v' || ch == '\f' || ch == ' ' ||
				ch == '\r' || ch == '\n' || ch == '/' || ch == '#') {
				switch (ch) {
				case '\n':
					parser->coord.line++;
					parser->reader.line++;
					parser->reader.lineHead = ++parser->reader.cursur;
					break;
				case '/':
					if (parser->reader.cursur[1] != '/' && parser->reader.cursur[1] != '*')
						return;
					parser->reader.cursur++;
					if (*(parser->reader.cursur) == '/') {
						parser->reader.cursur++;
						while (*(parser->reader.cursur) != '\n' &&
							!(*(parser->reader.cursur) == END_OF_FILE && parser->hasReachedBufferSize())) {
							parser->reader.cursur++;
						}
					}
					else {
						parser->reader.cursur++;
						while (parser->reader.cursur[1] != '*' && parser->reader.cursur[1] != '/') {
							if (*(parser->reader.cursur) != '\n') {
								parser->coord.line++;
								parser->reader.line++;
								parser->reader.lineHead = parser->reader.cursur+1;
							}
							else if ((*(parser->reader.cursur) == END_OF_FILE && parser->hasReachedBufferSize()) ||
								(*(parser->reader.cursur + 1) == END_OF_FILE && parser->bufferMargin()<=1)) {
								parser->error(parser->reader.filename, parser->reader.line, "Comment block is not closed.");
								return;
							}
							parser->reader.cursur++;
						}
						parser->reader.cursur+=2;
					}
					break;

				default:
					parser->reader.cursur++;
					break;
				}
				ch = *(parser->reader.cursur);
			}
		}

		static int scanEscapeChar(SasModelParser *parser) {
			int v = 0;
			parser->reader.cursur++;
			switch (*(parser->reader.cursur)++) {
			case 'a':
				return '\a';
			case 'b':
				return '\b';
			case 'f':
				return '\f';
			case 'n':
				return '\n';
			case 'r':
				return '\r';
			case 't':
				return '\t';
			case 'v':
				return '\v';
			case '\'':
			case '"':
			case '\\':
			case '\?':
				return *(parser->reader.cursur - 1);
			case 'x':
				if (!isHexDigit(*(parser->reader.cursur))) {
					parser->error(parser->reader.filename, parser->reader.line, "Expect Hex digit.");
					return 'x';
				}
				v = 0;
				while (isHexDigit(*(parser->reader.cursur))) {
					if (isDigit(*(parser->reader.cursur))) {
						v = (v << 4) + *(parser->reader.cursur) - '0';
					}
					else {
						v = (v << 4) + toUpper(*(parser->reader.cursur)) - 'A'+10;
					}
					parser->reader.cursur++;
				}
				return v;
			case '0': case '1': case '2': case '3':
			case '4': case '5': case '6': case '7':
				v = *(parser->reader.cursur - 1) - '0';
				if (isOctDigit(*(parser->reader.cursur))) {
					v = (v << 3) + *(parser->reader.cursur)++ - '0';
					if (isOctDigit(*(parser->reader.cursur))) {
						v = (v << 3) + *(parser->reader.cursur)++ - '0';
					}
				}
				return v;
			default:
				return *(parser->reader.cursur);
			}
		}

		static int findKeyWord(char* str, int len) {
			keyword *p = NULL;
			int index = 0;
			if (*str != '_') {
				index = toUpper(*str) - 'A' + 1;
			}
			p = keywords[index];
			while (p->name) {
				if (p->len == len && strncmp(str, p->name, len) == 0) {
					break;
				}
				p++;
			}
			return p->tok;
		}

		static int findDefaultFunctions(char* str, int len) {
			func *p = defaultFunctions;
			while (p->name) {
				if (p->len == len && strncmp(str, p->name, len) == 0) {
					break;
				}
				p++;
			}
			return p->tok;
		}

		static shared_ptr<IdPearl> lookupIdentifier(SasModelParser *parser,const char* iden, int len) {
			return parser->currentModel->lookupIdentifier(iden, len);
		}

		static int scanIdentifier(SasModelParser *parser) {
			uchar* start = parser->reader.cursur;
			parser->reader.cursur++;
			while (isLetterOrDigit(*(parser->reader.cursur))) {
				parser->reader.cursur++;
			}
			parser->currentId.reset();
			int tok = findKeyWord((char*)start, (int)(parser->reader.cursur - start));
			if (tok == TK_ID) {
				tok = findDefaultFunctions((char*)start, (int)(parser->reader.cursur - start));
				parser->currentId=lookupIdentifier(parser, (char*)start, (int)(parser->reader.cursur - start));
				if (tok != TK_ID) {
					parser->currentId->idTy = ID_FUNC;
				}
				parser->value.p = (char*)parser->currentId->str;
			}
			if (tok == TK_ID) {
				// record ID value
				parser->currentId=lookupIdentifier(parser,(char*)start, (int)(parser->reader.cursur - start));
				parser->value.p = (char*)parser->currentId->str;
			}
			return tok;
		}

		static int scanIntLiteral(SasModelParser *parser,uchar* start, int len) {
			uchar *p = start;
			uchar *endp = start + len;
			unsigned int i[2] = { 0,0 };
			int tok = TK_INTCONST;
			int d = 0;
			int carry0 = 0, carry1 = 0;
			int overflow = 0;

			while (p != endp) {
				d = *p - '0';
				unsigned int t1, t2;
				// number * 10 = number * 8 + number * 2 = (number << 3) + (number << 1)
				carry0 = high3Bit(i[0]) + high1Bit(i[0]);
				carry1 = high3Bit(i[1]) + high1Bit(i[1]);
				t1 = i[0] << 3;
				t2 = i[0] << 1;
				if (t1 > UINT_MAX - t2)	// In maths:  t1 + t2 > UINT_MAX
				{
					carry0++;
				}
				i[0] = t1 + t2;
				t1 = i[1] << 3;
				t2 = i[1] << 1;
				if (t1 > UINT_MAX - t2)
				{
					carry1++;
				}
				i[1] = t1 + t2;

				if (i[0] > UINT_MAX - d)	// for decimal, i[0] + d maybe greater than UINT_MAX
				{
					carry0 += i[0] - (UINT_MAX - d);
				}
				if (carry1 || (i[1] > UINT_MAX - carry0))
				{
					overflow = 1;
				}
				i[0] += d;
				i[1] += carry0;
				p++;
			}
			if (overflow || i[1] != 0) {
				parser->warn(parser->reader.filename, parser->reader.line, "Integer literal is too big.");
			}
			parser->value.i[1] = 0;
			parser->value.i[0] = i[0];
			tok = TK_INTCONST;
			if (*(parser->reader.cursur) == 'U' || *(parser->reader.cursur) == 'u')	{
				parser->reader.cursur++;
				parser->warn(parser->reader.filename, parser->reader.line, "Unsigned int are treated as int.");
			}
			if (*(parser->reader.cursur) == 'L' || *(parser->reader.cursur) == 'l') {
				parser->reader.cursur++;
				if (*(parser->reader.cursur) == 'L' || *(parser->reader.cursur) == 'l') {
					parser->reader.cursur++;
				}
				parser->warn(parser->reader.filename, parser->reader.line, "Long int are treated as int.");
			}
			return tok;
		}

		static int scanFloatLiteral(SasModelParser *parser, uchar* start) {
			double d;
			if (*(parser->reader.cursur) == '.' ) {
				parser->reader.cursur++;
				while (isDigit(*(parser->reader.cursur))) {
					parser->reader.cursur++;
				}
			}
			if (*(parser->reader.cursur) == 'E' || *(parser->reader.cursur) == 'e') {
				parser->reader.cursur++;
				if (*(parser->reader.cursur) == '+' || *(parser->reader.cursur) == '-') {
					parser->reader.cursur++;
				}
				if (!isDigit(*(parser->reader.cursur))) {
					parser->error(parser->reader.filename, parser->reader.line, "Expect exponent in float.");
				}
				else {
					while (isDigit(*(parser->reader.cursur))) {
						parser->reader.cursur++;
					}
				}
			}
			d = strtod((char*)start, NULL);
			parser->value.d = d;
			if (*(parser->reader.cursur) == 'F' || *(parser->reader.cursur) == 'f') {
				parser->reader.cursur++;
				parser->warn(parser->reader.filename, parser->reader.line, "Floats are treated as double.");
			}else if (*(parser->reader.cursur) == 'L' || *(parser->reader.cursur) == 'l') {
				parser->reader.cursur++;
				parser->warn(parser->reader.filename, parser->reader.line, "Long double are treated as double.");
			}
			return TK_DOUBLECONST;
		}

		static int scanNumericalLiteral(SasModelParser *parser) {
			uchar *start = parser->reader.cursur;
			int base = 10;
			if (*(parser->reader.cursur) == '.') {
				return scanFloatLiteral(parser, start);
			}
			if (*(parser->reader.cursur) == '0' && (parser->reader.cursur[1] == 'x' || parser->reader.cursur[1] == 'X')) {
				parser->reader.cursur+=2;
				start = parser->reader.cursur;
				base = 16;
				if (!isHexDigit(*(parser->reader.cursur))) {
					parser->error(parser->reader.filename, parser->reader.line, "Expect HEX digit.");
					parser->value.i[0] = 0;
					return TK_INTCONST;
				}
				while (isHexDigit(*(parser->reader.cursur))) {
					parser->reader.cursur++;
				}
			}
			else if (*(parser->reader.cursur) == '0') {
				parser->reader.cursur++;
				while (isOctDigit(*(parser->reader.cursur))) {
					parser->reader.cursur++;
				}
			}
			else {
				parser->reader.cursur++;
				while (isDigit(*(parser->reader.cursur))) {
					parser->reader.cursur++;
				}
			}

			if (base == 16 || (parser->reader.cursur[0] != '.'&&parser->reader.cursur[0] != 'E'&&parser->reader.cursur[0] != 'e')) {
				return scanIntLiteral(parser, start, (int)(parser->reader.cursur - start));
			}
			else {
				return scanFloatLiteral(parser, start);
			}
		}

		static int scanStringLiteral(SasModelParser *parser) {
			char tmp[MAX_STR_LEN + 1];
			tmp[MAX_STR_LEN] = 0;
			int len = 0;
			char ch=0;

			if (*(parser->reader.cursur) == 'L') {
				parser->reader.cursur++;
				parser->warn(parser->reader.filename, parser->reader.line, "Wide string not supported and will be stored as string.");
			}
			parser->reader.cursur++;
			while (*(parser->reader.cursur) != '"') {
				if (*(parser->reader.cursur) == END_OF_FILE && parser->hasReachedBufferSize()) {
					break;
				}
				if (*(parser->reader.cursur) == '\n') {
					break;
				}
				if (*(parser->reader.cursur) == '\\') {
					ch = scanEscapeChar(parser);
				}
				else {
					ch = *(parser->reader.cursur);
					parser->reader.cursur++;
				}
				if (len < MAX_STR_LEN) {
					tmp[len++] = ch;
				}
			}
			if (*(parser->reader.cursur) != '"') {
				parser->error(parser->reader.filename, parser->reader.line, "Expect end of string \".");
			}
			else {
				parser->reader.cursur++;
			}
			GlobalPool::getInstance().strTable.push_back(make_shared<StrPearl>((uchar*)tmp, len));
			parser->value.p = GlobalPool::getInstance().strTable.back()->str;

			return TK_STRING;
		}

		static int scanEqual(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) != '=') {
				return TK_EQN;
			}
			else {
				parser->reader.cursur++;
				return TK_EQUAL;
			}
		}

		static int scanAnd(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) == '&') {
				parser->reader.cursur++;
				return TK_AND;
			}
			else {
				return TK_BITAND;
			}
		}

		static int scanOr(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) == '|') {
				parser->reader.cursur++;
				return TK_OR;
			}
			else {
				return TK_BITOR;
			}
		}

		static int scanGreat(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) == '=') {
				parser->reader.cursur++;
				return TK_GE;
			}
			else {
				return TK_GREAT;
			}
		}

		static int scanLess(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) == '=') {
				parser->reader.cursur++;
				return TK_LE;
			}else if (*(parser->reader.cursur) == '>') {
				parser->reader.cursur++;
				return TK_UNEQ;
			}else {
				return TK_LESS;
			}
		}

		static int scanColon(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) == '=') {
				parser->reader.cursur++;
				return TK_ASSIGN;
			}
			else {
				return TK_COLON;
			}
		}

		static int scanDot(SasModelParser *parser) {
			parser->reader.cursur++;
			if (*(parser->reader.cursur) == '+') {
				parser->reader.cursur++;
				return TK_DOTADD;
			}
			else if (*(parser->reader.cursur) == '-') {
				parser->reader.cursur++;
				return TK_DOTSUB;
			}
			else if (*(parser->reader.cursur) == '*') {
				parser->reader.cursur++;
				return TK_DOTMUL;
			}
			else if (*(parser->reader.cursur) == '/') {
				parser->reader.cursur++;
				return TK_DOTDIV;
			}
			else if (*(parser->reader.cursur) == '^') {
				parser->reader.cursur++;
				return TK_DOTPOW;
			}
			else {
				return TK_DOT;
			}
		}

#define SINGLE_CHAR_SCANNER(t) \
static int scan##t(SasModelParser *parser)       \
{                              \
    parser->reader.cursur++;   \
    return TK_##t;             \
}

		SINGLE_CHAR_SCANNER(ADD)
		SINGLE_CHAR_SCANNER(SUB)
		SINGLE_CHAR_SCANNER(MUL)
		SINGLE_CHAR_SCANNER(DIV)
		SINGLE_CHAR_SCANNER(COMMA)
		SINGLE_CHAR_SCANNER(LPAREN)
		SINGLE_CHAR_SCANNER(RPAREN)
		SINGLE_CHAR_SCANNER(LBRACKET)
		SINGLE_CHAR_SCANNER(RBRACKET)
		SINGLE_CHAR_SCANNER(LBRACE)
		SINGLE_CHAR_SCANNER(RBRACE)
		SINGLE_CHAR_SCANNER(SEMICOLON)
		SINGLE_CHAR_SCANNER(QUESTION)
		SINGLE_CHAR_SCANNER(EXCLAMATION)
		SINGLE_CHAR_SCANNER(CARET)
		SINGLE_CHAR_SCANNER(PERCENT)
		SINGLE_CHAR_SCANNER(POUND)
		SINGLE_CHAR_SCANNER(DOLLAR)
		SINGLE_CHAR_SCANNER(AT)

		int getNextToken(SasModelParser *parser);
			
		static int scanBadChar(SasModelParser *parser) {
			parser->error(parser->reader.filename, parser->reader.line, "Illegal char: \\x%x.", *(parser->reader.cursur));
			parser->reader.cursur++;
			return getNextToken(parser);
		}

		static int scanEOF(SasModelParser *parser) {
			if (*(parser->reader.cursur) == END_OF_FILE && parser->hasReachedBufferSize()) {
				return TK_EOF;
			}
			else {
				return scanBadChar(parser);
			}
		}
			
		void setupLexer(scanner* scanners) {
			int i;
			for (i = 0; i < END_OF_FILE + 1; i++) {
				if (isLetter(i)) {
					scanners[i] = scanIdentifier;
				}
				else if (isDigit(i)) {
					scanners[i] = scanNumericalLiteral;
				}
				else {
					scanners[i] = scanBadChar;
				}
				scanners[END_OF_FILE] = scanEOF;
				scanners['"'] = scanStringLiteral;
				scanners['&'] = scanAnd;
				scanners['|'] = scanOr;
				scanners['<'] = scanLess;
				scanners['>'] = scanGreat;
				scanners['='] = scanEqual;
				scanners['.'] = scanDot;
				scanners[':'] = scanColon;

				scanners['+'] = scanADD;
				scanners['-'] = scanSUB;
				scanners['*'] = scanMUL;
				scanners['/'] = scanDIV;
				scanners[','] = scanCOMMA;
				scanners['('] = scanLPAREN;
				scanners[')'] = scanRPAREN;
				scanners['['] = scanLBRACKET;
				scanners[']'] = scanRBRACKET;
				scanners['{'] = scanLBRACE;
				scanners['}'] = scanRBRACE;
				scanners[';'] = scanSEMICOLON;
				scanners['?'] = scanQUESTION;
				scanners['!'] = scanEXCLAMATION;
				scanners['^'] = scanCARET;
				scanners['%'] = scanPERCENT;
				scanners['#'] = scanPOUND;
				scanners['$'] = scanDOLLAR;
				scanners['@'] = scanAT;
			}
		}

		int getNextToken(SasModelParser *parser) {
			int tok;
			parser->prevCoord = parser->coord;
			skipWhiteSpace(parser);
			parser->coord.line = parser->reader.line;
			parser->coord.col = (int)(parser->reader.cursur - parser->reader.lineHead + 1);
			tok = (*Lexer::getInstance().scanners[*(parser->reader.cursur)])(parser);
			parser->currentToken = tok;
			return tok;
		}

		void expect(SasModelParser *parser, int tok) {
			if (parser->currentToken == tok) {
				getNextToken(parser);
				return;
			}
			if (tok == TK_SEMICOLON && parser->coord.line - parser->prevCoord.line == 1) {
				parser->error("Expect ;.");
			}
			else {
				parser->error("Expect %s.", TokenStrings[tok]);
			}
		}
	}
}