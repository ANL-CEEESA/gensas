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

#include "SasExpr.h"
#include "SasLexico.h"
#include <iostream>

using namespace std;

namespace che {
	namespace core {

		static TokenOp tokenOps[]{
#define TOKENOP(tok, bop, uop) {bop, uop},
#include "sas/tokenop.txt"
			0,
#undef  TOKENOP
		};

		static int tokenPrec[]{
#define OPINFO(op, prec, name, func, opcode) prec,
#include "sas/opinfo.txt"
#undef OPINFO
		};

		static const char *opNames[]{
#define OPINFO(op, prec, name, func, opcode) name,
#include "sas/opinfo.txt"
			NULL,
#undef OPINFO
		};

		void parseSingleSasModel(SasModelParser* parser) {
			getNextToken(parser);
			if (parser->currentToken != TK_K_MDL) {
				parser->error("actual token %s: expect \"model\".",Lexer::getInstance().tokenStrings[parser->currentToken]);
			}
			getNextToken(parser);
			if (parser->currentToken != TK_ID) {
				parser->error("actual token %s: expect identifier.", Lexer::getInstance().tokenStrings[parser->currentToken]);
				parser->currentId=parser->currentModel->lookupIdentifier("Unknown", 7);
			}
			parser->currentModel->modelInfo = parser->currentId;
			parser->currentId->idTy = ID_MODEL;

			getNextToken(parser);
			if (parser->currentToken == TK_STRING) {
				parser->currentModel->modelInfo->setAnnotation((uchar*)parser->value.p);
				getNextToken(parser);
			}

			// parse variable decl portion
			while (parser->currentToken == TK_K_REAL || parser->currentToken == TK_K_PARA) {
				IdType idTy = ID_VAR;
				if (parser->currentToken == TK_K_PARA) {//parameter
					idTy = ID_PARAM;
					getNextToken(parser);
					expect(parser, TK_K_REAL);
				}
				else {
					getNextToken(parser);
				}
				if (parser->currentToken != TK_ID) {
					parser->error("actual token %s: expect identifier.", Lexer::getInstance().tokenStrings[parser->currentToken]);
				}
				shared_ptr<IdPearl> var = parser->currentId;
				getNextToken(parser);
				if (parser->currentToken == TK_LPAREN) {
					getNextToken(parser);
					expect(parser, TK_K_UNIT);
					expect(parser, TK_EQN);
					if (parser->currentToken == TK_STRING) {
						var->setUnit((uchar*)parser->value.p);
						getNextToken(parser);
					}
					else {
						parser->error("actual token %s: expect string.", Lexer::getInstance().tokenStrings[parser->currentToken]);
					}

					expect(parser, TK_RPAREN);
				}

				if (idTy == ID_PARAM) {
					expect(parser, TK_EQN);
					//expect a value
					int sign = 1;
					while (parser->currentToken == TK_ADD || parser->currentToken == TK_SUB) {
						sign *= (parser->currentToken == TK_ADD ? 1 : -1);
						getNextToken(parser);
					}
					if (parser->currentToken == TK_DOUBLECONST ||
						parser->currentToken == TK_FLOATCONST ||
						parser->currentToken == TK_REALCONST) {
						var->value.d = (double)sign*parser->value.d;
						var->initVal = true;
					}
					else if (parser->currentToken == TK_INTCONST || parser->currentToken == TK_LONGCONST) {
						var->value.d = (double)(sign*parser->value.i[0]);
						var->initVal = true;
					}
					else {
						parser->error("Expect parameter value.");
						var->value.d = NAN;
						var->initVal = true;
					}
					getNextToken(parser);
				}
				if (parser->currentToken == TK_STRING) {
					var->setAnnotation((uchar*)parser->value.p);
					getNextToken(parser);
				}
				var->idTy = idTy;
				var->ty = TY_REAL;
				expect(parser, TK_SEMICOLON);
			}

			//initial block
			if (parser->currentToken != TK_K_INIT) {
				parser->error("Expect keyword initial.");
			}
			getNextToken(parser);
			if (parser->currentToken != TK_K_EQN) {
				parser->error("Expect keyword equation.");
			}
			getNextToken(parser);
			while (parser->currentToken == TK_ID|| parser->currentToken == TK_F_DER) {
				shared_ptr<IdPearl> var = parser->currentId;
				bool initDer = false;
				if (parser->currentToken == TK_ID&&var->idTy == ID_VAR && var->ty == TY_REAL) {

				}
				else if(var->idTy==ID_FUNC&& parser->currentToken == TK_F_DER){
					getNextToken(parser);
					expect(parser, TK_LPAREN);
					if (parser->currentToken == TK_ID && var->idTy == ID_VAR && var->ty == TY_REAL) {
						initDer = true;
						getNextToken(parser);
						expect(parser, TK_RPAREN);
					}
					else {
						parser->error("der(): expect id inside ().");
					}
				}
				getNextToken(parser);
				expect(parser, TK_EQN);
				//expect a value
				int sign = 1;
				while (parser->currentToken == TK_ADD || parser->currentToken == TK_SUB) {
					sign *= (parser->currentToken == TK_ADD ? 1 : -1);
					getNextToken(parser);
				}
				if (parser->currentToken == TK_DOUBLECONST ||
					parser->currentToken == TK_FLOATCONST ||
					parser->currentToken == TK_REALCONST) {
					if (initDer) {
						var->derVal.d = (double)sign*parser->value.d;
						var->initDer = true;
					}
					else {
						var->value.d = (double)sign*parser->value.d;
						var->initVal = true;
					}
				}
				else if (parser->currentToken == TK_INTCONST || parser->currentToken == TK_LONGCONST) {
					if (initDer) {
						var->derVal.d = (double)(sign*parser->value.i[0]);
						var->initDer = true;
					}
					else {
						var->value.d = (double)(sign*parser->value.i[0]);
						var->initVal = true;
					}
				}
				else {
					parser->error("Expect variable initial value.");
					if (initDer) {
						var->derVal.d = NAN;
						var->initDer = true;
					}
					else {
						var->value.d = NAN;
						var->initVal = true;
					}
				}
				getNextToken(parser);
				expect(parser, TK_SEMICOLON);
			}

			//equation block
			if (parser->currentToken != TK_K_EQN) {
				parser->error("Expect keyword equation.");
			}
			getNextToken(parser);
			while (parser->currentToken != TK_K_END && parser->currentToken != TK_EOF) {
				uchar* cur = parser->reader.cursur;
				shared_ptr<AstTree> tree = make_shared<AstTree>();
				tree->pHead = parseExpr(parser);
				parser->currentModel->modelTrees.push_back(tree);
				expect(parser, TK_SEMICOLON);
				if (parser->reader.cursur == cur) {
					getNextToken(parser);
				}
			}

			if (parser->currentToken == TK_K_END) {
				getNextToken(parser);

				getNextToken(parser);
				expect(parser, TK_SEMICOLON);
			}
			else {
				parser->error("Missing end statement of model.");
			}
		}
	}
}