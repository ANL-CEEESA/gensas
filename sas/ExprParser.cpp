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
#include "sas/SasExpr.h"
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

		extern int getNextToken(SasModelParser *parser);
		extern void expect(SasModelParser *parser, int tok);

		static AstNode* parseAssignmentExpr(SasModelParser *parser);
		AstNode* parseExpr(SasModelParser *parser);
		static AstNode* parseConditionalExpr(SasModelParser *parser);
		static AstNode* parseBinaryExpr(SasModelParser *parser, AstNode *curExpr, int prec);

		static AstNode* badDefaultToken() {
			AstNode* expr = new AstNode();
			TokenValue v;
			v.i[1] = 0;
			v.i[0] = 0;
			expr->ty = TY_INT;
			expr->op = OP_CONST;
			expr->nodeKind = NK_Expression;
			expr->value = v;
			return expr;
		}

		static AstNode* parserPrimaryExpr(SasModelParser * parser) {
			AstNode* expr;

			switch (parser->currentToken) {
			case TK_ID:
				expr = new AstNode();
				expr->op = OP_ID;
				expr->value = parser->value;
				getNextToken(parser);
				return expr;

			case TK_INTCONST:
				expr = new AstNode();
				expr->ty = TY_INT;
				expr->op = OP_CONST;
				expr->value= parser->value;
				getNextToken(parser);
				return expr;
			case TK_FLOATCONST:
			case TK_DOUBLECONST:
				expr = new AstNode();
				expr->ty = TY_REAL;
				expr->op = OP_CONST;
				expr->value = parser->value;
				getNextToken(parser);
				return expr;
			case TK_STRING:
				expr = new AstNode();
				expr->ty = TY_STR;
				expr->op = OP_STR;
				expr->value = parser->value;
				getNextToken(parser);
				return expr;
			case TK_LPAREN:
				getNextToken(parser);
				expr = parseExpr(parser);
				expect(parser, TK_RPAREN);
				return expr;
			default:
				if (parser->currentToken >= TK_F_DER && parser->currentToken < NUM_TK) {
					expr = new AstNode();
					expr->ty = TY_FUNCTION;
					expr->op = OP_ID;
					expr->value = parser->value;
					getNextToken(parser);
					return expr;
				}
				else {
					parser->error("[primExpr]Expect identifier, string, constant or (*).");
					expr = badDefaultToken();
					return expr;
				}
			}
		}

		static AstNode* parsePostfixExpr(SasModelParser *parser) {
			AstNode *expr, *p;
			expr = parserPrimaryExpr(parser);

			while (true) {
				switch (parser->currentToken) {
					// brackets [] for array to be added in the future
				case TK_LPAREN:
					p = new AstNode();
					p->op = OP_CALL;
					p->subs[0] = expr;
					getNextToken(parser);
					if (parser->currentToken != TK_RPAREN) {
						p->subs[1] = parseAssignmentExpr(parser);
						AstNode *tail = p->subs[1];
						while (parser->currentToken == TK_COMMA) {
							getNextToken(parser);
							tail->next = parseAssignmentExpr(parser);
							tail = tail->next;
						}
					}
					expect(parser, TK_RPAREN);
					expr = p;
					break;
				case TK_DOT:
					p = new AstNode();
					p->op = OP_MEMBER;
					p->subs[0] = expr;
					getNextToken(parser);
					if (parser->currentToken != TK_ID) {
						parser->error("Expect identifier as member.");
					}else{
						p->value = parser->value;
						getNextToken(parser);
					}
					expr = p;
					break;
				default:
					return expr;
				}
			}
		}

		static AstNode* parseUnaryExpr(SasModelParser *parser,int prec) {
			AstNode *expr;
			int op, tkPrec;
			switch (parser->currentToken) {
			case TK_K_NOT:
			case TK_ADD:
			case TK_SUB:
				op = tokenOps[parser->currentToken - TK_ADD].uop;
				tkPrec = tokenPrec[op];
				if (prec <= tkPrec) { // prec should increase as parsing goes
					expr = new AstNode();
					expr->op = static_cast<Op>(op);
					getNextToken(parser);
					expr->subs[0] = parseUnaryExpr(parser, tkPrec+1);
				}
				else {// If prec gets lower, report ERROR
					parser->error("Precedence of operator \"%s\" is lower than the previous operators. This operator is ignored.",opNames[op]);
					expr= parseUnaryExpr(parser, prec);
				}
				break;
			default:
				// parse conditional expr, which includes unary/binary expr
				expr = parsePostfixExpr(parser);				
			}
			op = OP_NONE;
			if (parser->currentToken >= TK_ADD && parser->currentToken <= TK_K_NOT && (op = tokenOps[parser->currentToken - TK_ADD].bop) != OP_NONE) {
				return parseBinaryExpr(parser, expr, prec);
			}
			else {
				// not a binary operator, means that binary expression has ended.
				return expr;
			}
		}

		static AstNode* parseBinaryExpr(SasModelParser *parser,AstNode *curExpr, int prec) {
			AstNode *expr=curExpr;
			int tkPrec,op;
			while (parser->currentToken >= TK_ADD && parser->currentToken <= TK_K_NOT&&(op=tokenOps[parser->currentToken - TK_ADD].bop) != OP_NONE && (tkPrec = tokenPrec[op])>=prec) {
				expr = new AstNode();
				expr->op = static_cast<Op>(tokenOps[parser->currentToken - TK_ADD].bop);
				expr->subs[0] = curExpr;
				getNextToken(parser);
				expr->subs[1] = parseUnaryExpr(parser, tkPrec + 1);
				curExpr = expr;
			}
			return expr;
		}

		static AstNode* parseConditionalExpr(SasModelParser *parser) {
			return parseUnaryExpr(parser, tokenPrec[OP_OR]);
		}

		static AstNode* parseAssignmentExpr(SasModelParser *parser) {
			AstNode *expr = parseConditionalExpr(parser);
			if (parser->currentToken == TK_EQN || parser->currentToken == TK_ASSIGN) {
				AstNode *assnExpr = new AstNode();
				assnExpr->op = static_cast<Op>(tokenOps[parser->currentToken - TK_ADD].bop);
				assnExpr->subs[0] = expr;
				getNextToken(parser);
				assnExpr->subs[1] = parseConditionalExpr(parser);

				return assnExpr;
			}
			return expr;
		}

		AstNode* parseExpr(SasModelParser *parser) {
			AstNode *expr, *commaExpr;
			expr = parseAssignmentExpr(parser);
			while (parser->currentToken == TK_COMMA) {
				commaExpr = new AstNode();
				commaExpr->op = OP_COMMA;
				commaExpr->subs[0] = expr;
				getNextToken(parser);
				commaExpr->subs[1] = parseAssignmentExpr(parser);
				expr = commaExpr;
			}
			return expr;
		}
	}
}