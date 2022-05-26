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
#ifndef TOKEN
#error "YOU MUST DEFINE MACRO M_TOKEN FIRST"
#endif

TOKEN(TK_KW_function, "function")
TOKEN(TK_KW_end, "end")
TOKEN(TK_KW_RETURN, "return")
TOKEN(TK_SQBRKT_L, "[")
TOKEN(TK_SQBRKT_R, "]")
TOKEN(TK_BRKT_L, "(")
TOKEN(TK_BRKT_R, ")")
TOKEN(TK_CURBRKT_L, "{")
TOKEN(TK_CURBRKT_R, "}")
TOKEN(TK_TRIDOT, "...")
TOKEN(TK_COMMA, ",")
TOKEN(TK_SEMCOLON, ";")
TOKEN(TK_BKSLASH, "\\")
TOKEN(TK_EQUAL, "=")
TOKEN(TK_DOT, ".")
TOKEN(TK_SPACE, " ")
TOKEN(TK_PERCENT, "%")
TOKEN(TK_EMPTY, "")
TOKEN(TK_NEWLINE, "\n")
TOKEN(TK_CR, "\r")
TOKEN(TK_END, "EOF")