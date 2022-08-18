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
#include "SasModel.h"
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdarg>
#include <iomanip>

using namespace std;

namespace che {
	namespace core {
		IdPearl::IdPearl(const unsigned char* s, int sz, const unsigned char* ann, const unsigned char* un, Type ty, IdType idTy) {
			setName(s, sz);
			setAnnotation(ann);
			setUnit(un);

			this->ty = ty;
			this->idTy = idTy;
			this->value.d = NAN;
			this->derVal.d = NAN;
		}

		IdPearl::IdPearl(const unsigned char* s, int sz) {
			setName(s, sz);
		}

		IdPearl::IdPearl(const IdPearl& idp) {
			setName(idp.str, idp.size);
			setAnnotation(idp.annotation);
			setUnit(idp.unit);
			this->ty = idp.ty;
			this->idTy = idp.idTy;
			this->value = idp.value;
			this->initVal = idp.initVal;
			this->initDer = idp.initDer;
		}

		void IdPearl::setName(const unsigned char* s, int sz) {
			resetName();
			if (sz > 0) {
				size = sz;
				str = new unsigned char[sz + 1];
				for (int i = 0; i < sz; i++) {
					str[i] = s[i];
				}
				str[sz] = 0;
			}
		}

		void IdPearl::setAnnotation(const unsigned char* ann) {
			resetAnnotation();
			if (ann != NULL) {
				int annLen = strlen((const char*)ann);

				annotation = new unsigned char[annLen + 1];
				for (int i = 0; i < annLen; i++) {
					annotation[i] = ann[i];
				}
				annotation[annLen] = 0;
			}
		}

		void IdPearl::setUnit(const unsigned char* un) {
			resetUnit();
			if (un != NULL) {
				int unLen = strlen((const char*)un);

				unit = new unsigned char[unLen + 1];
				for (int i = 0; i < unLen; i++) {
					unit[i] = un[i];
				}
				unit[unLen] = 0;
			}
		}

		IdPearl::~IdPearl() {
			resetName();
			resetAnnotation();
			resetUnit();
			pair = NULL;
		}

		void IdPearl::resetName() {
			if (size > 0) {
				delete[] str;
				str = NULL;
				size = 0;
			}
		}

		void IdPearl::resetAnnotation() {
			if (annotation != NULL) {
				delete[] annotation;
				annotation = NULL;
			}
		}

		void IdPearl::resetUnit() {
			if (unit != NULL) {
				delete[] unit;
				unit = NULL;
			}
		}

		SasModel::SasModel() :modelTrees(), idTable(), idNecklaces(MAX_ID_PEARL), subModels() {
			errCount = 0;
			warnCount = 0;
		}

		SasModel::~SasModel() {
			if (modelInfo != nullptr) {
				modelInfo.reset();
				modelInfo = nullptr;
			}
			//Need to destroy tables and trees.
			for (int i = 0; i < idNecklaces.size(); i++) {
				idNecklaces[i].clear();
			}
			idNecklaces.clear();
			idTable.clear();
			modelTrees.clear();
			subModels.clear();
		}

		void SasModel::printAllModelTrees() {
			list<shared_ptr<AstTree>>::iterator it;
			for (it = modelTrees.begin(); it != modelTrees.end(); it++) {
				(*it)->printTree(cout);
			}
		}

		//ID_NONE, ID_PARAM, ID_VAR, ID_IVAR, ID_TMPVAR, ID_MODEL, ID_FUNC
		static const char* idTypeName[] = {
			"None","Param","Var","iVar","tmpVar","Model","Function",
		};

		void SasModel::printIdTable() {
			cout << "+-"<<string(MAX_ID_LEN, '-') << "-+-" 
				<< string(15, '-') << "-+-" 
				<< string(15, '-') << "-+-" << string(15, '-') << "-+-" << string(5, '-') << "-+-"
				<< string(MAX_ID_LEN, '-') << "-+-" 
				<< string(MAX_STR_LEN, '-') << "-+" << endl;
			cout << "| " << left << setw(MAX_ID_LEN) << "Name" << " | " <<
				left << setw(15) << "Type" << " | " <<
				left << setw(15) << "Init Value" << " | " << left << setw(15) << "Init Der" << " | " <<
				left << setw(5) << "Index" << " | " <<
				left << setw(MAX_ID_LEN) << "Unit" << " | " <<
				left << setw(MAX_STR_LEN) << "Annotation" << " |" << endl;
			cout << "+-" << string(MAX_ID_LEN, '-') << "-+-"
				<< string(15, '-') << "-+-"
				<< string(15, '-') << "-+-" << string(15, '-') << "-+-" << string(5, '-') << "-+-"
				<< string(MAX_ID_LEN, '-') << "-+-"
				<< string(MAX_STR_LEN, '-') << "-+" << endl;
			list<shared_ptr<IdPearl>>::iterator it;
			for (it = idTable.begin(); it != idTable.end(); it++) {
				cout << "| " << left << setw(MAX_ID_LEN) << ((*it)->str == NULL ? "" : (char*)(*it)->str) << " | " <<
					left << setw(15) << idTypeName[(*it)->idTy] << " | ";
				if ((*it)->idTy == ID_MODEL || (*it)->idTy == ID_FUNC || (*it)->idTy == ID_NONE) {
					cout << left << setw(15) << "" << " | ";
				}
				else if((*it)->initVal){
					cout << left << setw(15) << ((*it)->ty == TY_INT ? (*it)->value.i : (*it)->value.d) << " | ";
				}
				else {
					cout << left << setw(15) << "uninitialized" << " | ";
				}

				if ((*it)->idTy == ID_MODEL || (*it)->idTy == ID_FUNC || (*it)->idTy == ID_NONE) {
					cout << left << setw(15) << "" << " | ";
				}
				else if ((*it)->initDer) {
					cout << left << setw(15) << ((*it)->ty == TY_INT ? (*it)->derVal.i : (*it)->derVal.d) << " | ";
				}
				else {
					cout << left << setw(15) << "uninitialized" << " | ";
				}

				if ((*it)->idTy == ID_VAR || (*it)->idTy == ID_IVAR || (*it)->idTy == ID_TMPVAR ) {
					cout << left << setw(5) << (*it)->index << " | ";
				}
				else {
					cout << left << setw(5) << "" << " | ";
				}
				
				cout << left << setw(MAX_ID_LEN) << ((*it)->unit == NULL ? "" : (char*)(*it)->unit) << " | " <<
					left << setw(MAX_STR_LEN) << ((*it)->annotation == NULL ? "" : (char*)(*it)->annotation) << " |" << endl;
			}
			cout << "+-" << string(MAX_ID_LEN, '-') << "-+-"
				<< string(15, '-') << "-+-"
				<< string(15, '-') << "-+-" << string(15, '-') << "-+-" << string(5, '-') << "-+-"
				<< string(MAX_ID_LEN, '-') << "-+-"
				<< string(MAX_STR_LEN, '-') << "-+" << endl;
		}

		void SasModel::error(const char* formatMsg, ...) {
			va_list ap;
			char buffer[MAX_PRINT_LEN + 1];
			va_start(ap, formatMsg);
			vsprintf(buffer, formatMsg, ap);
			cerr << "(Model "<<modelInfo->str<<") error: " << buffer << endl;
			errCount++;
			va_end(ap);
		}

		void SasModel::warn(const char* formatMsg, ...) {
			va_list ap;
			char buffer[MAX_PRINT_LEN + 1];
			va_start(ap, formatMsg);
			vsprintf(buffer, formatMsg, ap);
			cerr << "(Model " << modelInfo->str << ") warn: " << buffer << endl;
			warnCount++;
			va_end(ap);
		}

		shared_ptr<IdPearl> SasModel::lookupIdFromModel(const char* iden, int len) {
			int h = SasCompUtil::ELFHashStr(iden, len)&MAX_ID_PEARL;
			list<shared_ptr<IdPearl>>::iterator it;
			for (it = this->idNecklaces[h].begin(); it != this->idNecklaces[h].end(); it++) {
				if ((*it)->size == len && strncmp(iden, (char*)(*it)->str, len) == 0) {
					return (*it);
				}
			}
			return nullptr;
		}

		shared_ptr<IdPearl> SasModel::injectNewId(const char* iden, int len) {
			int h = SasCompUtil::ELFHashStr(iden, len)&MAX_ID_PEARL;
			shared_ptr<IdPearl> pId = make_shared<IdPearl>((const unsigned char*)iden, len);
			this->idTable.push_back(pId);
			this->idNecklaces[h].push_back(pId);
			return pId;
		}

		shared_ptr<IdPearl> SasModel::lookupIdentifier(const char* iden, int len) {
			int h = SasCompUtil::ELFHashStr(iden, len)&MAX_ID_PEARL;
			shared_ptr<IdPearl> pId = lookupIdFromModel(iden, len);
			if (pId != nullptr) {
				return pId;
			}
			return injectNewId(iden, len);
		}


		GlobalPool& GlobalPool::getInstance() {
			static GlobalPool instance;
			return instance;
		}

		shared_ptr<SasModel> GlobalPool::createNewTopModel() {
			shared_ptr<SasModel> pModel = make_shared<SasModel>();
			topModels.push_back(pModel);
			return pModel;
		}

		GlobalPool::~GlobalPool() {
			strTable.clear();
			topModels.clear();
		}

		GlobalPool::GlobalPool() :strTable(), topModels() {}

		shared_ptr<IdPearl> lookupIdFromModel(shared_ptr<SasModel>& model, const char* iden, int len);

		shared_ptr<SasModel> GlobalPool::findModel(list<shared_ptr<SasModel>>& modelTree, AstNode* p) {
			return nullptr;
		}

	}
}