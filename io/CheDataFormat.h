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
#pragma once
#ifndef _Che_CheDataFormat_H_
#define _Che_CheDataFormat_H_

#define C_IDX(x) ((x)-1)
#define M_IDX(x) ((x)+1)

#include <string>
#include "util/SafeArmadillo.h"
#include <map>

using namespace std;
using namespace arma;

//This data format follows PSAT (Matlab toolbox) data format
//And defined as the standard data format for Che.

namespace che {
	namespace io {
		namespace chedata {}
		namespace psat = chedata;

		namespace chedata {
			class CheComponent {
			public:
				static int counter;

				CheComponent() {
					regNewId();
				}

				CheComponent(const CheComponent &c, bool reg = false) {
					id = c.id;
					regNewId(reg);
				}

				int getId() const {
					return id;
				}

				void regNewId(bool reg = true) {
					if (reg) id = counter++;
				}
			private:
				int id;
			};

			class Bus :public CheComponent {
			public:
				int busNumber;
				double baseV;
				double vMag;
				double vAng;
				int areaNumber;
				int regionNumber;
				string busName;

				Bus() {}

				Bus(const Bus& bus, bool regNew = false) :CheComponent(bus, regNew) {
					this->busNumber = bus.busNumber;
					this->baseV = bus.baseV;
					this->vMag = bus.vMag;
					this->vAng = bus.vAng;
					this->areaNumber = bus.areaNumber;
					this->regionNumber = bus.regionNumber;
					this->busName = string(bus.busName);
				}
			};

			class SW :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double vMag;
				double vAng;
				double qMax;
				double qMin;
				double vMax;
				double vMin;
				double pg0;
				double lossParticipation;
				unsigned char isRefBus;
				unsigned char status;

				SW() {}

				SW(const SW &sw, bool regNew = false) :CheComponent(sw, regNew) {
					this->busNumber = sw.busNumber;
					this->baseS = sw.baseS;
					this->baseV = sw.baseV;
					this->vMag = sw.vMag;
					this->vAng = sw.vAng;
					this->qMax = sw.qMax;
					this->qMin = sw.qMin;
					this->vMax = sw.vMax;
					this->vMin = sw.vMin;
					this->pg0 = sw.pg0;
					this->lossParticipation = sw.lossParticipation;
					this->isRefBus = sw.isRefBus;
					this->status = sw.status;
				}
			};

			class PV :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double P;
				double vMag;
				double qMax;
				double qMin;
				double vMax;
				double vMin;
				double lossParticipation;
				unsigned char status;

				PV() {}

				PV(const PV &pv, bool regNew = false) :CheComponent(pv, regNew) {
					this->busNumber = pv.busNumber;
					this->baseS = pv.baseS;
					this->baseV = pv.baseV;
					this->P = pv.P;
					this->vMag = pv.vMag;
					this->qMax = pv.qMax;
					this->qMin = pv.qMin;
					this->vMax = pv.vMax;
					this->vMin = pv.vMin;
					this->lossParticipation = pv.lossParticipation;
					this->status = pv.status;
				}
			};

			class PQ :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double P;
				double Q;
				double vMax;
				double vMin;
				unsigned char allowConversion;
				unsigned char status;

				PQ() {}

				PQ(const PQ &pq, bool regNew = false) :CheComponent(pq, regNew) {
					this->busNumber = pq.busNumber;
					this->baseS = pq.baseS;
					this->baseV = pq.baseV;
					this->P = pq.P;
					this->Q = pq.Q;
					this->vMax = pq.vMax;
					this->vMin = pq.vMin;
					this->allowConversion = pq.allowConversion;
					this->status = pq.status;
				}
			};

			class Shunt :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double baseF;
				double g;
				double b;
				unsigned char status;

				Shunt() {}

				Shunt(const Shunt &shunt, bool regNew = false) :CheComponent(shunt, regNew) {
					this->busNumber = shunt.busNumber;
					this->baseS = shunt.baseS;
					this->baseV = shunt.baseV;
					this->baseF = shunt.baseF;
					this->g = shunt.g;
					this->b = shunt.b;
					this->status = shunt.status;
				}
			};

			class Line :public CheComponent {
			public:
				int fromBus;
				int toBus;
				double baseS;
				double baseV;
				double baseF;
				double len;
				double kT;
				double r;
				double x;
				double b;
				double k;
				double ang;
				double iMax;
				double pMax;
				double sMax;
				unsigned char status;

				Line() {}

				Line(const Line &line, bool regNew = false) :CheComponent(line, regNew) {
					this->fromBus = line.fromBus;
					this->toBus = line.toBus;
					this->baseS = line.baseS;
					this->baseV = line.baseV;
					this->baseF = line.baseF;
					this->len = line.len;
					this->kT = line.kT;
					this->r = line.r;
					this->x = line.x;
					this->b = line.b;
					this->k = line.k;
					this->ang = line.ang;
					this->iMax = line.iMax;
					this->pMax = line.pMax;
					this->sMax = line.sMax;
					this->status = line.status;
				}
			};

			class Pl :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double baseF;
				double g;
				double Ip;
				double P;
				double b;
				double Iq;
				double Q;
				unsigned char initAfterPF;
				unsigned char status;

				Pl() {}

				Pl(const Pl &pl, bool regNew = false) :CheComponent(pl, regNew) {
					this->busNumber = pl.busNumber;
					this->baseS = pl.baseS;
					this->baseV = pl.baseV;
					this->baseF = pl.baseF;
					this->g = pl.g;
					this->Ip = pl.Ip;
					this->P = pl.P;
					this->b = pl.b;
					this->Iq = pl.Iq;
					this->Q = pl.Q;
					this->initAfterPF = pl.initAfterPF;
					this->status = pl.status;
				}
			};

			class Syn :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double baseF;
				double model;
				double xl;
				double ra;
				double xd;
				double xd1;
				double xd2;
				double Td01;
				double Td02;
				double xq;
				double xq1;
				double xq2;
				double Tq01;
				double Tq02;
				double M;
				double D;
				double Ko;
				double Kp;
				double gammaP;
				double gammaQ;
				double TAA;
				double sat1;
				double sat2;
				int nCOI;
				unsigned char status;

				Syn() {}

				Syn(const Syn &syn, bool regNew = false) :CheComponent(syn, regNew) {
					this->busNumber = syn.busNumber;
					this->baseS = syn.baseS;
					this->baseV = syn.baseV;
					this->baseF = syn.baseF;
					this->model = syn.model;
					this->xl = syn.xl;
					this->ra = syn.ra;
					this->xd = syn.xd;
					this->xd1 = syn.xd1;
					this->xd2 = syn.xd2;
					this->Td01 = syn.Td01;
					this->Td02 = syn.Td02;
					this->xq = syn.xq;
					this->xq1 = syn.xq1;
					this->xq2 = syn.xq2;
					this->Tq01 = syn.Tq01;
					this->Tq02 = syn.Tq02;
					this->M = syn.M;
					this->D = syn.D;
					this->Ko = syn.Ko;
					this->Kp = syn.Kp;
					this->gammaP = syn.gammaP;
					this->gammaQ = syn.gammaQ;
					this->TAA = syn.TAA;
					this->sat1 = syn.sat1;
					this->sat2 = syn.sat2;
					this->nCOI = syn.nCOI;
					this->status = syn.status;
				}
			};

			class Ind :public CheComponent {
			public:
				int busNumber;
				double baseS;
				double baseV;
				double baseF;
				double model;
				unsigned char startCtrl;
				double rs;
				double xs;
				double rr1;
				double xr1;
				double rr2;
				double xr2;
				double xm;
				double Hm;
				double Ta;
				double Tb;
				double Tc;
				double tup;
				unsigned char allowBrake;
				unsigned char status;

				Ind() {}

				Ind(const Ind &ind, bool regNew = false) :CheComponent(ind, regNew) {
					this->busNumber = ind.busNumber;
					this->baseS = ind.baseS;
					this->baseV = ind.baseV;
					this->baseF = ind.baseF;
					this->model = ind.model;
					this->startCtrl = ind.startCtrl;
					this->rs = ind.rs;
					this->xs = ind.xs;
					this->rr1 = ind.rr1;
					this->xr1 = ind.xr1;
					this->rr2 = ind.rr2;
					this->xr2 = ind.xr2;
					this->xm = ind.xm;
					this->Hm = ind.Hm;
					this->Ta = ind.Ta;
					this->Tb = ind.Tb;
					this->Tc = ind.Tc;
					this->tup = ind.tup;
					this->allowBrake = ind.allowBrake;
					this->status = ind.status;
				}
			};

			class Tg :public CheComponent {
			public:
				struct Tg1Data {
					double wref0;
					double R;
					double Tmax;
					double Tmin;
					double Ts;
					double Tc;
					double T3;
					double T4;
					double T5;
				};

				struct Tg2Data {
					double wref0;
					double R;
					double Tmax;
					double Tmin;
					double T2;
					double T1;
				};

				int synNumber;
				unsigned char tgType;
				union {
					Tg1Data tg1;
					Tg2Data tg2;
				} tgData;
				unsigned char status;

				Tg() {}

				Tg(const Tg &tg, bool regNew = false) :CheComponent(tg, regNew) {
					this->synNumber = tg.synNumber;
					this->tgType = tg.tgType;
					if (tg.tgType == 1) {
						this->tgData.tg1.wref0 = tg.tgData.tg1.wref0;
						this->tgData.tg1.R = tg.tgData.tg1.R;
						this->tgData.tg1.Tmax = tg.tgData.tg1.Tmax;
						this->tgData.tg1.Tmin = tg.tgData.tg1.Tmin;
						this->tgData.tg1.Ts = tg.tgData.tg1.Ts;
						this->tgData.tg1.Tc = tg.tgData.tg1.Tc;
						this->tgData.tg1.T3 = tg.tgData.tg1.T3;
						this->tgData.tg1.T4 = tg.tgData.tg1.T4;
						this->tgData.tg1.T5 = tg.tgData.tg1.T5;
					}
					else if (tg.tgType == 2) {
						this->tgData.tg2.wref0 = tg.tgData.tg2.wref0;
						this->tgData.tg2.R = tg.tgData.tg2.R;
						this->tgData.tg2.Tmax = tg.tgData.tg2.Tmax;
						this->tgData.tg2.Tmin = tg.tgData.tg2.Tmin;
						this->tgData.tg2.T2 = tg.tgData.tg2.T2;
						this->tgData.tg2.T1 = tg.tgData.tg2.T1;
					}
					this->status = tg.status;
				}
			};

			class Exc :public CheComponent {
			public:
				struct Exc1Data {
					double vMax;
					double vMin;
					double mu0;
					double T1;
					double T2;
					double T3;
					double T4;
					double Te;
					double Tr;
					double Ae;
					double Be;
				};
				struct Exc2Data {
					double vMax;
					double vMin;
					double Ka;
					double Ta;
					double Kf;
					double Tf;
					double Te;
					double Tr;
					double Ae;
					double Be;
				};
				struct Exc3Data {
					double vMax;
					double vMin;
					double mu0;
					double T2;
					double T1;
					double vf0;
					double V0;
					double Te;
					double Tr;
				};

				int synNumber;
				unsigned char excType;
				union {
					Exc1Data exc1;
					Exc2Data exc2;
					Exc3Data exc3;
				} excData;
				unsigned char status;

				Exc() {}

				Exc(const Exc &exc, bool regNew = false) :CheComponent(exc, regNew) {
					this->synNumber = exc.synNumber;
					this->excType = exc.excType;
					this->status = exc.status;
					if (exc.excType == 1) {
						this->excData.exc1.vMax = exc.excData.exc1.vMax;
						this->excData.exc1.vMin = exc.excData.exc1.vMin;
						this->excData.exc1.mu0 = exc.excData.exc1.mu0;
						this->excData.exc1.T1 = exc.excData.exc1.T1;
						this->excData.exc1.T2 = exc.excData.exc1.T2;
						this->excData.exc1.T3 = exc.excData.exc1.T3;
						this->excData.exc1.T4 = exc.excData.exc1.T4;
						this->excData.exc1.Te = exc.excData.exc1.Te;
						this->excData.exc1.Tr = exc.excData.exc1.Tr;
						this->excData.exc1.Ae = exc.excData.exc1.Ae;
						this->excData.exc1.Be = exc.excData.exc1.Be;
					}
					else if (exc.excType == 2) {
						this->excData.exc2.vMax = exc.excData.exc2.vMax;
						this->excData.exc2.vMin = exc.excData.exc2.vMin;
						this->excData.exc2.Ka = exc.excData.exc2.Ka;
						this->excData.exc2.Ta = exc.excData.exc2.Ta;
						this->excData.exc2.Kf = exc.excData.exc2.Kf;
						this->excData.exc2.Tf = exc.excData.exc2.Tf;
						this->excData.exc2.Te = exc.excData.exc2.Te;
						this->excData.exc2.Tr = exc.excData.exc2.Tr;
						this->excData.exc2.Ae = exc.excData.exc2.Ae;
						this->excData.exc2.Be = exc.excData.exc2.Be;
					}
					else if (exc.excType == 3) {
						this->excData.exc3.vMax = exc.excData.exc3.vMax;
						this->excData.exc3.vMin = exc.excData.exc3.vMin;
						this->excData.exc3.mu0 = exc.excData.exc3.mu0;
						this->excData.exc3.T2 = exc.excData.exc3.T2;
						this->excData.exc3.T1 = exc.excData.exc3.T1;
						this->excData.exc3.vf0 = exc.excData.exc3.vf0;
						this->excData.exc3.V0 = exc.excData.exc3.V0;
						this->excData.exc3.Te = exc.excData.exc3.Te;
						this->excData.exc3.Tr = exc.excData.exc3.Tr;
					}
				}
			};

			// The definition of the entire dataset
			class PsatDataSet {
			public:
				int nBus;
				int nSw;
				int nPv;
				int nPq;
				int nShunt;
				int nLine;
				int nPl;
				int nSyn;
				int nInd;
				int nTg;
				int nExc;

				Bus *buses;
				SW *sws;
				PV *pvs;
				PQ *pqs;
				Shunt *shunts;
				Line *lines;
				Pl *pls;
				Syn *syns;
				Ind *inds;
				Tg *tgs;
				Exc *excs;

				map<int, int> newToOld;
				map<int, int> oldToNew;

				bool isFormatted;

			public:
				PsatDataSet() {

					nBus = -1;
					nSw = -1;
					nPv = -1;
					nPq = -1;
					nShunt = -1;
					nLine = -1;
					nPl = -1;
					nSyn = -1;
					nInd = -1;
					nTg = -1;
					nExc = -1;

					newToOld = map<int, int>();
					oldToNew = map<int, int>();

					isFormatted = false;
				}

				PsatDataSet(const PsatDataSet &psat, bool regNew = false) {
					nBus = -1;
					nSw = -1;
					nPv = -1;
					nPq = -1;
					nShunt = -1;
					nLine = -1;
					nPl = -1;
					nSyn = -1;
					nInd = -1;
					nTg = -1;
					nExc = -1;

					newToOld = map<int, int>();
					oldToNew = map<int, int>();

					isFormatted = false;
					copyMembers(psat, regNew);
				}

				PsatDataSet& operator=(const PsatDataSet &psat) {
					nBus = -1;
					nSw = -1;
					nPv = -1;
					nPq = -1;
					nShunt = -1;
					nLine = -1;
					nPl = -1;
					nSyn = -1;
					nInd = -1;
					nTg = -1;
					nExc = -1;

					newToOld = map<int, int>();
					oldToNew = map<int, int>();

					isFormatted = false;
					copyMembers(psat, false);
					return *this;
				}

				void renumberBuses(map<int, int>& oldToNew) {
					map<int, int> busMap;
					for (int i = 0; i < nBus; i++) {
						busMap.insert(pair<int, int>(buses[i].busNumber, i + 1));
					}
					for (int i = 0; i < nSw; i++) {
						sws[i].busNumber = busMap[sws[i].busNumber];
					}
					for (int i = 0; i < nPv; i++) {
						pvs[i].busNumber = busMap[pvs[i].busNumber];
					}
					for (int i = 0; i < nPq; i++) {
						pqs[i].busNumber = busMap[pqs[i].busNumber];
					}
					for (int i = 0; i < nShunt; i++) {
						shunts[i].busNumber = busMap[shunts[i].busNumber];
					}
					for (int i = 0; i < nLine; i++) {
						lines[i].fromBus = busMap[lines[i].fromBus];
						lines[i].toBus = busMap[lines[i].toBus];
					}
					for (int i = 0; i < nPl; i++) {
						pls[i].busNumber = busMap[pls[i].busNumber];
					}
					for (int i = 0; i < nSyn; i++) {
						syns[i].busNumber = busMap[syns[i].busNumber];
					}
					for (int i = 0; i < nInd; i++) {
						inds[i].busNumber = busMap[inds[i].busNumber];
					}
					if (oldToNew.empty()) {
						this->oldToNew = busMap;
					}
					else {
						map<int, int>tempMap;
						map<int, int>::iterator itr;
						for (itr = oldToNew.begin(); itr != oldToNew.end(); ++itr) {
							tempMap.insert(pair<int, int>(itr->first, busMap[itr->second]));
						}
						this->oldToNew.clear();
						this->oldToNew = tempMap;
					}
					newToOld.clear();
					map<int, int>::iterator itr;
					for (itr = this->oldToNew.begin(); itr != this->oldToNew.end(); ++itr) {
						newToOld.insert(pair<int, int>(itr->second, itr->first));
					}
					isFormatted = true;
				}

				void renumberBuses() {
					renumberBuses(this->oldToNew);
				}

				virtual ~PsatDataSet() {
					freeSpace();
				}

				void reset() {
					freeSpace();

					nBus = -1;
					nSw = -1;
					nPv = -1;
					nPq = -1;
					nShunt = -1;
					nLine = -1;
					nPl = -1;
					nSyn = -1;
					nInd = -1;
					nTg = -1;
					nExc = -1;

					newToOld = map<int, int>();
					oldToNew = map<int, int>();

					isFormatted = false;
				}

#define GET_FIELD_VEC(type,data,field,size) type get_##data##_##field##_vec() const{\
type x(size>0?size:0);\
for(int i=0;i<size;i++)x(i)=data[i].field;\
return x;\
}

#define GET_UNION_SUBFIELD_VEC(type,data,un,field,subf,size) type get_##data##_##un##_##field##_##subf##_vec() const{\
type x(size>0?size:0);\
for(int i=0;i<size;i++)x(i)=data[i].un.field.subf;\
return x;\
}
				// Bus
				GET_FIELD_VEC(uvec, buses, busNumber, nBus);
				GET_FIELD_VEC(vec, buses, baseV, nBus);
				GET_FIELD_VEC(vec, buses, vMag, nBus);
				GET_FIELD_VEC(vec, buses, vAng, nBus);
				GET_FIELD_VEC(uvec, buses, areaNumber, nBus);
				GET_FIELD_VEC(uvec, buses, regionNumber, nBus);
				// SW
				GET_FIELD_VEC(uvec, sws, busNumber, nSw);
				GET_FIELD_VEC(vec, sws, baseS, nSw);
				GET_FIELD_VEC(vec, sws, baseV, nSw);
				GET_FIELD_VEC(vec, sws, vMag, nSw);
				GET_FIELD_VEC(vec, sws, vAng, nSw);
				GET_FIELD_VEC(vec, sws, qMax, nSw);
				GET_FIELD_VEC(vec, sws, qMin, nSw);
				GET_FIELD_VEC(vec, sws, vMax, nSw);
				GET_FIELD_VEC(vec, sws, vMin, nSw);
				GET_FIELD_VEC(vec, sws, pg0, nSw);
				GET_FIELD_VEC(vec, sws, lossParticipation, nSw);
				GET_FIELD_VEC(uvec, sws, isRefBus, nSw);
				GET_FIELD_VEC(uvec, sws, status, nSw);
				//PV
				GET_FIELD_VEC(uvec, pvs, busNumber, nPv);
				GET_FIELD_VEC(vec, pvs, baseS, nPv);
				GET_FIELD_VEC(vec, pvs, baseV, nPv);
				GET_FIELD_VEC(vec, pvs, P, nPv);
				GET_FIELD_VEC(vec, pvs, vMag, nPv);
				GET_FIELD_VEC(vec, pvs, qMax, nPv);
				GET_FIELD_VEC(vec, pvs, qMin, nPv);
				GET_FIELD_VEC(vec, pvs, vMax, nPv);
				GET_FIELD_VEC(vec, pvs, vMin, nPv);
				GET_FIELD_VEC(vec, pvs, lossParticipation, nPv);
				GET_FIELD_VEC(uvec, pvs, status, nPv);
				//PQ
				GET_FIELD_VEC(uvec, pqs, busNumber, nPq);
				GET_FIELD_VEC(vec, pqs, baseS, nPq);
				GET_FIELD_VEC(vec, pqs, baseV, nPq);
				GET_FIELD_VEC(vec, pqs, P, nPq);
				GET_FIELD_VEC(vec, pqs, Q, nPq);
				GET_FIELD_VEC(vec, pqs, vMax, nPq);
				GET_FIELD_VEC(vec, pqs, vMin, nPq);
				GET_FIELD_VEC(uvec, pqs, allowConversion, nPq);
				GET_FIELD_VEC(uvec, pqs, status, nPq);
				//Shunt
				GET_FIELD_VEC(uvec, shunts, busNumber, nShunt);
				GET_FIELD_VEC(vec, shunts, baseS, nShunt);
				GET_FIELD_VEC(vec, shunts, baseV, nShunt);
				GET_FIELD_VEC(vec, shunts, baseF, nShunt);
				GET_FIELD_VEC(vec, shunts, g, nShunt);
				GET_FIELD_VEC(vec, shunts, b, nShunt);
				GET_FIELD_VEC(uvec, shunts, status, nShunt);
				//Line
				GET_FIELD_VEC(uvec, lines, fromBus, nLine);
				GET_FIELD_VEC(uvec, lines, toBus, nLine);
				GET_FIELD_VEC(vec, lines, baseS, nLine);
				GET_FIELD_VEC(vec, lines, baseV, nLine);
				GET_FIELD_VEC(vec, lines, baseF, nLine);
				GET_FIELD_VEC(vec, lines, len, nLine);
				GET_FIELD_VEC(vec, lines, kT, nLine);
				GET_FIELD_VEC(vec, lines, r, nLine);
				GET_FIELD_VEC(vec, lines, x, nLine);
				GET_FIELD_VEC(vec, lines, b, nLine);
				GET_FIELD_VEC(vec, lines, k, nLine);
				GET_FIELD_VEC(vec, lines, ang, nLine);
				GET_FIELD_VEC(vec, lines, iMax, nLine);
				GET_FIELD_VEC(vec, lines, pMax, nLine);
				GET_FIELD_VEC(vec, lines, sMax, nLine);
				GET_FIELD_VEC(uvec, lines, status, nLine);
				//Pl
				GET_FIELD_VEC(uvec, pls, busNumber, nPl);
				GET_FIELD_VEC(vec, pls, baseS, nPl);
				GET_FIELD_VEC(vec, pls, baseV, nPl);
				GET_FIELD_VEC(vec, pls, baseF, nPl);
				GET_FIELD_VEC(vec, pls, g, nPl);
				GET_FIELD_VEC(vec, pls, Ip, nPl);
				GET_FIELD_VEC(vec, pls, P, nPl);
				GET_FIELD_VEC(vec, pls, b, nPl);
				GET_FIELD_VEC(vec, pls, Iq, nPl);
				GET_FIELD_VEC(vec, pls, Q, nPl);
				GET_FIELD_VEC(uvec, pls, initAfterPF, nPl);
				GET_FIELD_VEC(uvec, pls, status, nPl);
				//Syn
				GET_FIELD_VEC(uvec, syns, busNumber, nSyn);
				GET_FIELD_VEC(vec, syns, baseS, nSyn);
				GET_FIELD_VEC(vec, syns, baseV, nSyn);
				GET_FIELD_VEC(vec, syns, baseF, nSyn);
				GET_FIELD_VEC(vec, syns, model, nSyn);
				GET_FIELD_VEC(vec, syns, xl, nSyn);
				GET_FIELD_VEC(vec, syns, ra, nSyn);
				GET_FIELD_VEC(vec, syns, xd, nSyn);
				GET_FIELD_VEC(vec, syns, xd1, nSyn);
				GET_FIELD_VEC(vec, syns, xd2, nSyn);
				GET_FIELD_VEC(vec, syns, Td01, nSyn);
				GET_FIELD_VEC(vec, syns, Td02, nSyn);
				GET_FIELD_VEC(vec, syns, xq, nSyn);
				GET_FIELD_VEC(vec, syns, xq1, nSyn);
				GET_FIELD_VEC(vec, syns, xq2, nSyn);
				GET_FIELD_VEC(vec, syns, Tq01, nSyn);
				GET_FIELD_VEC(vec, syns, Tq02, nSyn);
				GET_FIELD_VEC(vec, syns, M, nSyn);
				GET_FIELD_VEC(vec, syns, D, nSyn);
				GET_FIELD_VEC(vec, syns, Ko, nSyn);
				GET_FIELD_VEC(vec, syns, Kp, nSyn);
				GET_FIELD_VEC(vec, syns, gammaP, nSyn);
				GET_FIELD_VEC(vec, syns, gammaQ, nSyn);
				GET_FIELD_VEC(vec, syns, TAA, nSyn);
				GET_FIELD_VEC(vec, syns, sat1, nSyn);
				GET_FIELD_VEC(vec, syns, sat2, nSyn);
				GET_FIELD_VEC(uvec, syns, nCOI, nSyn);
				GET_FIELD_VEC(uvec, syns, status, nSyn);
				//Ind
				GET_FIELD_VEC(uvec, inds, busNumber, nInd);
				GET_FIELD_VEC(vec, inds, baseS, nInd);
				GET_FIELD_VEC(vec, inds, baseV, nInd);
				GET_FIELD_VEC(vec, inds, baseF, nInd);
				GET_FIELD_VEC(vec, inds, model, nInd);
				GET_FIELD_VEC(uvec, inds, startCtrl, nInd);
				GET_FIELD_VEC(vec, inds, rs, nInd);
				GET_FIELD_VEC(vec, inds, xs, nInd);
				GET_FIELD_VEC(vec, inds, rr1, nInd);
				GET_FIELD_VEC(vec, inds, xr1, nInd);
				GET_FIELD_VEC(vec, inds, rr2, nInd);
				GET_FIELD_VEC(vec, inds, xr2, nInd);
				GET_FIELD_VEC(vec, inds, xm, nInd);
				GET_FIELD_VEC(vec, inds, Hm, nInd);
				GET_FIELD_VEC(vec, inds, Ta, nInd);
				GET_FIELD_VEC(vec, inds, Tb, nInd);
				GET_FIELD_VEC(vec, inds, Tc, nInd);
				GET_FIELD_VEC(vec, inds, tup, nInd);
				GET_FIELD_VEC(uvec, inds, allowBrake, nInd);
				GET_FIELD_VEC(uvec, inds, status, nInd);
				//Tg
				GET_FIELD_VEC(uvec, tgs, synNumber, nTg);
				GET_FIELD_VEC(uvec, tgs, tgType, nTg);
				GET_FIELD_VEC(uvec, tgs, status, nTg);
				//Tg - Tg1
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, wref0, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, R, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, Tmax, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, Tmin, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, Ts, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, Tc, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, T3, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, T4, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg1, T5, nTg);
				// Tg - Tg2
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg2, wref0, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg2, R, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg2, Tmax, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg2, Tmin, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg2, T2, nTg);
				GET_UNION_SUBFIELD_VEC(vec, tgs, tgData, tg2, T1, nTg);
				//Exc
				GET_FIELD_VEC(uvec, excs, synNumber, nExc);
				GET_FIELD_VEC(uvec, excs, excType, nExc);
				GET_FIELD_VEC(uvec, excs, status, nExc);
				// Exc - Exc1
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, vMax, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, vMin, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, mu0, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, T1, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, T2, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, T3, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, T4, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, Te, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, Tr, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, Ae, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc1, Be, nExc);
				// Exc - Exc2
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, vMax, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, vMin, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Ka, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Ta, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Kf, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Tf, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Te, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Tr, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Ae, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc2, Be, nExc);
				// Exc - Exc3
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, vMax, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, vMin, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, mu0, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, T2, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, T1, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, vf0, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, V0, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, Te, nExc);
				GET_UNION_SUBFIELD_VEC(vec, excs, excData, exc3, Tr, nExc);

			private:
				void freeSpace() {
					if (nBus > 0) {
						delete[] buses;
					}
					if (nSw > 0) {
						delete[] sws;
					}
					if (nPv > 0) {
						delete[] pvs;
					}
					if (nPq > 0) {
						delete[] pqs;
					}
					if (nShunt > 0) {
						delete[] shunts;
					}
					if (nLine > 0) {
						delete[] lines;
					}
					if (nPl > 0) {
						delete[] pls;
					}
					if (nSyn > 0) {
						delete[] syns;
					}
					if (nInd > 0) {
						delete[] inds;
					}
					if (nTg > 0) {
						delete[] tgs;
					}
					if (nExc > 0) {
						delete[] excs;
					}
				}

				void copyMembers(const PsatDataSet& psat, bool regNew) {
					reset();
					this->nBus = psat.nBus;
					if (this->nBus > 0) {
						this->buses = new Bus[this->nBus];
						for (int i = 0; i < this->nBus; i++) this->buses[i] = Bus(psat.buses[i], regNew);
					}
					this->nSw = psat.nSw;
					if (this->nSw > 0) {
						this->sws = new SW[this->nSw];
						for (int i = 0; i < this->nSw; i++) this->sws[i] = SW(psat.sws[i], regNew);
					}
					this->nPv = psat.nPv;
					if (this->nPv > 0) {
						this->pvs = new PV[this->nPv];
						for (int i = 0; i < this->nPv; i++) this->pvs[i] = PV(psat.pvs[i], regNew);
					}
					this->nPq = psat.nPq;
					if (this->nPq > 0) {
						this->pqs = new PQ[this->nPq];
						for (int i = 0; i < this->nPq; i++) this->pqs[i] = PQ(psat.pqs[i], regNew);
					}
					this->nShunt = psat.nShunt;
					if (this->nShunt > 0) {
						this->shunts = new Shunt[this->nShunt];
						for (int i = 0; i < this->nShunt; i++) this->shunts[i] = Shunt(psat.shunts[i], regNew);
					}
					this->nLine = psat.nLine;
					if (this->nLine > 0) {
						this->lines = new Line[this->nLine];
						for (int i = 0; i < this->nLine; i++) this->lines[i] = Line(psat.lines[i], regNew);
					}
					this->nPl = psat.nPl;
					if (this->nPl > 0) {
						this->pls = new Pl[this->nPl];
						for (int i = 0; i < this->nPl; i++) this->pls[i] = Pl(psat.pls[i], regNew);
					}
					this->nSyn = psat.nSyn;
					if (this->nSyn > 0) {
						this->syns = new Syn[this->nSyn];
						for (int i = 0; i < this->nSyn; i++) this->syns[i] = Syn(psat.syns[i], regNew);
					}
					this->nInd = psat.nInd;
					if (this->nInd > 0) {
						this->inds = new Ind[this->nInd];
						for (int i = 0; i < this->nInd; i++) this->inds[i] = Ind(psat.inds[i], regNew);
					}
					this->nTg = psat.nTg;
					if (this->nTg > 0) {
						this->tgs = new Tg[this->nTg];
						for (int i = 0; i < this->nTg; i++) this->tgs[i] = Tg(psat.tgs[i], regNew);
					}
					this->nExc = psat.nExc;
					if (this->nExc > 0) {
						this->excs = new Exc[this->nExc];
						for (int i = 0; i < this->nExc; i++) this->excs[i] = Exc(psat.excs[i], regNew);
					}

					this->newToOld = psat.newToOld;
					this->oldToNew = psat.oldToNew;
					this->isFormatted = psat.isFormatted;
				}
			};
		}
	}
}

#endif