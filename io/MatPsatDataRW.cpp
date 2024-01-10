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
#include "io/MatPsatDataRW.h"
#include "matio.h"
#include <iostream>
#include <cmath>

using namespace std;

#define IND2(ir,ic,nr) ((ic)*(nr)+ir)

namespace che {
	namespace io {
		namespace chedata {

			int MatPsatReader::parse(const char *filePath, PsatDataSet *psatData) {
				this->psatData = psatData;
				return parse(filePath);
			}

			int MatPsatReader::parse(const char *filePath) {
				mat_t *matfp;
				matvar_t *matvar;
				double* matdata;

				matfp = Mat_Open(filePath, MAT_ACC_RDONLY);
				if (NULL == matfp) {
					cerr << "Error opening MAT file \"<<filePath<<\"!" << endl;
					return CHE_IO_FAIL;
				}

				psatData->reset(); // Drop all existing data and get ready to obtain new data.

#define NEXT_C matdata[IND2(row, col++, matvar->dims[0])]

				// Read bus data
				matvar = Mat_VarRead(matfp, "bus");
				if (matvar == NULL) { cerr << "Does not contain bus data." << endl; return CHE_IO_FAIL; }
				if (matvar->rank != 2) { cerr << "Bus data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
				if (matvar->dims[1] != 6) { cerr << "Expect 6 columns in bus data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

				psatData->nBus = matvar->dims[0];
				matdata = (double*)matvar->data;
				if (psatData->nBus > 0) {
					psatData->buses = new Bus[psatData->nBus];
					for (int row = 0; row < psatData->nBus; row++) {
						int col = 0;
						psatData->buses[row].busNumber = (int)NEXT_C;
						psatData->buses[row].baseV = NEXT_C;
						psatData->buses[row].vMag = NEXT_C;
						psatData->buses[row].vAng = NEXT_C;
						psatData->buses[row].areaNumber = (int)NEXT_C;
						psatData->buses[row].regionNumber = (int)NEXT_C;
					}
				}
				Mat_VarFree(matvar);

				// Read SW data
				matvar = Mat_VarRead(matfp, "sw");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "SW data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 13) { cerr << "Expect 13 columns in SW data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nSw = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nSw > 0) {
						psatData->sws = new SW[psatData->nSw];
						for (int row = 0; row < psatData->nSw; row++) {
							int col = 0;
							psatData->sws[row].busNumber = (int)NEXT_C;
							psatData->sws[row].baseS = NEXT_C;
							psatData->sws[row].baseV = NEXT_C;
							psatData->sws[row].vMag = NEXT_C;
							psatData->sws[row].vAng = NEXT_C;
							psatData->sws[row].qMax = NEXT_C;
							psatData->sws[row].qMin = NEXT_C;
							psatData->sws[row].vMax = NEXT_C;
							psatData->sws[row].vMin = NEXT_C;
							psatData->sws[row].pg0 = NEXT_C;
							psatData->sws[row].lossParticipation = NEXT_C;
							psatData->sws[row].isRefBus = (unsigned char)NEXT_C;
							psatData->sws[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read PV data
				matvar = Mat_VarRead(matfp, "pv");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "PV data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 11) { cerr << "Expect 11 columns in PV data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nPv = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nPv > 0) {
						psatData->pvs = new PV[psatData->nPv];
						for (int row = 0; row < psatData->nPv; row++) {
							int col = 0;
							psatData->pvs[row].busNumber = (int)NEXT_C;
							psatData->pvs[row].baseS = NEXT_C;
							psatData->pvs[row].baseV = NEXT_C;
							psatData->pvs[row].P = NEXT_C;
							psatData->pvs[row].vMag = NEXT_C;
							psatData->pvs[row].qMax = NEXT_C;
							psatData->pvs[row].qMin = NEXT_C;
							psatData->pvs[row].vMax = NEXT_C;
							psatData->pvs[row].vMin = NEXT_C;
							psatData->pvs[row].lossParticipation = NEXT_C;
							psatData->pvs[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read PQ data
				matvar = Mat_VarRead(matfp, "pq");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "PQ data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 9) { cerr << "Expect 9 columns in PQ data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nPq = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nPq > 0) {
						psatData->pqs = new PQ[psatData->nPq];
						for (int row = 0; row < psatData->nPq; row++) {
							int col = 0;
							psatData->pqs[row].busNumber = (int)NEXT_C;
							psatData->pqs[row].baseS = NEXT_C;
							psatData->pqs[row].baseV = NEXT_C;
							psatData->pqs[row].P = NEXT_C;
							psatData->pqs[row].Q = NEXT_C;
							psatData->pqs[row].vMax = NEXT_C;
							psatData->pqs[row].vMin = NEXT_C;
							psatData->pqs[row].allowConversion = (unsigned char)NEXT_C;
							psatData->pqs[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read Shunt data
				matvar = Mat_VarRead(matfp, "shunt");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Shunt data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 7) { cerr << "Expect 7 columns in Shunt data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nShunt = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nShunt > 0) {
						psatData->shunts = new Shunt[psatData->nShunt];
						for (int row = 0; row < psatData->nShunt; row++) {
							int col = 0;
							psatData->shunts[row].busNumber = (int)NEXT_C;
							psatData->shunts[row].baseS = NEXT_C;
							psatData->shunts[row].baseV = NEXT_C;
							psatData->shunts[row].baseF = NEXT_C;
							psatData->shunts[row].g = NEXT_C;
							psatData->shunts[row].b = NEXT_C;
							psatData->shunts[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read Line data
				matvar = Mat_VarRead(matfp, "line");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Line data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 16) { cerr << "Expect 16 columns in Line data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nLine = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nLine > 0) {
						psatData->lines = new Line[psatData->nLine];
						for (int row = 0; row < psatData->nLine; row++) {
							int col = 0;
							psatData->lines[row].fromBus = (int)NEXT_C;
							psatData->lines[row].toBus = (int)NEXT_C;
							psatData->lines[row].baseS = NEXT_C;
							psatData->lines[row].baseV = NEXT_C;
							psatData->lines[row].baseF = NEXT_C;
							psatData->lines[row].len = NEXT_C;
							psatData->lines[row].kT = NEXT_C;
							psatData->lines[row].r = NEXT_C;
							psatData->lines[row].x = NEXT_C;
							psatData->lines[row].b = NEXT_C;
							psatData->lines[row].k = NEXT_C;
							psatData->lines[row].ang = NEXT_C;
							psatData->lines[row].iMax = NEXT_C;
							psatData->lines[row].pMax = NEXT_C;
							psatData->lines[row].sMax = NEXT_C;
							psatData->lines[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read ZIP data
				matvar = Mat_VarRead(matfp, "zip");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Zip data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 12) { cerr << "Expect 12 columns in Zip data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nPl = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nPl > 0) {
						psatData->pls = new Pl[psatData->nPl];
						for (int row = 0; row < psatData->nPl; row++) {
							int col = 0;
							psatData->pls[row].busNumber = (int)NEXT_C;
							psatData->pls[row].baseS = NEXT_C;
							psatData->pls[row].baseV = NEXT_C;
							psatData->pls[row].baseF = NEXT_C;
							psatData->pls[row].g = NEXT_C;
							psatData->pls[row].Ip = NEXT_C;
							psatData->pls[row].P = NEXT_C;
							psatData->pls[row].b = NEXT_C;
							psatData->pls[row].Iq = NEXT_C;
							psatData->pls[row].Q = NEXT_C;
							psatData->pls[row].initAfterPF = (unsigned char)NEXT_C;
							psatData->pls[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read Syn data
				matvar = Mat_VarRead(matfp, "syn");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Syn data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] < 19) { cerr << "Expect >=19 columns in Syn data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nSyn = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nSyn > 0) {
						psatData->syns = new Syn[psatData->nSyn];
						for (int row = 0; row < psatData->nSyn; row++) {
							int col = 0;
							psatData->syns[row].busNumber = (int)NEXT_C;
							psatData->syns[row].baseS = NEXT_C;
							psatData->syns[row].baseV = NEXT_C;
							psatData->syns[row].baseF = NEXT_C;
							psatData->syns[row].model = NEXT_C;
							psatData->syns[row].xl = NEXT_C;
							psatData->syns[row].ra = NEXT_C;
							psatData->syns[row].xd = NEXT_C;
							psatData->syns[row].xd1 = NEXT_C;
							psatData->syns[row].xd2 = NEXT_C;
							psatData->syns[row].Td01 = NEXT_C;
							psatData->syns[row].Td02 = NEXT_C;
							psatData->syns[row].xq = NEXT_C;
							psatData->syns[row].xq1 = NEXT_C;
							psatData->syns[row].xq2 = NEXT_C;
							psatData->syns[row].Tq01 = NEXT_C;
							psatData->syns[row].Tq02 = NEXT_C;
							psatData->syns[row].M = NEXT_C;
							psatData->syns[row].D = NEXT_C;
							psatData->syns[row].Ko = (matvar->dims[1] >= 20) ? NEXT_C : NAN;
							psatData->syns[row].Kp = (matvar->dims[1] >= 21) ? NEXT_C : NAN;
							psatData->syns[row].gammaP = (matvar->dims[1] >= 22) ? NEXT_C : NAN;
							psatData->syns[row].gammaQ = (matvar->dims[1] >= 23) ? NEXT_C : NAN;
							psatData->syns[row].TAA = (matvar->dims[1] >= 24) ? NEXT_C : NAN;
							psatData->syns[row].sat1 = (matvar->dims[1] >= 25) ? NEXT_C : NAN;
							psatData->syns[row].sat2 = (matvar->dims[1] >= 26) ? NEXT_C : NAN;
							psatData->syns[row].nCOI = (matvar->dims[1] >= 27) ? (int)NEXT_C : NAN;
							psatData->syns[row].status = (matvar->dims[1] >= 28) ? (unsigned char)NEXT_C : 1;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read Ind data
				matvar = Mat_VarRead(matfp, "ind");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Ind data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 20) { cerr << "Expect 20 columns in Ind data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nInd = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nInd > 0) {
						psatData->inds = new Ind[psatData->nInd];
						for (int row = 0; row < psatData->nInd; row++) {
							int col = 0;
							psatData->inds[row].busNumber = (int)NEXT_C;
							psatData->inds[row].baseS = NEXT_C;
							psatData->inds[row].baseV = NEXT_C;
							psatData->inds[row].baseF = NEXT_C;
							psatData->inds[row].model = NEXT_C;
							psatData->inds[row].startCtrl = (unsigned char)NEXT_C;
							psatData->inds[row].rs = NEXT_C;
							psatData->inds[row].xs = NEXT_C;
							psatData->inds[row].rr1 = NEXT_C;
							psatData->inds[row].xr1 = NEXT_C;
							psatData->inds[row].rr2 = NEXT_C;
							psatData->inds[row].xr2 = NEXT_C;
							psatData->inds[row].xm = NEXT_C;
							psatData->inds[row].Hm = NEXT_C;
							psatData->inds[row].Ta = NEXT_C;
							psatData->inds[row].Tb = NEXT_C;
							psatData->inds[row].Tc = NEXT_C;
							psatData->inds[row].tup = NEXT_C;
							psatData->inds[row].allowBrake = (unsigned char)NEXT_C;
							psatData->inds[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read Tg data
				matvar = Mat_VarRead(matfp, "tg");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Tg data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] <8 ) { cerr << "Expect at least 8 columns in Tg data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nTg = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nTg > 0) {
						psatData->tgs = new Tg[psatData->nTg];
						for (int row = 0; row < psatData->nTg; row++) {
							int col = 0;
							psatData->tgs[row].synNumber = (int)NEXT_C;
							psatData->tgs[row].tgType = (unsigned char)NEXT_C;
							if (psatData->tgs[row].tgType == 1) {
								if (matvar->dims[1] < 11) { cerr << "Expect at least 11 columns in Tg1 data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
								psatData->tgs[row].tgData.tg1.wref0 = NEXT_C;
								psatData->tgs[row].tgData.tg1.R = NEXT_C;
								psatData->tgs[row].tgData.tg1.Tmax = NEXT_C;
								psatData->tgs[row].tgData.tg1.Tmin = NEXT_C;
								psatData->tgs[row].tgData.tg1.Ts = NEXT_C;
								psatData->tgs[row].tgData.tg1.Tc = NEXT_C;
								psatData->tgs[row].tgData.tg1.T3 = NEXT_C;
								psatData->tgs[row].tgData.tg1.T4 = NEXT_C;
								psatData->tgs[row].tgData.tg1.T5 = NEXT_C;
							}
							else if (psatData->tgs[row].tgType == 2) {
								psatData->tgs[row].tgData.tg2.wref0 = NEXT_C;
								psatData->tgs[row].tgData.tg2.R = NEXT_C;
								psatData->tgs[row].tgData.tg2.Tmax = NEXT_C;
								psatData->tgs[row].tgData.tg2.Tmin = NEXT_C;
								psatData->tgs[row].tgData.tg2.T2 = NEXT_C;
								psatData->tgs[row].tgData.tg2.T1 = NEXT_C;
								col += 3;
							}
							else {
								cerr << "Unknown Tg type." << endl; return CHE_IO_FAIL;
							}
							psatData->tgs[row].status = (matvar->dims[1] >= 12) ? (unsigned char)NEXT_C : 1;
						}
					}
					Mat_VarFree(matvar);
				}

				// Read Esc data
				matvar = Mat_VarRead(matfp, "exc");
				if (matvar != NULL) {
					if (matvar->rank != 2) { cerr << "Exc data rank not 2." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }
					if (matvar->dims[1] != 14) { cerr << "Expect 14 columns in Exc data." << endl; Mat_VarFree(matvar); return CHE_IO_FAIL; }

					psatData->nExc = matvar->dims[0];
					matdata = (double*)matvar->data;
					if (psatData->nExc > 0) {
						psatData->excs = new Exc[psatData->nExc];
						for (int row = 0; row < psatData->nExc; row++) {
							int col = 0;
							psatData->excs[row].synNumber = (int)NEXT_C;
							psatData->excs[row].excType = (unsigned char)NEXT_C;
							if (psatData->excs[row].excType == 1) {
								psatData->excs[row].excData.exc1.vMax = NEXT_C;
								psatData->excs[row].excData.exc1.vMin = NEXT_C;
								psatData->excs[row].excData.exc1.mu0 = NEXT_C;
								psatData->excs[row].excData.exc1.T1 = NEXT_C;
								psatData->excs[row].excData.exc1.T2 = NEXT_C;
								psatData->excs[row].excData.exc1.T3 = NEXT_C;
								psatData->excs[row].excData.exc1.T4 = NEXT_C;
								psatData->excs[row].excData.exc1.Te = NEXT_C;
								psatData->excs[row].excData.exc1.Tr = NEXT_C;
								psatData->excs[row].excData.exc1.Ae = NEXT_C;
								psatData->excs[row].excData.exc1.Be = NEXT_C;
							}
							else if (psatData->excs[row].excType == 2) {
								psatData->excs[row].excData.exc2.vMax = NEXT_C;
								psatData->excs[row].excData.exc2.vMin = NEXT_C;
								psatData->excs[row].excData.exc2.Ka = NEXT_C;
								psatData->excs[row].excData.exc2.Ta = NEXT_C;
								psatData->excs[row].excData.exc2.Kf = NEXT_C;
								psatData->excs[row].excData.exc2.Tf = NEXT_C;
								col++;
								psatData->excs[row].excData.exc2.Te = NEXT_C;
								psatData->excs[row].excData.exc2.Tr = NEXT_C;
								psatData->excs[row].excData.exc2.Ae = NEXT_C;
								psatData->excs[row].excData.exc2.Be = NEXT_C;
							}
							else if (psatData->excs[row].excType == 3) {
								psatData->excs[row].excData.exc3.vMax = NEXT_C;
								psatData->excs[row].excData.exc3.vMin = NEXT_C;
								psatData->excs[row].excData.exc3.mu0 = NEXT_C;
								psatData->excs[row].excData.exc3.T2 = NEXT_C;
								psatData->excs[row].excData.exc3.T1 = NEXT_C;
								psatData->excs[row].excData.exc3.vf0 = NEXT_C;
								psatData->excs[row].excData.exc3.V0 = NEXT_C;
								psatData->excs[row].excData.exc3.Te = NEXT_C;
								psatData->excs[row].excData.exc3.Tr = NEXT_C;
							}
							else {
								cerr << "Unknown Tg type." << endl; return CHE_IO_FAIL;
							}
							psatData->excs[row].status = (unsigned char)NEXT_C;
						}
					}
					Mat_VarFree(matvar);
				}

				Mat_Close(matfp);

				return CHE_IO_SUCCESS;
			}
		}
	}
}