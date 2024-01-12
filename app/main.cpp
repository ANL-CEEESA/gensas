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
#include "util/AbstractCheCalculator.h"
#include "pf/ChePFCalculator.h"
#include "io/MatPsatDataRW.h"
#include "util/SafeArmadillo.h"
#include "util/CheCompUtil.h"
#include "nvwa/pctimer.h"

#include "sas/SasInput.h"
#include "sas/SasLexico.h"
#include "sas/SasExpr.h"
#include "sas/SasComputation.h"
#include <memory>
// #include <memory>
// #include "debug_new.h"

using namespace arma;
using namespace che::core;
using namespace che::util;
using namespace nvwa;

#ifdef __DIR__
#undef __DIR__
#endif

#ifdef WIN32
#define __DIR__ std::string(__FILE__).substr(0,std::string(__FILE__).find_last_of("\\"))
#elif defined linux
#define __DIR__ (std::getenv("PWD"))
#else
#define __DIR__ ""
#endif

#ifdef WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>

// #define new DEBUG_NEW

// using namespace nvwa;

std::string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

int main(int argc, char **argv)
{

	string compMode = argv[1];
	cout << "compMode=" << compMode << endl;

	if (compMode == "-g") {

		bool verbose = false;
		bool fileMode = false;
		bool writeJson = false;
		SasComputationOptions options;
		options.alphaTol = 1e-3;
		options.errorTol = 1e-5;
		options.segment = 1.0;
		double outStep = 0.01;

		const double MAX_TIME = 1000.0;
		const double MIN_TIME = 0.0;
		const double MAX_SEGMENT = 10.0;
		const double MIN_SEGMENT = 1e-4;
		const double MAX_ALPHA_IN_SEG = 1;
		const double MIN_ALPHA_IN_SEG = 1e-9;
		const double MAX_ERRTOL = 1e-2;
		const double MIN_ERRTOL = 1e-10;
		const double MAX_OUTSTEP = 1;
		const double MIN_OUTSTEP = 1e-3;
		string input = ""; 
		string output = "";
		string jsonOutput = "";
		//	parsing command line arguments
		for (int iArg = 2; iArg < argc; iArg++)
		{
			string arg = argv[iArg];
			if (arg == "--mode" || arg == "-m")
			{
				if (++iArg < argc)
				{
					string subArg = argv[iArg];
					if (subArg == "file")
					{
						fileMode = true;
					}
					else if (subArg == "string")
					{
						fileMode = false;
					}
					else
					{
						cerr << "Mode name [file|string] should be specified after --mode or -m. Using string as default." << endl;
					}
				}
				else
				{
					cerr << "Mode name [file|string] should be specified after --mode or -m. Using string as default." << endl;
				}
			}
			else if (arg == "--input" || arg == "-i")
			{
				if (++iArg < argc)
				{
					input = argv[iArg];
				}
				else
				{
					cerr << "Input should be specified after --input or -i. Using empty string as default." << endl;
					input = "";
				}
			}
			else if (arg == "--output" || arg == "-o")
			{
				if (++iArg < argc)
				{
					output = argv[iArg];
				}
				else
				{
					cerr << "Output should be specified after --output or -o. Using empty string as default." << endl;
					output = "";
				}
			}
			else if (arg == "--json" || arg == "-j")
			{
				if (++iArg < argc)
				{
					jsonOutput = argv[iArg];
					writeJson = true;
				}
				else
				{
					cerr << "JSON output should be specified after --json or -j. Using empty string as default." << endl;
					jsonOutput = "";
				}
			}
			else if (arg == "--time" || arg == "-t")
			{
				if (++iArg < argc)
				{
					options.endTime = stod(argv[iArg]);
				}
				else
				{
					cerr << "Simulation length (s) should be specified after --time or -t. Using default value " << options.endTime << "." << endl;
				}
			}
			else if (arg == "--level" || arg == "-l")
			{
				if (++iArg < argc)
				{
					options.nLvl = stoi(argv[iArg]);
				}
				else
				{
					cerr << "SAS order --level or -l. Using default value " << options.nLvl << "." << endl;
				}
			}
			else if (arg == "--atol" || arg == "-a")
			{
				if (++iArg < argc)
				{
					options.alphaTol = stod(argv[iArg]);
				}
				else
				{
					cerr << "Alpha tolerance should be specified after --atol or -a. Using default value " << options.alphaTol << "." << endl;
				}
			}
			else if (arg == "--segment" || arg == "-s")
			{
				if (++iArg < argc)
				{
					options.segment = stod(argv[iArg]);
				}
				else
				{
					cerr << "Segment length should be specified after --segment or -s. Using default value " << options.segment << "." << endl;
				}
			}
			else if (arg == "--etol" || arg == "-e")
			{
				if (++iArg < argc)
				{
					options.errorTol = stod(argv[iArg]);
				}
				else
				{
					cerr << "Error tolerance should be specified after --etol or -e. Using default value " << options.errorTol << "." << endl;
				}
			}
			else if (arg == "--step" || arg == "-k")
			{
				if (++iArg < argc)
				{
					outStep = stod(argv[iArg]);
				}
				else
				{
					cerr << "Output step should be specified after --step or -k. Using default value " << outStep << "." << endl;
				}
			}
			else if (arg == "--verbose" || arg == "-v")
			{
				verbose = true;
			}
		}

		if (options.endTime > MAX_TIME) {
			cout << "Simulation run length are limited up to " << MAX_TIME << " s." << endl;
			options.endTime = MAX_TIME;
		}

		if (options.endTime < MIN_TIME) {
			cout << "Simulation run length must be at least " << MIN_TIME << " s." << endl;
			options.endTime = MIN_TIME;
		}

		if (options.segment > MAX_SEGMENT) {
			cout << "Simulation segment length are limited up to " << MAX_SEGMENT << " s." << endl;
			options.segment = MAX_SEGMENT;
		}

		if (options.segment < MIN_SEGMENT) {
			cout << "Simulation segment length must be at least " << MIN_SEGMENT << " s." << endl;
			options.segment = MIN_SEGMENT;
		}

		if (options.alphaTol > MAX_ALPHA_IN_SEG*options.segment) {
			cout << "Alpha tolerance are limited up to " << MAX_ALPHA_IN_SEG * options.segment << " s." << endl;
			options.alphaTol = MAX_ALPHA_IN_SEG * options.segment;
		}

		if (options.alphaTol < MIN_ALPHA_IN_SEG*options.segment) {
			cout << "Alpha tolerance must be at least " << MIN_ALPHA_IN_SEG * options.segment << " s." << endl;
			options.alphaTol = MIN_ALPHA_IN_SEG * options.segment;
		}

		if (options.errorTol > MAX_ERRTOL) {
			cout << "Simulation errorTol are limited up to " << MAX_ERRTOL << " s." << endl;
			options.errorTol = MAX_ERRTOL;
		}

		if (options.errorTol < MIN_ERRTOL) {
			cout << "Simulation errorTol must be at least " << MIN_ERRTOL << " s." << endl;
			options.errorTol = MIN_ERRTOL;
		}

		if (outStep > MAX_OUTSTEP) {
			cout << "Simulation output step is limited up to " << MAX_OUTSTEP << " s." << endl;
			outStep = MAX_OUTSTEP;
		}

		if (outStep < MIN_OUTSTEP) {
			cout << "Simulation output step must be at least " << MIN_OUTSTEP << " s." << endl;
			outStep = MIN_OUTSTEP;
		}

		// end parsing command line arguments, begin main program
		pctimer_t stTime = pctimer();

		SasModelParser parser;
		if (fileMode)
		{
			parser.parseFromFile(input.c_str());
		}
		else
		{
			if (input.length() >= 2) {
				cout << input.substr(1, input.length() - 2).c_str() << endl;
				parser.parseString(input.substr(1, input.length() - 2).c_str());
			}
			else {
				cout << input.c_str() << endl;
				parser.parseString(input.c_str());
			}
		}
		if (parser.reader.buffer == NULL || parser.reader.size <= 1) {
			cerr << "Error: reader fails." << endl;
			return -1;
		}

		parser.currentModel = GlobalPool::getInstance().createNewTopModel();
		parseSingleSasModel(&parser);
		if (parser.errCount > 0)
		{
			cerr << "Error: Parsing model has errors." << endl;
			return -1;
		}

		SasComputationModel compModel;
		compModel.generateCompModel(GlobalPool::getInstance());
		compModel.preprocessModels();
		// check if all the ids are initialized, need to be repacked in the future.
		list<shared_ptr<SasModel>>::iterator modelIt;
		int nonInitCnt = 0;
		for (modelIt = compModel.topModels.begin(); modelIt != compModel.topModels.end(); modelIt++)
		{
			list<shared_ptr<IdPearl>>::iterator idIt;
			for (idIt = (*modelIt)->idTable.begin(); idIt != (*modelIt)->idTable.end(); idIt++) {
				IdType idt = (*idIt)->idTy;
				if (idt == ID_VAR || idt == ID_TMPVAR || idt == ID_IVAR) {
					if (!(*idIt)->initVal) {
						nonInitCnt++;
					}
				}
			}
		}
		if (nonInitCnt > 0) {
			cerr << "Error: There are " << nonInitCnt << " uninitialized variables." << endl;
		}
		compModel.generateDAEs();

		if (verbose)
		{
			cout << "Original model:" << endl;
			parser.currentModel->printAllModelTrees();
			cout << "Copied model:" << endl;

			list<shared_ptr<SasModel>>::iterator modelIt;
			for (modelIt = compModel.topModels.begin(); modelIt != compModel.topModels.end(); modelIt++)
			{
				(*modelIt)->printAllModelTrees();
				(*modelIt)->printIdTable();
			}
		}

		// Check if DAEs are valid.
		bool isValidDAE = true;
		if (compModel.nX != compModel.nDE) {
			cerr << "Error: Number of state variable (" << compModel.nX << ") not equal number of diff. eqs. (" << compModel.nDE << ")." << endl;
			isValidDAE = false;
		}
		if (compModel.nY != compModel.nAE) {
			cerr << "Error: Number of algebraic variable (" << compModel.nY << ") not equal number of algeb. eqs. (" << compModel.nAE << ")." << endl;
			isValidDAE = false;
		}
		if (!isValidDAE) {
			cerr << "Error: DAE is not well defined." << endl;
			return -1;
		}
		SasSolutionSet *sasSolLink = compModel.solve(options);

		pctimer_t endTime = pctimer();

		cout << "Computation time: " << endTime - stTime << " s." << endl;

		if(writeJson){
			sasSolLink->writeJSONFile(jsonOutput.c_str(), outStep);
			cout << "Json output written to file "<<jsonOutput<<"."<<endl;
		}

		delete sasSolLink;

		return 0;

	}
	else if (compMode == "-p") {

		string filePath = "";
		string outputPath = "detault.mat";
		int nlvl = 15;
		double segment = 1.0;
		double alphaTol = 1e-4;
		double diffTol = 1e-6;
		double diffTolMax = 1e-2;
		int repeat = 10;
		for (int iArg = 2; iArg < argc; iArg++)
		{
			string arg = argv[iArg];
			if (arg == "--file" || arg == "-f")
			{
				if (++iArg < argc)
				{
					filePath = argv[iArg];
				}
				else
				{
					cerr << "File name should be specified after --file or -f." << endl;
				}
			}
			else if (arg == "--level" || arg == "-l")
			{
				if (++iArg < argc)
				{
					nlvl = stoi(argv[iArg]);
				}
				else
				{
					cerr << "nlvl should be specified after --level or -l. Using nlvl=" << nlvl << " as default." << endl;
				}
			}
			else if (arg == "--repeat" || arg == "-r")
			{
				if (++iArg < argc)
				{
					repeat = stoi(argv[iArg]);
				}
				else
				{
					cerr << "repeat should be specified after --repeat or -r. Using repeat=" << repeat << " as default." << endl;
				}
			}
			else if (arg == "--segment" || arg == "-s")
			{
				if (++iArg < argc)
				{
					segment = stod(argv[iArg]);
				}
				else
				{
					cerr << "segment should be specified after --segment or -s. Using segment=" << segment << " as default." << endl;
				}
			}
			else if (arg == "--alphatol" || arg == "-a")
			{
				if (++iArg < argc)
				{
					alphaTol = stod(argv[iArg]);
				}
				else
				{
					cerr << "alphaTol should be specified after --alphatol or -a. Using alphaTol=" << alphaTol << " as default." << endl;
				}
			}
			else if (arg == "--difftol" || arg == "-d")
			{
				if (++iArg < argc)
				{
					diffTol = stod(argv[iArg]);
				}
				else
				{
					cerr << "diffTol should be specified after --difftol or -d. Using diffTol=" << diffTol << " as default." << endl;
				}
			}
			else if (arg == "--output" || arg == "-o")
			{
				if (++iArg < argc)
				{
					string subArg = argv[iArg];
					outputPath = GetCurrentWorkingDir() + "/" + subArg;
				}
				else
				{
					cerr << "Output file name should be specified after --file or -f." << endl;
				}
			}
			else if (arg == "--difftolmax" || arg == "-e")
			{
				if (++iArg < argc)
				{
					diffTolMax = stod(argv[iArg]);
				}
				else
				{
					cerr << "diffTolMax should be specified after --difftolmax or -e. Using diffTolMax=" << diffTolMax << " as default." << endl;
				}
			}
		}

		if (nlvl < 3) {
			nlvl = 3;
		}
		if (nlvl > 100) {
			nlvl = 100;
		}
		if (repeat > 100) {
			repeat = 100;
		}
		if (segment < 0.01) {
			segment = 0.01;
		}
		if (segment > 1.0) {
			segment = 1.0;
		}
		if (alphaTol < 1e-8) {
			alphaTol = 1e-8;
		}
		if (alphaTol > 1e-2) {
			alphaTol = 1e-2;
		}
		if (diffTol < 1e-10) {
			diffTol = 1e-10;
		}
		if (diffTol > 1e-2) {
			diffTol = 1e-2;
		}
		if (diffTolMax < 1e-4) {
			diffTolMax = 1e-4;
		}
		if (diffTolMax < 100.0*diffTol) {
			diffTolMax = 100.0*diffTol;
		}

		chedata::MatPsatReader matReader;
		chedata::PsatDataSet psatData;
		//  string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_dcase2383wp_mod2_ind_zip3.mat";
		//string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_014_ind_zip1.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_ei_458_100.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_dcase2383wp_mod2_zip9x.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_70k.mat";
		int flag = matReader.parse(filePath.c_str(), &psatData);
		psatData.renumberBuses();
		uvec islands = CheCompUtil::searchIslands(psatData);
		// islands.print("islands");

		CheCompOptions compOpt(nlvl, 1.0, alphaTol, segment, diffTol, diffTolMax);

		pctimer_t totalTime = 0.;

		for (int i = 0; i < repeat; i++) {
			AbstractCheCalculator* pCalculator = new ChePfCalculator(psatData, compOpt, islands);
			pctimer_t stTime = pctimer();
			int pfFlag = pCalculator->calc();
			pctimer_t endTime = pctimer();
			pCalculator->writeMatFile(outputPath.c_str());
			delete pCalculator;
			totalTime += endTime - stTime;

			cout << "Computation No." << i + 1 << endl;
		}

		cout << "Computation time: " << totalTime << " s." << endl;

		return 0;

	}
	else {

		cerr << "The first arg should either be -g (general SAS) or -p (power flow)." << endl;
		return 0;
	}

}