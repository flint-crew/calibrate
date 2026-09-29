#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <stdexcept>
#include "calibrator.h"
#include "calibrationmethod.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout
			<< "Usage: calibrate [-p <phases.txt> <gains.txt>] [-refmode <0|1|2>] [-minuv <min uvw dist in m>] [-maxuv <max uvw dist in m>] [-startscan <scan>] [-endscan <scan>] [-a <min-accuracy> <stop-accuracy>] [-i <niter>] [-j <threads>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-t timesteps] [-ch <channels per solution>] [-interval <start timestep> <end timestep>] [-absmem <memory in GB>] [-datacolumn <name>] [-weightcolumn <name>] [-quiet] <measurementset.ms> <solutions.bin>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n"
			<< "The algorithm is described by Offringa et al. (2016), MNRAS 458, 1057, doi:10.1093/mnras/stw310; please cite it when using this program.\n"
            << "refmode=0 process all baselines; =1 only include baselines to reference antenna; =2 exclude baselines to reference antenna.\n"
            << "rhs: accepted for compatibility, but has no effect.\n"
            << "ch: solve one solution per block of this many channels (default 1: every channel).\n"
            << "interval: only use timesteps start to end-1 of the measurement set (counting from 0).\n"
            << "absmem: memory to plan for in GB, instead of the machine's physical memory.\n"
            << "weightcolumn: WEIGHT_SPECTRUM, WEIGHT, SIGMA_SPECTRUM or SIGMA (SIGMA columns are used as 1/sigma^2). Default: WEIGHT_SPECTRUM if it has values, otherwise WEIGHT.\n";
	} else {
		int argi = 1;
		// bool saveCrossTermsPlotFile = false, saveFaradayPlotFiles = false;
		bool
			savePlotFiles = false, //beamOnSource = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false, doQuiet = false;
		std::string plotPhaseFile, plotGainFile, plotFaradayFile, crossTermsPlotFile, rhsSolutionFile, modelFile;
		size_t niter = CalibrationMethod::DefaultNIter(), solutionInterval = 0, startScan = -1, endScan = -1, refMode=0;
		size_t solutionChannels = 1, intervalStart = 0, intervalEnd = 0;
		bool hasInterval = false;
		double absMem = 0.0;
		std::string dataColumnName = "DATA", weightColumnName;
		double
			minAccuracy = CalibrationMethod::DefaultMinAccuracy(),
			stopAccuracy = CalibrationMethod::DefaultStoppingAccuracy(),
			minUVW = -1.0,
			maxUVW = -1.0;
		size_t threadCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
		while(argv[argi][0] == '-')
		{
			std::string param(&argv[argi][1]);
			if(param == "p")
			{
				savePlotFiles = true;
				plotPhaseFile = argv[argi+1];
				plotGainFile = argv[argi+2];
				argi += 3;
			}
			else if(param == "datacolumn")
			{
				dataColumnName = argv[argi+1];
				argi += 2;
			}
			else if(param == "weightcolumn")
			{
				weightColumnName = argv[argi+1];
				argi += 2;
			}
			else if(param == "i")
			{
				niter = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "j")
			{
				threadCount = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "a")
			{
				minAccuracy = atof(argv[argi+1]);
				stopAccuracy = atof(argv[argi+2]);
				argi += 3;
			}
			else if(param == "m")
			{
				modelFile = argv[argi + 1];
                argi +=2 ;
                std::cout << modelFile << "\n";    
			}
			else if(param == "t")
			{
				solutionInterval = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "ch")
			{
				solutionChannels = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "interval")
			{
				hasInterval = true;
				intervalStart = atoi(argv[argi+1]);
				intervalEnd = atoi(argv[argi+2]);
				argi += 3;
			}
			else if(param == "absmem")
			{
				absMem = atof(argv[argi+1]);
				argi += 2;
			}
			else if(param == "startscan")
			{
				startScan = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "refmode")
			{
				refMode = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "endscan")
			{
				endScan = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "minuv")
			{
				minUVW = atof(argv[argi+1]);
				argi += 2;
			}
			else if(param == "maxuv")
			{
				maxUVW = atof(argv[argi+1]);
				argi += 2;
			}
			else if(param == "scalar")
			{
				onlyScalar = true;
				++argi;
			}
			else if(param == "diag")
			{
				onlyDiag = true;
				++argi;
			}
			else if(param == "rhs")
			{
				rhsSolutionFile = argv[argi+1];
				argi += 2;
			}
			else if(param == "rotation")
			{
				onlyRotation = true;
				argi++;
			}
			else if(param == "quiet")
			{
				doQuiet = true;
				++argi;
			}
			else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
		}
		
		if(argc <= argi + 1) throw std::runtime_error("Incorrect parameters");

		std::cout << "Min Accuracy: " << minAccuracy << "; Stopping Accuracy: " << stopAccuracy << "\n";
		
		const char *msName = argv[argi];
		const char *outName = argv[argi+1];
		casacore::MeasurementSet ms(msName);
		
		Calibrator calibrator(ms, threadCount);
		calibrator.SetNIter(niter);
		calibrator.SetAccuracy(minAccuracy, stopAccuracy);
		calibrator.SetModelFilename(modelFile);
		calibrator.SetSolutionInterval(solutionInterval);
		calibrator.SetSolutionChannels(solutionChannels);
		calibrator.SetAbsMem(absMem);
		if(hasInterval)
			calibrator.SetInterval(intervalStart, intervalEnd);
		calibrator.SetStartScan(startScan);
		calibrator.SetEndScan(endScan);
		calibrator.SetMinUVW(minUVW);
		calibrator.SetMaxUVW(maxUVW);
		calibrator.SetRefMode(refMode);
		calibrator.SetOnlyScalar(onlyScalar);
		calibrator.SetOnlyDiag(onlyDiag);
		calibrator.SetRHSSolutionFile(rhsSolutionFile);
		calibrator.SetOnlyRotation(onlyRotation);
		calibrator.SetSolutionOutputFilename(outName);
		calibrator.SetVerbose(!doQuiet);
		calibrator.SetSavePlotFiles(savePlotFiles);
		calibrator.SetDataColumnName(dataColumnName);
		calibrator.SetWeightColumnName(weightColumnName);
		if(savePlotFiles)
		{
			calibrator.SetPlotFilenames(plotPhaseFile, plotGainFile);
		}
		calibrator.Perform();
	}
}
