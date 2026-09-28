#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include "solutionfile.h"

#include "model/model.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <boost/thread/mutex.hpp>

#include <queue>

class Calibrator
{
public:
	Calibrator(casacore::MeasurementSet& ms, size_t threadCount);
	
	void Perform();

//	void SetBeamOnSource(bool beamOnSource)
//	{
//		_beamOnSource = beamOnSource;
//	}
	
//	void SetApplyBeam(bool applyBeam)
//	{
//		_applyBeam = applyBeam;
//	}
	
	void SetOnlyScalar(bool onlyScalar)
	{
		_onlyScalar = onlyScalar;
	}
	
	void SetOnlyDiag(bool onlyDiag)
	{
		_onlyDiag = onlyDiag;
	}
	
	void SetOnlyRotation(bool onlyRotation)
	{
		_onlyRotation = onlyRotation;
	}
	
	void SetModelFilename(const std::string& modelFilename)
	{
		_modelFilename = modelFilename;
	}
	
	void SetSolutionOutputFilename(const std::string& solutionOutputFilename)
	{
		_solutionFilename = solutionOutputFilename;
	}

	void SetModel(const Model& model)
	{
		_model = model;
	}

	void SetDataColumnName(const std::string& dataColumnName)
	{
		_dataColumnName = dataColumnName;
	}
	
	/**
	 * Column to read the weights from. An empty name selects automatically:
	 * WEIGHT_SPECTRUM, then WEIGHT, then SIGMA_SPECTRUM, then SIGMA. The
	 * SIGMA columns are converted to weights as 1/sigma^2.
	 */
	void SetWeightColumnName(const std::string& weightColumnName)
	{
		_weightColumnName = weightColumnName;
	}
	
	void SetRHSSolutionFile(const std::string& rhsSolutionFile)
	{
		_rhsSolutionFilename = rhsSolutionFile;
	}
	
	void SetNIter(size_t nIter)
	{
		_nIter = nIter;
	}
	
	void SetAccuracy(double minAccuracy, double stoppingAccuracy)
	{
		_minAccuracy = minAccuracy;
		_stoppingAccuracy = stoppingAccuracy;
	}
	
	void SetMinUVW(double minUVW)
	{
		_minUVW = minUVW;
	}
	
	void SetMaxUVW(double maxUVW)
	{
		_maxUVW = maxUVW;
	}
	
	void SetRefMode(size_t refMode)
	{
		_refMode = refMode;
	}
	
	void SetSolutionInterval(size_t solutionInterval)
	{
		_solutionInterval = solutionInterval;
	}
	
	/** Solve one solution per block of this many channels (1 = every channel). */
	void SetSolutionChannels(size_t solutionChannels)
	{
		_solutionChannels = solutionChannels == 0 ? 1 : solutionChannels;
	}
	
	/** Memory to plan for, in GB. 0 uses the physical memory of the machine. */
	void SetAbsMem(double absMem)
	{
		_absMem = absMem;
	}
	
	/** Only use timesteps [start, end) of the measurement set. */
	void SetInterval(size_t start, size_t end)
	{
		_hasInterval = true;
		_intervalStart = start;
		_intervalEnd = end;
	}
	
	void SetStartScan(int startScan)
	{
		_startScan = startScan;
	}
	
	void SetEndScan(int endScan)
	{
		_endScan = endScan;
	}
	
	void SetVerbose(bool verbose)
	{
		_verbose = verbose;
	}
	
	void SetSavePlotFiles(bool savePlotFiles)
	{
		_savePlotFiles = savePlotFiles;
	}
	
	void SetPlotFilenames(const std::string& phasePlotFilename, const std::string& gainPlotFilename)
	{
		_phasePlotFilename = phasePlotFilename;
		_gainPlotFilename = gainPlotFilename;
	}
	
	SolutionFile& GetSolutionFile() {
		return _solutionFile;
	}
private:
	struct ThreadData
	{
		ThreadData() { }
		
		boost::mutex *mutex;
		std::queue<size_t> *tasks;
		std::vector<class CalibrationMethod*> *calMethods;
        size_t index;
	};
	
	void threadFunction(ThreadData data);
	
	/** Selects which column the weights are read from, see SetWeightColumnName(). */
	casacore::MSMainEnums::PredefinedColumns selectWeightColumn() const;

	casacore::MeasurementSet _ms;
	SolutionFile _solutionFile;
	std::string _modelFilename, _solutionFilename, _rhsSolutionFilename, _dataColumnName, _weightColumnName;
	Model _model;
	double _minAccuracy, _stoppingAccuracy;
	size_t _nIter, _solutionInterval, _solutionChannels;
    int _startScan, _endScan;
    size_t _threadCount, _refMode;
	bool _onlyScalar, _onlyDiag, _onlyRotation;
	double _minUVW, _maxUVW;
	bool _savePlotFiles, _saveFaradayPlotFiles, _saveCrossTermsPlotFile, _verbose;
	double _absMem;
	bool _hasInterval;
	size_t _intervalStart, _intervalEnd;
	std::string _phasePlotFilename, _gainPlotFilename, _faradayPlotFilename, _crossTermsPlotFilename;
};

#endif
