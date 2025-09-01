#ifndef MS_PREDICTER_H
#define MS_PREDICTER_H

#include <ms/MeasurementSets/MeasurementSet.h>

#include "lane.h"
#include "predicter.h"
#include "banddata.h"

#include <boost/thread/thread.hpp>

#include "model/model.h"

#include <complex>
#include <memory>

class MSPredicter
{
public:
	struct RowData
	{
		RowData() : modelData(0)
		{	}
		
		std::complex<double> *modelData;
		size_t rowIndex, timeIndex, a1, a2;
		double u, v, w;
	};
	
	explicit MSPredicter(casacore::MeasurementSet &ms, size_t threadCount) :
		_ms(ms),
		_startChannel(0),
		_endChannel(0),
		_useModelColumn(true),
		_laneSize(64),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_startRow(0),
		_endRow(ms.nrow()),
        _startScan(-1),
        _endScan(-1),
    	_threadCount(threadCount)
	{ }

	MSPredicter(casacore::MeasurementSet &ms, size_t threadCount, const Model &model) :
		_ms(ms),
		_startChannel(0),
		_endChannel(0),
		_useModelColumn(false),
		_model(model),
		_laneSize(64),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_startRow(0),
		_endRow(ms.nrow()),
        _startScan(-1),
        _endScan(-1),
    	_threadCount(threadCount)
	{ }

	~MSPredicter();
	
	void Start(bool reportSources = false);
	
	bool GetNextRow(RowData& data)
	{
		return _outputLane.read(data);
	}
	void FinishRow(RowData& data)
	{
		_availableBufferLane.write(data);
	}
	
	boost::mutex &IOMutex() { return _mutex; }
	
	void SetStartRow(size_t startRow) { _startRow = startRow; }
	void SetEndRow(size_t endRow) { _endRow = endRow; }
	void SetStartScan(int startScan) { _startScan = startScan; }
	void SetEndScan(int endScan) { _endScan = endScan; }
private:
	void ReadThreadFunc();
	void PredictThreadFunc();
	void clearBuffers();
	
	casacore::MeasurementSet &_ms;
	size_t _channelCount;
	size_t _startChannel, _endChannel;
	bool _useModelColumn;
	
	Model _model;
	boost::mutex _mutex;
	
	const size_t _laneSize;
	ao::lane<RowData> _workLane, _outputLane, _availableBufferLane;
	
	std::unique_ptr<boost::thread> _readThread;
	std::vector<boost::thread> _workThreadGroup;
	std::unique_ptr<Predicter> _predicter;
	std::vector<std::complex<double>*> _buffers;
	std::unique_ptr<BandData> _bandData;
	std::string _solutionFile;
	size_t _startRow, _endRow;
    int _startScan, _endScan;
    size_t _threadCount;
};

#endif
