#include "mspredicter.h"

#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

MSPredicter::~MSPredicter()
{
//    if(_readThread != 0)
//        _readThread->join();
    clearBuffers();
}

void MSPredicter::clearBuffers()
{
    for(std::vector<std::complex<double>*>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
        delete[] *i;
}

void MSPredicter::Start(bool reportSources)
{
    boost::mutex::scoped_lock lock(_mutex);
    if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
    
    casacore::ROArrayColumn<casacore::Complex> dataColumn(_ms, _ms.columnName(casacore::MSMainEnums::DATA));
    casacore::IPosition dataShape = dataColumn.shape(0);
    unsigned polarizationCount = dataShape[0];
    if(polarizationCount != 4)
        throw std::runtime_error("Expecting MS with 4 polarizations");
    
    _bandData.reset(new BandData(_ms.spectralWindow()));
    if(_endChannel == 0)
    {
        _startChannel = 0;
        _endChannel = _bandData->ChannelCount();
    }
    _channelCount = _bandData->ChannelCount();
    
    // By setting the time beforehand, we don't waste time calculating a time step we don't need.
    if(!_useModelColumn)
    {
        casacore::ROScalarColumn<int> fieldid(_ms, _ms.columnName(casacore::MSMainEnums::FIELD_ID));
        int field_id = fieldid(0);
        std::cout << "Using FIELD_ID = " << fieldid(0) << "\n";
        casacore::MSField fieldTable = _ms.field();
        casacore::ROArrayColumn<double> phaseDirColumn(fieldTable, fieldTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
//        if(phaseDirColumn.nrow() != 1)
//            throw std::runtime_error("Field table nrow != 1");
        casacore::Array<double> phaseDir = phaseDirColumn(field_id);
        casacore::Array<double>::const_iterator phaseDirIter = phaseDir.begin();
        long double phaseCentreRA = *phaseDirIter; ++phaseDirIter;
        long double phaseCentreDec = *phaseDirIter;

        casacore::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
        casacore::MEpoch startTime = timeColumn(_startRow);
        _predicter.reset(new Predicter(phaseCentreRA, phaseCentreDec, _bandData->LowestFrequency(), _bandData->HighestFrequency(), _bandData->ChannelCount()));
        _predicter->Initialize(_model);
        if(reportSources)
            _predicter->ReportSources(_model);
    }
    
    // Create buffers
    if(_buffers.empty())
    {
        _buffers.resize(_laneSize);
        for(size_t i=0; i!=_laneSize; ++i)
        {
            RowData rowData;
            _buffers[i] = new std::complex<double>[_channelCount*4];
            rowData.modelData = _buffers[i];
            _availableBufferLane.write(rowData);
        }
    }

    // Start all threads
    if(!_useModelColumn)
        _workThreadGroup.clear();
    _readThread.reset(new boost::thread(&MSPredicter::ReadThreadFunc, this));
}
    
void MSPredicter::ReadThreadFunc()
{
    size_t actualThreadCount = _threadCount;
    if(!_useModelColumn)
    {
        if(_model.SourceCount() == 0)
            actualThreadCount = 1;
        for(size_t i=0; i!=actualThreadCount; ++i)
        {
            _workThreadGroup.emplace_back(&MSPredicter::PredictThreadFunc, this);
        }
    }
    
    boost::mutex::scoped_lock lock(_mutex);

    casacore::ROScalarColumn<int> ant1Column(_ms, _ms.columnName(casacore::MSMainEnums::ANTENNA1));
    casacore::ROScalarColumn<int> ant2Column(_ms, _ms.columnName(casacore::MSMainEnums::ANTENNA2));
    casacore::ROScalarColumn<int> scanColumn(_ms, _ms.columnName(casacore::MSMainEnums::SCAN_NUMBER));
    casacore::ROArrayColumn<double> uvwColumn(_ms, _ms.columnName(casacore::MSMainEnums::UVW));
    casacore::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
    std::unique_ptr<casacore::ROArrayColumn<casacore::Complex> > modelColumn;
    
    RowData rowData;
    
    casacore::Array<casacore::Complex> modelData;
    if(_useModelColumn)
    {
        modelColumn.reset(new casacore::ROArrayColumn<casacore::Complex>(_ms, _ms.columnName(casacore::MSMainEnums::MODEL_DATA)));
        modelData = casacore::Array<casacore::Complex>(modelColumn->shape(_startRow));
    }
    
    size_t timeIndex = 0;
    casacore::MEpoch previousTime = timeColumn(_startRow);
    for(size_t rowIndex=_startRow; rowIndex!=_endRow; ++rowIndex)
    {
        size_t a1 = ant1Column(rowIndex);
        size_t a2 = ant2Column(rowIndex);
        int sn = scanColumn(rowIndex);
        casacore::MEpoch time = timeColumn(rowIndex);
         // Only calibrate selected scans (if selected)
        if(_startScan != -1 && sn < _startScan)
            continue;
        if(_endScan != -1 && sn > _endScan)
            continue;
//        std::cout << "Processing " << a1 << "-" << a2 << " scan=" << sn << "\n";
        if((a1 != a2))
        {
            casacore::Array<double> uvwArray = uvwColumn(rowIndex);
            casacore::Array<double>::const_contiter uvwI = uvwArray.cbegin();
            double u = *uvwI; ++uvwI;
            double v = *uvwI; ++uvwI;
            double w = *uvwI;
            if(_useModelColumn)
                modelColumn->get(rowIndex, modelData);
            lock.unlock();
            
            if(time.getValue() != previousTime.getValue())
            {
                ++timeIndex;
                previousTime = time;
            }
            
            _availableBufferLane.read(rowData);
            rowData.u = u;
            rowData.v = v;
            rowData.w = w;
            rowData.rowIndex = rowIndex;
            rowData.a1 = a1;
            rowData.a2 = a2;
            rowData.timeIndex = timeIndex;
            if(_useModelColumn)
            {
                std::complex<double> *outptr = rowData.modelData;
                casacore::Complex* inptr = modelData.cbegin();
                for(size_t ch=0; ch!=_channelCount*4; ++ch)
                {
                    *outptr = *inptr;
                    ++outptr; ++inptr;
                }
                _outputLane.write(rowData);
            }
            else {
                _workLane.write(rowData);
            }
            lock.lock();
        }
    }
    
    lock.unlock();
    if(!_useModelColumn)
    {
        _workLane.write_end();
        for(boost::thread& t : _workThreadGroup)
            t.join();
    }
    _outputLane.write_end();
}

void MSPredicter::PredictThreadFunc()
{
    RowData rowData;
    while(_workLane.read(rowData))
    {
        std::complex<double> *valIter = rowData.modelData + 4*_startChannel;
        for(size_t ch=_startChannel; ch!=_endChannel; ++ch)
        {
            double lambda = _bandData->ChannelWavelength(ch);
            _predicter->Predict4(valIter, _model, rowData.u/lambda, rowData.v/lambda, rowData.w/lambda, ch, rowData.a1, rowData.a2);
            valIter += 4;
        }
        _outputLane.write(rowData);
    }
}
