#include "calibrator.h"

#include "calibrationmethod.h"
#include "banddata.h"
#include "matrix2x2.h"
#include "progressbar.h"

#include "mspredicter.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <boost/thread/thread.hpp>

#include <fstream>
#include <stdexcept>
#include <memory>
#include <queue>
#include <complex>

Calibrator::Calibrator(casacore::MeasurementSet& ms, size_t threadCount) :
    _ms(ms),
    _dataColumnName("DATA"),
    _minAccuracy(CalibrationMethod::DefaultMinAccuracy()),
    _stoppingAccuracy(CalibrationMethod::DefaultStoppingAccuracy()),
    _nIter(1000),
    _solutionInterval(0),
    _startScan(-1),
    _endScan(-1),
    _threadCount(threadCount),
    _refMode(0),
    _onlyScalar(false),
    _onlyDiag(false),
    _onlyRotation(false),
    _minUVW(-1.0),
    _maxUVW(-1.0),
    _savePlotFiles(false),
    _saveFaradayPlotFiles(false),
    _saveCrossTermsPlotFile(false),
    _verbose(false)
{
    std::cout << "Calibrator::Calibrator minAccuracy:" << _minAccuracy << "; stoppingAccuracy:" << _stoppingAccuracy << "\n";
}

void Calibrator::Perform()
{
    if(_verbose)
        std::cout << "Reading meta data... " << std::flush;
    
    /**
        * Read some meta data from the measurement set
        */
    casacore::MSAntenna aTable = _ms.antenna();
    size_t antennaCount = aTable.nrow();
    
    BandData bandData(_ms.spectralWindow());
    size_t channelCount = bandData.ChannelCount();
    if(channelCount == 0) throw std::runtime_error("No channels in set");
    if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
    
    typedef float num_t;
    typedef std::complex<num_t> complex_t;
    casacore::ROScalarColumn<double> timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
    casacore::ROArrayColumn<complex_t> dataColumn(_ms, _dataColumnName);
    casacore::ROArrayColumn<bool> flagColumn(_ms, _ms.columnName(casacore::MSMainEnums::FLAG));
    casacore::ROScalarColumn<int> scanColumn(_ms, _ms.columnName(casacore::MSMainEnums::SCAN_NUMBER)); // Check for scan number
    
    // Use SIGMA_SPECTRUM if available, otherwise use SIGMA
    casacore::MSMainEnums::PredefinedColumns weight_column = casacore::MSMainEnums::SIGMA_SPECTRUM;
    if(!_ms.tableDesc().isColumn("SIGMA_SPECTRUM"))
    {
        std::cout << "No SIGMA_SPECTRUM column found, using SIGMA instead" << std::endl;
        weight_column = casacore::MSMainEnums::SIGMA;
    }
    else
    {
        std::cout << "Using SIGMA_SPECTRUM column for weights" << std::endl;
    }
    casacore::ROArrayColumn<float> weightColumn(_ms, _ms.columnName(weight_column));
    
    casacore::IPosition weightShape = weightColumn.shape(0);
    casacore::IPosition dataShape = dataColumn.shape(0);
    
    unsigned polarizationCount = dataShape[0];

    if(polarizationCount != 4)
        throw std::runtime_error("Pol count in MS != 4");

    if(_verbose)
        std::cout << "DONE\nCounting timesteps... " << std::flush;
    double time = -1.0;
    std::vector<size_t> timestepRows;
    for(size_t rowIndex=0;rowIndex!=_ms.nrow();++rowIndex)
    {
        if((_startScan == -1 || scanColumn(rowIndex) >= _startScan) && (_endScan == -1 || scanColumn(rowIndex) <= _endScan))
        if(timeColumn(rowIndex) != time)
        {
            timestepRows.push_back(rowIndex);
            time = timeColumn(rowIndex);
        }
    }
    size_t timestepCount = timestepRows.size();
    timestepRows.push_back(_ms.nrow());
    size_t intervalCount = (_solutionInterval!=0) ? (timestepCount + _solutionInterval - 1) / _solutionInterval : 1;
    if(_verbose)
     std::cout << "DONE (" << timestepCount << " timesteps, " << intervalCount << " intervals)\n";

    if(!_modelFilename.empty()) {
        if(_verbose)
            std::cout << "Reading model... " << std::flush;
        _model = Model(_modelFilename.c_str());
        if(_verbose)
            std::cout << "DONE\n";
    }

    _solutionFile.SetAntennaCount(antennaCount);
    _solutionFile.SetChannelCount(channelCount);
    _solutionFile.SetIntervalCount(intervalCount);
    _solutionFile.SetPolarizationCount(4);
    if(_solutionFilename.empty())
        _solutionFile.OpenInMemory();
    else
        _solutionFile.OpenForWriting(_solutionFilename.c_str());

    long int pageCount = sysconf(_SC_PHYS_PAGES);
    long int pageSize = sysconf(_SC_PAGE_SIZE);
    int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
    double memSizeInGB = (double) memSize / (1024.0*1024.0*1024.0);
    size_t nBaselines = antennaCount * (antennaCount-1) / 2;
    
    for(size_t intervalIndex=0; intervalIndex!=intervalCount; ++intervalIndex)
    {
        std::cout << " >>> INTERVAL " << (intervalIndex+1) << " / " << intervalCount << " <<<\n";
        size_t intervalTimestepStart = (intervalIndex*timestepCount) / intervalCount;
        size_t intervalTimestepEnd = ((intervalIndex+1)*timestepCount) / intervalCount;
        size_t intervalRowStart = timestepRows[intervalTimestepStart];
        size_t intervalRowEnd = timestepRows[intervalTimestepEnd];
        size_t timestepsInInterval = intervalTimestepEnd - intervalTimestepStart;
            
        size_t samplesPerChannel = nBaselines * timestepsInInterval * 4;
        // 2 for complex data, 2 for complex model, 1 for weights
        double memPerChannel = samplesPerChannel * 5 * sizeof(double);
        if(_verbose)
        {
            std::cout << "Will use " << _threadCount << " cores.\n";
            std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";
            std::cout << "One channel takes " << round(memPerChannel*10.0/(1024*1024))/10.0 << " MB of mem.\n";
        }
        size_t channelsPerPass = memSize / memPerChannel;
        if(channelsPerPass > channelCount)
            channelsPerPass = channelCount;
        if(channelsPerPass == 0) {
            if(_verbose)
                std::cout << "WARNING: NOT ENOUGH MEMORY FOR EVEN ONE CHANNEL, expect very bad performance.\n";
            channelsPerPass = 1;
        }
        size_t passCount = (channelCount + channelsPerPass - 1) / channelsPerPass;
        if(_verbose)
            std::cout << "Number of channels that fit in memory: " << channelsPerPass << " (" << passCount << " passes)\n";
        
        for(size_t pass=0; pass!=passCount; ++pass) {
            size_t startChannel = (channelCount * pass) / passCount;
            size_t endChannel = (channelCount * (pass+1)) / passCount;
            size_t partChannelCount = endChannel - startChannel;
            std::cout << "pass = " << pass << "; start=" << startChannel << "; end=" << endChannel << "; count=" << partChannelCount << "\n";
            BandData partBandData(bandData, startChannel, endChannel);

            std::vector<CalibrationMethod*> calMethods(partChannelCount);
            for(size_t ch=0; ch!=partChannelCount; ++ch)
            {
                calMethods[ch] = new CalibrationMethod(1, antennaCount, timestepsInInterval);
                calMethods[ch]->SetOnlySolveScalar(_onlyScalar);
                calMethods[ch]->SetOnlySolveDiag(_onlyDiag);
                calMethods[ch]->SetOnlySolveRotation(_onlyRotation);
            }
            std::unique_ptr<MSPredicter> predicter;
			std::unique_ptr<ProgressBar> progress;
            if(_model.Empty()) {
                if(_verbose)
                    progress.reset(new ProgressBar("Reading data and model column"));
                predicter.reset(new MSPredicter(_ms, _threadCount));
            }
            else {
                predicter.reset(new MSPredicter(_ms, _threadCount, _model));
                if(_verbose)
                    progress.reset(new ProgressBar("Reading data & predicting model"));
            }
            predicter->SetStartRow(intervalRowStart);
            predicter->SetEndRow(intervalRowEnd);
            predicter->SetStartScan(_startScan);
            predicter->SetEndScan(_endScan);

            std::vector<std::complex<double> > modelValues(4 * channelCount);
            casacore::Array<complex_t> data(dataShape);
            casacore::Array<float> weights(dataShape);
            casacore::Array<bool> flags(dataShape);
            time = timeColumn(intervalRowStart);
            size_t selectedCount = 0, notSelected = 0, previousTime = 0;
            MSPredicter::RowData rowData;

            predicter->Start(_verbose);
            while(predicter->GetNextRow(rowData))
            {
                size_t rowIndex = rowData.rowIndex;
                // Cross correlation?
                size_t antenna1 = rowData.a1, antenna2 = rowData.a2;
                if(previousTime < rowData.timeIndex)
                {
                    previousTime = rowData.timeIndex;
                    if(_verbose)
                        std::cout << '.' << std::flush;
                }
                // Check reference antenna usage mode. If _refMode == 1 only include baselines to reference antenna 0.
                bool useRow = true;
                if(_refMode == 1)
                    if((antenna1 != 0) && (antenna2 != 0))
                        useRow = false;
                // If _refMode == 2 only include baselines not including the reference antenna 0.
                if(_refMode == 2)
                    if((antenna1 == 0) || (antenna2 == 0))
                        useRow = false;
                if(antenna1 == antenna2)
                    useRow = false;
                if(useRow == true)
                {
                    boost::mutex::scoped_lock lock(predicter->IOMutex());
                    dataColumn.get(rowIndex, data);
                    // If using SIGMA_SPECTRUM then have a weight per channel
                    if(weight_column == casacore::MSMainEnums::SIGMA_SPECTRUM)
                        weightColumn.get(rowIndex, weights);
                    else
                    {
                        // If using SIGMA then have to copy weight across for each channel.
                        casacore::Array<float> vweights(weightShape);
                        weightColumn.get(rowIndex, vweights);
                        float *weightscPtr = weights.cbegin();
                        float *vweightsPtr = vweights.cbegin();
                        for(size_t ch = 0; ch!=partChannelCount; ++ch)
                        {
                            size_t chIndex = (ch + startChannel) * 4;
                            for(size_t p=0; p!=4; ++p)
                            {
                                weightscPtr[chIndex+p] = vweightsPtr[p];
                            }
                        }
                    }
                    flagColumn.get(rowIndex, flags);
                    lock.unlock();

                    std::complex<float> *dataPtr = data.cbegin();
                    float *weightsPtr = weights.cbegin();
                    bool *flagPtr = flags.cbegin();

                    double u = rowData.u;
                    double v = rowData.v;
                    double w = rowData.w;

                    bool selected = true;
                    if((_minUVW >= 0.0 && (u*u + v*v + w*w < _minUVW*_minUVW))||(_maxUVW>=0.0 && (u*u + v*v + w*w > _maxUVW*_maxUVW)))
                        selected = false;
                    if(selected)
                        selectedCount++;
                    else
                        notSelected++;

                    for(size_t ch = 0; ch!=partChannelCount; ++ch)
                    {
                        size_t chIndex = (ch + startChannel) * 4;
                        for(size_t p=0; p!=4; ++p)
                        {
                            modelValues[chIndex+p] = rowData.modelData[chIndex+p];
                            if(flagPtr[chIndex+p] || !selected)
                            {
//                                if(weight_column == casa::MSMainEnums::SIGMA_SPECTRUM)
//                                weightsPtr[chIndex+p] = 0.0;
                                // If any one of the instrumental polarisations is flagged then flag everything.
                                weightsPtr[chIndex] = 0.0;
                                weightsPtr[chIndex+1] = 0.0;
                                weightsPtr[chIndex+2] = 0.0;
                                weightsPtr[chIndex+3] = 0.0;
                            }
                        }
                        calMethods[ch]->AddData(&dataPtr[chIndex], &weightsPtr[chIndex], &modelValues[chIndex], antenna1, antenna2, rowData.timeIndex);
                    }
                }
                
                predicter->FinishRow(rowData);
            }
            if(_verbose)
                std::cout << "DONE (" << selectedCount<< "/" << (selectedCount+notSelected) << " rows selected)\nCalibrating...\n";
        
            std::queue<size_t> tasks;
            for(size_t ch=0; ch!=partChannelCount; ++ch)
                tasks.push(ch);
            boost::thread_group threadGroup;
            boost::mutex mutex;
            for(size_t i=0; i!=_threadCount; ++i)
            {
                ThreadData threadData;
                threadData.mutex = &mutex;
                threadData.tasks = &tasks;
                threadData.calMethods = &calMethods;
                threadData.index = i;
                threadGroup.add_thread(new boost::thread(&Calibrator::threadFunction, this, threadData));
            }
            std::cout << "Created thread(s)\n";
            threadGroup.join_all();
            std::cout << "All threads finished\n";

            // Save solutions
            for(size_t ant=0; ant!=antennaCount; ++ant)
            {
                for(size_t ch=0; ch!=partChannelCount; ++ch)
                {
                    std::complex<double> val[4];
                    for(size_t p=0; p!=4; ++p)
                        val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
                    Matrix2x2::Invert(val);
                    
                    for(size_t p=0; p!=4; ++p)
                    {
                        _solutionFile.WriteSolution(val[p], intervalIndex, ant, ch+startChannel, p);
                    }
                }
            }
            
            if(_savePlotFiles)
            {
                std::ofstream phasePlotStream(_phasePlotFilename.c_str()), gainPlotStream(_gainPlotFilename.c_str());
                phasePlotStream << antennaCount << ' ' << partChannelCount << " 4\n";
                gainPlotStream << antennaCount << ' ' << partChannelCount << " 4\n";
                
                for(size_t ch=0; ch!=partChannelCount; ++ch)
                {
                    phasePlotStream << (ch+startChannel) << '\t';
                    gainPlotStream << (ch+startChannel) << '\t';
                    
                    for(size_t p=0; p!=4; ++p)
                    {
                        for(size_t ant=0; ant!=antennaCount; ++ant)
                        {
                            std::complex<double> val[4];
                            for(size_t p2=0; p2!=4; ++p2)
                                val[p2] = calMethods[ch]->JonesSolution(ant, 0, p2);
                            Matrix2x2::Invert(val);
                    
                            double s1, s2;
                            Matrix2x2::SingularValues(val, s1, s2);
                            switch(p)
                            {
                            case 0: gainPlotStream << '\t' << s1; break;
                            case 1: case 2: gainPlotStream << '\t' << 0.0; break;
                            case 3: gainPlotStream << '\t' << s2;
                            }
                            phasePlotStream << '\t' << std::arg(val[p]);
                        }
                    }
                    phasePlotStream << '\n';
                    gainPlotStream << '\n';
                }
            }
            
            if(_saveFaradayPlotFiles)
            {
                std::ofstream faradayPlotStream(_faradayPlotFilename.c_str());
                
                for(size_t ch=0; ch!=partChannelCount; ++ch)
                {
                    faradayPlotStream << (ch+startChannel) << '\t';
                    
                    for(size_t ant=0; ant!=antennaCount; ++ant)
                    {
                        std::complex<double> val[4];
                        for(size_t p=0; p!=4; ++p)
                            val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
                
                        faradayPlotStream << '\t' << -Matrix2x2::RotationAngle(val);
                    }
                    faradayPlotStream << '\n';
                }
            }
            
            if(_saveCrossTermsPlotFile)
            {
                std::ofstream crossTermPlotStream(_crossTermsPlotFilename.c_str());
                
                for(size_t ch=0; ch!=partChannelCount; ++ch)
                {
                    crossTermPlotStream << (ch+startChannel) << '\t';
                    
                    for(size_t ant=0; ant!=antennaCount; ++ant)
                    {
                        std::complex<double> val[4];
                        for(size_t p=0; p!=4; ++p)
                            val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
                        Matrix2x2::Invert(val);
                        double totalPower = std::abs(val[0]) + std::abs(val[1]) + std::abs(val[2]) + std::abs(val[3]);
                        crossTermPlotStream << '\t' << (std::abs(val[1]) + std::abs(val[2]))*100.0/totalPower;
                    }
                    crossTermPlotStream << '\n';
                }
            }
            
            for(size_t ch=0; ch!=partChannelCount; ++ch)
                delete calMethods[ch];
        }
    }

}

void Calibrator::threadFunction(ThreadData data)
{
    boost::mutex::scoped_lock lock(*data.mutex);
//    size_t lastSuccessfulChannel = data.tasks->front();
    while(!data.tasks->empty()) {
        size_t taskIndex = data.tasks->front();
        data.tasks->pop();
        lock.unlock();
// *** len067 : Took out this section as it seems to cause subsequent solutions to fail in some cases. ***
//        if(lastSuccessfulChannel != taskIndex)
//        {
//            std::cout << "Thread " << data.index << " initSolutions last=" << lastSuccessfulChannel << " taskIndex=" << taskIndex << "\n";
//            (*(data.calMethods))[taskIndex]->InitSolutions(*(*(data.calMethods))[lastSuccessfulChannel]);
//        }
        size_t iters = _nIter;
        double limit = _stoppingAccuracy;
        (*(data.calMethods))[taskIndex]->Execute(limit, iters);
        if((iters >= _nIter || !std::isfinite(limit)) && !(*(data.calMethods))[taskIndex]->OnlySolveRotation())
        {
            std::cout << "Recalculating channel " << taskIndex << " (accuracy=" << limit << ").\n";
            (*(data.calMethods))[taskIndex]->InitSolutionsToUnity();
            iters = _nIter;
            limit = _stoppingAccuracy;
            (*(data.calMethods))[taskIndex]->Execute(limit, iters);

            if((iters >= _nIter && limit > _minAccuracy) || !std::isfinite(limit))
            {
                std::cout << "Channel " << taskIndex << " did not converge (accuracy=" << limit << "), setting gains to NaN.\n";
                (*(data.calMethods))[taskIndex]->InitSolutionsToNaN();
            }
            else {
                if(iters >= _nIter && limit > _stoppingAccuracy)
                {
                    std::cout << "Channel " << taskIndex << " converged (accuracy=" << limit << ") but did not reach stopping accuracy.\n";
                }
//                lastSuccessfulChannel = taskIndex;
            }
        }
        else {
//            lastSuccessfulChannel = taskIndex;
        }
        lock.lock();
//        if(taskIndex>200)
//        {
//            if(_verbose)
//                std::cout << "Current value of Jones matrix for ant 0, ch " << taskIndex << ":\n"
//            << CalibrationMethod::MatrixToString(& (*(data.calMethods))[taskIndex]->JonesSolution(1, 0, 0));
//        }
    
        if(_verbose)
            std::cout << "Thread " << data.index << " finished calibrating channel " << taskIndex << " in " << iters << " iterations, precision=" << limit << ".\n";
    }
}

