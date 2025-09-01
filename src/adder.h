#ifndef ADDER_H
#define ADDER_H

#include "banddata.h"
#include "predicter.h"
#include "mspredicter.h"
#include "progressbar.h"

#include "model/model.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <cmath>
#include <fstream>

class Adder
{
public:
    Adder(size_t threadCount) :
        _addNoise(false), _cmode('s'),
        _noiseSigma(1.0), _threadCount(threadCount), _dataColumn("DATA")
    {
    }
    
    void SetMode(char mode) { _cmode = mode; }
    void SetAddNoise(bool addNoise) { _addNoise = addNoise; }
    void SetNoiseSigma(double noiseSigma) { _noiseSigma = noiseSigma; }
    void SetDataColumn(const std::string& dataColumn) { _dataColumn = dataColumn; }
    
    void Add(casacore::MeasurementSet& ms, const Model& model)
    {
        /**
         * Read some meta data from the measurement set
         */
        BandData bandData(ms.spectralWindow());
        size_t channelCount = bandData.ChannelCount();
        
        typedef float num_t;
        typedef std::complex<num_t> complex_t;
        casacore::ArrayColumn<complex_t> dataColumn(ms, _dataColumn);
        
        casacore::IPosition dataShape = dataColumn.shape(0);
        unsigned polarizationCount = dataShape[0];
        
        MSPredicter predicter(ms, _threadCount, model);
        predicter.Start(true);
        
        if(_addNoise)
            std::cout << "Adding noise of " << _noiseSigma << " Jy.\n";
        
        /**
         * Subtract
         */
        std::ostringstream taskDesc;
        switch(_cmode)
        {
            case 's':
                taskDesc << "Subtracting model with " << model.SourceCount() << " sources";
                break;
            case 'a':
                taskDesc << "Adding model with " << model.SourceCount() << " sources";
                break;
            case 'c':
                taskDesc << "Copying model with " << model.SourceCount() << " sources";
                break;
            case 'z':
                taskDesc << "Zeroing";
                break;
        }
        ProgressBar progress(taskDesc.str());
        
        casacore::Array<complex_t> data(dataShape);
        MSPredicter::RowData rowData;
        while(predicter.GetNextRow(rowData))
        {
            size_t rowIndex = rowData.rowIndex;
            
            boost::mutex::scoped_lock lock(predicter.IOMutex());
            progress.SetProgress(rowIndex, ms.nrow());
            dataColumn.get(rowIndex, data);
            lock.unlock();
            
            casacore::Array<complex_t>::iterator dataPtr = data.begin();
            std::complex<double> *modelDataPtr = rowData.modelData;
            for(size_t ch=0; ch!=channelCount; ++ch)
            {
                for(size_t p=0; p!=polarizationCount; ++p)
                {
                    std::complex<double> predicted;
                    predicted = *modelDataPtr;
                    if(_addNoise)
                        addGausNoise(predicted, _noiseSigma);
                    switch(_cmode)
                    {
                        case 's':
                            *dataPtr -= predicted;
                            break;
                        case 'a':
                            *dataPtr += predicted;
                            break;
                        case 'c':
                            *dataPtr = predicted;
                            break;
                        case 'z':
                            *dataPtr -= *dataPtr;
                            break;
                    }
                    ++dataPtr;
                    ++modelDataPtr;
                }
            }
            
            lock.lock();
            dataColumn.put(rowIndex, data);
            lock.unlock();
            
            predicter.FinishRow(rowData);
        }
    }
private:
    template<typename T>
    void addGausNoise(std::complex<T> &value, double sigma)
    {
        long double x1, x2, w;

        do {
            long double r1 = (long double) rand() / (long double) RAND_MAX; 
            long double r2 = (long double) rand() / (long double) RAND_MAX; 
            x1 = 2.0 * r1 - 1.0;
            x2 = 2.0 * r2 - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = std::sqrt( (-2.0 * std::log( w ) ) / w ) * sigma;
        value += std::complex<T>(x1 * w, x2 * w);
    }
    
    bool _addNoise;
    char _cmode;
    double _noiseSigma;
    size_t _threadCount;
    std::string _dataColumn;
};

#endif
