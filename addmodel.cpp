#include <iostream>
#include <stdexcept>
#include <math.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrColDesc.h>

#include "model/model.h"

#include "adder.h"

using namespace casacore;

typedef std::complex<float> complex_t;
    
// -m <a|s|c|z> = add, subtract, copy, zero
int main(int argc, char **argv)
{
    if(argc < 3)
    {
        std::cout << "Usage: addmodel [-usemodelcol] [-datacolumn <COLUMN>] [-m <a|s|c|z>] [-n <σ>] <model> <ms>\n"
            "Modify visibilities using a model. If -usemodelcol is specified then the MODEL_DATA column is used as the source model otherwise the specified component model file is used. The modification to use is defined with the mode switch(-m) where a=add model to visibilities (default), s=subtract model from visibilities, c=copy model to visibilities, z=zero visibilities.\n";
    } else {
        bool addNoise = false, useModelCol = false;
        double noiseSigma = 1.0;
        char cMode = 's'; // Default to model subtraction
        size_t argi = 1;
        size_t threadCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
        std::string dataColumnName = "DATA";
        
        while(argv[argi][0] == '-')
        {
            if(strcmp(argv[argi], "-m") == 0) { argi++; cMode = argv[argi][0]; }
            else if(strcmp(argv[argi], "-n") == 0) { addNoise=true; ++argi; noiseSigma = atof(argv[argi]); }
            else if(strcmp(argv[argi], "-datacolumn") == 0) { ++argi; dataColumnName=argv[argi]; }
            else if(strcmp(argv[argi], "-usemodelcol") == 0) { useModelCol=true; }
            else throw std::runtime_error("Invalid param");
            ++argi;
        }
        
        std::cout << "Opening measurement set... " << std::flush;
        MeasurementSet ms(argv[argi+1], Table::Update);
        std::cout << "DONE\n";
        
        // Get access to the DATA column
        casacore::ArrayColumn<complex_t> dataColumn(ms, "DATA");
        
        // For non-DATA columns, check to see if they exist
        if(dataColumnName.compare("DATA") != 0)
        {
            if(!ms.tableDesc().isColumn(dataColumnName))
            {
                std::cout << "Adding column '" << dataColumnName << "'... " << std::flush;
                casacore::IPosition shape = dataColumn.shape(0);
                casacore::ArrayColumnDesc<casacore::Complex> columnDesc(dataColumnName, shape);
                try {
                    ms.addColumn(columnDesc, "StandardStMan", true, true);
                } catch(std::exception& e)
                {
                    ms.addColumn(columnDesc, "StandardStMan", false, true);
                }
                std::cout << "DONE\n";
            }
        }
//        copyColumn.reset(new casacore::ArrayColumn<complex_t>(ms, data_column_name));
//        outputColumn = &*copyColumn;
        
        switch(cMode)
        {
            case 's':
                std::cout << "Subtracting model from " << dataColumnName << "\n";
                break;
            case 'a':
                std::cout << "Adding model to " << dataColumnName << "\n";
                break;
            case 'c':
                std::cout << "Copying model to " << dataColumnName << "\n";
                break;
            case 'z':
                std::cout << "Zeroing " << dataColumnName << "\n";
                break;
            default:
                throw std::runtime_error("Invalid mode specified");
        }

        if(useModelCol)
        {
            ProgressBar progress("Using MODEL_DATA column");
            casacore::ArrayColumn<casacore::Complex> dataColumn(ms, dataColumnName);
            casacore::ROArrayColumn<casacore::Complex> modelColumn(ms, casacore::MeasurementSet::columnName(casacore::MeasurementSet::MODEL_DATA));
            casacore::Array<casacore::Complex>
                dataArr(dataColumn.shape(0)), modelArr(modelColumn.shape(0));
            for(size_t i=0; i!=ms.nrow(); ++i)
            {
                dataColumn.get(i, dataArr);
                modelColumn.get(i, modelArr);
                
                casacore::Array<casacore::Complex>::contiter d=dataArr.cbegin();
                for(casacore::Array<casacore::Complex>::const_contiter m=modelArr.cbegin(); m!=modelArr.cend(); ++m, ++d)
                {
                    switch(cMode)
                    {
                        case 's':
                            *d -= *m;
                            break;
                        case 'a':
                            *d += *m;
                            break;
                        case 'c':
                            *d = *m;
                            break;
                        case 'z':
                            *d = (*m - *m);
                            break;
                    }
                }
                dataColumn.put(i, dataArr);
                progress.SetProgress(i+1, ms.nrow());
            }
        }
        else {
            std::cout << "Reading model... " << std::flush;
            Model model(argv[argi]);
            std::cout << "DONE\n";
        
            Adder adder(threadCount);
            adder.SetMode(cMode);
            adder.SetAddNoise(addNoise);
            adder.SetNoiseSigma(noiseSigma);
            adder.SetDataColumn(dataColumnName);
            adder.Add(ms, model);
        }
    }
}
