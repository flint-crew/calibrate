#include "solutionfile.h"
#include "solutionapplier.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "stdint.h"

using namespace casacore;

int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std::cout << "Usage: applysolutions [-datacolumn <name>] [-gflag <solutions-flag-file.txt>] [-startscan <scan>] [-endscan <scan>] [-interval <start timestep> <end timestep>] [-copy/-nocopy] [-s xx xy yx yy] <ms> <gains-bin-file>\n"
				"Will apply the found solution matrices.\n"
				"Options:\n"
				"  -copy/-nocopy Don't(/do) alter the original DATA column but store the corrected data in the CORRECTED_DATA (this is std CASA behaviour)\n"
				"    default: -copy\n"
				"  -interval Only correct timesteps start to end-1 (counting from 0), as calibrate -interval\n"
				"Solutions for blocks of channels (calibrate -ch) are applied to every channel of their block.\n";
    } else {
		size_t argi = 1;
		int startScan = -1, endScan = -1;
		size_t intervalStart = 0, intervalEnd = 0;
		bool hasInterval = false;
		double xx=0.0, xy=0.0, yx=0.0, yy=0.0;
		bool preset = false, copyData = true;
		std::string dataColumnName = "DATA";
		while(argv[argi][0] == '-')
		{
			std::string p(&argv[argi][1]);
			if(p == "s")
			{
				xx = atof(argv[argi+1]);
				xy = atof(argv[argi+2]);
				yx = atof(argv[argi+3]);
				yy = atof(argv[argi+4]);
				argi += 4;
				preset = true;
			}
			else if(p == "startscan")
			{
				startScan = atoi(argv[argi+1]);
				++argi;
			}
			else if(p == "endscan")
			{
				endScan = atoi(argv[argi+1]);
				++argi;
			}
			else if(p == "interval")
			{
				hasInterval = true;
				intervalStart = atoi(argv[argi+1]);
				intervalEnd = atoi(argv[argi+2]);
				argi += 2;
			}
			else if(p == "datacolumn")
			{
				dataColumnName = argv[argi+1];
				++argi;
			}
			else if(p == "copy")
			{
				copyData = true;
			}
			else if(p == "nocopy")
			{
				copyData = false;
			}
			else throw std::runtime_error("What?");
			++argi;
		}
		MeasurementSet ms(argv[argi], Table::Update);
		SolutionFile solutionFile;
		solutionFile.OpenForReading(argv[argi+1]);
		
		SolutionApplier applier;
		if(preset)
			applier.SetPresets(xx, xy, yx, yy);
		applier.SetInputColumn(dataColumnName);
        applier.SetStartScan(startScan);
        applier.SetEndScan(endScan);
		if(hasInterval)
			applier.SetInterval(intervalStart, intervalEnd);
		if(copyData)
			applier.SetOutputColumn(casacore::MeasurementSet::columnName(casacore::MSMainEnums::CORRECTED_DATA));
		else
			applier.SetOutputColumn(dataColumnName);
		applier.Apply(ms, solutionFile);
	}
}
