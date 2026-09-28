#ifndef ROW_ORDER_H
#define ROW_ORDER_H

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ScalarColumn.h>

#include <sstream>
#include <stdexcept>

/**
 * How the cross-correlation rows of a measurement set order their antennas.
 * Normal: ANTENNA1 < ANTENNA2 (e.g. ASKAP). Reversed: ANTENNA1 > ANTENNA2 on
 * every row (e.g. SKA-Low).
 *
 * The solutions of a reversed measurement set are stored as the complex
 * conjugate of the physical solutions, which is what calibrate has always
 * written for such data, so that existing solution files stay valid.
 */
enum class RowOrder { Normal, Reversed };

/**
 * Determine the row order of a measurement set. Autocorrelations are ignored.
 * A measurement set that mixes both orders is not supported and throws.
 */
inline RowOrder GetRowOrder(const casacore::MeasurementSet& ms)
{
	casacore::ROScalarColumn<int> ant1Column(ms, ms.columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> ant2Column(ms, ms.columnName(casacore::MSMainEnums::ANTENNA2));
	const casacore::Vector<int> ant1 = ant1Column.getColumn();
	const casacore::Vector<int> ant2 = ant2Column.getColumn();
	size_t normalRows = 0, reversedRows = 0;
	for(size_t row=0; row!=ant1.size(); ++row)
	{
		if(ant1[row] < ant2[row])
			++normalRows;
		else if(ant1[row] > ant2[row])
			++reversedRows;
	}
	if(normalRows != 0 && reversedRows != 0)
	{
		std::ostringstream s;
		s << "Measurement set mixes row orders: " << normalRows << " rows have ANTENNA1 < ANTENNA2 and "
			<< reversedRows << " rows have ANTENNA1 > ANTENNA2. Only one order per measurement set is supported.";
		throw std::runtime_error(s.str());
	}
	return reversedRows != 0 ? RowOrder::Reversed : RowOrder::Normal;
}

#endif
