# Clean up
rm -fr ChAvg_Beam0.calfield.ms ChAvg_Beam0.calfield.fits calfield.bin ChAvg_Beam0.cal1934.ms ChAvg_Beam0.cal1934.fits cal1934.bin calfield_bandpass.png cal1934_bandpass.png

# Calibrate against field
cp -r ChAvg_Beam0.ms ChAvg_Beam0.calfield.ms
./addfield.py ChAvg_Beam0.calfield.ms
./calibrate ChAvg_Beam0.calfield.ms calfield.bin
./applysolutions ChAvg_Beam0.calfield.ms calfield.bin
casapy -c exportfield.py
./plotbp.py calfield.bin calfield_bandpass.png

# Calibrate against B1934-638 only
cp -r ChAvg_Beam0.ms ChAvg_Beam0.cal1934.ms
./add1934.py ChAvg_Beam0.cal1934.ms
./calibrate ChAvg_Beam0.cal1934.ms cal1934.bin
./applysolutions ChAvg_Beam0.cal1934.ms cal1934.bin
casapy -c export1934.py
./plotbp.py cal1934.bin cal1934_bandpass.png





rm -fr ChAvg_Beam0.test* cal_0.*
cp -r ChAvg_Beam0.ms ChAvg_Beam0.test.ms
./flagcal.py ChAvg_Beam0.test.ms 0
./addfield.py ChAvg_Beam0.test.ms 0
./calibrate -minuv 50.0 -startscan 0 -endscan 0 ChAvg_Beam0.test.ms cal_0.bin
./applysolutions -startscan 0 -endscan 0 ChAvg_Beam0.test.ms cal_0.bin
./plotbp.py cal_0.bin cal_0.png
open cal_0.png
casapy -c exporttest.py

rm -fr ChAvg_Beam1.test* cal_1.*
cp -r ChAvg_Beam1.ms ChAvg_Beam1.test.ms
./flagcal.py ChAvg_Beam1.test.ms 1
./addfield.py ChAvg_Beam1.test.ms 1
./calibrate -minuv 50.0 -startscan 1 -endscan 1 ChAvg_Beam1.test.ms cal_1.bin
./applysolutions -startscan 1 -endscan 1 ChAvg_Beam1.test.ms cal_1.bin
./plotbp.py cal_1.bin cal_1.png
open cal_1.png

rm -fr ChAvg_Beam15.test* cal_15* 
cp -r ChAvg_Beam15.ms ChAvg_Beam15.test.ms
./flagcal.py ChAvg_Beam15.test.ms 15
./addfield.py ChAvg_Beam15.test.ms 15
./calibrate -minuv 50.0 -startscan 15 -endscan 15 ChAvg_Beam15.test.ms cal_15.bin
./applysolutions -startscan 15 -endscan 15 ChAvg_Beam15.test.ms cal_15.bin
./plotbp.py cal_15.bin cal_15.png
open cal_15.png

fits in=ChAvg_Beam0.test.fits op=uvin out=ChAvg_Beam0.test.mir
fits in=ChAvg_Beam0.test.mir op=uvout out=ChAvg_Beam0.testf.fits
splitc ChAvg_Beam0.testf.fits ChAvg_Beam0.testfs.fits

fits in=ChAvg_Beam1.test.fits op=uvin out=ChAvg_Beam1.test.mir
fits in=ChAvg_Beam1.test.mir op=uvout out=ChAvg_Beam1.testf.fits
splitc ChAvg_Beam1.testf.fits ChAvg_Beam1.testfs.fits

fits in=ChAvg_Beam15.test.fits op=uvin out=ChAvg_Beam15.test.mir
fits in=ChAvg_Beam15.test.mir op=uvout out=ChAvg_Beam15.testf.fits
splitc ChAvg_Beam15.testf.fits ChAvg_Beam15.testfs.fits

difmap
obs ChAvg_Beam0.testfs.fits
select xx,1,196,201,288
mapsize 4096
rmodel var.mod
mapl


Iteration 04: Reduced Chi-squared=6.4256975  Degrees of Freedom=1628738
! Flux (Jy) Radius (mas)  Theta (deg)  Major (mas)  Axial ratio   Phi (deg) T \
! Freq (Hz)     SpecIndex
   14.1210v     21.4446v     68.9284v     0.00000      1.00000     0.00000  0 \
 8.87730e+08    0.448269v


# recipe to do a field selfcal for holography
./applysolutions ChAvg_Beam0.ms cal_0.bin
./addfield.py ChAvg_Beam0.ms
./applysolutions ChAvg_Beam15.ms cal_15.bin
./addfield.py ChAvg_Beam15.ms
./calibrate -datacolumn CORRECTED_DATA -minuv 50.0 -startscan 0 -endscan 0 ChAvg_Beam0.ms holo_000.bin
./plotbp.py holo_000.bin holo_000.png


