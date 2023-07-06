
g++ -c predicter.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c mspredicter.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c calibrator.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c calibrationmethod.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c polynomialfitter.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c model/model.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c nlplfitter.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include
g++ -c progressbar.cpp -Wall -std=c++11 -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include

g++ calibrate.cpp calibrationmethod.o predicter.o progressbar.o model.o nlplfitter.o polynomialfitter.o  mspredicter.o calibrator.o -Wall -std=c++11 -DNDEBUG -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include -L/usr/lib -lgsl -lgslcblas -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -lboost_thread -lboost_system -o calibrate -Wl,-rpath,/usr/lib

g++ addmodel.cpp calibrationmethod.o predicter.o progressbar.o model.o nlplfitter.o polynomialfitter.o  mspredicter.o calibrator.o -Wall -std=c++11 -DNDEBUG -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include -L/usr/lib -lgsl -lgslcblas -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -lboost_thread -lboost_system -o addmodel -Wl,-rpath,/usr/lib

g++ applysolutions.cpp -Wall -std=c++11 -DNDEBUG -msse2 -O3 -march=native -I/usr/include/casacore -I/usr/include -L/usr/lib -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -o applysolutions -Wl,-rpath,/usr/lib

cp calibrate /usr/local/bin
cp applysolutions /usr/local/bin
cp addmodel /usr/local/bin
