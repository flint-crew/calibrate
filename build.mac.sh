export CONDA_PATH=$HOME/miniforge3/envs/casacore
set -e
g++ -c predicter.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c mspredicter.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c calibrator.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c calibrationmethod.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c polynomialfitter.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c model/model.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c nlplfitter.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include
g++ -c progressbar.cpp -Wall -std=c++11 -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include

g++ calibrate.cpp calibrationmethod.o predicter.o progressbar.o model.o nlplfitter.o polynomialfitter.o  mspredicter.o calibrator.o -Wall -std=c++11 -DNDEBUG -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include -L$CONDA_PATH/lib -lgsl -lgslcblas -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -lboost_thread -lboost_system -o calibrate -Wl,-rpath,$CONDA_PATH/lib

g++ addmodel.cpp calibrationmethod.o predicter.o progressbar.o model.o nlplfitter.o polynomialfitter.o  mspredicter.o calibrator.o -Wall -std=c++11 -DNDEBUG -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include -L$CONDA_PATH/lib -lgsl -lgslcblas -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -lboost_thread -lboost_system -o addmodel -Wl,-rpath,$CONDA_PATH/lib

g++ applysolutions.cpp -Wall -std=c++11 -DNDEBUG -O3 -march=native -I$CONDA_PATH/include/casacore -I$CONDA_PATH/include -L$CONDA_PATH/lib -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -o applysolutions -Wl,-rpath,$CONDA_PATH/lib

cp calibrate ~/bin
cp applysolutions ~/bin
cp addmodel ~/bin
