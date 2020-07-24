g++ -std=c++11 -c -O4 *.cpp 
g++ -std=c++11 -O4 Source.cpp com.o Individual.o Intervention.o linpack.o Mosquito.o Params.o Population.o randlib.o Simulation.o -o a.out
./a.out model_parameters.txt domestic_mosquitoes.txt occupational_mosquitoes.txt intervention_coverage.txt model_output.txt
