
main: *.cpp
	g++ -fopenmp -O2 -I /home/dirk/adv_fem/eigen/eigen-3.4.0 -I /home/dirk/adv_fem/eigen/json main.cpp read_json.cpp wrte_lvtk.cpp dmat.cpp smat.cpp scg.cpp  -o main
