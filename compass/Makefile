CPP_FLAGS=$(shell apfelxx-config --cppflags) $(shell lhapdf-config -cxx) $(shell apfel-config --cppflags)
LD_FLAGS=$(shell apfelxx-config --ldflags) $(shell lhapdf-config -ldflags) $(shell apfel-config --ldflags) -lm -lyaml-cpp


all:
	g++ -o test main.cpp $(CPP_FLAGS) $(LD_FLAGS) -O3

clean:
	rm -rf test
