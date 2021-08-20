CPP_FLAGS=$(shell apfel-config --cppflags) $(shell lhapdf-config -cxx)
LD_FLAGS=$(shell apfel-config --ldflags) $(shell lhapdf-config -ldflags)

all:
	g++ -o test main.cpp -lyaml-cpp -lm $(CPP_FLAGS) $(LD_FLAGS) -O3

benchmark:
	g++ -o benchmark benchmark.cpp -lyaml-cpp -lm $(CPP_FLAGS) $(LD_FLAGS)	
clean:
	rm -rf test benchmark
