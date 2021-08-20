LDFLAGS=-lm -lyaml-cpp

all: reweight chi2

reweight:
	g++ -o reweight main.cpp -O3 $(LDFLAGS)

chi2:
	g++ -o chi2 analysis.cpp -O3 $(LDFLAGS)

clean:
	rm -f testff_p/*
	rm -f testff_m/*
	rm -f reweight
	rm -f chi2
