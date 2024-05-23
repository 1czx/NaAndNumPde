SRC = $(wildcard ./*.cpp)
JSONSRC = $(wildcard ./json_src/lib_json/*.cpp)

all:
	g++ -o ../bin/test $(SRC) $(JSONSRC) -I include -std=c++20 -llapacke -llapack

.PHONY: run
run:
	./test;\

.PHONY: plot
plot: plot.py
	python3 plot.py;
