all: ./main.cpp ./include/function.h ./include/splines.h ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp
	g++ -o test main.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include;

.PHONY: run
run:
	./test;\

.PHONY: plot
plot: plot.py
	python3 plot.py;