all: ./main1D.cpp ./main2D.cpp ./accuracy.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp
	g++ -o main1D main1D.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include -std=c++20 -llapacke -llapack;\
	g++ -o main2D main2D.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include -std=c++20 -llapacke -llapack;\
	g++ -o accuracy accuracy.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include -std=c++20 -llapacke -llapack;

.PHONY: run
run: 
	./main1D;\
	./main2D;\
	./accuracy

