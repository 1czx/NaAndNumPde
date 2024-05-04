all: ./main1.cpp ./main2.cpp  ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp
	g++ -o test1 main1.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include -Ofast;
	g++ -o test2 main2.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include -Ofast;
	g++ -o test3 main3.cpp ./json_src/lib_json/json_reader.cpp ./json_src/lib_json/json_value.cpp ./json_src/lib_json/json_writer.cpp -I include -Ofast;

.PHONY: run
run: ./test1 ./test2 ./test3
	./test1 & ./test2 & ./test3;\

