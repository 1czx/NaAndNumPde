all:
	make -C src;\
	mkdir bin;\
	mv ./src/test1 ./bin;\
	mv ./src/test2 ./bin;\
	mv ./src/test3 ./bin;\

.PHONY: run
run: all
	cp ./src/Makefile ./bin/Makefile;\
	make run -C bin;\
	rm ./bin/Makefile;\

.PHONY: clean
clean:
	rm -f ./result/*;\
	rm -f ./fig/*.jpg;\
	rm -rf ./bin/
