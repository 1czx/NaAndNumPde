all:
	make -C src;\
	mkdir bin;\
	mv ./src/main1D ./bin;\
	mv ./src/main2D ./bin;\
	mv ./src/accuracy ./bin;\

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
