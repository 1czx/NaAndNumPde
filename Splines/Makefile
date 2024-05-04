all:
	make -C src;\
	mkdir bin;\
	mv ./src/test ./bin;\
	cp ./src/parameter.json ./bin/parameter.json;\

.PHONY: run
run:
	cp ./src/Makefile ./bin/Makefile;\
	make run -C bin;\
	rm ./bin/Makefile;\
	make plot;\

.PHONY:plot
plot: ./src/plot.py
	make plot -C src;\

.PHONY: clean
clean:
	rm -f ./result/*;\
	rm -f ./fig/*.jpg;\
	rm -rf ./bin/
