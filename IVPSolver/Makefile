all:
	make -C src;\

.PHONY: run
run:
	make;\
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
	rm -f ./bin/*;\
	rm -f ./doc/report.pdf;\

.PHONY: story
story:
	make story -C doc;\