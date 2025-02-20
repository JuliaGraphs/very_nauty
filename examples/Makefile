OPT=-g
OPT=-O3
INC=-I..
LIB=-L.. -lvn_graph

all: example_01 example_02 example_03

example_02: example_02.c
	gcc -Wall ${OPT} ${INC} -fomit-frame-pointer example_02.c -o example_02 ${LIB} -lm
	@echo Run example_02 like this "./example_02 [n [r [ngraphs]]]"
	
example_01: example_01.c
	gcc -Wall ${OPT} ${INC} -fomit-frame-pointer example_01.c -o example_01 ${LIB} -lm
	
%: %.c
	gcc -Wall ${OPT}  ${INC} -fomit-frame-pointer $< -o $@ ${LIB} -lm

clean:
	rm -f *o a.out *dat *dot core
	./execlean.sh
