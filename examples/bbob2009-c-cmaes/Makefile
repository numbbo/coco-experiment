CFLAGS ?= -std=c99 -D_XOPEN_SOURCE=600 -D_POSIX_C_SOURCE=200809L -O0 -g -march=native -ffp-contract=off
CXXFLAGS ?= -std=c++17 -O0 -g -march=native -ffp-contract=off -fno-exceptions
OBJS := run_bbob.o coco.o cmaes.o
LIBS := -lm

.PHONEY: all run clean

all: run_bbob

run: run_bbob
	./run_bbob

clean:
	rm -fR run_bbob run_bbob.o coco.o cmaes.o c-cmaes* actparcmaes.par 

run_bbob: ${OBJS}
	${CXX} ${LDFLAGS} ${OBJS} ${LIBS} -o run_bbob

run_bbob.o: run_bbob.c coco.h cmaes_interface.h 

coco.o: coco.cpp coco.h
	${CXX} -c ${CXXFLAGS} -o coco.o coco.cpp

cmaes.o: cmaes.c cmaes_interface.h cmaes.h
