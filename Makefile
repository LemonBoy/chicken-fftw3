CSC ?= csc

all: fftw.so

fftw.so: fftw.scm Makefile
	$(CSC) -k -C -g3 -d3 -s -J fftw.scm -lfftw3 -lfftw3f
