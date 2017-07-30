CSC ?= csc

all: fftw.so

fftw.so: fftw.scm
	$(CSC) -s -j fftw fftw.scm -lfftw3
