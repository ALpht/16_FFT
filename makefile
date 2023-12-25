fft : fft.o
	g++ -o fft fft.o
fft.o : fft.cpp
	g++ -c fft.cpp
clean:
	rm *.o fft
