all: 
	g++ -mavx512f -O3  .\main.cpp .\mandelbrot.cpp .\benchmark.cpp

