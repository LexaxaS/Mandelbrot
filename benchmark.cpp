#include "C:\Programming\TX\TXlib.h"
#include <stdlib.h>
#include "mandelbrot.hpp"

static size_t _FuncBenchmark(DrawFunc_t drawFunc, scr_t screen, RGBQUAD* colors);

error_t Benchmark(scr_t screen, RGBQUAD* colors)
    {
    size_t  NaiveTime = _FuncBenchmark(DrawFuncs[0], screen, colors);
    size_t ArraysTime = _FuncBenchmark(DrawFuncs[1], screen, colors);
    size_t AVX512Time = _FuncBenchmark(DrawFuncs[2], screen, colors);

    printf(" Naive / AVX512 = %.4lg\n", (double) NaiveTime / (double) AVX512Time);
    printf("Arrays / AVX512 = %.4lg\n", (double) ArraysTime / (double) AVX512Time);
    printf(" Naive / Arrays = %.4lg\n", (double) NaiveTime / (double) ArraysTime);
    printf("in %d runs\n", BenchmarkRuns);

    return NO_ERROR;
    }

static size_t _FuncBenchmark(DrawFunc_t drawFunc, scr_t screen, RGBQUAD* colors)
    {
    size_t start = __rdtsc();

    for (size_t i = 0; i < BenchmarkRuns; i++)
        drawFunc(screen, colors, 0, 0, 1.f);
    
    size_t end   = __rdtsc();

    return end - start;
    }
