#ifndef MANDELBROT_HPP
#define MANDELBROT_HPP

#include "C:\Users\1\OneDrive\Рабочий стол\TX\TXlib.h"

const int scr_width = 800;
const int scr_height = 800;

const int PAO = 4;                  // Pixels at once
const size_t nMax = 256;
const double r2Max = 100.f;

typedef RGBQUAD (&scr_t) [scr_width][scr_height];
typedef size_t error_t;
typedef error_t (*DrawFunc_t)(scr_t, RGBQUAD*, double, double, double);

error_t MandelbrotLoop(scr_t screen, RGBQUAD* colors);

const size_t BenchmarkRuns = 500;
error_t Benchmark(scr_t screen, RGBQUAD* colors);

error_t MandelbrotNaive(scr_t screen, RGBQUAD* colors, double xShift, double yShift, double scale);
error_t MandelbrotArrays(scr_t screen, RGBQUAD* colors, double xShift, double yShift, double scale);
error_t MandelbrotAVX512(scr_t screen, RGBQUAD* colors, double xShift, double yShift, double scale);

const DrawFunc_t DrawFuncs[] = {MandelbrotNaive, MandelbrotArrays, MandelbrotAVX512};

RGBQUAD* MakeColors();
scr_t CreateTxWindow();


#endif