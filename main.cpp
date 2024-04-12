#include "C:\Programming\TX\TXlib.h"
#include <stdlib.h>
#include "mandelbrot.hpp"

int main()
    {
    scr_t screen = CreateTxWindow();
    RGBQUAD* colors = MakeColors();

    // Benchmark(screen, colors);
    MandelbrotLoop(screen, colors);
    
    txEnd();
    return 0;
    }