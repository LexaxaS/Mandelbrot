#include "C:\Programming\TX\TXlib.h"
#include <emmintrin.h>
#include <stdlib.h>
#include "mandelbrot.hpp"

static const __m512d _01234567 = _mm512_set_pd (7.f, 6.f, 5.f, 4.f, 3.f, 2.f, 1.f, 0.f);
static const __m512d R2MAX_512 = _mm512_set1_pd (r2Max); 

const double dx = 1 / 800.f;
const double dy = 1 / 800.f;

scr_t CreateTxWindow()
    {
    txCreateWindow (scr_width, scr_height);
    Win32::_fpreset();
    txBegin();

    scr_t scr = (scr_t) *txVideoMemory();
    for (int iy = 0; iy < 800; iy++)
        for (int ix = 0; ix < 800; ix++)
            scr[iy][ix] = {0, 0, 0};

    return scr;
    }

RGBQUAD* MakeColors()
    {
    RGBQUAD* colors = (RGBQUAD*) calloc (256, sizeof(*colors));
    for (int N = 0; N < 256; N++) 
        {
        double I = sqrtf( sqrtf( (double) N / (double) nMax )) * 255.f;

        BYTE b_color = (BYTE) I;
        RGBQUAD color = {(BYTE) (255 - b_color), (BYTE) (b_color), (BYTE) (b_color)};
        colors[N] = color;
        }
    return colors;
    }

error_t MandelbrotLoop(scr_t screen, RGBQUAD* colors)
    {
    double xShift = 0;
    double yShift = 0;
    double scale  = 1.f;

    int DrawFuncNum = 0;
    while (true)
        {
        if (GetAsyncKeyState (VK_ESCAPE)) break;

        if (txGetAsyncKeyState ('W'))  yShift += dy * 50 / scale;
        if (txGetAsyncKeyState ('S'))  yShift -= dy * 50 / scale;
        if (txGetAsyncKeyState ('D'))  xShift += dx * 50 / scale;
        if (txGetAsyncKeyState ('A'))  xShift -= dx * 50 / scale;

        if (txGetAsyncKeyState ('E'))  scale  *=  1.2;
        if (txGetAsyncKeyState ('Q'))  scale  /=  1.2;

        if (txGetAsyncKeyState ('1'))  DrawFuncNum = 0;
        if (txGetAsyncKeyState ('2'))  DrawFuncNum = 1;
        if (txGetAsyncKeyState ('3'))  DrawFuncNum = 2;

        DrawFuncs[DrawFuncNum](screen, colors, xShift, yShift, scale);

        printf("\t\rfps=[%.0lf]          ", txGetFPS());
        txSleep();
        }

    return NO_ERROR;
    }

error_t MandelbrotNaive(scr_t screen, RGBQUAD* colors, double xShift, double yShift, double scale)
    {
    for (size_t iy = 0; iy < scr_height; iy++)
        {

        double y0 = (((double) iy - 400.f) * dy) / scale + yShift;

        for (size_t ix = 0; ix < scr_width; ix++)
            {
            double x0 = - 1.325f +  ( ((double) ix - 400.f) * dx) / scale + xShift;

            double x = x0;
            double y = y0;

            int N = 0;
            for (size_t n = 0; n < nMax; n++)
                {
                double x2 = x * x;
                double y2 = y * y;
                double xy = x * y;
                
                double r2 = x2 + y2;
                
                int cmp = 0; if (r2 < r2Max) cmp = 1; else break;

                N += cmp;

                x = x2 - y2 + x0;
                y = xy + xy + y0;

                }
            screen[iy][ix] = colors[N];
            }
        }

    return NO_ERROR;
    }


error_t MandelbrotArrays(scr_t screen, RGBQUAD* colors, double xShift, double yShift, double scale)
    {
    for (size_t iy = 0; iy < scr_height; iy++)
        {
        double y0 = (((double) iy - 400.f) * dy) / scale + yShift;
        double YO[PAO] = {}; for (int i = 0; i < PAO; i++) YO[i] = y0;

        for (size_t ix = 0; ix < scr_width; ix += PAO)
            {
            double x0 = - 1.325f +  ( ((double) ix - 400.f) * dx) / scale + xShift;
            double XO[PAO] = {}; for (int i = 0; i < PAO; i++) XO[i] = x0 + i * dx / scale;

            double  X[PAO] = {}; for (int i = 0; i < PAO; i++)  X[i] = XO[i];
            double  Y[PAO] = {}; for (int i = 0; i < PAO; i++)  Y[i] = YO[i];

            int    N[PAO] = {}; for (int i = 0; i < PAO; i++)  N[i] = 0;
            for (size_t n = 0; n < nMax; n++)
                {
                double x2[PAO] = {}; for (int i = 0; i < PAO; i++) x2[i] =  X[i] * X[i];
                double y2[PAO] = {}; for (int i = 0; i < PAO; i++) y2[i] =  Y[i] * Y[i];
                double xy[PAO] = {}; for (int i = 0; i < PAO; i++) xy[i] =  X[i] * Y[i];
                
                double r2[PAO] = {}; for (int i = 0; i < PAO; i++) r2[i] = x2[i] + y2[i];
                
                int flag = 0;

                int cmp[PAO]  = {}; for (int i = 0; i < PAO; i++) if (r2[i] <= r2Max) {cmp[i] = 1; flag = 1;}

                if (!flag) break;

                for (int i = 0; i < PAO; i++) N[i] =  N[i] + cmp[i];

                for (int i = 0; i < PAO; i++) X[i] = x2[i] - y2[i] + XO[i];
                for (int i = 0; i < PAO; i++) Y[i] = xy[i] + xy[i] + YO[i];
                }

            for (int i = 0; i < PAO; i++) 
                {
                screen[iy][ix + i] = colors[N[i]];
                }
            }
        }
    return NO_ERROR;
    }


error_t MandelbrotAVX512(scr_t screen, RGBQUAD* colors, double xShift, double yShift, double scale)
    {
    const __m512d DX = _mm512_mul_pd (_mm512_set1_pd (dx / scale), _01234567);
    
    for (size_t iy = 0; iy < scr_height; iy++)
        {
        double y0 = (((double) iy - 400.f) * dy) / scale + yShift;
        const __m512d Y0 = _mm512_set1_pd(y0);

        for (size_t ix = 0; ix < scr_width; ix += 8)
            {
            double x0 = - 1.325f +  ( ((double) ix - 400.f) * dx) / scale + xShift;
            const __m512d X0 = _mm512_add_pd (_mm512_set1_pd(x0), DX);

            __m512d X = X0;
            __m512d Y = Y0;

            __m512i N = _mm512_set1_epi32(0);

            for (size_t n = 0; n < nMax; n++)
                {
                __m512d X2 = _mm512_mul_pd(X, X);
                __m512d Y2 = _mm512_mul_pd(Y, Y);
                __m512d XY = _mm512_mul_pd(X, Y);

                __mmask8 cmp = _mm512_cmplt_pd_mask (_mm512_add_pd(X2, Y2), R2MAX_512);
                N = _mm512_add_epi32 (N, _mm512_maskz_expand_epi32 (cmp, _mm512_set1_epi32(1))); 

                if (cmp == 0) break;

                X = _mm512_add_pd(_mm512_sub_pd(X2, Y2), X0);
                Y = _mm512_add_pd(_mm512_add_pd(XY, XY), Y0);
                }

            int row[8] = {};
            _mm512_storeu_epi64(row, N);

            for (int i = 0; i < 8; i++)
                screen[iy][ix + i] = *(colors + row[i]);
            }
        }

    return NO_ERROR;
    }

