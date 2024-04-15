![](img/mandelbrot.png)

## Железо

Процессор: Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz   1.50 GHz

Компилятор: g++ 13.2.0

ОС: Windows 10 22H2

Этот процессор поддерживает AVX-512 инструкции

## Работа

Я написал фрактал Мандельброта на С, используя библиотеку TXLib(http://ded32.net.ru/) для графики.
Для ускорения программы я написал тот же код, используя массивы, что позволяет процессору на О3 векторизовать програму.
Тем не менее, если векторизовать программу самостоятельно, при помощи интринсиков, результат будет значительно лучше.

## Управление

W A S D - Вверх вправо вниз влево

E Q     - Приблизить отдалить

1       - Наивная реализация

2       - Реализация на массивах

3       - AVX-512

# Результаты

Я провел тест, в котором каждая реализация проходит 100 циклов обработки без отрисовывания на О3.

* $\frac{\text{Naive}}{\text{Arrays}} = 2.704$
* $\frac{\text{Naive}}{\text{AVX}} = 6.209$
* $\frac{\text{Arrays}}{\text{AVX}} = 2.296$

Ускорение в 6.2 раза, впечатляет!

## Build
```
make
```
```
.\a.exe
```
