#include "complex_own.h"

complex_double mul_by_factor(complex_double x, double a)
{
    complex_double y = {
        .re = x.re * a,
        .im = x.im * a
    };
    return y;
}

complex_double mul(complex_double x1, complex_double x2)
{
    complex_double y = {
        .re = x1.re * x2.re - x1.im * x2.im,
        .im = x1.re * x2.im + x1.im * x2.re
    };
    return y;
}

complex_double sub(complex_double x1, complex_double x2)
{
    complex_double y = {
        .re = x1.re - x2.re,
        .im = x1.im - x2.im
    };
    return y;
}

complex_double add(complex_double x1, complex_double x2)
{
    complex_double y = {
        .re = x1.re + x2.re,
        .im = x1.im + x2.im
    };
    return y;
}
