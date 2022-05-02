#ifndef _COMPLEX_H_
#define _COMPLEX_H_

typedef struct 
{
    double re;
    double im;
} complex_double;

complex_double mul_by_factor(complex_double x, double a);
complex_double mul(complex_double x1, complex_double x2);
complex_double sub(complex_double x1, complex_double x2);
complex_double add(complex_double x1, complex_double x2);

#endif
