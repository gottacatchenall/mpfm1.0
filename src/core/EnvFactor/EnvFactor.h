
#ifndef ENV_FACTOR_H
#define ENV_FACTOR_H

#include "include.h"

class EnvFactor{
    private:
        int id;
        double** map;
    public:
        EnvFactor(int size, double h_val);
        double get_env_factor_value(double x, double y);
        double** generate_fractal(int size, double H_VAL);
        double f4(double delta, double a, double b, double c, double d);
        double f3(double delta, double a, double b, double c);
};

#endif
