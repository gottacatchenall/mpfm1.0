#include "EnvFactor.h"

EnvFactor::EnvFactor(int size, double h_val){
    this->map = generate_fractal(size, h_val);
}

double EnvFactor::get_env_factor_value(double x, double y){
    int res = params["ENV_FACTOR_RESOLUTION"];
    double side_len = params["SIDE_LENGTH"];
    int x_cell = int(double(res)*(x/side_len));
    int y_cell = int(double(res)*(y/side_len));
    return this->map[x_cell][y_cell];
}

double EnvFactor::get_cell_value(int x, int y){
    return this->map[x][y];
}

double** EnvFactor::generate_fractal(int size, double H_VAL){
    int n = size+1;
    double** tmp = new double*[n];
    for (int i = 0; i < n; i++){
        tmp[i] = new double[n];
    }

    double delta = 1;
    int N = size;

    int maxlevel = log2(size);

    double H = H_VAL;

    int D = N;
    int d = N/2;

    for (int i = 0; i <= N; i++){
        for (int j = 0; j <= N; j++){
            tmp[i][j] = 0.0;
        }
    }

    // Init corner values
    tmp[0][0] = delta*std_normal(ef_generator);
    tmp[0][N] = delta*std_normal(ef_generator);
    tmp[N][0] = delta*std_normal(ef_generator);
    tmp[N][N] = delta*std_normal(ef_generator);

    for (int stage = 0; stage < maxlevel; stage++){
        delta = delta * double(pow(double(0.5), double(0.5*H)));

        // interpolate and offset points
        for (int x = d; x <= N-d; x+=D){
            for (int y = d; y <= N-d; y+=D){
                tmp[x][y] = this->f4(delta, tmp[x+d][y+d], tmp[x+d][y-d], tmp[x-d][y+d],tmp[x-d][y-d]);
            }
        }

        delta = delta * double(pow(double(0.5), double(0.5*H)));

        // boundary grid points
        for (int x = d; x <= N-d; x+=D){
            tmp[x][0] = this->f3(delta, tmp[x+d][0], tmp[x-d][0], tmp[x][d]);
            tmp[x][N] = this->f3(delta, tmp[x+d][N], tmp[x-d][N], tmp[x][N-d]);
            tmp[0][x] = this->f3(delta, tmp[0][x+d], tmp[0][x-d], tmp[d][x]);
            tmp[N][x] = this->f3(delta, tmp[N][x+d], tmp[N][x-d], tmp[N-d][x]);
        }

        // interpolate and offset interior grid points
        for (int x = d; x <= N-d; x+=D){
            for (int y = D; y <= N-d; y+=D){
                tmp[x][y] = this->f4(delta, tmp[x][y+d], tmp[x][y-d], tmp[x+d][y], tmp[x-d][y]);
            }
        }

        for (int x = D; x <= N-d; x+=D){
            for (int y = d; y <= N-d; y+=D){
                tmp[x][y] = this->f4(delta, tmp[x][y+d], tmp[x][y-d], tmp[x+d][y], tmp[x-d][y]);
            }
        }

        D = D/2;
        d = d/2;
    }

    double** grid = new double*[n];
    for (int i = 0; i < n; i++){
        grid[i] = new double[n];
    }


    double min = tmp[0][0];
    double max = tmp[0][0];


    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (tmp[i][j] > max){
                max = tmp[i][j];
            }
            if (tmp[i][j] < min){
                min = tmp[i][j];
            }
        }
    }

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            grid[i][j] = double(max - tmp[i][j])/(max-min);
        }
    }


    for (int i = 0; i < n; i++){
        delete tmp[i];
    }
    delete[] tmp;

    return grid;
}

double EnvFactor::f4(double delta, double a, double b, double c, double d){
    double sum = a+b+c+d;
    double avg = sum/4.0;
    double val =  avg + double(delta) * double(std_normal(ef_generator));
    return val;
}

double EnvFactor::f3(double delta, double a, double b, double c){
    double sum = a+b+c;
    double avg = sum/3.0;
    double val =  avg + double(delta) * double(std_normal(ef_generator));
    return val;
}
