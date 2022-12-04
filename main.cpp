#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;

float h, x0,y0,b;
ifstream input;
FILE *output;

float dydx(float x, float y){
    return (3*y-x*x)/x;
}
float u(float x){
    return 3*x*x*x + x*x;
}

void method_rk2(float h, float x, float y, float b) {
    float delta = 0;
    output = fopen("method_rk2_output","w");
    while ( x <= b + h  ) {
        fprintf(output, "%.4f | %.4f | %.4f\n", x, y, delta);
        x = x + h;
        y = y + h * dydx(x + h / 2, y + (h / 2) * dydx(x, y));
        delta = abs(u(x) - y);
    }
    fclose(output);
}


void method_rk3(float h, float x, float y, float b) {
    float delta = 0;
    float k1, k2, k3;
    output = fopen("method_rk3_output","w");
    while ( x <= b + h  ) {
        fprintf(output, "%.4f | %.4f | %.4f\n", x, y, delta);
        x = x + h;
        k1 = h * dydx(x,y);
        k2 = h * dydx(x + h/2, y + k1/2);
        k3 = h * dydx(x + h, y - k1 + 2*k2);
        y = y + (k1 + 4*k2 + k3)/6;
        delta = abs(u(x) - y);
    }
    fclose(output);
}


void method_rk4(float h, float x, float y, float b) {
    float delta = 0;
    float k1, k2, k3, k4;
    output = fopen("method_rk4_output","w");
    while ( x <= b + h  ) {
        fprintf(output, "%.4f | %.4f | %.4f\n", x, y, delta);
        x = x + h;
        k1 = h * dydx(x,y);
        k2 = h * dydx(x + h/2, y + k1/2);
        k3 = h * dydx(x + h/2, y + k2/2);
        k4 = h * dydx(x + h, y + k3);
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
        delta = abs(u(x) - y);
    }
    fclose(output);
}

int main() {

    input.open("input");
    input >> h >> x0 >> y0 >> b;
    input.close();

    method_rk2(h, x0, y0, b);
    method_rk3(h, x0, y0, b);
    method_rk4(h, x0, y0, b);


    return 0;
}
