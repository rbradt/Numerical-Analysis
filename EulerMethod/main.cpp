// EulerMethod.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "Eigen/Core"

using Eigen::Matrix;
using Eigen::Vector;

template <typename Y>
class Euler_Method {
private:
    struct coord {
        Y y;
        double x;
        coord(Y _y, double _x) : y{ _y }, x{ _x } {};
        coord(const coord* other) : y{ other->y }, x{ other->x } {};
    };
    Y (*f)(Y, double);
    coord init;
    double dx, a, b;
    std::vector<coord*> pts;

    bool inDomain(double x) {
        return (dx > 0) ? x + dx < b + 0.0000001 : x + dx > a - 0.0000001;
    }

    void clear() {
        for (coord* v : pts)
            delete v;
        pts.clear();
    }

public:
    Euler_Method(Y(*function)(Y, double), double a, double b, double increment) : f{ function }, init{new coord(1,1)}, a{ a }, b{ b }, dx{ increment } {}

    void function(Y(*function)(Y, double)) {

    }

    void initial_condition(Y y_at_X0, double x0) {

    }

    void domain(double lower_bound, double upper_bound) {

    }

    void increment(double increment) {

    }

    void euler_method(Y y0, double x0) {
        if (a <= x0 && x0 <= b) {
            if (!pts.empty())
                clear();
            pts.push_back(new coord(y0, x0));
            while(inDomain(pts.back()->x)) {
                coord* c = new coord(pts.back());
                c->y = dx * f(c->y, c->x) + c->y;
                c->x += dx;
                pts.push_back(c);
            }
        }
    }

    void print() const {
        for (coord* v : pts) {
            std::cout << "(" << v->x << ", " << v->y << ")" << std::endl;
        }
    }

    ~Euler_Method() {
        clear();
    }
};

double func(double y, double x) {
    return 2*y*y + 2;
}

Matrix<double,2,1> func2(Matrix<double, 2, 1> y, double x) {
    y(1, 1) += y(1, 1);
    return y;
}

int main()
{
    Euler_Method<double> a(func, -0.1, 0.1, -0.001);
    a.euler_method(1, 0);
    a.print();
}