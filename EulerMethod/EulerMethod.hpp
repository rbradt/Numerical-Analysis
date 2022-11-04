#ifndef __EULERMETHOD_HPP__
#define __EULERMETHOD_HPP__

#include <iostream>
#include <vector>
#include "../Resources/Eigen/Core"
#include "../Resources/MathPlotLib/matplotlibcpp.h"
#include "FirstOrderIVP.hpp"

namespace plt = matplotlibcpp;

using Eigen::Matrix;
using Eigen::Vector;

template <typename Y>
class Euler_Method {
private:
    IVP<Y>* f;
    double dx, a, b;
    std::vector<coord<Y>*> pts;

    bool inDomain(double x) {
        return (dx > 0) ? x + dx < b + 0.0000001 : x + dx > a - 0.0000001;
    }

    void clear() {
        for (coord<Y>* v : pts)
            delete v;
        pts.clear();
    }

public:
    Euler_Method(IVP<Y>* f, double a, double b, double increment) : f{ f }, a{ a }, b{ b }, dx{ increment } {}

    void function(IVP<Y>* f) {
        f = f;
    }

    void domain(double lower_bound, double upper_bound) {
        a = lower_bound;
        b = upper_bound;
    }

    void increment(double increment) {
        dx = increment;
    }

    void compute() {
        if (a <= f->get_ic().x && f->get_ic().x <= b) {
            if (!pts.empty())
                clear();
            pts.push_back(new coord(f->get_ic()));
            while (inDomain(pts.back()->x)) {
                coord<Y>* c = new coord(pts.back());
                c->y = dx * (*f)(c->y, c->x) + c->y;
                c->x += dx;
                pts.push_back(c);
                f->parse(c);
            }
        }
    }

    ~Euler_Method() {
        clear();
    }
};

#endif