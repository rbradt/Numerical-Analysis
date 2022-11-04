#include "EulerMethod.hpp"


// One Dimensional Ligand Receptor Initial Value Problem
class LR1D_IVP : public IVP<Matrix<double, 2, 1>> {
private:
    double ah1, bh2, v, h1, h2;             // parameters
    coord<Matrix<double, 2, 1>>* init;      // initial condition
    std::vector<double> u_n, v_n, t;        // solution
    Euler_Method<Matrix<double, 2, 1>>* s;  // initial value problem solver

    // f(x) = x^h1/(a^h1 + x^h1)
    double f(double x) {
        double xh1 = pow(x, h1);
        return xh1 / (ah1 + xh1);
    }

    // g(x) = b^h2/(b^h2 + x^h2)
    double g(double x) { return bh2 / (bh2 + pow(x, h2)); }

public:
    LR1D_IVP(double a, double b, double v, double h1, double h2) : ah1{ pow(a, h1) }, bh2{ pow(b, h2) }, v{ v }, h1{ h1 }, h2{ h2 }, init{ nullptr }, s{ nullptr } {}

    // set the various parameters of the initial value problem
    void set_parameters(double a, double b, double _v, double _h1, double _h2) {
        ah1 = pow(a, h1);
        bh2 = pow(b, h2);
        v = _v;
        h1 = _h1;
        h2 = _h2;
    }

    // set the initial condition
    void set_ic(Matrix<double, 2, 1> x0, double t0) {
        if(init != nullptr)
            delete init;
        init = new coord<Matrix<double, 2, 1>>(x0, t0);
        
        u_n.clear();
        v_n.clear();
        t.clear();
        parse(init);
    }

    // set the interval and increment
    void set_interval(double a, double b, double increment) {
        if (s == nullptr)
            s = new Euler_Method<Matrix<double, 2, 1>>(this, a, b, increment);
        else {
            s->domain(a, b);
            s->increment(increment);
        }
    }

    // retrieve the initial condition
    coord<Matrix<double, 2, 1>>& get_ic() {
        return *init;
    }

    // parse a single solution point
    void parse(coord<Matrix<double, 2, 1>>* data) {
        u_n.push_back(data->y(0, 0));
        v_n.push_back(data->y(1, 0));
        t.push_back(data->x);
    }

    // compute and display the full solution
    void display_solution() {
        s->compute();
        plt::named_plot("u_n(t)", t, u_n);
        plt::named_plot("v_n(t)", t, v_n);
        plt::xlabel("Time");
        plt::legend();
        plt::grid(1);
        plt::show();
    }

    // initial value problem
    // x' = [u_n', v_n'] = [f(v_n) - u_n, v(g(u_n) - v_n)]
    Matrix<double, 2, 1> operator() (Matrix<double, 2, 1> x, double t) {
        x(0, 0) = f(x(1,0)) - x(0, 0);
        x(1, 0) = v * (g(x(0,0)) - x(1, 0));
        return x;
    }

    ~LR1D_IVP() {
        delete init;
    }
};

// g++ EulerMethod/main.cpp -I C:/Python310/include -I C:/python310/lib/site-packages/numpy/core/include -L C:\Python310\libs -lpython310
// g++ EulerMethod/main.cpp -I <drive>:/<path>/Python310/include -I <drive>:/<path>/python310/lib/site-packages/numpy/core/include -L <drive>:/<path>\Python310\libs -lpython310

int main() {
    LR1D_IVP* f = new LR1D_IVP(0.5,0.5,1,2,2);
    Matrix<double, 2, 1> y0(0,0);
    f->set_ic(y0, 0);
    f->set_interval(0, 100, 0.01);
    f->display_solution();

    delete f;
}