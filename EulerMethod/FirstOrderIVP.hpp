#ifndef __FIRSTORDERIVP_HPP__
#define __FIRSTORDERIVP_HPP__

template <typename Y>
struct coord {
    Y y;
    double x;
    coord(Y _y, double _x) : y{ _y }, x{ _x } {};
    coord(const coord* other) : y{ other->y }, x{ other->x } {};
    coord(const coord& other) : y{ other.y }, x{ other.x } {};
};

template<typename Y>
class IVP {
public:
	virtual Y operator() (Y, double) = 0; // Function of the first order initial value problem
    virtual void set_ic(Y, double) = 0; // Set the initial condition
    virtual coord<Y>& get_ic() = 0; // Retrieve the initial condition
    virtual void parse(coord<Y>*) = 0;  // Parse a single data point
    virtual void display_solution() = 0; // Display solution to the initial value problem
};

#endif // __FIRSTORDERIVP_HPP__
