#include <vector>
#include <algorithm>
#include <cpp11.hpp>

[[cpp11::register]]
double nth_element(std::vector<double> copy, int n) {
    auto it = copy.begin() + n;
    std::nth_element(copy.begin(), it, copy.end());
    return *it;
}
