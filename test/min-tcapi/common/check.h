#pragma once
#include "tcapi/tcapi.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <limits>
#include <type_traits>
inline int checks = 0;
inline void check(bool condition, const char* message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
template<class Exception, class F> void throws(F&& f) {
    bool caught = false;
    try { f(); } catch (const Exception&) { caught = true; }
    check(caught, "expected exception was not thrown");
}
template<class E> E value(double r, double i = 0) {
    if constexpr (std::is_arithmetic_v<E>) return static_cast<E>(r);
    else return E(static_cast<typename E::value_type>(r), static_cast<typename E::value_type>(i));
}
#define RUN_FOUR_TYPES(test) \
int main() { try { \
    test<float>(); test<double>(); test<std::complex<float>>(); test<std::complex<double>>(); \
    std::cout << "PASS " << checks << " checks\n"; return 0; \
} catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; } }
