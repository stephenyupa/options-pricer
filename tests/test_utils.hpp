#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace testutils {

inline int& failure_count() {
    static int count = 0;
    return count;
}

inline void check(bool condition, const std::string& description) {
    if (!condition) {
        std::cerr << "  FAIL: " << description << "\n";
        ++failure_count();
    }
}

inline void check_near(double actual, double expected, double tol,
                        const std::string& description) {
    const bool ok = std::fabs(actual - expected) <= tol;
    if (!ok) {
        std::cerr << "  FAIL: " << description << " (actual=" << actual
                   << ", expected=" << expected << ", tol=" << tol << ")\n";
        ++failure_count();
    }
}

inline int report(const std::string& suite_name) {
    if (failure_count() == 0) {
        std::cout << suite_name << ": all checks passed\n";
        return EXIT_SUCCESS;
    }
    std::cerr << suite_name << ": " << failure_count() << " check(s) failed\n";
    return EXIT_FAILURE;
}

} // namespace testutils
