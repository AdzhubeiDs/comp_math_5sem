#include <iostream>
#include <array>
#include <fstream>
#include <cmath>
#include <iomanip>


const int N = 3;

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
private:
    std::array<xType, N> xPoints;
    std::array<yType, N> yValues;
    std::array<yType, N> dividedDifferences;

    void computeDividedDifferences() {
        for (unsigned int i = 0; i < N; ++i) {
            dividedDifferences[i] = yValues[i];
        }
        for (unsigned int j = 1; j < N; ++j) {
            for (unsigned int i = N - 1; i >= j; --i) {
                dividedDifferences[i] = (dividedDifferences[i] - dividedDifferences[i - 1]) / (xPoints[i] - xPoints[i - j]);
            }
        }
    }

public:
    NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values) noexcept
        : xPoints(points), yValues(values) {
        computeDividedDifferences();
    }

    yType interpolate(const xType& x) const noexcept {
        yType result = dividedDifferences[0];
        yType term = 1;

        for (unsigned int i = 1; i < N; ++i) {
            term *= (x - xPoints[i - 1]);
            result += dividedDifferences[i] * term;
        }

        return result;
    }
};

double f(double x) {
    return exp(x);
}


void evaluateInterpolation( const std::array<double, 2>& interval, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    double a = interval[0];
    double b = interval[1];
    std::array<double, N> xPoints;
    std::array<double, N> yValues;

    for (int i = 0; i < N; ++i) {
        xPoints[i] = a + (b - a) * i / (N - 1);
        yValues[i] = f(xPoints[i]);
    }

    NewtonInterpolator<double, double, N> interpolator(xPoints, yValues);

    double maxError = 0.0;

    for (int i = 0; i < 1000; ++i) {
        double x = a + (b - a) * i / 999;
        double error = std::abs(interpolator.interpolate(x) - f(x));

        file << error << ' ' << x << "\n";
    }
    file.close();
}

int main() {
    const std::array<std::array<double, 2>, 6> intervals = {{{0, 2}, {0, 1}, {0, 0.5}, {0, 0.25}, {0, 0.125}, {0, 0.0625}}};


    for (const auto& interval : intervals) {
        evaluateInterpolation( interval, "data_" + std::to_string(N) + ".txt");
    }


    return 0;
}
