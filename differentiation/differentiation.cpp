#include <iostream>
#include <vector>
#include <array>
#include <type_traits>
#include <typeinfo>
#include <cmath>
#include <fstream>
#include <limits>

typedef long double Float;
// typedef float Float;
typedef /*__int128*/ long long Integer;
typedef std::vector<std::vector<Float>> Matrix;
typedef std::vector<Float> Vector;


template<size_t N, size_t M>
class Task{
public:
    Task() : eps(std::numeric_limits<Float>::epsilon()) {
        for (size_t i = 0; i < N; ++i) {
            NArray[i] = i + 3;
        }

        hArray[0] = 1;
        for (size_t i = 1; i < M; ++i) {
            hArray[i] = hArray[i - 1] / static_cast<Float>(10);
        }

        precalcFactorials();
    }

    void Solve() {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                solveForN(NArray[i], hArray[j]);
            }
        }
    }

private:
    const Float eps;
    Integer x = 1;
    std::array<Float, M> hArray;
    std::array<Integer, N> NArray;

    void solveForN(Integer n, Float h) {
        Vector solution = calculateSLHE(n, h);

        Float sigma = std::exp(static_cast<Float>(1));
        Float err = 0;
        for (size_t i = 0; i < n; ++i) {
            // sigma -= solution[i] * std::exp(x + i * h);
            err += std::exp(i * h) * solution[i];
        }
        sigma -= std::exp(1) * err;
        sigma += (eps * std::exp(1)) / h;

        writeInFile(std::to_string(n), n, h, sigma);
    }




    Vector calculateSLHE(Integer n, Float h) {
        Matrix c(n, Vector(n));
        Vector b(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == 0) {
                    c[i][j] = 1;
                } else {
                    if (j != 0) {
                        Float t1 = binPow(j * h, i);
                        Float t2 = factorial(i);
                        c[i][j] = binPow(j * h, i) / factorial(i);
                    }
                }
            }
            b[i] = static_cast<Integer>(i == 1);
        }
        return kramer(c, b);
    }






    Float calculateDeterminant(const Matrix& matrix) {
        int n = matrix.size();
        if (n == 1) {
            return matrix[0][0];
        }
        if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        Float determinant = 0;
        for (size_t k = 0; k < n; ++k) {
            Matrix subMatrix(n - 1, Vector(n - 1));
            for (size_t i = 1; i < n; ++i) {
                int colIndex = 0;
                for (size_t j = 0; j < n; ++j) {
                    if (j == k) {
                        continue;
                    }
                    subMatrix[i - 1][colIndex] = matrix[i][j];
                    ++colIndex;
                }
            }

            determinant += (k % 2 == 0 ? 1 : -1) * matrix[0][k] * calculateDeterminant(subMatrix);
        }
        return determinant;
    }

    Vector kramer(const Matrix& c, const Vector& b) {
        int n = c.size();
        Float det = calculateDeterminant(c);
        Vector solutions(n);
        for (size_t i = 0; i < n; ++i) {
            Matrix tempMatrix = c;
            for (size_t j = 0; j < n; ++j) {
                tempMatrix[j][i] = b[j];
            }
            solutions[i] = calculateDeterminant(tempMatrix) / det;
        }

        return solutions;
    }



    void writeInFile(std::string fileName, Integer n, Float h, Float error) {
        std::ofstream file(fileName, std::ios::app);
        file << h << "," << error << std::endl;
        file.close();
    }





    Float binPow(Float x, Integer p) {
        if (!p) {
            return static_cast<Float>(1);
        }
        if (p % 2) {
            return binPow(x, p - 1) * x;
        }
        Float r = binPow(x, p / 2);
        return r * r;
    }
    std::vector<Integer> precalcFactorials() {
        std::vector<Integer> memo(std::max(N, M));
        memo[0] = 1;
        for (size_t i = 1; i < memo.size(); i++) {
            memo[i] = memo[i - 1] * i;
        }
        return memo;
    }
    Float factorial(size_t n) {
        static std::vector<Integer> memo = precalcFactorials();
        return static_cast<Float>(memo[n]);
    }
};


int main(){
    // testWrite();
    const size_t n = 3;
    const size_t m = 20;
    Task<n, m> t;
    t.Solve();
    return 0;
}
