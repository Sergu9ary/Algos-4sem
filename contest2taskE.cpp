#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numbers>

const double PI = acos(-1);

size_t next_power_of_two(size_t n) {
  size_t power = 1;
  while (power < n) {
    power *= 2;
  }
  return power;
}

void FFT(std::vector<std::complex<long double>>& values, bool invert) {
  size_t n = values.size();
  int lg_n = std::floor(std::log2(n));

  std::vector<int> rev(n);
  rev[0] = 0;
  int oldest = -1;

  for (size_t i = 1; i < n; ++i) {
    if ((i & (i - 1)) == 0) {
      ++oldest;
    }
    rev[i] = rev[i ^ (1 << oldest)] | (1 << (lg_n - oldest - 1));
  }

  for (size_t i = 0; i < n; ++i) {
    if (i < rev[i]) {
      std::swap(values[i], values[rev[i]]);
    }
  }

  for (size_t len = 2; len <= n; len <<= 1) {
    long double ang = 2 * PI / len * (invert ? -1 : 1);
    std::complex<long double> root(cosl(ang), sinl(ang));
    for (size_t i = 0; i < n; i += len) {
      std::complex<long double> cur_root(1);
      for (size_t j = 0; j < len / 2; ++j) {
        std::complex<long double> left = values[i + j];
        std::complex<long double> right = values[i + j + len / 2] * cur_root;
        values[i + j] = left + right;
        values[i + j + len / 2] = left - right;
        cur_root *= root;
      }
    }
  }

  if (invert) {
    for (auto& val : values) {
      val /= n;
    }
  }
}

int Multiply(const std::vector<int>& lhs, const std::vector<int>& rhs, std::vector<int>& result) {
  size_t n = next_power_of_two(lhs.size() + rhs.size() - 1);

  std::vector<std::complex<long double>> fur_lhs(lhs.begin(), lhs.end());
  std::vector<std::complex<long double>> fur_rhs(rhs.begin(), rhs.end());

  fur_lhs.resize(n);
  fur_rhs.resize(n);

  FFT(fur_lhs, false);
  FFT(fur_rhs, false);

  for (size_t i = 0; i < n; ++i) {
    fur_lhs[i] *= fur_rhs[i];
  }

  FFT(fur_lhs, true);

  result.resize(n);
  int max_degree = 0;

  for (size_t i = 0; i < n; ++i) {
    result[i] = std::lround(fur_lhs[i].real());
  }
  for (size_t i = 0; i < n; ++i) {
    if (result[n - 1 - i] != 0) {
      max_degree = n - 1 - i;
      break;
    }
  }

  return max_degree;
}

int main() {
  int degree1;
  std::cin >> degree1;

  std::vector<int> first(degree1 + 1);
  for (int i = 0; i <= degree1; ++i) {
    std::cin >> first[degree1 - i];
  }

  int degree2;
  std::cin >> degree2;

  std::vector<int> second(degree2 + 1);
  for (int i = 0; i <= degree2; ++i) {
    std::cin >> second[degree2 - i];
  }

  std::vector<int> result;

  int max_degree = Multiply(first, second, result);

  std::cout << max_degree << ' ';
  for (int i = max_degree; i >= 0; --i) {
    std::cout << result[i] << " ";
  }

  std::cout << std::endl;

  return 0;
}
