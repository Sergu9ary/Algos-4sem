#include <iostream>
#include <vector>

std::string preprocess(const std::string& s) {
  std::string processed = "#";
  for (char c : s) {
    processed += c;
    processed += "#";
  }
  return processed;
}

size_t countPalindromes(const std::string& s) {
  std::string processed = preprocess(s);
  size_t n = processed.size();
  std::vector<size_t> palindromeCount(n, 0);
  size_t center = 0;
  size_t radius = 0;

  for (size_t i = 0; i < n; ++i) {
    size_t mirror = 2 * center - i;
    if (i < center + radius) {
      palindromeCount[i] = std::min(radius - (i - center), palindromeCount[mirror]);
    }

    while (i - palindromeCount[i] - 1 >= 0 && i + palindromeCount[i] + 1 < n && processed[i - palindromeCount[i] - 1] == processed[i + palindromeCount[i] + 1]) {
      ++palindromeCount[i];
    }
    if (i + palindromeCount[i] > center + radius) {
      center = i;
      radius = palindromeCount[i];
    }
  }
  size_t totalPalindromes = 0;
  for (size_t count : palindromeCount) {
    totalPalindromes += count / 2;
  }
  return totalPalindromes;
}

int main() {
  std::string s;
  std::cin >> s;
  std::cout << countPalindromes(s) << "\n";
  return 0;
}
