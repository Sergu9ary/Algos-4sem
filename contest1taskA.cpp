#include <iostream>
#include <vector>
#include <string>

std::vector<int> pFunction(const std::string& str) {
    std::vector<int> prefix(str.size(), 0);
    for (size_t i = 1; i < str.size(); ++i) {
        int k = prefix[i - 1];
        while (k > 0 && str[i] != str[k]) {
            k = prefix[k - 1];
        }
        if (str[i] == str[k]) {
            prefix[i] = k + 1;
        }
    }
    return prefix;
}

std::vector<int> KMP(const std::string& text, const std::string& pattern) {
    std::vector<int> entry;
    std::string concatenated = pattern + "#" + text;
    std::vector<int> prefix = pFunction(concatenated);
    int patternLength = pattern.size();
    for (int i = 0; i < static_cast<int>(prefix.size()) - patternLength - 1; ++i) {
        if (prefix[i + patternLength + 1] == patternLength) {
            entry.push_back(i - patternLength + 1);
        }
    }
    return entry;
}

int main() {
    std::string text;
    std::string pattern;
    std::cin >> text >> pattern;
    std::vector<int> entry = KMP(text, pattern);
    for (int i = 0; i < static_cast<int>(entry.size()); ++i) {
        std::cout << entry[i] << '\n';
    }
    return 0;
}
