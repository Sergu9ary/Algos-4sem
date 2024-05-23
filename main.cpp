#include <iostream>
#include <vector>
#include <unordered_set>

using namespace std;

struct Dot {
  int x, y;
  Dot& operator+=(const Dot& other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  Dot& operator-=(const Dot& other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  bool operator==(const Dot& other) const {
    return x == other.x && y == other.y;
  }
};

struct Cmp {
  size_t operator()(const Dot& dot) const {
    return hash<int>()(dot.x) ^ hash<int>()(dot.y);
  }
};

int Solution(const vector<Dot>& points) {
  unordered_set<Dot, Cmp> s(points.begin(), points.end());
  int count = 0;

  for (const auto& i : points) {
    for (const auto& j : points) {
      if (i == j) continue;
      Dot tmp = j;
      tmp -= i;
      Dot normal = {tmp.y, -tmp.x};

      Dot d = i;
      d += normal;
      Dot dd = j;
      dd += normal;

      if (s.count(d) && s.count(dd)) {
        count++;
      }
    }
  }
  return count / 4;
}

int main() {
  int n;
  cin >> n;
  vector<Dot> points(n);

  for (int i = 0; i < n; ++i) {
    cin >> points[i].x >> points[i].y;
  }

  cout << Solution(points) << endl;

  return 0;
}
