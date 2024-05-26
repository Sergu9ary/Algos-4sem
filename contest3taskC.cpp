#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

const long double kDelta = 0.000001;
const long double kPi = 3.1415926535897932384626433832795;
const long double kAnglePi = 180;
const long double kMult = 0.5;
const long double kPrec = 1;
const long double kMax = 10e6;
const long double kRand1 = -156;
const long double kRand2 = 22.4;

struct Point;
class Line;

Point Intersection(const Line& l1, const Line& l2);
Line Ortogonal(const Point& pp, const Line& ll);
long double Dist(const Point& p1, const Point& p2);

long double RootCount(const long double& aa, const long double& bb) {
  long double res = std::pow(aa * aa - bb * bb, kMult);
  return res;
}

bool Compare(const long double& aa, const long double& bb) {
  return (aa - bb <= kDelta);
}

bool Equal(const long double& aa, const long double& bb) {
  return (std::abs(aa - bb) <= kDelta);
}

struct Point {
  long double x;
  long double y;
  Point() : x(0), y(0) {}
  Point(long double x_0, long double y_0) : x(x_0), y(y_0) {}
  Point(std::pair<long double, long double> pa2r)
      : x(pa2r.first), y(pa2r.second) {}
  bool operator==(const Point& pp) const {
    return (std::abs(x - pp.x) < kDelta && std::abs(y - pp.y) < kDelta);
  }
  bool operator!=(const Point& pp) const { return !(*this == pp); }
  Point& operator-=(const Point& pp);
  Point& operator+=(const Point& pp);
  void Rotate(const Point& pp, long double angle);
  void Reflect(const Point& pp);
  void Reflect(const Line& ll);
  void Scale(const Point& pp, long double kk);
  long double Dist(const Point& p1, const Point& p2);
  long double Module() const;
};

Point& Point::operator-=(const Point& pp) {
  x -= pp.x;
  y -= pp.y;
  return *this;
}

Point& Point::operator+=(const Point& pp) {
  x += pp.x;
  y += pp.y;
  return *this;
}

Point operator-(const Point& p1, const Point& p2) {
  Point result = p1;
  result -= p2;
  return result;
}

Point operator+(const Point& p1, const Point& p2) {
  Point result = p1;
  result += p2;
  return result;
}

void Point::Rotate(const Point& pp, long double angle) {
  x -= pp.x;
  y -= pp.y;
  long double temp = x;
  x = x * std::cos(kPi * angle / kAnglePi) -
      y * std::sin(kPi * angle / kAnglePi);
  y = y * std::cos(kPi * angle / kAnglePi) +
      temp * std::sin(kPi * angle / kAnglePi);
  x += pp.x;
  y += pp.y;
}

void Point::Reflect(const Point& pp) {
  x -= pp.x;
  y -= pp.y;
  x *= -1;
  y *= -1;
  x += pp.x;
  y += pp.y;
}

void Point::Scale(const Point& pp, long double kk) {
  x -= pp.x;
  y -= pp.y;
  x *= kk;
  y *= kk;
  x += pp.x;
  y += pp.y;
}

long double Point::Module() const { return std::pow(x * x + y * y, kMult); }
long double Dist(const Point& p1, const Point& p2) {
  return std::pow(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2), kMult);
}

long double DotProduct(const Point& p1, const Point& p2) {
  return p1.x * p2.x + p1.y * p2.y;
}

long double CrossProduct(const Point& p1, const Point& p2) {
  return p1.x * p2.y - p1.y * p2.x;
}
std::pair<long double, long double> SplitRelation(const Point& p1,
                                                  const Point& p2,
                                                  long double rel) {
  return std::make_pair((p1.x + rel * p2.x) / (rel + 1),
                        (p1.y + rel * p2.y) / (rel + 1));
}

class Line {
private:
  long double aa_;
  long double bb_;
  long double cc_;

public:
  Line() : aa_(1), bb_(1), cc_(0) {}
  Line(const Point& p1, const Point& p2) {
    aa_ = p2.y - p1.y;
    bb_ = p1.x - p2.x;
    cc_ = -(aa_ * p1.x + bb_ * p1.y);
  }
  Line(const Point& pp, long double kk) : aa_(kk), bb_(-1) {
    cc_ = pp.y - aa_ * pp.x;
  }
  Line(long double kk, long double shift) {
    *this = Line(Point(0, shift), Point(1, kk + shift));
  }
  bool operator==(const Line& ll) const {
    return (Equal(std::abs(aa_ / ll.aa_ - bb_ / ll.bb_), 0) &&
            Equal(std::abs(aa_ / ll.aa_ - cc_ / ll.cc_), 0));
  }
  bool operator!=(const Line& ll) const { return !(*this == ll); }
  std::string GetKoef() const {
    return std::to_string(aa_) + " " +
           std::to_string(bb_) + " " +
           std::to_string(cc_);
  }
  long double GetA() const { return aa_; }
  long double GetB() const { return bb_; }
  long double GetC() const { return cc_; }
  Point DirectVect() const;
};

Line Ortogonal(const Point& pp, const Line& ll) {
  return Line(pp, Point(ll.GetA() + pp.x, ll.GetB() + pp.y));
}

Point Intersection(const Line& l1, const Line& l2) {
  long double dd = l1.GetA() * l2.GetB() - l1.GetB() * l2.GetA();
  long double d1 = l1.GetC() * l2.GetB() - l1.GetB() * l2.GetC();
  long double d2 = l1.GetA() * l2.GetC() - l1.GetC() * l2.GetA();
  long double x = (d1 / dd) * (-1);
  long double y = (d2 / dd) * (-1);
  return Point(x, y);
}

void Point::Reflect(const Line& ll) {
  Point intr = Intersection(ll, Ortogonal(*this, ll));
  x += 2 * (intr.x - x);
  y += 2 * (intr.y - y);
}

Point Line::DirectVect() const {
  long double x = -GetB();
  long double y = GetA();
  long double ab = Dist(Point(0, 0), Point(x, y));
  x /= ab;
  y /= ab;
  return Point(x, y);
}

long double Angle(const Point& p1, const Point& p2, const Point& center) {
  return std::abs(std::atan2(CrossProduct(p1 - center, p2 - center),
                             DotProduct(p1 - center, p2 - center)));
}

long double DistToLine(const Line& ll, const Point& pp) {
  Line l1 = Ortogonal(pp, ll);
// std::cout << l1.GetKoef() << " getkoef " << "\n"; 
  Point intr = Intersection(ll, l1);
  return Dist(intr, pp);
}

long double DistToRay(const Point& st, const Point& en, const Point& pp) {
  Line ll(st, en);
  Line l1 = Ortogonal(pp, ll);
  Point intr = Intersection(ll, l1);
  if (Compare(0, DotProduct(en - st, intr - st))) {
    return Dist(pp, intr);
  }
  return Dist(pp, st);
}

long double DistToSegment(const Point& st, const Point& en, const Point& pp) {
  Line ll(st, en);
  Line l1 = Ortogonal(pp, ll);
  Point intr = Intersection(ll, l1);
  if (Compare(DotProduct(intr - en, intr - st), 0)) {
    return Dist(pp, intr);
  }
  return std::min(Dist(pp, en), Dist(pp, st));
}

long double DistBetwSegments(const Point& a1, const Point& a2, const Point& b1,
                             const Point& b2) {
  Line l1(a1, a2);
  Line l2(b1, b2);
  Point intr = Intersection(l1, l2);
  if (Compare(DotProduct(intr - a1, intr - a2), 0) &&
      Compare(DotProduct(intr - b1, intr - b2), 0)) {
    return 0;
  }
  long double res = kMax;
  res = std::min(res, DistToSegment(a1, a2, b1));
  res = std::min(res, DistToSegment(a1, a2, b2));
  res = std::min(res, DistToSegment(b1, b2, a1));
  res = std::min(res, DistToSegment(b1, b2, a2));
  return res;
}

class Polygon {
  std::vector<Point> vertices_ = {};

public:
  Polygon(std::vector<Point> vp) : vertices_(vp) {}
  template <typename... Point>
  Polygon(Point... args) : vertices_({args...}) {}
  ~Polygon() = default;
  std::vector<Point>& GetVertices();
  int VerticesCount() const { return vertices_.size(); }
  long double Perimeter() const {
    long double perimeter = 0;
    for (int i = 0; i < VerticesCount(); ++i) {
      perimeter += Dist(Getv(i), Getv(i + 1));
    }
    return perimeter;
  }

  long double Area() const {
    long double area = 0;
    for (int i = 1; i < VerticesCount() - 1; ++i) {
      area += CrossProduct(vertices_[i] - vertices_[0],
                           vertices_[i + 1] - vertices_[0]);
    }
    area = std::abs(area) / 2;
    return area;
  }

  void Rotate(const Point& pp, double angle) {
    for (int i = 0; i < VerticesCount(); ++i) {
      vertices_[i].Rotate(pp, angle);
    }
  }

  void Reflect(const Point& pp) {
    for (int i = 0; i < VerticesCount(); ++i) {
      vertices_[i].Reflect(pp);
    }
  }

  void Reflect(const Line& ll) {
    for (int i = 0; i < VerticesCount(); ++i) {
      vertices_[i].Reflect(ll);
    }
  }

  void Scale(const Point& pp, long double kk) {
    for (int i = 0; i < VerticesCount(); ++i) {
      vertices_[i].Scale(pp, kk);
    }
  }

  bool ContainsPoint(const Point& pp) const {
    for (auto to : vertices_) {
      if (to == pp) {
        return true;
      }
    }
    for (int i = 0; i < VerticesCount(); ++i) {
      if (Equal(CrossProduct(vertices_[i] - pp,
                             pp - vertices_[(i + 1) % vertices_.size()]),
                0) &&
          DotProduct(vertices_[i] - pp,
                     pp -
                     vertices_[(i + 1) % vertices_.size()]) > 0) {
        return true;
      }
    }
    Point anotherp = pp + Point(kRand1, kRand2);
    Line ll = Line(pp, anotherp);
    int intersections = 0;
    for (int i = 0; i < VerticesCount(); ++i) {
      Point inter = Intersection(
          ll, Line(vertices_[i], vertices_[(i + 1) % vertices_.size()]));
      if (DotProduct(vertices_[i] - inter,
                     vertices_[(i + 1) % vertices_.size()] - inter) <= 0 &&
          DotProduct(inter - pp, anotherp - pp) >= 0) {
        ++intersections;
      }
    }
    return (intersections % 2 != 0);
  }
  Point Getv(int ind) const;
  bool operator==(const Polygon& pol) const;
  bool IsConvex() const;
  Point& operator[](int ind) { return vertices_[ind]; }
  const Point& operator[](int ind) const { return vertices_[ind]; }
};

Point Polygon::Getv(int ind) const {
  int cc = (ind + VerticesCount()) % VerticesCount();
  return vertices_[cc];
}

bool Polygon::operator==(const Polygon& pol) const {
  if (VerticesCount() != pol.VerticesCount()) {
    return false;
  }

  for (int shift = 0; shift < VerticesCount(); ++shift) {
    bool is_true1 = true;
    bool is_true2 = true;
    for (int i = 0; i < VerticesCount(); ++i) {
      if (Getv(i) != pol.Getv(i + shift)) {
        is_true1 = false;
      }
      if (Getv(i) != pol.Getv(shift - i)) {
        is_true2 = false;
      }
      if (!is_true1 && !is_true2) {
        break;
      }
    }
    if (is_true1 || is_true2) {
      return true;
    }
  }
  return false;
}

std::vector<Point>& Polygon::GetVertices() { return vertices_; }

bool Polygon::IsConvex() const {
  long double cross = CrossProduct(
      vertices_[VerticesCount() - 1] - vertices_[VerticesCount() - 2],
      vertices_[0] - vertices_[VerticesCount() - 1]);
  long double cross1 =
      CrossProduct(vertices_[0] - vertices_[VerticesCount() - 1],
                   vertices_[1] - vertices_[0]);
  if (cross1 * cross < 0) {
    return false;
  }
  for (int i = 1; i < VerticesCount() - 1; ++i) {
    long double cross2 = CrossProduct(vertices_[i] - vertices_[i - 1],
                                      vertices_[i + 1] - vertices_[i]);
    if (cross * cross2 < 0) {
      return false;
    }
  }
  return true;
}

class AngleSort {
  Point base_;

public:
  AngleSort(Point pp) : base_(pp) {}
  bool operator()(const Point& aa, const Point& bb) {
    if (bb == base_) {
      return false;
    }
    if (aa == base_) {
      return true;
    }
    if (CrossProduct(bb - base_, aa - base_) >= 0) {
      if (CrossProduct(bb - base_, aa - base_) == 0) {
        return (Dist(aa, base_) < Dist(bb, base_));
      }
      return true;
    }
    return false;
  }
};

int FindInd(const Polygon& conv_pol) {
  int index = 0;
  Point p0 = conv_pol.Getv(0);
  for (int i = 1; i < conv_pol.VerticesCount(); ++i) {
    if (conv_pol.Getv(i).x == p0.x && conv_pol.Getv(i).y < p0.y) {
      index = i;
      p0 = conv_pol.Getv(i);
    }
    if (conv_pol.Getv(i).x < p0.x) {
      index = i;
      p0 = conv_pol.Getv(i);
    }
  }
  return index;
}

Polygon BuildConvexHull(const Polygon& pol) {
  Polygon conv_pol(pol);
  int index = 0;
  index = FindInd(conv_pol);
  std::swap(conv_pol[index], conv_pol[0]);
  AngleSort cmp(conv_pol.Getv(0));
  std::sort(conv_pol.GetVertices().begin() + 1, conv_pol.GetVertices().end(),
            cmp);
  int ind = 0;
  std::vector<Point> hull_base;
  hull_base.push_back(conv_pol[0]);
  conv_pol.GetVertices().push_back(conv_pol[0]);
  while (conv_pol[ind] == conv_pol[0]) {
    ++ind;
  }
  hull_base.push_back(conv_pol[ind]);
  ++ind;
  Point p1 = hull_base[hull_base.size() - 1];
  Point p2 = hull_base[hull_base.size() - 2];
  for (; ind < conv_pol.VerticesCount(); ++ind) {
    while (hull_base.size() >= 2 && ind < conv_pol.VerticesCount() &&
           CrossProduct(p1 - p2, conv_pol[ind] - p1) >= 0) {
      if (p1 == conv_pol[ind]) {
        break;
      }
      if (CrossProduct(p1 - p2, conv_pol[ind] - p1) == 0) {
        if (Dist(p2, p1) < Dist(p2, conv_pol[ind])) {
          hull_base.pop_back();
        }
      } else {
        hull_base.pop_back();
      }
      p1 = hull_base[hull_base.size() - 1];
      p2 = hull_base[hull_base.size() - 2];
    }
    hull_base.push_back(conv_pol[ind]);
    p1 = hull_base[hull_base.size() - 1];
    p2 = hull_base[hull_base.size() - 2];
  }
  hull_base.pop_back();
  Polygon pl = Polygon(hull_base);
  return pl;
}

int
main() {
  std::cin.tie(NULL);
  int nn;
  std::cin >> nn;
  long double aa;
  long double bb;
  std::vector<Point> pts;
  for (int i = 0; i < nn; ++i) {
    std::cin >> aa >> bb;
    pts.push_back(Point(aa, bb));
  }
  Polygon pol(pts);
  std::cout << "\n";
  Polygon hull(BuildConvexHull(pol));
  std::cout.setf(std::ios::dec);
  std::cout.precision(kPrec);
  std::cout << hull.VerticesCount() << "\n";
  for (int i = 0; i < hull.VerticesCount(); ++i) {
    std::cout << std::to_string(int(hull[i].x)) << " "
    << std::to_string(int(hull[i].y)) << "\n";
  }
  std::cout << std::to_string(hull.Area());
}
