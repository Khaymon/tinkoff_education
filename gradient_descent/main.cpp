#include <iostream>
#include <tuple>

const double EPS = 1e-9;

double Power(double x, int power) {
  if (power == 0) {
    return 1;
  }
  return x * Power(x, power - 1);
}

struct Point {
  double x, y;

  Point(double x_, double y_): x(x_), y(y_) {}

  Point operator*(double value) const {
    return Point(this->x * value, this->y * value);
  }

  Point operator+(Point rhs) const {
    return Point(this->x + rhs.x, this->y + rhs.y);
  }

  Point operator-() const {
    return Point(-this->x, -this->y);
  }

  Point operator+=(Point rhs) {
    this->x += rhs.x;
    this->y += rhs.y;
    return *this;
  }
  double SquaredDistance(Point &rhs) {
    return Power(this->x - rhs.x, 2) + Power(this->y - rhs.y, 2);
  }
};

double Function(Point point) {
	return Power(point.x, 6) * 6 + Power(point.x, 4) * Power(point.y, 2) * 2 + Power(point.x, 2) * 10 +
      point.x * point.y * 6 + Power(point.y, 2) * 10 - point.x * 6 + 4;
}

double FuctionXDerivative(Point point) {
  return Power(point.x, 5) * 36 + Power(point.x, 3) * Power(point.y, 2) * 8 + point.x * 20 + point.y * 6 - 6;
}

double FunctionYDerivative(Point point) {
  return Power(point.x, 4) * point.y * 4 + point.x * 6 + point.y * 20;
}

// Step doubling method for the solution localization
std::tuple<double, double> FindPivotPoints(Point point, Point anti_gradient) {
  double step = 0.1;
  double left_point = 0;
  double right_point = left_point + step;

  while (Function(point + anti_gradient * left_point) > Function(point + anti_gradient * right_point)) {
    left_point = right_point;
    step *= 2;
    right_point += step;
  }

  return std::tie(left_point, right_point);
}

double DihotomieMethod(Point current_point, Point anti_gradient, double left_point, double right_point) {
  while (right_point - left_point >= EPS) {
    double m1 = left_point + (right_point - left_point) / 3;
    double m2 = left_point + (right_point - left_point) * 2 / 3;
    if (Function(current_point + anti_gradient * m1) >= Function(current_point + anti_gradient * m2)) {
      left_point = m1;
    } else {
      right_point = m2;
    }
  }

  return (left_point + right_point) / 2;
}

Point GradientDescent(Point current_point) {
  Point previous_point = current_point;
  do {
    Point anti_gradient(-FuctionXDerivative(current_point), -FunctionYDerivative(current_point));
    double left_point = 0;
    double right_point = 0;
    std::tie(left_point, right_point) = FindPivotPoints(current_point, anti_gradient);
    double optimum_lambda = DihotomieMethod(current_point, anti_gradient, left_point, right_point);
    previous_point = current_point;
    current_point += anti_gradient * optimum_lambda;
  } while (current_point.SquaredDistance(previous_point) >= EPS);

  return current_point;
}

int main() {
  Point optimum(0, 0);

  optimum = GradientDescent(optimum);
  std::cout << Function(optimum) << std::endl;

  return 0;
}