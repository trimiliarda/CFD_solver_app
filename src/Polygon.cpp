#include "Polygon.h"

using namespace std;

Polygon::Polygon() : n(0) { p = new point_t[0]; }

Polygon::Polygon(int n1) : n(n1) { p = new point_t[n]; }

Polygon::Polygon(point_t *p1, int n1) : n(n1) {
  p = new point_t[n];
  for (int i = 0; i < n; i++)
    p[i] = p1[i];
}

Polygon::~Polygon() { delete[] p; }

point_t Polygon::Center() {
  int m = n + 1;
  point_t *pm = new point_t[m];
  for (int i = 0; i < n; i++) {
    pm[i].x = p[i].x;
    pm[i].y = p[i].y;
  }
  pm[n].x = p[0].x;
  pm[n].y = p[0].y;

  double A = 0., Cx = 0., Cy = 0.;
  for (int i = 0; i < n; i++) {
    A += pm[i].x * pm[i + 1].y - pm[i + 1].x * pm[i].y;
    Cx += (pm[i].x + pm[i + 1].x) *
          (pm[i].x * pm[i + 1].y - pm[i + 1].x * pm[i].y);
    Cy += (pm[i].y + pm[i + 1].y) *
          (pm[i].x * pm[i + 1].y - pm[i + 1].x * pm[i].y);
  }
  A *= 3.;
  point_t answer;
  answer.x = Cx / A;
  answer.y = Cy / A;

  return answer;
}

double Polygon::Square() {
  int m = n + 1;
  point_t *pm = new point_t[m];
  for (int i = 0; i < n; i++) {
    pm[i].x = p[i].x;
    pm[i].y = p[i].y;
  }
  pm[n].x = p[0].x;
  pm[n].y = p[0].y;

  double A = 0.;
  for (int i = 0; i < n; i++) {
    A += pm[i].x * pm[i + 1].y - pm[i + 1].x * pm[i].y;
  }
  A /= 2.;

  return abs(A);
}

void Polygon::load_data(string filename) {
  int n1;
  point_t *p1;
  ifstream reading(filename);

  if (reading) {
    cout << "Open file: " << filename << endl;

    reading >> n;
    p = new point_t[n];
    for (int i = 0; i < n; i++) {
      reading >> p[i].x >> p[i].y;
    }

  } else {
    cout << "Error: file \"" << filename << "\" can\'t open" << endl;
  }
}
