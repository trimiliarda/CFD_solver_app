#pragma once
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef struct Point {
  double x, y;
} point_t;

class Polygon {
private:
  point_t *p;
  int n;

public:
  Polygon();
  Polygon(int n1);
  Polygon(point_t *p1, int n1);
  ~Polygon();

  void load_data(string filename);
  point_t Center();
  double Square();

  int get_n() { return n; }
  point_t get_p(int i) { return p[i]; }

  void set_p(int i, point_t pnt) { p[i] = pnt; }
};
