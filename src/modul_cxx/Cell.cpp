#include "Cell.h"

Cell::Cell()
    : S(0), nFaces(0), nNodes(0), faces(0), nodes(0), type_f(0), cells(0) {
  c.x = 0;
  c.y = 0;
  cout << " i was born";
}

Cell::~Cell() {}

void Cell::set_c(point_t c_) {
  c.x = c_.x;
  c.y = c_.y;
}

void Cell::set_nFaces(int nf) {
  nFaces = nf;
  faces = new int[nFaces];
}

void Cell::set_nNodes(int nn) {
  nNodes = nn;
  nodes = new int[nNodes];
}
