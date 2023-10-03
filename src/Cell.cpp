#include "Cell.h"

Cell::Cell()
    : S(0), nFaces(0), nNodes(0), faces(0), nodes(0), fType(0), cells(0) {
  c.x = 0;
  c.y = 0;
}

Cell::~Cell() {}

void Cell::set_c(point_t c_) {
  c.x = c_.x;
  c.y = c_.y;
}

void Cell::set_nFaces(int nf) {
  nFaces = nf;
  faces = new int[nFaces];
  fType = new int[nFaces];
  cells = new int[nFaces];
}

void Cell::set_nNodes(int nn) {
  nNodes = nn;
  nodes = new int[nNodes];
}

void Cell::set_Nodes(int *nodes_, int nn) {
  nNodes = nn;
  nodes = new int[nNodes];
  for (int i = 0; i < nNodes; i++) {
    nodes[i] = nodes_[i];
  }
}

void Cell::show(int i) {
  string f = "../data/cells.txt";
  ofstream record(f, ios::out | ios::app);

  if (record) {
    record << "\t\t\tcell index = " << i << endl;
    record << "center = " << c.x << ", " << c.y << endl;
    record << "square = " << S << endl;
    record << "number of faces = " << nFaces << endl;
    record << "face indexes and types:" << endl;
    for (int j = 0; j < nFaces; j++) {
      record << "\tind = " << faces[j] << ", type = " << fType[j] << endl;
    }

    record << "number of nodes = " << nNodes << endl;
    record << "node indexes:" << endl;
    for (int j = 0; j < nNodes; j++) {
      record << "\tind = " << nodes[j] << endl;
    }

    record << "neighbour indexes: " << endl;
    for (int j = 0; j < nFaces; j++) {
      record << "\tind = " << cells[j] << endl;
    }
    record << "****************************************************" << endl;

  } else {
    cout << "Problem with file: " << f << endl;
  }
  record.close();
}
