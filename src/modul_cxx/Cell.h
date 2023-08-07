#pragma once
#include "Polygon.h"


using namespace std;

class Cell {
private:
  point_t c; //  координаты ц.т.
  double S;  //  площадь

  int nFaces; //  число граней
  int nNodes; //  число  узлов

  int *faces; //  номера граней
  int *nodes; //  номера узлов

  int *fType; //  тип граней: внутр \ граничная

  int *cells; //  номкра соседних ячеек

public:
  Cell();
  ~Cell();

  point_t get_c() { return c; }
  double get_S() { return S; }
  int get_nFaces() { return nFaces; }
  int get_nNodes() { return nNodes; }
  int get_Node(int i) { return i < nNodes ? nodes[i] : -1; }
  int get_Face(int i) { return i < nFaces ? faces[i] : -1; }

  void set_c(point_t c_);
  void set_S(double S_) { S = S_; }
  void set_nFaces(int nf);
  void set_nNodes(int nn);
  void set_Nodes(int *nodes_, int nn);
  void set_Face(int iFace, int fIndex) { faces[iFace] = fIndex; }
  void set_fType(int iFace, int new_fType) { fType[iFace] = new_fType; }
  void set_cells(int iFace, int cell) { cells[iFace] = cell; }

  void show(int i);
};
