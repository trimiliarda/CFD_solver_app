#include "Polygon.h"

#define INNER 0
#define BOUND 1

using namespace std;

class Cell {
private:
  point_t c; //  координаты ц.т.
  double S;  //  площадь

  int nFaces; //  число граней
  int nNodes; //  число  узлов

  int *faces; //  номера граней
  int *nodes; //  номера узлов

  int *type_f; //  тип граней: внутр \ граничная

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
};
