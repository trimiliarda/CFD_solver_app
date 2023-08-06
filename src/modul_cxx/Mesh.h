#include "Cell.h"

typedef struct Face {
  int nodes[2];
  int cr;
  int cl;
  bool is_bound;
  double length;
  point_t f_centr;
} face_t;

class Mesh {
private:
  /* data */
public:
  int Nx, Ny;
  int nNodes, nCells, nFaces;

  //  массивы
  point_t *nodes;
  face_t * faces;

  Mesh();
  ~Mesh();

  void load_struct_mesh(string filename);
  void create_cells(Cell *(&cells));
  void create_faces();
};
