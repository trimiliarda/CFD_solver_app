#pragma once
#include "Cell.h"

#define INNER -1
#define BOUND_Left 1
#define BOUND_Bottom 2
#define BOUND_Right 3
#define BOUND_Up 4

typedef struct Face {
  int nodes[2];
  int cr;
  int cl;
  bool is_bound;
  double length;
  point_t f_centr;
  int zone;   // номер зоня для граничных условий

  void show(int k) {
  cout << "\t\tface index = " << k << endl;
  cout << "\tfaces nodes[0] = " << nodes[0] << endl;
  cout << "\tfaces nodes[1] = " << nodes[1] << endl;
  cout << "cr = " << cr << ", cl = " << cl << endl;
  cout << "is_bound = " << is_bound << endl;
  cout << "length = " << length << endl;
  cout << "f_centr: x = " << f_centr.x << ", y = " << f_centr.y << endl;
  if (zone == INNER)
    cout << "\t\tzone: INNER = " << zone << endl;
  else if (zone == BOUND_Left)
    cout << "\t\tzone: BOUND_Left = " << zone << endl;
  else if (zone == BOUND_Right)
    cout << "\t\tzone: BOUND_Right = " << zone << endl;
  else if (zone == BOUND_Bottom)
    cout << "\t\tzone: BOUND_Bottom = " << zone << endl;
  else if (zone == BOUND_Up)
    cout << "\t\tzone: BOUND_Up = " << zone << endl;
  cout << "***************************************************" << endl;
  }
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
  void cell_funcs(Cell *(&cells));
};

void show(face_t f);
