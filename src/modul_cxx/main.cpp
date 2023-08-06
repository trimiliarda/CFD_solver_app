// #include "Polygon.h"
// #include "Cell.h"
#include "Mesh.h"

// #define und <<","<<
// #define end <<std::endl

using namespace std;

int main() {
  Polygon pl;
  pl.load_data("../data/point.txt");

  // for (int i = 0; i < pl.get_n(); i++) {
  //   point_t p = pl.get_p(i);
  //   cout << p.x << " " << p.y << endl;
  // }

  Mesh mesh;
  mesh.load_struct_mesh("../data/grid.txt");

  int nCells = mesh.nCells;
  cout << "nCells = " << nCells << endl;

  Cell *cells = new Cell[nCells];
  mesh.create_cells(cells);

  for (int i = 0; i < nCells; i++) {
    cells[i].show(i);
  }

  return 0;
}
