// #include "Polygon.h"
#include "Cell.h"

// #define und <<","<<
// #define end <<std::endl

using namespace std;

int main() {
  Polygon pl;

  pl.DataEntry("../data/point.txt");

  for (int i = 0; i < pl.get_n(); i++) {
    point_t p = pl.get_p(i);
    cout << p.x << " " << p.y << endl;
  }

  Cell cell;
  cell.set_c(pl.Center());
  cell.set_S(pl.Square());

  cout << cell.get_c().x << " " << cell.get_c().y << endl;

  cout << cell.get_S() << endl;

  int n = pl.get_n();
  cell.set_nFaces(n);
  cout << cell.get_Face(3);
  cell.set_nNodes(n);
  cout << cell.get_Node(3);
  return 0;
}
