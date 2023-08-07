#include "Polygon.h"
#include "Cell.h"
#include "Mesh.h"
#include "Functions.h"

// #define und <<","<<
// #define end <<std::endl
#define ItMAX 100
#define EPS 1.e-6

using namespace std;

int main() {
  // Polygon pl;
  // pl.load_data("../data/point.txt");

  // for (int i = 0; i < pl.get_n(); i++) {
  //   point_t p = pl.get_p(i);
  //   cout << p.x << " " << p.y << endl;
  // }

  gen_mesh(100, 64);

  Mesh mesh;
  mesh.load_struct_mesh("../data/grid.txt");
  // mesh.load_struct_mesh("../data/mesh.dat");

  int nCells = mesh.nCells;
  cout << "nCells = " << nCells << endl;

  Cell *cells = new Cell[nCells];

  mesh.create_cells(cells);

  mesh.create_faces();

  mesh.cell_funcs(cells);

  remove("../data/cells.txt");
  remove("../data/faces.txt");

  for (int i = 0; i < nCells; i++) {
    cells[i].show(i);
  }
  for (int i = 0; i < mesh.nFaces; i++) {
    show(mesh.faces[i]);
  }

  param_t * p = new param_t[nCells];
  changes * du = new changes[nCells];

  // initialization
  int Nm = 1;
  Init(p, nCells, Nm);

  int It = 0;
  double res = 1.;

  for (int i = 0; i < nCells; i++) {
    du[i].dU = new double[Nm];
  }

  while (It++ < ItMAX && res > EPS) {
    // cout << "Iter = " << It << ", res = " << res << endl;
    
    for (int i = 0; i < nCells; i++) {
      for (int j = 0; j < Nm; j++) {
        p[i].U[j] = p[i].U1[j];
        du[i].dU[j] = 0.;

      }
    }
    double dt = 0.01;
    //  приращение за счёт вязкости 
    Viscous(p, du, mesh, cells, dt);

    for (int i = 0; i < nCells; i++) {
      for (int j = 0; j < Nm; j++) {
        p[i].U1[j] = p[i].U[j] + du[i].dU[j];
      }
    }

    get_params(p, nCells, Nm);

    res = 0.;
    for (int i = 0; i < nCells; i++) {
      double res_ = abs(du[i].dU[0] / p[i].U1[0]);
      res = res < res_ ? res_ : res;
    }

  }
  
  cout << "Iter = " << --It << ", res = " << res << endl;

  int Nx = mesh.Nx;
  int Ny = mesh.Ny;
  
  plot(p, cells, Nx, Ny, nCells);

  return 0;
}
