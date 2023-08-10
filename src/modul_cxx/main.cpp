#include "Cell.h"
#include "Functions.h"
#include "Mesh.h"
#include "Polygon.h"

// #define und <<","<<
// #define end <<std::endl
#define ItMAX 50000
#define EPS 1.e-6

using namespace std;

int main() {
  // gen_mesh(64, 50);

  Mesh mesh;
  // mesh.load_struct_mesh("../data/grid.txt");
  mesh.load_struct_mesh("../data/mesh.dat");

  int nCells = mesh.nCells;
  cout << "nCells = " << nCells << endl;

  Cell *cells = new Cell[nCells];

  mesh.create_cells(cells);

  mesh.create_faces();

  mesh.cell_funcs(cells);

  // remove("../data/cells.txt");
  // remove("../data/faces.txt");

  // for (int i = 0; i < nCells; i++) {
  //   cells[i].show(i);
  // }

  param_t *p = new param_t[nCells];
  changes *du = new changes[nCells];

  //
  // initialization
  //
  int Nm = 1;
  Init(p, nCells, Nm);

  mesh.set_zones();

  // int nZone = 4;
  // zone_t * z = new zone_t[nZone];
  set_gran(mesh);

  Yw(mesh, cells, nCells);

  // mesh.set_zones();

  int It = 0;
  double res = 1.;

  for (int i = 0; i < nCells; i++) {
    du[i].dU = new double[Nm];
  }
  //////////////////////////////////////////////////////////////////////////////////
  while (It++ < ItMAX && res > EPS) {
    // cout << "Iter = " << It << ", res = " << res << endl;

    for (int i = 0; i < nCells; i++) {
      for (int j = 0; j < Nm; j++) {
        p[i].U[j] = p[i].U1[j];
        du[i].dU[j] = 0.;
      }
    }

    Gradients(cells, mesh, gr, p, Nm);

    double dt = 7.e+2;


    //  приращение за счёт невязких потоков
    // Convect(p, du, mesh, cells, It, dt);

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
