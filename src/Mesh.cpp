#include <algorithm>
#include <math.h>

#include "Mesh.h"

using namespace std;

Mesh::Mesh() : Nx(0), Ny(0), nNodes(0), nCells(0), nFaces(0), nodes(0) {}

Mesh::~Mesh() {}

void Mesh::load_struct_mesh(string filename) {
  int z, k = 0;

  ifstream reading(filename);

  if (reading) {
    cout << "Open file: " << filename << endl;

    reading >> Nx >> Ny >> z >> z;
    nNodes = Nx * Ny;
    nFaces = Nx * (Ny - 1) + Ny * (Nx - 1);
    nCells = (Nx - 1) * (Ny - 1);

    nodes = new point_t[nNodes];

    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        reading >> z >> z >> nodes[k].x >> nodes[k].y;
        ++k;
      }
    }

  } else {
    cout << "Error: file \"" << filename << "\" can\'t open" << endl;
  }
  reading.close();
}

void Mesh::create_cells(Cell *(&cells)) {
  for (int i = 0; i < Nx - 1; i++) {
    for (int j = 0; j < Ny - 1; j++) {
      int nc = (Ny - 1) * i + j;

      // cout << "nc = " << nc << endl;

      int nNodes = 4; //  число узлов вокруг ячейки
      int *nodes = new int[nNodes];

      //  элемент под номером 0
      int n2 = Ny * i + j;
      nodes[0] = n2;

      //  элемент под номером 1
      n2 = Ny * (i + 1) + j;
      nodes[1] = n2;

      //  элемент под номером 2
      n2 = Ny * (i + 1) + (j + 1);
      nodes[2] = n2;

      //  элемент под номером 3
      n2 = Ny * i + (j + 1);
      nodes[3] = n2;

      cells[nc].set_Nodes(nodes, nNodes);
    }
  }
}

void Mesh::create_faces() {
  faces = new Face[nFaces];
  int k = 0;
  //  вертикальный грани
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny - 1; j++) {
      int n1 = Ny * i + j;
      int n2 = Ny * i + j + 1;
      faces[k].nodes[0] = n1;
      faces[k].nodes[1] = n2;

      if (i == 0) {
        faces[k].is_bound = true;
        faces[k].cr = -1;
        faces[k].cl = (Ny - 1) * i + j;
        faces[k].zone = BOUND_Left;
      } else if (i == Nx - 1) {
        faces[k].is_bound = true;
        faces[k].cr = (Ny - 1) * (i - 1) + j;
        faces[k].cl = -1;
        faces[k].zone = BOUND_Right;
      } else {
        faces[k].is_bound = false;
        faces[k].cr = (Ny - 1) * (i - 1) + j;
        faces[k].cl = (Ny - 1) * i + j;
        faces[k].zone = INNER;
      }
      faces[k].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
      faces[k].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

      double dx = (nodes[n1].x - nodes[n2].x);
      double dy = (nodes[n1].y - nodes[n2].y);

      faces[k].length = sqrt(dx * dx + dy * dy);

      // faces[k].show(k);
      ++k;
    }
  }
  //  горизонтальные грани
  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx - 1; i++) {
      int n1 = Ny * i + j;
      int n2 = Ny * (i + 1) + j;
      faces[k].nodes[0] = n1;
      faces[k].nodes[1] = n2;

      if (j == 0) {
        faces[k].is_bound = true;
        faces[k].cr = (Ny - 1) * i + j;
        faces[k].cl = -1;
        faces[k].zone = BOUND_Bottom;
      } else if (j == Ny - 1) {
        faces[k].is_bound = true;
        faces[k].cr = -1;
        faces[k].cl = (Ny - 1) * i + j - 1;
        faces[k].zone = BOUND_Up;
      } else {
        faces[k].is_bound = false;
        faces[k].cr = (Ny - 1) * i + j;
        faces[k].cl = (Ny - 1) * i + j - 1;
        faces[k].zone = INNER;
      }
      faces[k].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
      faces[k].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

      double dx = (nodes[n1].x - nodes[n2].x);
      double dy = (nodes[n1].y - nodes[n2].y);

      faces[k].length = sqrt(dx * dx + dy * dy);

      // faces[k].show(k);
      ++k;
    }
  }
}

void Mesh::grad_coeffs(Cell *(&cells)) {
  for (int i = 0; i < nCells; i++) {
    int nFaces = cells[i].get_nFaces();

    //  создаём массивы весовых коэффициентов wk и векторов ck
    cells[i].wk = new double[nFaces];
    cells[i].ck = new vector_t[nFaces];

    int i_dim = 2;

    for (int k = 0; k < nFaces; k++) {
      cells[i].ck[k].cx = new double[i_dim];
    }

    point_t xc = cells[i].get_c();

    double *dx = new double[nFaces];
    double *dy = new double[nFaces];

    double axx = 0., axy = 0., ayy = 0.;

    for (int k = 0; k < nFaces; k++) {
      int nf = cells[i].get_Face(k);
      if (!faces[nf].is_bound) {
        int nc = cells[i].get_Cell(k);
        point_t xk = cells[i].get_c();

        dx[k] = xk.x - xc.x;
        dy[k] = xk.y - xc.y;

      } else {

        point_t xk = faces[nf].f_centr;

        dx[k] = xk.x - xc.x;
        dy[k] = xk.y - xc.y;
      }

      double wk = 1. / sqrt(dx[k] * dx[k] + dy[k] * dy[k]);
      cells[i].wk[k] = wk;

      axx += wk * dx[k] * dx[k];
      axy += wk * dx[k] * dy[k];
      ayy += wk * dy[k] * dy[k];
    }

    double det = axx * ayy - axy * axy;
    double Mxx, Mxy, Myy;

    Mxx = ayy / det;
    Mxy = -axy / det;
    Myy = axx / det;

    for (int k = 0; k < nFaces; k++) {
      cells[i].ck[k].cx[0] = Mxx * dx[k] + Mxy * dy[k];
      cells[i].ck[k].cx[1] = Mxy * dx[k] + Myy * dy[k];
    }
  }
}

void Mesh::set_zones() {
  nZone = 4;
  zones = new zone_t[nZone];
}

void Mesh::cell_funcs(Cell *(&cells)) {
  for (int i = 0; i < nCells; i++) {
    // Cell c(cells[i]);
    int n = cells[i].get_nNodes(); // 4
    Polygon pl(n);
    for (int k = 0; k < n; k++) {
      pl.set_p(k, nodes[cells[i].get_Node(k)]);
    }
    cells[i].set_c(pl.Center());
    cells[i].set_S(pl.Square());
    cells[i].set_nFaces(n);

    int f_index;

    f_index = min(cells[i].get_Node(0), cells[i].get_Node(1)); // 0
    f_index = Nx * (Ny - 1) + f_index / Ny + Ny * (f_index % Ny);
    cells[i].set_Face(0, f_index);
    cells[i].set_fType(0, faces[f_index].is_bound);
    cells[i].set_cells(0, faces[f_index].cl == i ? faces[f_index].cr
                                                 : faces[f_index].cl);

    f_index = min(cells[i].get_Node(1), cells[i].get_Node(2)); // 1
    f_index -= f_index / Ny;
    cells[i].set_Face(1, f_index);
    cells[i].set_fType(1, faces[f_index].is_bound);
    cells[i].set_cells(1, faces[f_index].cl == i ? faces[f_index].cr
                                                 : faces[f_index].cl);

    f_index = min(cells[i].get_Node(2), cells[i].get_Node(3)); // 3
    f_index = Nx * (Ny - 1) + f_index / Ny + Ny * (f_index % Ny);
    cells[i].set_Face(2, f_index);
    cells[i].set_fType(2, faces[f_index].is_bound);
    cells[i].set_cells(2, faces[f_index].cl == i ? faces[f_index].cr
                                                 : faces[f_index].cl);

    f_index = min(cells[i].get_Node(0), cells[i].get_Node(3)); // 0
    f_index -= f_index / Ny;
    cells[i].set_Face(3, f_index);
    cells[i].set_fType(3, faces[f_index].is_bound);
    cells[i].set_cells(3, faces[f_index].cl == i ? faces[f_index].cr
                                                 : faces[f_index].cl);
  }
}
