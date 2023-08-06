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
    // k = 9;
    // cout << nodes[k].x << nodes[k].y << endl;

  } else {
    cout << "Error: file \"" << filename << "\" can\'t open" << endl;
  }
  reading.close();
}

void Mesh::create_cells(Cell *(&cells)) {
  // cout << Nx << endl;
  // cout << Ny << endl;
  for (int i = 0; i < Nx - 1; i++) {
    for (int j = 0; j < Ny - 1; j++) {
      int nc = (Ny - 1) * i + j;

      cout << "nc = " << nc << endl;

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
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny - 1; j++) {
      int n1 = Ny * i + j;
      int n2 = Ny * i + j + 1;
      faces[k].nodes[0] = n1;
      faces[k].nodes[0] = n2;
      
      if (i == 0) {
        faces[k].is_bound = true;
        faces[k].cr = -1;
        faces[k].cr = (Ny - 1) * i + j;
      } else if (i == Nx - 1) {
        faces[k].is_bound = true;
        faces[k].cr = (Ny - 1) * (i - 1) + j;
        faces[k].cr = -1;
      } else {
        faces[k].is_bound = false;
        faces[k].cr = -1;
        faces[k].cr = (Ny - 1) * i + j;
      }
      faces[k].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
      faces[k].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

      double dx = (nodes[n1].x - nodes[n2].x);
      double dy = (nodes[n1].y - nodes[n2].y);

      faces[k].length = sqrt(dx * dx + dy * dy);

      // faces[k].show();
      ++k;
    }
  }
}
