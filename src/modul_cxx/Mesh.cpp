#include <math.h> 
#include <algorithm>

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

      faces[k].show(k);
      ++k;
    }
  }
  //  горизонтальные грани
  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx - 1; i++)  {
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



void  Mesh::cell_funcs(Cell *(&cells)) {
  for (int i = 0; i < nCells; i++) {
    // Cell c(cells[i]);
    int n = cells[i].get_nNodes(); // 4
    Polygon pl(n);
    for (int k = 0; k < n; k++) {
      pl.set_p(k, nodes[cells[i].get_Node(k)]);
    }
    cells[i].set_c(pl.Center());
    cells[i].set_S(pl.Square());
    cells[i].set_nFaces(4);

    vector<int> v = {cells[i].get_Node(0), cells[i].get_Node(1), cells[i].get_Node(2), cells[i].get_Node(3)};
    sort(v.begin(), v.end());
    cout << *v.begin() << " " << *--v.end() << endl;
    cells[i].set_Face(0, 0);
    cells[i].set_fType(0, 1);
    cells[i].set_Face(1, 1);
    cells[i].set_fType(1, 1);
    cells[i].set_Face(2, 2);
    cells[i].set_fType(2, 1);
    cells[i].set_Face(3, 3);
    cells[i].set_fType(3, 1);

    // for (int k = 0; k < nFaces; k++) {
      
    //   switch (faces[k].nodes[0])
    //   {
    //   case cells[i].get_Face(0):
    //     /* code */
    //     break;
      
    //   default:
    //     break;
    //   }
    //   if (faces[k].nodes[0] < faces[k].nodes[1] && cells[i].get_Face(0))
    // }
  }
}


void show(face_t f) {
  string fname = "../data/faces.txt";
  ofstream record(fname, ios::out | ios::app);

  if (record) {
    record << "\t\tfaces nodes[0] = " << f.nodes[0] << endl;
    record << "\t\tfaces nodes[1] = " << f.nodes[1] << endl;
    record << "cr = " << f.cr << ", cl = " << f.cl << endl;
    record << "is_bound = " << f.is_bound << endl;
    record << "length = " << f.length << endl;
    record << "f_centr: x = " << f.f_centr.x << ", y = " << f.f_centr.y << endl;
    
    record << "****************************************************" << endl;

  } else {
    cout << "Problem with file: " << fname << endl;
  }
  record.close();
}
