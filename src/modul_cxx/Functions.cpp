#include <math.h>

#include "Functions.h"

void Init(param_t *(&p), int nCells, int Nm) {
  double T0, Cp, la, P0, Gam, Gm, R, ro, mu;
  double U0;

  //  входный данные
  T0 = 293.;      //  K
  la = 2.5685e-2; // W/(m k)
  mu = 17.863e-6;

  U0 = 0.1;

  P0 = 101325.; // Pa
  Gam = 1.4;
  Gm = 28.97;

  // расчёт
  R = 8314.41 / Gm; //  J / (kg K)
  ro = P0 / (R * T0);
  Cp = Gam / (Gam - 1.) * R;

  cout << "Cp = " << Cp << endl;

  for (int i = 0; i < nCells; i++) {
    p[i].ro = ro;
    p[i].p = P0;
    p[i].u = U0;
    p[i].v = 0.;
    p[i].w = 0.;
    p[i].T = T0;

    p[i].Cp = Cp;
    p[i].la = la;
    p[i].mu = mu;
    p[i].Gam = Gam;
    p[i].Gm = Gm;

    p[i].Pr = p[i].mu * p[i].Cp / p[i].la;
    p[i].h = Cp * T0;

    double q2 = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

    p[i].H = p[i].h + q2;
    p[i].E = p[i].H - p[i].p / p[i].ro;

    p[i].e = p[i].h / p[i].Gam;

    p[i].U = new double[Nm];
    p[i].U1 = new double[Nm];
    p[i].V = new double[Nm];

    p[i].U1[0] = p[i].ro * p[i].E;
    p[i].V[0] = p[i].h; //энтальпия
  }
}

void Viscous(param_t *p, changes *du, Mesh mesh, Cell *cells, double dt) {
  int nFaces = mesh.nFaces;

  for (int i = 0; i < nFaces; i++) {
    Face face = mesh.faces[i];

    int cr = face.cr;
    int cl = face.cl;

    if (face.is_bound) {
      int c = max(cl, cr);

      //  координаты ц.т. ячейки
      point_t xc = cells[c].get_c();

      //  координаты ц.т. грани
      point_t fc = face.f_centr;

      //  расстояние до стенки
      double dx = xc.x - fc.x;
      double dy = xc.y - fc.y;
      double dl = sqrt(dx * dx + dy * dy);

      // длина грани и площадь ячейки с
      double len = face.length;
      double S = cells[c].get_S();

      int z = face.zone;




      double Tw;
      if (face.zone == BOUND_Bottom)
        Tw = 500;
      if (face.zone == BOUND_Up)
        Tw = 200;

      if (face.zone == BOUND_Up || face.zone == BOUND_Bottom) {
        double hw = p[c].Cp * Tw;

        // dh / dn
        double dh_dn = (p[c].h - hw) / dl;

        // mu / Pr
        double mu_Pr = (p[c].mu / p[c].Pr);

        //  Fv - поток через грань
        double Fv = mu_Pr * dh_dn;

        du[c].dU[0] += -Fv * len / S * dt;
      }

      if (face.zone == BOUND_Left) {
        //  Fv - поток через грань
        double Fv = 0.;
        // double Fv = 5.;

        du[c].dU[0] += -Fv * len / S * dt;
      }
      if (face.zone == BOUND_Right) {
        //  Fv - поток через грань
        double Fv = 0.;
        // double Fv = -5.;

        du[c].dU[0] += -Fv * len / S * dt;
      }

    } else {

      //  координаты ц.т. правой и левой ячеик
      point_t xr = cells[cr].get_c();
      point_t xl = cells[cl].get_c();

      //  расстояние между ячейками
      double dx = (xr.x - xl.x);
      double dy = (xr.y - xl.y);
      double dl = sqrt(dx * dx + dy * dy);

      // dh / dn
      double dh_dn = (p[cr].h - p[cl].h) / dl;

      //  среднее mu / Pr
      double mu_Pr = 0.5 * (p[cr].mu / p[cr].Pr + p[cl].mu / p[cl].Pr);

      double Fv = mu_Pr * dh_dn; // поток через грань

      //  приращение
      //  если p[cr].h > p[cl].h ==>> Fv > 0
      //  тепло втекает в левую ячейку

      double len = face.length;
      double Sr = cells[cr].get_S();
      double Sl = cells[cl].get_S();

      du[cr].dU[0] += -Fv * len / Sr * dt;
      du[cl].dU[0] += Fv * len / Sl * dt;
    }
  }
}

void Convect(param_t *p, changes *du, Mesh mesh, Cell *cells, int it,
             double dt) {
  int nFaces = mesh.nFaces;

  for (int i = 0; i < nFaces; i++) {
    Face face = mesh.faces[i];

    int cr = face.cr;
    int cl = face.cl;

    int n1 = face.nodes[0];
    int n2 = face.nodes[1];

    double dl = face.length;

    point_t x1 = mesh.nodes[n1];
    point_t x2 = mesh.nodes[n2];

    double nx = -(x2.y - x1.y) / dl;
    double ny = (x2.x - x1.x) / dl;

    //  Fc - поток через грань
    double Fc;

    if (face.is_bound) {
      int c = max(cl, cr);
      double S = cells[c].get_S();

      if (face.zone == BOUND_Left) {
        double U_inlet, V_inlet, T_inlet;
        U_inlet = p[c].u;
        V_inlet = p[c].v;
        T_inlet = 800.;
        // T_inlet = 0.;

        double un = U_inlet * nx + V_inlet * ny;

        double h = p[c].Cp * T_inlet;
        double H = h + 0.5 * (U_inlet * U_inlet + V_inlet * V_inlet);

        double E = h / p[c].Gam + 0.5 * (U_inlet * U_inlet + V_inlet * V_inlet);

        //  Fc - поток через грань
        Fc = p[c].ro * H * un;
        du[c].dU[0] += -Fc * dl / S * dt;
      }
      if (face.zone == BOUND_Right) {
        double un = p[c].u * nx + p[c].v * ny;

        //  Fc - поток через грань
        Fc = p[c].ro * p[c].H * un;
        du[c].dU[0] += Fc * dl / S * dt;
      }
      if (face.zone == BOUND_Bottom || face.zone == BOUND_Up) {
        du[c].dU[0] += 0.;
      }

    } else {
      //  средние значения
      double u_, v_, un_, H_, E_, ro_;

      u_ = 0.5 * (p[cr].u + p[cl].u);
      v_ = 0.5 * (p[cr].v + p[cl].v);
      un_ = u_ * nx + v_ * ny;
      H_ = 0.5 * (p[cr].H + p[cl].H);
      E_ = 0.5 * (p[cr].E + p[cl].E);
      ro_ = 0.5 * (p[cr].ro + p[cl].ro);

      double A = H_ / E_ * un_;
      double Apl, Amn;
      Apl = 0.5 * (A + abs(A));
      Amn = 0.5 * (A - abs(A));

      double UL = p[cl].U[0];
      double UR = p[cr].U[0];

      Fc = Apl * UL + Amn * UR;

      double Sr = cells[cr].get_S();
      double Sl = cells[cl].get_S();

      du[cr].dU[0] += +Fc * dl / Sr * dt;
      du[cl].dU[0] += -Fc * dl / Sl * dt;
    }
  }
}

void Yw(Mesh mesh, Cell *cells, int nCells) {
  int nFace = mesh.nFaces;
  for (int k = 0; k < nCells; k++) {
    double z1 = 1.e10;
    point_t E = cells[k].get_c();

    for (int i = 0; i < nFace; i++) {
      if (mesh.faces[i].is_bound) {
        int n1 = mesh.faces[i].nodes[0];
        int n2 = mesh.faces[i].nodes[1];

        point_t A = mesh.nodes[n1];
        point_t B = mesh.nodes[n2];

        double z2 = dist(A, B, E);

        if (z1 > z2)
          z1 = z2;
      }
    }
    cells[k].Yw = z1;
  }

  //     int Nx = mesh.Nx;
  //     int Ny = mesh.Ny;

  //     string fname = "../output/Yw_dist.csv";
  //     // remove(fname);
  //     ofstream record(fname, ios::out);

  //     if (record) {
  //         record << "x,y,z,Yw(m)" << endl;

  //         for (int i = 0; i < nCells; i++) {
  //             record << cells[i].get_c().x und cells[i].get_c().y und 0. und
  //             cells[i].Yw << endl;
  //         }
  //     } else {
  //         cout << "Problem with file: " << fname << endl;
  //     }
  //     record.close();
}

double dist(point_t A, point_t B, point_t E) {
  double AB[2], BE[2], AE[2];

  double ab_be, ab_ae;
  double x, x1, x2, y, y1, y2, mod;
  double dist;

  //  vector AB
  AB[0] = B.x - A.x;
  AB[1] = B.y - A.y;

  //  vector BE
  BE[0] = E.x - B.x;
  BE[1] = E.y - B.y;

  //  vector AE
  AB[0] = E.x - A.x;
  AB[1] = E.y - A.y;

  ab_be = (AB[0] * BE[0] + AB[1] * BE[1]);
  ab_ae = (AB[0] * AE[0] + AB[1] * AE[1]);

  if (ab_be > 0) {
    y = E.y - B.y;
    x = E.x - B.x;
    dist = sqrt(x * x + y * y);
  } else if (ab_ae < 0) {
    y = E.y - A.y;
    x = E.x - A.x;
    dist = sqrt(x * x + y * y);
  } else {
    x1 = AB[0];
    y1 = AB[1];
    x2 = AE[0];
    y2 = AE[1];
    mod = sqrt(x1 * x1 + y1 * y1);
    dist = abs(x1 * y2 - y1 * x2) / mod;
  }

  return dist;
}

void get_params(param_t *(&p), int nCells, int Nm) {
  for (int i = 0; i < nCells; i++) {

    p[i].E = p[i].U1[0] / p[i].ro;

    double q2 = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

    double e = p[i].E - q2;

    p[i].h = e * p[i].Gam;
    p[i].T = p[i].h / p[i].Cp;

    p[i].H = p[i].h + q2;

    //  газовая постоянная (Дж/кМол) / молекулярную массу
    double R = 8314.41 / p[i].Gm;
    p[i].p = p[i].ro * R * p[i].T;

    p[i].V[0] = p[i].h;
  }
}


void set_gran(Mesh& mesh) {
    zone_t * z = mesh.zones;

    mesh.zones[0].grantype = GranType::WALL;
    mesh.zones[1].grantype = GranType::WALL;
    mesh.zones[2].grantype = GranType::WALL;
    mesh.zones[3].grantype = GranType::WALL;

    int nZone = mesh.nZone;
    for (int i = 0; i < nZone; i++) {
        int tp = mesh.zones[i].grantype;
        if (tp == GranType::WALL) {
            mesh.zones[i].wall = new wall_t[1];
        }
        // if (tp == GranType::SUPERSONIC_INLET) {
        //     // mesh.zones[i].wall = new wall_t[1];
        // }
    }

    //  left
        mesh.zones[0].wall[0].vel = Slip::NO_slip;
        mesh.zones[0].wall[0].temp = Temp::qw;
        mesh.zones[0].wall[0].value = 0.;
    //  bottom
        mesh.zones[1].wall[0].vel = Slip::NO_slip;
        mesh.zones[1].wall[0].temp = Temp::Tw;
        mesh.zones[1].wall[0].value = 500.;
    //  right
        mesh.zones[2].wall[0].vel = Slip::NO_slip;
        mesh.zones[2].wall[0].temp = Temp::qw;
        mesh.zones[2].wall[0].value = 0.;
    //  up
        mesh.zones[3].wall[0].vel = Slip::NO_slip;
        mesh.zones[3].wall[0].temp = Temp::Tw;
        mesh.zones[3].wall[0].value = 200.;
}


void Gradients(Cell * cells, Mesh mesh, gradient_t gr, param_t *(&p), int Nm) {
    
}




void plot(param_t *(&p), Cell *cells, int Nx, int Ny, int nCells) {
  string fname = "../output/out_data.csv";
  ofstream record(fname, ios::out);

  if (record) {
    record << "x,y,z,T(K)" << endl;

    for (int i = 0; i < nCells; i++) {
      record << cells[i].get_c().x und cells[i].get_c().y und 0. und p[i].T
             << endl;
    }
  } else {
    cout << "Problem with file: " << fname << endl;
  }
  record.close();
}

void gen_mesh(int Nx, int Ny) {
#define I << "\t" <<
  string fname = "../data/mesh.dat";
  ofstream record(fname, ios::out);

  if (record) {
    record << Nx I Ny I 0 I 0 << endl;

    double step = 1.;
    for (double i = 0.; i < Nx; i += step) {
      for (double j = 0.; j < Ny; j += step) {
        record << 0 I 0 I i I j << endl;
      }
    }
  } else {
    cout << "Problem with file: " << fname << endl;
  }
  record.close();
}
