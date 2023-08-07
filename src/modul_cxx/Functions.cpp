#include <math.h>

#include "Functions.h"


void Init(param_t * (&p), int nCells, int Nm) {
    double T0, Cp, la, P0,  Gam, Gm, R, ro, mu;
    double U0;

//  входный данные
    T0 = 293.;      //  K
    la = 2.5685e-2; // W/(m k)
    mu = 17.863e-6;

    U0 = 0.1;

    P0 = 101325.;   // Pa
    Gam = 1.4;
    Gm = 28.97;

// расчёт
    R = 8314.41 / Gm;   //  J / (kg K)
    ro = P0 / (R * T0);
    Cp = Gam / (Gam - 1.) * R;

    cout << "Cp = " << Cp << endl;

    for (int i = 0; i < nCells ; i++) {
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
        p[i].V[0] = p[i].h; //энтальрия
         
    }

}


void Viscous(param_t* p, changes * du, Mesh mesh, Cell * cells, double dt) {
    int nFaces = mesh.nFaces;

    for (int i = 0; i < nFaces; i++ ) {
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

            double Tw;
            if (face.zone == BOUND_Bottom) Tw = 500;
            if (face.zone == BOUND_Up) Tw = 200;

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

            if (face.zone == BOUND_Left || face.zone == BOUND_Right) {
                //  Fv - поток через грань
                double Fv = 0.;
                
                du[c].dU[0] += -Fv * len / S * dt;
            }


        } else {
            
            //  координаты ц.т. правой и левой ячеик
            point_t xr = cells[cr].get_c();
            point_t xl = cells[cl].get_c();

            //  расстояние между ячейками
            double dx = xr.x - xl.x;
            double dy = xr.y - xl.y;
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
            du[cl].dU[0] +=  Fv * len / Sl * dt;
        }
    }
}



void get_params(param_t * (&p), int nCells, int Nm) 
{ 
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

void plot(param_t * (&p),  Cell * cells,int Nx,int Ny,int nCells) {
  string fname = "../output/out_data.csv";
  ofstream record(fname, ios::out);

  if (record) {
    record << "x,y,z,T(K)" << endl;

    for (int i = 0; i < nCells; i++) {
        record << cells[i].get_c().x und cells[i].get_c().y und 0. und p[i].T << endl;
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
    record << Nx I Ny I 0  I 0 << endl;

    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        record << 0 I 0 I i I j << endl;
      }
    }
  } else {
    cout << "Problem with file: " << fname << endl;
  }
  record.close();
}
