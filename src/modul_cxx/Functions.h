#pragma once
#include "Mesh.h"

#define und << "," <<

typedef struct Paremeters
{
    double ro, p, h, H, u, v, w, T, E, e;

    double mu, la, Pr;
    
    double Cp, Gam, Gm;

    double *U, *U1;
    double *V; //h ro u v
}param_t;

typedef struct Changes
{
    double * dU;
} changes;

void gen_mesh(int Nx, int Ny);
void Init(param_t * (&p), int nCells, int Nm);
void Viscous(param_t* p, changes * du, Mesh mesh, Cell * cells, double dt);
void get_params(param_t * (&p), int nCells, int Nm);
void plot(param_t * (&p),  Cell * cells,int Nx,int Ny,int nCells);
