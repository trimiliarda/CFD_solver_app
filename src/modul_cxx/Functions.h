#pragma once
#include "Mesh.h"

#define und << "," <<

typedef struct Paremeters {
  double ro, p, h, H, u, v, w, T, E, e;

  double mu, la, Pr;

  double Cp, Gam, Gm;

  double *U, *U1;
  double *V; // h ro u v
} param_t;

typedef struct Changes {
  double *dU;
} changes;

typedef struct Gradient
{
    vector_t * g;
} gradient_t;


void Init(param_t *(&p), int nCells, int Nm);
void set_gran(Mesh mesh);

void Yw(Mesh mesh, Cell *cells, int nCells);
void Viscous(param_t *(&p), changes *du, Mesh mesh, Cell *cells, double dt);
void Convect(param_t *(&p), changes *du, Mesh mesh, Cell *cells, int it,
             double dt);

void Gradients(Cell * cells, Mesh mesh, gradient_t gr, param_t *(&p), int Nm);

void get_params(param_t *(&p), int nCells, int Nm);

void gen_mesh(int Nx, int Ny);
void plot(param_t *(&p), Cell *cells, int Nx, int Ny, int nCells);
//  расстояние от т.Е до грани AB
double dist(point_t A, point_t B, point_t E);
