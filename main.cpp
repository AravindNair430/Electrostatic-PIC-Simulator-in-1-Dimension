//
//  main.cpp
//  sample1DGrid3
//
//  Created by Aravind Nair on 27/06/18.
//  Copyright © 2018 Aravind Nair. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;

double Lmin;
double Lmax;  //Domain of 1d solution Lmin<= x <= Lmax
int N;    //No. of electrons
double J;    //No. of grid points
double dx;

struct particle {
    double q;
    double m;
    double v[3];
    double x[3];
    void setPart() {
        q = m = 1;
        cout << "Enter the velocity of the particle : "; cin >> v[0] >> v[1] >> v[2];
        cout << "Enter the displacement of the particle : "; cin >> x[0] >> x[1] >> x[2];
    }
}p[100];

struct grid {
    double E[3];
    double B[3];
    double Emag;
    double Bmag;
}G[1000];

void setdx() { dx = (Lmax - Lmin) / (J - 1); }

void setGrid() {
    cout << "Enter the no. of grid points : "; cin >> J;
    setdx();
    for(int i = 0 ; i < 1000 ; i++) {
        if(i < J) {
            G[i].Emag = 10;
            G[i].E[0] = G[i].Emag * sin(dx * i);
            G[i].E[1] = G[i].E[2] = 0;
            cout << "\n" << G[i].E[0];
        }
        else G[i].E[0] = G[i].E[1] = G[i].E[2] = 0;
        G[i].Bmag = 0;
        G[i].B[0] = G[i].B[1] = G[i].B[2] = 0;
        //cout << "\n" << i; //G[i].E[0];
    }
}

/*void em(grid g, particle& p,double dt) {
    double vn;
    //cout << "\n" << g.E[0];
    for(int j = 0 ; j < 3 ; j++) {
        vn = p.v[j];
        p.v[j] += (p.q / p.m) * g.E[j] * dt;
        p.x[j] += ((vn + p.v[j])/2)*dt;
    }
    if(p.x[0] < Lmax) cout << "\n" << p.x[0];
}*/

double magSquare(double b[3]) {
    double bSqr = (b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]);
    return bSqr;
}

double crossProduct( double v[], double B[], int i) {
    int j = 1,k = 2;
    switch(i) {
        case 0 : j = 1, k = 2;
            break;
        case 1 : j = 2, k = 0;
            break;
        case 2 : j = 0, k = 1;
            break;
    }
    return (v[k]*B[j] - v[j] * B[k]);
}

void em(grid g, particle& p,double dt) {
    int i; double T[3], vprime[3], s[3], vplus[3];
    for(i = 0 ; i < 3 ; i++) p.v[i] -= 0.5 * (p.q / p.m) * g.E[i] * dt;
    for(i = 0 ; i < 3 ; i++) T[i] = 0.5 * (p.q / p.m) * g.B[i] * dt;
    for(i = 0 ; i < 3 ; i++) s[i] = (2 * T[i])/(1 + magSquare(T));
    for(i = 0 ; i < 3 ; i++) vprime[i] = p.v[i] + crossProduct(p.v, T, i);
    for(i = 0 ; i < 3 ; i++) vplus[i] = p.v[i] + crossProduct(vprime, s, i);
    for(i = 0 ; i < 3 ; i++) {
        p.v[i] = vplus[i] + 0.5 * (p.q / p.m) * g.E[i] * dt;
        p.x[i] += p.v[i] * dt;
    }
    /*if(t >= tout) {
        //cout << "\n" << p.x[0];
        tout += 0.01;
    }*/
    //if(p.x[0] < Lmax) cout << "\n" << p.x[0];
}

void linearInterpolation(particle& p,double dt) {
    grid gr;
    int a = int(p.x[0] / dx), b = a + 1, point;
    if(a < J) {
        double dista = p.x[0] - a, distb = p.x[0] - b, magdista = sqrt(dista * dista), magdistb = sqrt(distb * distb);
        if(magdista > magdistb) point = b;
        else point = a;
        //cout << "\n" << a; //G[a].E[0];// << "\tG[b].E[0] : " << G[b].E[0];
        gr.E[0] = G[a].E[0] * ((p.x[0] - (double(a) * dx))/dx) + G[b].E[0] * (((double(b) * dx) - p.x[0] )/dx);
        //cout<< "\n" << gr.E[0];
        em(gr, p, dt);
        //if(p.x[0] < Lmax) cout << "\n" << p.x[0];
        //return G[point];
    }
}

/*void eulerMethod(int N, particle& p, double tstop, double dt) {
    double t = 0; //vn;
    //grid g; //g = linearInterpolation(p);
    for(double i = 0; i < tstop; i += dt) {
        linearInterpolation(p, dt);
        //g = linearInterpolation(p);
        //cout << "\n" << g.E[0];
        //for(int j = 0 ; j < 3 ; j++) {
        //    vn = p.v[j];
        //    p.v[j] += (p.q / p.m) * g.E[j] * dt;
        //    p.x[j] += ((vn + p.v[j])/2)*dt;
        //} t++;
        //if(t == tstop) break;
    }
    cout << "\n";
}*/

void borisMethod(int N, particle p, double tstop, double dt) {
    //grid g;
    //double tout = 0.01, omega1 = 1,omega2 = 1, phi1 = 0, phi2 = 0,T[3], s[3], vprime[3], vplus[3], x0[3]; int i;
    //for(i = 0 ; i < 3 ; i++) x0[i] = p.x[i];
    for(double t = 0 ; t < tstop; t += dt) {
        linearInterpolation(p, dt);/* cout << "\n" <<g.E[0];
        for(i = 0 ; i < 3 ; i++) p.v[i] -= 0.5 * (p.q / p.m) * g.E[i] * dt;
        for(i = 0 ; i < 3 ; i++) T[i] = 0.5 * (p.q / p.m) * g.B[i] * dt;
        for(i = 0 ; i < 3 ; i++) s[i] = (2 * T[i])/(1 + magSquare(T));
        for(i = 0 ; i < 3 ; i++) vprime[i] = p.v[i] + crossProduct(p.v, T, i);
        for(i = 0 ; i < 3 ; i++) vplus[i] = p.v[i] + crossProduct(vprime, s, i);
        for(i = 0 ; i < 3 ; i++) {
            p.v[i] = vplus[i] + 0.5 * (p.q / p.m) * g.E[i] * dt;
            p.x[i] += p.v[i] * dt;
        }
        if(t >= tout) {
            //cout << "\n" << p.x[0];
            tout += 0.01;
        }
        if(p.x[0] > Lmax) break;*/
    }
}

int main(int argc, const char * argv[]) {
    double t; int n, i; N = 1;
    //cout << "Enter the no. of particles : "; cin >> N;
    for(i = 0 ; i < N ; i++) p[i].setPart();
    
    cout << "Enter the number of time steps : "; cin >> n;
    cout << "Enter the time upto which integration is to be performed : "; cin >> t;
    double dt = t/(n - 1);
    cout << "Enter the upper limit for the particle : "; cin >> Lmax;
    cout << "Enter the lower limit for the particle : "; cin >> Lmin;
    setGrid();
    for(i = 0 ; i < N ; i++)  {
        if(p[i].x[0] < Lmax) {
        //eulerMethod(N, p[i], t, dt);
            borisMethod(N, p[i], t, dt);
            cout << "\n" << p[i].x[0];
        }
    }
    return 0;
}
