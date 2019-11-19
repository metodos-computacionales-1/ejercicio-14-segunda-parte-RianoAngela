#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double derivada1(double t, double y, double x); 
double derivada2(double t, double v, double x);
void euler(double dt, string nombre);
void runge_kutta(double dt, string nombre);

const double K = 100;
const double M = 2;
const double LAMBDA = 1;

int main ()
{
    euler(0.01, "euler.dat");
    runge_kutta(0.01, "rk.dat");
    cout<<"Se espera que la solucion sea en terminos de seno y coseno, puesto que es un oscilador armonico ";
    return 0;
}

double derivada1(double t, double v, double x)
{
  return v;
}

double derivada2(double t, double v, double x)
{
  return (-K*x);
}

void euler(double DeltaT, string nombre)
{
    int pasos = 1000;
    double x[pasos],vx[pasos];

    ofstream outfile;
    outfile.open(nombre);
    
    //condiciones iniciales
    x[0] = 1;
    vx[0] = 0;
    
    outfile<< x[0]<<" "<<vx[0]<<" "<<endl;
    
    double t=0.0;
    for (int i=1;i<=20;i++)
    {
        x[i] = (vx[i-1]*DeltaT)+x[i-1];
        vx[i] = (x[i-1]*DeltaT)+vx[i-1];   
        outfile<< x[i]<<" "<<vx[i]<<" "<<endl;
        t=t+DeltaT;
    }
    outfile.close();
    
}

void runge_kutta(double dt, string nombre)
{
    double t, t_futuro;
    double x_presente, vx_presente;
    double x_futuro, vx_futuro;
    double k0_x, k0_vx;
    double k1_x, k1_vx;
    double k2_x, k2_vx;
    double k3_x, k3_vx;
    double k_x, k_vx;
    double dvx, dvx1;
  
    ofstream outfile;
    outfile.open(nombre);
    //condiciones iniciales
    x_presente = 1;
    vx_presente = 0;

    t=0.0;
    for (int i=0; i<=20;i++)
    {
        dvx = (x_presente*dt)+vx_presente;
        outfile<<x_futuro<<" "<<vx_futuro<<" "<<endl;
        
        k0_x = vx_presente;
        k0_vx = dvx;
        
        x_futuro = x_presente + ((0.5 * dt) * k0_x);
        vx_futuro = vx_presente + ((0.5 * dt) * k0_vx);

        k1_x = vx_futuro;
        k1_vx = dvx;
        
        x_futuro = x_presente + ((0.5 * dt) * k1_x);
        vx_futuro = vx_presente + ((0.5 * dt) * k1_vx);
        k2_x = vx_futuro;
        k2_vx = dvx;
        
       
        x_futuro = x_presente + (dt * k2_x);
        vx_futuro = vx_presente + (dt * k2_vx);
        k3_x = vx_futuro;
        k3_vx = dvx;
        
        
        k_x  = (k0_x/6.0) + (k1_x/3.0) + (k2_x/3.0) + (k3_x/6.0);
        k_vx = (k0_vx/6.0) + (k1_vx/3.0) + (k2_vx/3.0) + (k3_vx/6.0);
        
        x_futuro = x_presente + (dt * k_x);
        vx_futuro = vx_presente + (dt * k_vx);
        t = t + dt;
        
        x_presente = x_futuro;
        vx_presente = vx_futuro;
     
    }
    outfile.close();
}