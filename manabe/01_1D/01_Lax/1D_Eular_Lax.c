#include <stdio.h>
#include <stdlib.h>
#define Nx 1000
#define Nt 500000

/*  CFD parameters  */
const double lambda = 2e-3;       //  m
const double dx = 2e-4;           //  m
const double dt = 1e-7;           //  s
/*  Experimental parameters */
const double U_p = 600;         //  m/s
const double S_ave = 8e+7;          //  W/m^2
const double eta = 1;
const double Q_in = S_ave*eta/lambda;   //  W/m^3
/*  Atomospheric parameters */
const double R = 287.0;             //  J/kg/K
const double Gamma = 1.4;
const double T0 = 300.0;             //  K
const double p0 = 1013*100.0;       //  J/m^3
const double rho0 = p0/R/T0;
const double E0 = rho0*(0.5*U_p*U_p + R*T0/(Gamma-1));     //  J/m^3
/*  Current parameters  */
double rho[Nx];
double u[Nx];
double E[Nx];
double p[Nx];
/*  parameters after dt */
double rho_next[Nx];
double u_next[Nx];
double E_next[Nx];
double p_next[Nx];
/*  Boundary Conditions */
void setBC(void){
    rho_next[0] = rho_next[1];
    u_next[0] = u_next[1];
    E_next[0] = E_next[1];
    p_next[0] = p_next[1];

    rho_next[Nx-1] = rho0;
    u_next[Nx-1] = -U_p;
    E_next[Nx-1] = E0;
    p_next[Nx-1] = p0;
}

void setIC(void){
    int i;
    for(i=0;i<Nx;i++){
        rho[i]=rho0;
        u[i]=-U_p;
        E[i]=E0;
        p[i]=p0;
    }
}

void setNEXT(void){
    int i;
    for(i=1;i<Nx-1;i++){
        rho_next[i] = (rho[i-1]+rho[i+1])*0.5 - dt*(rho[i+1]*u[i+1]-rho[i-1]*u[i-1])/2.0/dx;
        u_next[i] = ((rho[i-1]*u[i-1]+rho[i+1]*u[i+1])*0.5 - dt*(rho[i+1]*u[i+1]*u[i+1]+p[i+1]-rho[i-1]*u[i-1]*u[i-1]-p[i-1])/2.0/dx) / (rho_next[i]);
        if(((Nx) - (int)(0.01 + lambda/dx))/2<=i&&i<((Nx) + (int)(0.01 + lambda/dx))/2){
            E_next[i] = (E[i-1]+E[i+1]) * 0.5 - dt*((E[i+1]+p[i+1])*u[i+1] - (E[i-1]+p[i-1])*u[i-1]) / 2.0 / dx + Q_in * dt;
        }else{
            E_next[i] = (E[i-1]+E[i+1]) * 0.5 - dt*((E[i+1]+p[i+1])*u[i+1] - (E[i-1]+p[i-1])*u[i-1]) / 2.0 / dx;
        }
        p_next[i] = (Gamma-1.0)*(E_next[i]-(rho_next[i]*u_next[i]*u_next[i])*0.5);
    }
}

void copyNEXT(void){
    int i;
    for(i=0;i<Nx;i++){
        rho[i] = rho_next[i];
        u[i] = u_next[i];
        E[i] = E_next[i];
        p[i] = p_next[i];
    }
}

int main(){
    FILE *fp1,*fp2,*fp3,*fp4;
    fp1 = fopen("1D_Eular_Lax_p.csv","w");
    fp2 = fopen("1D_Eular_Lax_u.csv","w");
    fp3 = fopen("1D_Eular_Lax_T.csv","w");
    fp4 = fopen("1D_Eular_Lax_rho.csv","w");
    int i,j;
    setIC();
    for(i=0;i<=Nt;i++){
        setNEXT();
        setBC();
        copyNEXT();
        if((i % 100) == 0){
            for(j=0;j<Nx;j++){
                fprintf(fp1,"%f,",p[j]);
                fprintf(fp2,"%f,",u[j]);
                fprintf(fp3,"%f,",(p[j]/rho[j]/R));
                fprintf(fp4,"%f,",rho[j]);
            }
            fprintf(fp1,"\n");
            fprintf(fp2,"\n");
            fprintf(fp3,"\n");
            fprintf(fp4,"\n");
            
            printf("%d steps done.\n",i);
        }
    }
    return 0;
}

