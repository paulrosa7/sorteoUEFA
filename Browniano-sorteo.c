#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define NormRANu (2.3283063671E-10F)
#define PUNTOSMAX  250000
#define BLOQUESMAX 500
#define MEDIDASMAX 100000

#define NormRANu (2.3283063671E-10F)
#define PUNTOSMAX  250000
#define BLOQUESMAX 500
#define MEDIDASMAX 100000
#define pi 3.141592
#define m   1
#define k   1
#define const_estocastica  2*mu*T

double parisi_rapuano ()
{
    unsigned int irr [256];
    unsigned char ind_ran=0, ig1=0, ig2=0, ig3=0;
    double numero_aleatorio;
    int i;

    for (i=0; i<256; i++)
        irr [i] = (rand () << 16) + rand ();

    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;
    irr [ind_ran] = irr [ig1] + irr [ig2];

    numero_aleatorio = (irr[ind_ran]^irr [ig3]);
    ind_ran++;

    return numero_aleatorio*NormRANu;
}

void box_muller (double *g1, double *g2, double varianza)
{
    double w1, w2;
    w1 = parisi_rapuano ();
    w2 = parisi_rapuano ();

    *g1 = -sqrt(-2*log(w1))*cos(2*pi*w2)*sqrt(varianza);
    *g2 = -sqrt(-2*log(w1))*sin(2*pi*w2)*sqrt(varianza);
}

void eulerMaruyama(double x[6],double y[6],double h,double mu, double T, double num_aleatorio){
    double random1, random2,random3, random4;

    box_muller (&random1, &random4, const_estocastica);
    box_muller (&random3, &random4, const_estocastica);
    box_muller (&random2, &random4, const_estocastica);

    y[0]=x[0]+h*x[1]/m;
    y[2]=x[2]+h*x[3]/m;
    y[4]=x[4]+h*x[5]/m;

    y[1]=x[1]+(-mu*x[1]/m)*h+sqrt(2*mu*T*h)*random1;
    y[3]=x[3]+(-mu*x[3]/m-m*9.8)*h+sqrt(2*mu*T*h)*random3;
    y[5]=x[5]+(-mu*x[5]/m)*h+sqrt(2*mu*T*h)*random2;

}
void min_max (int npuntos, double *puntos, double *minimo, double *maximo)
{
    int i;

    *minimo = *maximo = puntos [0];
    for (i=1; i<npuntos; i++)
    {
        if (puntos[i]<*minimo) *minimo=puntos[i];
        if (puntos[i]>*maximo) *maximo=puntos[i];
    }
}
void construye_histograma (int nbloques, int npuntos, double *puntos, double *hist, double *x)
{
    double minimo, maximo, delta;
    int i, aux;

    min_max (npuntos, puntos, &minimo, &maximo);
    delta = (maximo-minimo)/((double)(nbloques));
    printf ("minimo=%lf\tmaximo=%lf\tdelta=%lf\n", minimo, maximo, delta);

    for (i=0; i<npuntos; i++)
    {
        aux = (int)((puntos[i]-minimo)/delta);
        if (aux==nbloques) aux--;
        hist[aux] ++;
    }

    for (i=0; i<nbloques; i++)
    {
        hist[i]/=npuntos;
        x[i]=minimo+i*delta;
    }


}


void exporta_histograma (int nbloques, double *x, double *hist, char *nombre)
{
    int i;
    FILE *f = fopen (nombre, "w");

    for (i=0; i<nbloques; i++)
    {
        fprintf (f, "%lf\t%lf\n", x[i], hist[i]);
    }
    fclose (f);
}


int main(){

double t_total=50;
double h=0.01;
double mu=1.5;
double T=0.1;
double hist [BLOQUESMAX], hist_ejex[BLOQUESMAX];
double xaux, x[MEDIDASMAX], histx [BLOQUESMAX], x_ejex[BLOQUESMAX];
FILE *f,*p,*d;
double nbloques=10;

int num_pasos=(int)(t_total/h);
double altura=100;

f=fopen("estados.txt","w");
p=fopen("puntos.txt","a");
d=fopen("distancias.txt","w");

double estados[num_pasos+1][6];

// Condiciones iniciales
estados[0][0]=0;
estados[0][1]=0;
estados[0][2]=altura;
estados[0][3]=0;
estados[0][4]=0;
estados[0][5]=0;

double itmax=10000;
double x_aux,y_aux;

int seguir=1;
for(int j=0;j<=itmax;j++){
    seguir=1;
    for(int i=1;i<=num_pasos && seguir;i++){
        eulerMaruyama(estados[i-1],estados[i],h,mu,T,0.1);
        if (estados[i][2]<=0){
            x_aux=estados[i][0];
            y_aux=estados[i][4];
            fprintf(p,"%f %f\n",x_aux,y_aux);
            //dd[i]=sqrt(x_aux*x_aux+y_aux*y_aux);
            fprintf(d,"%f\n",sqrt(x_aux*x_aux+y_aux*y_aux));
            seguir=0;
        }
        fprintf(f,"%f %f %f \n",estados[i][0],estados[i][4],estados[i][2]);
    }
}



fclose(p);
fclose(d);
fclose(f);
//system("distancias.plt");
//system("distancias2.plot");


}
