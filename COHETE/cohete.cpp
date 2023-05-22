// PROGRAMA PARA RESOLVER EL PROBLEMA DE LOS TRES CUERPOS USANDO EL METODO DE RUNGEKUTTA4

// Algunas aproximaciones fisicas: suponemos los tres objetos giran en el mismo plano, suponemos que la luna
// gira a velocidad angular constantes. Asimismo, la tierra sera el origen del sistema de referencia.

#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

#define ML 0.07349e24     // Masa Luna
#define MT 5.9736e24      // Masa Tierra
#define v 11157.3         // Velocidad de lanzamiento del cohete. Aproximada a la velocidad de escape
#define G 6.67384e-11     // Constante de gravitacion universal
#define dtl 3.844e8       // Distancia en metros entre la Tierra y la Luna
#define w 2.6617e-6       // Velocidad angular de la Luna respecto a la Tierra
#define rt 6.37816e6      // Radio de la Tierra
#define rl 1.7374e6       // Radio de la Luna
#define h 60.0            // Paso en segundos
#define teta (M_PI) / 5.0 // Angulo de lanzamiento del cohete

// Algunas constantes necesarias para las funciones reescaladas

#define tri 7.018782618e-12 // G*MT/(dtl)^3
#define mu 0.01230246       // Cociente masa luna entre masa tierra

// Funciones para reescalar las magnitudes (inversa)
double resr(double r)
{
    double rr;
    rr = r * dtl;
    return rr;
}

double respphi(double pphi, double m) // m denota la masa del cohete
{
    double resphi;
    resphi = pphi * m * dtl * dtl;
    return resphi;
}

double respr(double pr, double m)
{
    double respr;
    respr = pr * m * dtl;
    return respr;
}

// Funciones para el calculo de los elementos de r', phi' y los momentos generalizados
// con las expresiones reescaladas
double rpun(double pr)
{
    double respr;
    respr = pr;
    return respr;
}

double phipun(double pphi, double r)
{
    double respr;
    respr = pphi / (r * r);
    return respr;
}

double prpun(double pphi, double phi, double r, double t)
{
    double s;
    s = pphi * pphi / pow(r, 3) - tri * (1.0 / pow(r, (2)) + mu / pow(sqrt(1 + r * r - 2 * r * cos(phi - w * t)), (3)) * (r - cos(phi - w * t)));
    return s;
}

double pphipun(double r, double phi, double t)
{
    double s;
    s = -tri * mu * r * sin(phi - w * t) / pow(sqrt(1 + r * r - 2 * r * cos(phi - w * t)), 3);
    return s;
}

int main()
{
    // Declaramos las variables del problema

    double r, pr, phiprima, pphi, phi, rluna[2], rtierra[2], rcohete[2];
    double auxr, auxpr, auxpphi; // Magnitudes auxiliares para reescalar
    double k[4][4];              // Valores de las constantes k_i para cada magnitud en cada instante
    double t, m, H;

    ofstream ficheropos, ficherokas, fichh, fichvec; // Ficheros sobre los que escribir los datos de las posiciones, k_i y H'

    int i, j;

    // Inicio las variables

    phi = (M_PI) / 4.0;
    pr = v / dtl * cos(teta - phi);
    pphi = rt * v / (dtl * dtl) * sin(teta - phi);
    r = rt / dtl;
    auxr = 0.0;
    auxpr = 0.0;
    auxpphi = 0.0;
    m = 777.0;

    // CONDICIONES INICIALES

    rtierra[0] = 0.0;
    rtierra[1] = 0.0; // La Tierra se queda quieta en el origen
    rluna[0] = 1.0;   // Valor ya reescalado
    rluna[1] = 0.0;   // La luna situada inicialmente en el eje OX
    rcohete[0] = rt / dtl * cos(phi);
    rcohete[1] = rt / dtl * sin(phi);

    // Meto los datos iniciales en el fichero
    ficheropos.open("posiciones.txt");
    ficheropos << rtierra[0] << ", " << rtierra[1] << endl;
    ficheropos << rluna[0] << ", " << rluna[1] << endl;
    ficheropos << rcohete[0] << ", " << rcohete[1] << endl;
    ficheropos << endl;

    ficherokas.open("kas.txt");
    fichh.open("hamiltoniano.txt");
    fichvec.open("tp.txt");

    // ALGORITMO DE RUNGEKUTTA 4

    for (t = 0.0; t <= 604800.0; t = t + h) // Tiempo hasta una semana en segundos
    {
        // Calculo las k_1 para cada variable evaluando los valores iniciales
        //______________________________________________________________________________
        //______________________________________________________________________________
        k[0][0] = h * (rpun(pr));
        k[0][1] = h * (phipun(pphi, r));
        k[0][2] = h * (prpun(pphi, phi, r, t));
        k[0][3] = h * (pphipun(r, phi, t));

        // Calculo los siguientes k_2 de las variables
        //______________________________________________________________________________
        //______________________________________________________________________________
        k[1][0] = h * (rpun(pr + (k[0][2] / 2.0)));
        k[1][1] = h * (phipun(pphi + (k[0][3] / 2.0), r + (k[0][0] / 2.0)));
        k[1][2] = h * (prpun(pphi + (k[0][3] / 2.0), phi + (k[0][1] / 2.0), r + (k[0][0] / 2.0), t + (h / 2.0)));
        k[1][3] = h * (pphipun(r + (k[0][0] / 2.0), phi + (k[0][1] / 2.0), t + (h / 2.0)));

        // Calculo los siguientes k_3 de las variables
        //______________________________________________________________________________
        //______________________________________________________________________________
        k[2][0] = h * (rpun(pr + (k[1][2] / 2.0)));
        k[2][1] = h * (phipun(pphi + (k[1][3] / 2.0), r + (k[1][0] / 2.0)));
        k[2][2] = h * (prpun(pphi + (k[1][3] / 2.0), phi + (k[1][1] / 2.0), r + (k[1][0] / 2.0), t + (h / 2.0)));
        k[2][3] = h * (pphipun(r + (k[1][0] / 2.0), phi + (k[1][1] / 2.0), t + (h / 2.0)));

        // Calculo los siguientes k_4 de las variables
        //______________________________________________________________________________
        //______________________________________________________________________________
        k[3][0] = h * (rpun(pr + (k[2][2])));
        k[3][1] = h * (phipun(pphi + (k[2][3]), r + (k[2][0])));
        k[3][2] = h * (prpun(pphi + (k[2][3]), phi + (k[2][1]), r + (k[2][0]), t + (h)));
        k[3][3] = h * (pphipun(r + (k[2][0]), phi + (k[2][1]), t + (h)));

        // Una vez tenemos los coeficientes estamos en condiciones de obtener el valor del vector
        // en el siguiente instante
        phi = phi + (k[0][1] + 2.0 * k[1][1] + 2.0 * k[2][1] + k[3][1]) / (6.0);

        pphi = pphi + (k[0][3] + 2.0 * k[1][3] + 2.0 * k[2][3] + k[3][3]) / (6.0);

        pr = pr + (k[0][2] + 2.0 * k[1][2] + 2.0 * k[2][2] + k[3][2]) / (6.0);

        r = r + (k[0][0] + 2.0 * k[1][0] + 2.0 * k[2][0] + k[3][0]) / (6.0); // Modulo del vector del cohete,
        rcohete[0] = r * cos(phi);                                           // a continuacion el vector posicion del cohete
        rcohete[1] = r * sin(phi);

        // Movimiento de la Luna
        rluna[0] = cos(w * t);
        rluna[1] = sin(w * t);

        // Escribo en el archivo los valores de las posiciones para la animacion
        ficheropos << rtierra[0] << ", " << rtierra[1] << endl;
        ficheropos << rluna[0] << ", " << rluna[1] << endl;
        ficheropos << rcohete[0] << ", " << rcohete[1] << endl;
        ficheropos << endl;

        // Meto en un archivo los valores de k_i
        for (j = 0; j < 4; j++)
        {
            ficherokas << "k" << j + 1 << "-  ";
            for (i = 0; i < 4; i++)
            {
                ficherokas << k[j][i] << "     ";
            }
            ficherokas << endl;
        }
        ficherokas << endl;

        // Calculamos el halmiltoniano
        // Reescalamos las magnitudes utilizadas en las variables axiliares
        auxr = resr(r);
        auxpr = respr(pr, m);
        auxpphi = respphi(pphi, m);

        H = (auxpr * auxpr) / (2.0) / m + (auxpphi * auxpphi) / (2.0) / m / auxr / auxr - G * m * MT / auxr - G * ML * m / sqrt(dtl * dtl + auxr * auxr - (2.0) * dtl * auxr * cos(phi - w * t));

        // Escribimos en el fichero los valores de H'
        fichh << H - w * auxpphi << endl;

        fichvec << r << "     " << phi << "     " << pr << "     " << pphi;

        fichvec << endl;
    }
    ficheropos.close();
    ficherokas.close();
    fichh.close();
    fichvec.close();
}