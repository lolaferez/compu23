// PROGRAMA PARA RESOLVER LA ECUACION DE ONDA DE SCHRODINGER QUE ES UNA EC. DIF. EN DERIVADAS PARCIALES
//_____________________________________________________________________________________________________

#include <iostream>
#include <cmath>
#include <fstream>
#include "complex.h" //Biblioteca externa para trabajar con numeros complejos

using namespace std;

#define N 100 // Numero de puntos en los que dividimos el espacio
#define nciclo (N / (4.0))
#define cte 11.0778 //Norma (HACER PREVIA)

// hacer tres graficas con diferentes alturas de potencial

int main()
{

    // Declaro las variables
    int j, t;
    double lamda, s, k0, V[N + 1];
    double norma;

    fcomplex phi[N + 1], alpha[N], A0[N], beta[N], b[N], X[N + 1];

    ofstream fich, fichnorma; // Para escribir los datos

    fich.open("fonda.txt");
    fichnorma.open("norma.txt");

    // Inicializo los valores de los parametros
    lamda = 0.1;
    norma = 0.0;

    // Calculamos los valores de s y de k0
    k0 = 2.0 * (M_PI) * nciclo / (1.0 * N);
    s = 1.0 / (4.0 * k0 * k0);

    // Inicializamos la barrera de potencial
    for (j = 0; j <= N; j++)
    {
        if (j >= (2.0 / 5.0 * N) && (j <= (3.0 / 5.0 * N)))
            V[j] = lamda * k0 * k0;
        else
            V[j] = 0.0;
    }

    // Ahora aplicamos las conidiciones de contorno a la funcion de onda
    phi[0] = Complex(0.0, 0.0);
    phi[N] = Complex(0.0, 0.0);
    // Y el valor de la funcion de onda en el resto de puntos sigue la siguiente expresion en el primer instante
    for (j = 1; j < N; j++)
    {
        phi[j] = (RCmul(exp((-(8.0) * (4 * j - N) * (4 * j - N) / (N * N))), (Complex(cos(k0 * j), sin(k0 * j)))));
        //Normalizo una vez conocida la norma
        phi[j]= RCmul(1.0 / sqrt(cte) , phi[j]);
        norma = norma + pow(Cabs(phi[j]), 2);
    }
    
    //Guardo los valores de la funcion de onda
    for (j = 0; j <= N; j++)
    {
        fich << (1.0) * j / (100.0) << ", " << Cabs(phi[j]) << ", " << V[j] << endl; // Nodo, valor de funcion de onda en el nodo y el pozo de potencial
    }
    // Guardo la primera norma
    fichnorma << norma << endl;
    fich << endl;

    // Obtenemos los valores para alpha
    // Primeramente los coeficientes A0

    for (j = 0; j < N; j++)
    {
        A0[j] = Complex((-2.0) - V[j], (2.0) / s);
    }

    // Dado que A- y A+  son 1 en cualquier punto j
    alpha[N - 1] = Complex(0.0, 0.0);
    for (j = N - 2; j >= 0; j--)
    {
        alpha[j] = Cdiv(Complex(-1.0, 0.0), (Cadd((A0[j + 1]), (alpha[j + 1]))));
    }

    for (t = 0; t <= 5000; t++) // tiempo discretizado. Fotogramas que mostrara la animacion 10000
    {

        norma=0.0;

        // Calculamos los terminos b
        for (j = 0; j < N; j++)
        {
            b[j] = RCmul(4.0 / s, Cmul(phi[j], Complex(0.0, 1.0)));
        }

        // Calculamos ahora los coeficientes beta
        beta[N - 1] = Complex(0.0, 0.0);
        for (j = N - 2; j >= 0; j--)
        {
            beta[j] = Cdiv((Csub(b[j + 1], beta[j + 1])), (Cadd((A0[j + 1]), (alpha[j + 1]))));
        }

        // Ahora las X
        X[0] = Complex(0.0, 0.0);
        X[N] = Complex(0.0, 0.0);
        for (j = 1; j < N; j++)
        {
            X[j] = Cadd(Cmul(alpha[j - 1], X[j - 1]), beta[j - 1]);
        }

        // Ahora podemos calcular la funcion de onda en el tiempo
        for (j = 0; j < N; j++)
        {
            phi[j] = Csub(X[j], phi[j]);
            norma = norma + pow(Cabs(phi[j]), 2);
        }
        
        // Los guardo en los ficheros
        for (j = 0; j <= N; j++)
        {
            fich << (1.0) * j / (100.0) << ", " << Cabs(phi[j]) << ", " << V[j] << endl; // Nodo, valor de funcion de onda en el nodo y el pozo de potencial
        }
        
        fich << endl;
        fichnorma << norma << endl;
    }

    fich.close();
    fichnorma.close();
}
