// MODELO DE RED NEURONAL USANDO EL ALGORITMO DE METROPOLIS
// MODELO DE HOPFIELD

#include <iostream>
#include <cmath>
#include <fstream>
#include "gsl_rng.h" //Libreria para generación de números aleatorios

gsl_rng *tau;

using namespace std;

#define N 150     // Dimension de la matriz cuadrada CAMBIAR A N=20 PARA 4)
#define T 0.0001 // Temperatura del sistema
#define monte 50 // Numero de pasos montecarlo a realizar
#define pat 1    // Cantidad de patrones para aprender
#define c 1      // Para comenzar la matriz aleatoria c=1; patron deformado c=2
#define ka 1     // Poner de valor cero para realizar apartado 4

int main()
{

    extern gsl_rng *tau;   // Puntero al estado del número aleatorio
    int semilla = 84102925; // Semilla del generador de números aleatorios

    ifstream patron1, patron2, patron3, patron4, patron5;         // Ficheros para lectura de los patrones iniciales
    ofstream fichero, solapa, solapa1, solapa2, solapa3, solapa4; // Fichero para la escritura de datos iterativos

    int i, j, i2, j2, k, q, n, m, cuenta, c1; // Multiples contadores y el valor a caracteristico de cada patron
    double t, suma, suma2, H, aux, valor, p, eps;
    int M[N][N], P[N][N][pat];                     // Matriz de la red neuronal y del patron
    double w[N][N], a[pat], tetha[N][N], sol[pat]; // Para la interaccion entre neuronas, el a de cada patron y el solapamiento

    tau = gsl_rng_alloc(gsl_rng_taus); // Inicializamos el puntero
    gsl_rng_set(tau, semilla);         // Inicializamos la semilla

    fichero.open("matame.txt");
    patron1.open("patron2.txt");
    patron2.open("patron3.txt");
    patron3.open("patron4.txt");
    patron4.open("FLOR.txt");
    patron5.open("CALF.txt");
    solapa.open("Solvar1defor.txt");
    solapa1.open("Solvar2defor.txt");
    solapa2.open("Solvar3defor.txt");
    solapa3.open("Solvar4defor.txt");
    solapa4.open("Solvar5defor.txt");

    // Llenamos la matriz con los datos del primer patron deseado
    if (ka == 0)
    { // Apartado 4, necesitamos matrices aleatorias
        for (k = 0; k < pat; k++)
        {
            for (i = 0; i <= (N - 1); i++)
            {
                for (j = 0; j <= (N - 1); j++)
                {
                    P[i][j][k] = gsl_rng_uniform_int(tau, 2); // Numero random entero o 0 o 1
                }
            }
        }
    }
    else
    {
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron1 >> P[i][j][0];
            }
        }

        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron2 >> P[i][j][1];
            }
        }

        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron3 >> P[i][j][2];
            }
        }

        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron4 >> P[i][j][3];
            }
        }

        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                patron5 >> P[i][j][4];
            }
        }
    }

    if (c == 2) // Llenamos la matriz para el caso de patron con ruido
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                M[i][j] = P[i][j][3]; // NO OLVIDARSE DE SELECIONAR EL PATRON INICIAL CON RUIDO

                t = gsl_rng_uniform_int(tau, 5); // para que el patron se parezca en 4/5 partes

                if ((t == 0) && (M[i][j] == 0))
                {
                    M[i][j] = 1;
                }
                else
                {
                    if ((t == 0) && (M[i][j] == 1))
                    {
                        M[i][j] = 0;
                    }
                }
            }
        }
    }
    else // Llenamos la matriz con valor de actividad neuronal 1 y 0 aleatoriamente
    {
        if (c == 1)
        {
            for (i = 0; i <= (N - 1); i++)
            {
                for (j = 0; j <= (N - 1); j++)
                {
                    M[i][j] = gsl_rng_uniform_int(tau, 2); // Numero random entero o 0 o 1
                }
            }
        }
    }

    // Guardamos la matriz inicial en el archivo
    for (i = 0; i < N; i++)
    { // en las primeras posiciones del archivo van los datos y dejamos un espacio en blanco entre instante e instante
        for (j = 0; j < N; j++)
        {
            if (j < (N - 1))
            {
                fichero << M[i][j] << ", ";
            }
            else
            {
                fichero << M[i][j] << endl;
            }
        }
    }
    // Cambio de linea para nueva iteracion
    fichero << endl;

    // Calculo el valor de a del patron k-esimo
    for (k = 0; k < pat; k++)
    {
        a[k] = 0.0;
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                a[k] += P[i][j][k];
            }
        }
        a[k] = a[k] / (N * N);
        cout << "Coeficiente a del patron " << k+1 << " es " << a[k] << endl;
    }

    // Calculo del segundo termino de la expresion del incremento de energia. Una matriz de tethas_ij
    for (i = 0; i < N; i++) // Para cada neurona del sistema
    {
        for (j = 0; j < N; j++)
        {
            tetha[i][j] = 0.0;

            for (i2 = 0; i2 < N; i2++) // Para calcular el sumatorio de interacciones w de una
            {                          // determinada neurona [i][j] con las demas
                for (j2 = 0; j2 < N; j2++)
                {
                    if ((i2 != i) || (j2 != j))
                    {
                        aux = 0.0;
                        for (k = 0; k < pat; k++) // Para cada patron
                        {
                            aux += (1.0) / (N * N) * (P[i2][j2][k] - a[k]) * (P[i][j][k] - a[k]); 
                        }

                        tetha[i][j] = tetha[i][j] + aux; // Aux representa cada w_ij_i2j2
                    }
                }
            }

            tetha[i][j] = tetha[i][j] / (2.0);
        }
    }

    cuenta = 0;
    c1 = 0;

    for (q = 1; q <= (monte * N * N); q++)
    {
        // Elegimos un punto a azar de la red
        n = gsl_rng_uniform_int(tau, N);
        m = gsl_rng_uniform_int(tau, N);

        // Calculamos la energía en esa determinada configuracion M[n][m] paso a paso
        // Comenzamos con el calculo del primer termino de la expresion para el incremento de energia

        suma = 0.0;

        for (i = 0; i < N; i++) // Para cada neurona
        {
            for (j = 0; j < N; j++)
            {
                w[i][j] = 0.0;

                if ((n != i) || (m != j)) // Excluye el caso en que la neurona autointeracciona con si misma
                {
                    for (k = 0; k < pat; k++)
                    {
                        w[i][j] += (1.0) / (N * N) * (P[n][m][k] - a[k]) * (P[i][j][k] - a[k]);
                    }
                    suma = suma + w[i][j] * M[i][j];
                }
            }
        }

        // Completo el segundo termino de la expresion del delta de H
        suma2 = tetha[n][m];

        // Ahora si calculamos el  incremento del hamiltoniano del sistema para una determinada configuracion
        if (M[n][m] == 0)
        {
            H = suma2 - suma;
        }
        else
        {
            H = suma - suma2;
        }

        // Una vez tenemos la energia de la configuracion, seguimos con el algoritmo de metropolis
        // Recordemos que no se considera el estado para temperatura 0, ya que esta no es accesible

        valor = exp(-H / T);
        if ((1.0) <= valor)
        {
            p = 1.0;
        }
        else
        {
            p = valor;
        }

        eps = gsl_rng_uniform(tau); // Numero aleatorio real entre 0 y 1

        if ((eps < p) && (M[n][m] == 0))
        {
            M[n][m] = 1;
        }
        else
        {
            if ((eps < p) && (M[n][m] == 1))
            {
                M[n][m] = 0;
            }
        }

        // Ahora guardo la matriz cada paso montecarlo y el solapamiento
        if (cuenta == (N * N))
        {
            // Calculamos el solapamiento entre los patrones y la matriz
            for (k = 0; k < pat; k++)
            {
                sol[k] = 0.0;

                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        sol[k] += (P[i][j][k] - a[k]) * (M[i][j] - a[k]);
                    }
                }
                sol[k] = sol[k] / (N * N * a[k] * (1 - a[k]));
            }

            // Guardo el valor del solapamiento para la determinada configuracion en un archivo
            solapa << sol[0] << endl;
            solapa1 << sol[1] << endl;
            solapa2 << sol[2] << endl;
            solapa3 << sol[3] << endl;
            solapa4 << sol[4] << endl;

            for (i = 0; i <= (N); i++)
            { // en las primeras posiciones del archivo van los datos y dejamos un
              // espacio en blanco entre instante e instante

                if (i < (N))
                {
                    for (j = 0; j <= (N - 1); j++)
                    {
                        if (j < (N - 1))
                        {
                            fichero << M[i][j] << ", ";
                        }
                        else
                            fichero << M[i][j] << endl;
                    }
                }
                else
                {
                    fichero << endl;
                }
            }

            cuenta = 0;
        }

        cuenta = cuenta + 1;

        // Apartado 4 estudiamos ahora en funcion del solapamiento al final del bucle, si lo recuerda exitosamente o no
        if (q == monte * N * N)
        {
            for (k = 0; k < pat; k++) // Solo para la parte 4 ultimo solapamiento
            {
                sol[k] = 0.0;

                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        sol[k] += (P[i][j][k] - a[k]) * (M[i][j] - a[k]);
                    }
                }
                sol[k] = sol[k] / (N * N * a[k] * (1 - a[k]));

                cout << sol[k] << endl;
            }

            //con esto guardado podemos ver en el fichero que conforme cambiamos pat al comienzo los solapamientos 
            //obtenidos dado un punto no hay convergencia a uno con eso se obtiene la fraccion maxima

        }
    }

    //cout << "El numero de patrones bien aprendidos es de: " << c1 << endl;
    //cout << "La fraccion de patrones que la red puede recordar es: " << c1 / (N * N);

    fichero.close();
    patron1.close();
    solapa.close();
    solapa1.close();
    solapa2.close();
    solapa3.close();
    solapa4.close();
    patron2.close();
    patron3.close();
    patron4.close();
    patron5.close();
}