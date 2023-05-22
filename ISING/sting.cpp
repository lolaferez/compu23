//________________________________________________________________________________________________//
//________________________________________________________________________________________________//
//____________________________ CODIGO APLICANDO EL MODELO ISING __________________________________//
//________________________________________________________________________________________________//
//________________________________________________________________________________________________//

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;



int main()
{

    srand(time(NULL));

    int i = 0, j = 0, N, num, n, m, q, c, cuenta; // contadores y valor de la matriz N
    int M[100][100];                         // Matriz de espines de NxN dimension
    float T, eps, p, E, valor;
    ofstream fichero;                  // Para escribir en el los datos de la matriz de espines a cada instante


    N=100;

    cout << "Si desea comenzar con un sistema ordenado, pulse 0, si lo desea desordenado pulse 1: ";
    cin >> c;

    cout << "Introduzca la temperatura T del sistema, entre 0 y 5: ";
    cin >> T;

    cuenta=0;

    //__________________________________________________________________________________________________//
    //______________________________________ALGORITMO DE METROPOLIS_____________________________________//
    //__________________________________________________________________________________________________//

    if (c == 0)
    {
        // Llenamos la matriz con valor de spin 1
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
                M[i][j] = 1;
        }
    }
    else
    {
        // Llenamos la matriz con valor de spin 1 y -1 aleatoriamente
        for (i = 0; i <= (N - 1); i++)
        {
          
            for (j = 0; j <= (N - 1); j++)
            {
                num = rand() % 2; // Numero random o 0 o 1 y con eso asignamos el valor de spin 1 y -1 a la matriz
                
                if (num == 0)
                {
                    M[i][j] = -1;
                }
                else
                {
                    M[i][j] = 1;
                }

            }
        }
    }

    

    fichero.open("spines.txt");

    //Muestro la primera matriz y la meto en el fichero
    cout << "LA MATRIZ INICAL ES:" << endl;
    for (i = 0; i < N ; i++)
    { // en las primeras posiciones del archivo van los datos y dejamos un espacio en blanco entre instante e instante
        for (j = 0; j < N; j++)
        {
            if (j < (N - 1))
            {
                cout << M[i][j] << ", ";
                fichero << M[i][j] << ", ";
            }
            else
            {
                cout << M[i][j] << endl;
                fichero << M[i][j] << endl;
            }
        }
    }


    fichero << endl;


    for (q = 1; q <= 2000*N*N; q++) //2000 pasos montecarlo
    {

        // Ahora elegimos un punto al azar de la red
        i = rand() % N; // Genera un numero aleatorio entre 0 y N-1, justo las posiciones del array
        j = rand() % N;

        // Calculamos el incremento de energia de ese punto, implementando ya las condiciones de contorno periodicas
        // encadenando condiciones if else
        if ((i != N - 1) && (j != N - 1) && (i != 0) && (j != 0))
        {
            E = (2.0) * (M[i][j]) * ((M[i + 1][j]) + (M[i - 1][j]) + (M[i][j + 1]) + (M[i][j - 1]));
        }
        else
        {
            if ((i == N - 1) && (j != N - 1) && (j != 0))
            {
                E = (2.0) * (M[i][j]) * ((M[0][j]) + (M[i - 1][j]) + (M[i][j + 1]) + (M[i][j - 1]));
            }
            else
            {
                if ((i != N - 1) && (j == N - 1) && (i != 0))
                {
                    E = (2.0) * (M[i][j]) * ((M[i + 1][j]) + (M[i - 1][j]) + (M[i][0]) + (M[i][j - 1]));
                }
                else
                {
                    if ((j != N - 1) && (i == 0) && (j != 0))
                    {
                        E = (2.0) * (M[i][j]) * ((M[i + 1][j]) + (M[N - 1][j]) + (M[i][j + 1]) + (M[i][j - 1]));
                    }
                    else
                    {
                        if ((i != N - 1) && (i != 0) && (j == 0))
                        {
                            E = (2.0) * (M[i][j]) * ((M[i + 1][j]) + (M[i - 1][j]) + (M[i][j + 1]) + (M[i][N - 1]));
                        }
                        else
                        {
                            if ((i == 0) && (j == 0))
                            {
                                E = (2.0) * (M[i][j]) * ((M[i + 1][j]) + (M[N - 1][j]) + (M[i][j + 1]) + (M[i][N - 1]));
                            }
                            else
                            {
                                if ((i == N - 1) && (j == N - 1))
                                {
                                    E = (2.0) * (M[i][j]) * ((M[0][j]) + (M[i - 1][j]) + (M[i][0]) + (M[i][j - 1]));
                                }
                                else
                                {
                                    if ((i == 0) && (j == N - 1))
                                    {
                                        E = (2.0) * (M[i][j]) * ((M[i + 1][j]) + (M[N - 1][j]) + (M[i][0]) + (M[i][j - 1]));
                                    }
                                    else
                                    {
                                        E = (2.0) * (M[i][j]) * ((M[0][j]) + (M[i - 1][j]) + (M[i][j + 1]) + (M[i][N - 1]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Evaluamos el valor de p teniendo en cuenta el caso T=0 [Este caso no se considera en la 
        // practica porque no es una temperatura accesible]
        if (T != 0)
        {
            valor = exp(-E / T);
            if ((1.0) <= valor)
            {
                p = 1.0;
            }
            else
            {
                p = valor;
            }
        }
        else // Al dividir entre cero en el calculo de la exponencial el resultado de ese minimo es cero, menor que 1
        {
            p = 0.0;
        } 

        // Ahora generamos otro numero aleatorio, esta vez uniforme entre 0 y 1
        //CUIDADO TIENE QUE SER DIVISION REAL Y ESTAS TRABAJANDO CON ENTEROS
        eps = (double)rand() /(double) RAND_MAX ;

        if (eps <= p)
        {
            M[i][j] = (-1) * (M[i][j]);
        }

        //Guardo con un paso MONTECARLO la matriz para la animacion
        if(cuenta==(N*N))
        {
            for (int n = 0; n <= (N); n++)
            { // en las primeras posiciones del archivo van los datos y dejamos un espacio en blanco entre instante e instante

                if(n<N)
                {
                    for (int m = 0; m <= (N - 1); m++)
                    {
                    if (m < (N - 1))
                    {
                    fichero << M[n][m] << ", ";
                    }
                    else
                        fichero << M[n][m] << endl;
                    }
                }
                else
                {
                    fichero << endl;
                }
                
            }
            cuenta=0;
        }
        
        cuenta=cuenta+1;
    }


    //Muestro la ultima matriz en pantalla
    cout << "LA ULTIMA MATRIZ ES:" << endl;

    for (i = 0; i < N ; i++)
    { // en las primeras posiciones del archivo van los datos y dejamos un espacio en blanco entre instante e instante
        for (j = 0; j < N; j++)
        {
            if (j < (N - 1))
            {
                cout << M[i][j] << ", ";
            }
            else
            {
                cout << M[i][j] << endl;
            }
        }
    }



    fichero.close();
    return 0;
}