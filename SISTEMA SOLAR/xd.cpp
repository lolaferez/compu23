#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

// Definimos fuera del main las constantes porque son valores que no van a cambiar y ya
// sabe lo que tiene que ocupar de espacio.

#define MS 1.989e30   // Masa del Sol
#define G 6.67384e-11 // Constante de gravitacion universal
#define c 1.4953e11   // Distancia en metros entre la Tierra y el Sol

// Funciones para reescalar las magnitudes de forma que la constante de gravitacion universal G valga 1.

double tiempo(double t)
{
    double tiemp;
    tiemp = sqrt(G * MS / (c * c * c)) * t;
    return tiemp;
}

double tiempoinver(double t)
{
    double tiemp;
    tiemp = t * 58.1; // Asi devuelve el resultado en dias
    return tiemp;
}

double masa(double m)
{
    double mas;
    mas = m / MS;
    return mas;
}

double posicion(double r)
{
    double pos;
    pos = (r) / c;
    return pos;
}

double velocidad(double v) // Para reescalarla, dividimos la posicion reescalada entre el tiempo reescalado
{                        // y lo expresamos en funcion de la velocidad.
    double vel;
    vel = v / sqrt(G * MS / (c));
    return vel;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// NOTA: no se va a considerar Pluton como planeta para agilizar la ejecucion del programa.
// NOTA2: Cada planeta se aproxima a una masa puntual, despreciando su radio ya que es mucho menor que su
// distancia al sol.

////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{

    // Definimos los vectores con las magnitudes propias de cada planeta (en orden desde mas cercanos
    // a mas lejanos al sol en el Sistema Solar).

    double m[9], r[9][2], v[9][2], a[9][2], w[9][2]; // Denotan la masa, posicion, velocidad y aceleracion y w respectivamente
                                                    // Los vectores posicion, velocidad y aceleracion son matrices pues tienen componente x e y en el plano.
    double rmasa[9], rpos[9][2], rvel[9][2]; // Los vectores reescalados
    double h, t, k;
    double cinetica[9], potencial[9], energia[9], periodo[9];

    ifstream fichm, fichr, fichv;                             // Ficheros de lectura con los datos iniciales
    ofstream ficheroposicion, ficheroenergia, ficheroperiodo; // Ficheros sobre los que escribir los datos

    int i, j;

    i = 0;
    j = 0;
    k = 0;
    t = 0.0;

    fichm.open("masas.txt");
    // Vamos a tomar los datos del .txt y llenar el array m definido anteriormente.
    while (!fichm.eof())
    { // Para comprobar que el fichero contiene datos.

        for (i = 0; i <= 8; i++)
            fichm >> m[i];
    }

    fichm.close();

    fichr.open("posiciones.txt");
    while (!fichr.eof())
    { // Para comprobar que el fichero contiene datos.
        // Vamos a tomar los datos del .txt y llenar el array p definido anteriormente.

        for (i = 0; i <= 8; i++)
        {
            for (j = 0; j <= 1; j++)
                fichr >> r[i][j];
        }
    }

    fichr.close();

    fichv.open("velocidades.txt");
    while (!fichv.eof())
    { // Para comprobar que el fichero contiene datos.
        // Vamos a tomar los datos del .txt y llenar el array v definido anteriormente.

        for (i = 0; i <= 8; i++)
        {
            for (j = 0; j <= 1; j++)
                fichv >> v[i][j];
        }
    }
    fichv.close();

    // Una vez tenemos los vectores definidos, los reescalamos usando las funciones previamente
    // definidas, posicion a posicion del vector

    // Llenamos los nuevos vectores reescalados
    for (i = 0; i <= 8; i++)
        rmasa[i] = masa(m[i]);

    for (i = 0; i <= 8; i++)
    {
        rpos[i][0] = posicion(r[i][0]);
        // Solo se reescala la coordenada x pues la y inicialmente es nula.
        rpos[i][1] = 0.0;
    }

    for (i = 0; i <= 8; i++)
    {
        rvel[i][1] = velocidad(v[i][1]);
        // Igual que con la posicion, solo se reescala la coordenada y
        // pues la x inicialmente es nula.
        rvel[i][0] = 0.0;
    }

    // Muestra en la terminal los valores reescalados de las condiciones iniciales
    cout << "MASAS REESCALADAS:" << endl;
    for (i = 0; i < 9; i++)
        cout << rmasa[i] << endl;

    cout << endl
         << "POSICIONES REESCALADAS:" << endl;
    for (i = 0; i < 9; i++)
    {
        cout << rpos[i][0] << "  " << rpos[i][1] << endl;
    }

    cout << endl
         << "VELOCIDADES REESCALADAS:" << endl;
    for (i = 0; i < 9; i++)
    {
        cout << rvel[i][0] << "  " << rvel[i][1] << endl;
    }

    // Comenzamos evaluando la aceleracion con los datos iniciales
    for (i = 0; i <= 8; i++)
    {
        // Inicio a cero el array
        a[i][0] = 0.0;
        a[i][1] = 0.0;
    }

    cout << endl << "ACELERACIONES INICIALES: " << endl;
    for (i = 0; i <= 8; i++)
    {
        for (j = 0; j < 9; j++)
        {
            if (j != i)
            {
                a[i][0] = a[i][0] - rmasa[j] * (rpos[i][0] - rpos[j][0]) / (pow(sqrt((rpos[i][0] - rpos[j][0]) * (rpos[i][0] - rpos[j][0]) + (rpos[i][1] - rpos[j][1]) * (rpos[i][1] - rpos[j][1])), 3));
                a[i][1] = a[i][1] - rmasa[j] * (rpos[i][1] - rpos[j][1]) / (pow(sqrt((rpos[i][0] - rpos[j][0]) * (rpos[i][0] - rpos[j][0]) + (rpos[i][1] - rpos[j][1]) * (rpos[i][1] - rpos[j][1])), 3));
            }
        }
        cout << a[i][0] << "  " << a[i][1] << endl;
    }

    cout << endl
         << "Introduzca el paso h para las iteraciones: ";
    cin >> h;

    // Incializo el vector de periodos en cero
    for (i = 0; i < 9; i++)
    {
        periodo[i] = 0.0;
    }

    // Ahora procedemos a emplear el algoritmo de Verlet
    // Comenzamos abriendo los ficheros donde queremos escribir los datos de las iteraciones

    ficheroposicion.open("posicionesitera.txt");
    ficheroenergia.open("energia.txt");
    ficheroperiodo.open("periodo.txt");

    for (k = 0; k <= 100000; k++)
    {

        // A continuacion, obtenemos los valores de posicion para el instante siguiente usando un desarrollo en serie de Taylor truncado
        // y los w.

        for (i = 0; i <= 8; i++)
        {
            rpos[i][0] = (rpos[i][0]) + h * (rvel[i][0]) + (a[i][0]) * (h * h * 0.5);
            rpos[i][1] = (rpos[i][1]) + h * (rvel[i][1]) + (a[i][1]) * (h * h * 0.5);
        }

        // Los w los inicializamos a cero asi como los vectores de energia
        for (i = 0; i < 9; i++)
        {
            energia[i] = 0.0;
            potencial[i] = 0.0;
            cinetica[i] = 0.0;
            for (j = 0; j < 2; j++)
                w[i][j] = 0.0;
        }

        for (i = 0; i <= 8; i++)
        {
            w[i][0] = (rvel[i][0]) + h * (a[i][0]) * 0.5;
            w[i][1] = (rvel[i][1]) + h * (a[i][1]) * 0.5;
        }

        // Reevaluamos la aceleracion en las nuevas posiciones
        for (i = 0; i <= 8; i++)
        {
            a[i][0] = 0.0;
            a[i][1] = 0.0;
            for (j = 0; j < 9; j++)
            {
                if (j != i)
                {
                    a[i][0] = a[i][0] - rmasa[j] * (rpos[i][0] - rpos[j][0]) / (pow(sqrt((rpos[i][0] - rpos[j][0]) * (rpos[i][0] - rpos[j][0]) + (rpos[i][1] - rpos[j][1]) * (rpos[i][1] - rpos[j][1])), 3));
                    a[i][1] = a[i][1] - rmasa[j] * (rpos[i][1] - rpos[j][1]) / (pow(sqrt((rpos[i][0] - rpos[j][0]) * (rpos[i][0] - rpos[j][0]) + (rpos[i][1] - rpos[j][1]) * (rpos[i][1] - rpos[j][1])), 3));
                }
            }
        }

        // Calculamos las velocidades en el instante siguiente
        for (i = 0; i <= 8; i++)
        {
            rvel[i][0] = w[i][0] + h * 0.5 * a[i][0];
            rvel[i][1] = w[i][1] + h * 0.5 * a[i][1];
        }

        // Calculo de la energía para cada planeta en el instante dado y posteriormente su promedio (es lo que
        // queremos ver que se mantiene constante)
        for (i = 0; i < 9; i++)
        {
            cinetica[i] = rmasa[i] * 0.5 * (rvel[i][0] * rvel[i][0] + rvel[i][1] * rvel[i][1]);
            for (j = 0; j < 9; j++)
            {
                if (j != i)
                {
                    potencial[i] = potencial[i] - rmasa[i] * rmasa[j] / (sqrt((rpos[i][0] - rpos[j][0]) * (rpos[i][0] - rpos[j][0]) + (rpos[i][1] - rpos[j][1]) * (rpos[i][1] - rpos[j][1])));
                }
            }
            energia[i] = cinetica[i] + potencial[i];
        }

        if (k == 0)
        {
            ficheroperiodo << "El periodo en dias de cada planeta en el Sistema solar es: " << endl;
        }

        // Calculo el periodo de cada planeta, menos el sol pero lo calculo solo una vez
        for (i = 1; i < 9; i++)
        {
            if (periodo[i] == 0.0)
            {
                // Lo calculamos cogiendo el instante en que la coordenada y pasa de ser negativa a positiva. De paso lo metemos en el fichero
                if (((rpos[i][1]) < 0.0) && ((rpos[i][1]) + h * (rvel[i][1]) + (a[i][1]) * (h * h * 0.5)) > 0.0)
                {
                    periodo[i] = t;
                    ficheroperiodo << tiempoinver(periodo[i]) << endl;
                }
            }
        }

        // Escribimos las posiciones, el promedio de la energia energia y periodo obtenidos en el fichero deseado
        for (i = 0; i <= 9; i++)
        {
            if (i < 9)
                ficheroenergia << energia[i] << endl;
            else
                ficheroenergia << endl;
        }

        for (i = 0; i <= 9; i++)
        { // en las primeras posiciones del archivo van los datos y dejamos un espacio en blanco entre instante e instante
            if (i < 9)
                ficheroposicion << rpos[i][0] << ", " << rpos[i][1] << endl;
            else
                ficheroposicion << endl;
        }

        t = t + h;
    }

    ficheroenergia.close();
    ficheroperiodo.close();
    ficheroposicion.close();

    return 0;
}