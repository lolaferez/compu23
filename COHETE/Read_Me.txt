Mediante el uso del algoritmo de runge-kutta 4 se ha resuelto un sistema de ecuaciones diferenciales de primer orden que caracteriza el hamiltoniano y las ecuaciones de movimiento de una nave que despega a la Luna desde la Tierra.

No se ha realizado el problema para un paso h variable. 
Se ha considerado como tiempo máximo de la simulación una semana y se ha trabajado en unidades del SI (segundos y metros). La masa de cohete para el cálculo de H' se ha elegido 777kg. 

Se adjuntan archivos .txt con los valores de Hamiltoniano prima, los k_i para la resolución del sistema de ecuaciones diferenciales, las posiciones de Tierra, Luna y cohete para las animaciones y finalmente en "tp.txt" los valores de cada magnitud del sistema de ecuaciones en cada instante.

Se adjuntan dos vídeos de la trayectoria, uno de ellos ampliado.
Así como una gráfica resultado de plotear el valor de H', mostrando que este se mantiene constante a lo largo de toda la trayectoria salvo al interactuar con el campo gravitatorio de la Luna que perturba los valores, pero finalmente se estabilizan en el valor inicial.

*El uso de la unidad segundos ha podido influir en el error para el cálculo de H'.