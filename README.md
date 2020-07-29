# Simplex
## mSimplexFaseII_eq
Implementación de la Fase II para un problema de programación lineal del tipo
```
min c⊤x
sujeto a Ax=bx≥0.
```
Para esta función puede suponer que conoce un conjunto de índices B de una SBF. Es importante mantener los conjuntos de índices básicas B y no básicas N. Las entradas y salidas de la función se describen a continuación:
In : 
- A mxn matrix
- b vector columna con  m  renglones
- c vector columna con  n  renglones
- jB vector de indices de la SBF inicial.

Out: 
- xo SFB óptima del problema
- zo valor óptimo del problema
- ban 
  - 0  si se encontro una solución óptima
  - 1  si la función objectivo no es acotada.
- iter es el número de iteraciones (cambios de variables basicas) que hizo el método
- B vector de indices de la SBF optima

## mSimplex_leq
Usa la función mSimplexFaseII_eq para resolver problemasde tipo:
```
min c⊤x
sujeto a Ax≤bx≥0.
```
Esa función debe introducir variables de holgura y tal vez una variable artificial (depende deb).Eso lo debe hacer antes de pasar el problema a la función mSimplexFaseII_eq.
In :
- A mxn matrix
- b vector columna con  m  renglones
- c vector columna con  n  renglones

Out: 
- xo SFB óptima del problema
- zo valor óptimo del problema
- ban 
  - -1  si el conjunto factible es vacio
  - 0  si se encontro una solución óptima
  - 1  si la función objectivo no es acotada.
- iter es el número de iteraciones (cambios de variables basicas) que hizo el método
- B vector de indices de la SBF optima
