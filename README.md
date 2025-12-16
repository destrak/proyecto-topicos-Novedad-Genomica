# proyecto-topicos-Novedad-Genomica

Este repositorio contiene la implementación y experimentación para el cálculo de
novedad genómica utilizando Minimizers y HyperLogLog (HLL). Incluye el código fuente,
procedimiento de compilación, ejemplos de ejecución y un conjunto de datos de prueba,
tal como se solicita en el enunciado.

---

## Requisitos

- Sistema operativo Linux
- Compilador con soporte C++17 (g++ >= 9)
- Python 3 (opcional, para gráficos)

---

## Contenido del repositorio

- minimizers.cpp : generación de minimizers desde genomas FASTA/FNA
- hll.cpp        : construcción de sketches HyperLogLog (.hll)
- catalogo.cpp   : construcción del catálogo genómico R
- compararhll.cpp: cálculo de novedad genómica y exportación a CSV
- grafico.py     : generación de gráficos desde resultados
- Dataset de prueba en formato FASTA/FNA

---

## Compilación

Desde la carpeta raíz del repositorio:

g++ -O2 -std=c++17 minimizers.cpp  -o minimizers
g++ -O2 -std=c++17 hll.cpp         -o hll
g++ -O2 -std=c++17 catalogo.cpp    -o catalogo
g++ -O2 -std=c++17 compararhll.cpp -o compararhll

---

## Ejecución

### 1) Generación de minimizers

./minimizers

Parámetros solicitados:
- k: tamaño del k-mer (ej. 21 o 31)
- w: tamaño de la ventana
- carpeta de entrada con FASTA/FNA
- carpeta de salida

---

### 2) Construcción de sketches HLL

./hll

Parámetros solicitados:
- carpeta con minimizers
- carpeta de salida para sketches

---

### 3) Construcción del catálogo genómico

./catalogo

Parámetros solicitados:
- carpeta con sketches .hll
- nombre del archivo del catálogo
- número de genomas a seleccionar (N)

Salidas:
- catalogo_R.hll
- catalogo_R.hll_usados.txt
- catalogo_R.hll_no_usados.txt

---

### 4) Cálculo de novedad genómica

./compararhll catalogo_R.hll carpeta_hll catalogo_R.hll_usados.txt resultados.csv

---

## Definición de novedad genómica

|S \ R| ≈ |R ∪ S| − |R|
ρ = |S \ R| / |S|

---

## Observaciones

- Se utiliza p = 14 (16384 registros HLL)
- Error relativo esperado cercano al 1 %
- Memoria constante por genoma (~16 KB)

