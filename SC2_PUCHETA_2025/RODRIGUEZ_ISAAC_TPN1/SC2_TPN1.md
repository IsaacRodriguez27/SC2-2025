# Trabajo Práctico N°1 - Sistemas de Control II

## Tabla de contenidos
1. [Objetivo](#objetivo)
2. [Caso 1: Circuito RLC](#caso-1-circuito-rlc-2-variables-de-estado)
    - [1. Simulación con entrada tipo escalón](#1-simulación-con-entrada-tipo-escalón)
    - [2. Estimación de R, L y C a partir de datos experimentales](#2-estimación-de-r-l-y-c-a-partir-de-datos-experimentales)
    - [3. Validación con curva de corriente medida](#3-validación-con-curva-de-corriente-medida)
3. [Caso 2: Motor de Corriente Continua (3 variables de estado)](#caso-2-motor-de-corriente-continua-3-variables-de-estado)
4. [Código en Matlab](#código-en-matlab)
5. [Conclusiones](#conclusiones)

---


## Objetivo

El objetivo del presente trabajo es modelar y analizar diferentes sistemas dinámicos (lineales y no lineales) utilizando simulaciones, herramientas de integración, métodos de identificación como Chen, y controladores PID.

---

## Caso 1: Circuito RLC 
![
](<Imagenes/circuito RLC.png>)
![
](<Imagenes/Variables estado.png>)

### 1. Simulación con entrada tipo escalón

Al simular podemos ver el comportamiento de tension y corriente en el capacitor.

![
](Imagenes/v_e_i_RLC.png)

El circuito RLC se comporta así porque está subamortiguado: la combinación de inductancia y capacitancia genera oscilaciones cada vez que cambia la entrada tipo escalon. Como la resistencia no es suficientemente alta para amortiguarlas rápidamente, la corriente y la tensión en el capacitor muestran un comportamiento ondulante con resonancia antes de estabilizarse.

### 2. Estimación de R, L y C a partir de datos experimentales

Se utilizan los datos del archivo `Curvas_Medidas_RLC_2025.xls`, específicamente la hoja 1, que contiene los valores temporales de la respuesta del sistema ante una entrada escalón, considerando como salida la **tensión en el capacitor**.  
Se aplica el **método de la respuesta al escalón(CHEN)** para estimar los parámetros \( R, L, C \).  
> Aclaración: se recomienda seleccionar los tres puntos característicos (p.ej., valor inicial, punto de inflexión y valor final) para aplicar el **método de Chen**.

![
](Imagenes/Tension_excel.png)

Podemos ver el comportamiento de la tension del capacitor en el tiempo segun los datos obtenidos en el metodo experimental.

---

### 3. Validación con curva de corriente medida

Una vez obtenidos los parámetros estimados, se vuelve a simular el circuito y se superpone la curva simulada con la **curva medida de la corriente** (disponible a partir de 0.05 s), comparando ambas gráficas.

> Se busca verificar si los valores hallados representan adecuadamente el comportamiento dinámico real del sistema.

---

## Caso 2: Motor de Corriente Continua (3 variables de estado)

📌 *(Este apartado se desarrolla luego de terminar el Caso 1, ver ítems [4] a [8] del enunciado)*

---

## Código en Matlab

```matlab
% Fragmentos clave de simulación y análisis

