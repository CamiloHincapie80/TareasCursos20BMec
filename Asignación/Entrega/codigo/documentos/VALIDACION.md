# Validación y límites del resultado

## Verificaciones ejecutadas

1. Pruebas de física con unittest: Hohmann analítico frente a integración independiente de dos cuerpos; conservación de energía y momento angular; hipérbola integrada en ambas direcciones frente a desviación asintótica; resta exacta de aceleraciones; integración P/N en distintos marcos con tolerancia estricta; lector etiquetado/CSV, duplicados, incompatibilidad temporal y derivada Hermite. Las pruebas se ejecutan con el comando indicado en el README. La prueba integrada necesitó reducir tolerancia y paso: no se relajó su umbral para ocultar el fallo.
2. DOP853 con tolerancias 1e−9, 1e−11, 1e−12 y pasos 1 día, 6 h, 3 h; Radau como método alternativo. Sobre junio–agosto de 1979, frente al cálculo más estricto, los residuos DOP853 son 0.519 y 0.0172 km; Radau 0.0181 km. Esta convergencia local no certifica sola el error en 1981.
3. Prueba global de cambio de coordenadas: para f=1 los residuos P−N al 1981-10-01 pasan de 7045 km (rtol 1e−11, paso 1 d) a 23.2 km (1e−12, 12 h), 19.5 km (1e−13, 6 h) y 6.30 km (2.3e−14, 3 h). La tolerancia local no es una cota del error global. Los saltos de R'' y la sensibilidad del encuentro obligan a controlar pasos. No se afirma igualdad a precisión de máquina de las trayectorias integradas.
4. En K jovicéntrico, rango relativo de energía 2.33e−10 y de |h| 9.37e−11. En N varían físicamente: la energía pasa 29.437→29.191→29.341 km²/s² (entrada/periapsis/salida), y |h| 1.2643e7→1.4605e7→1.2737e7 km²/s. No se busca conservarlos en el problema perturbado.
5. Continuidad de todos los empalmes registrada en `eventos_*.json`, además de segmentos y marcos. El estado cambia de representación, no de valor físico. K que no entra en Saturno no se fuerza a hacerlo.
6. Excluir muestras alternas al construir Hermite da errores máximos de posición de 5.02 km para Voyager horario cerca de J y 0.000961 km para J horario. Son errores del interpolante de paso doble en muestras omitidas, no incertidumbres de los datos ni cotas garantizadas para el interpolante completo.
7. Hermite diario de Voyager contra estados horarios no utilizados: hasta 17995 km en 8–12 de julio de 1979. El problema no se arregla reduciendo rtol. Los centros planetarios diarios son mucho más suaves, pero al reducir a dos días la propagación junio–agosto cambia 939 km. Omitir Saturno en ese tramo cambia 3555 km: experimento distinto del error de cambio de marco.

## Referencia histórica

Voyager Horizons horario en J da aproximadamente 721017 km, 20.230 km/s, 1979-07-09 22:29:54 TDB, e osculante 1.3291 y δ osculante 97.595°. NASA NSSDCA da 10.11×71398=721833.8 km y 22:29:51 UTC. Al convertir TDB a UTC hay una diferencia aproximada de −47 s y −817 km. La precisión de la trayectoria de diseño no permite interpretar esa diferencia como error de medición.

N, iniciado el 1978-11-01 00 TDB, da rp_J=721796 km y 22:23:14 TDB: diferencia de unos 779 km respecto al interpolante horario. El ángulo entre velocidades **finitas** en la SOI N es 96.317°, distinto de δ osculante en periapsis 97.250°. Las energías solares de entrada/salida son −121.450 y +30.705 km²/s²; la ganancia entre esas fechas es 152.155 km²/s² y no debe confundirse con una diferencia en una fecha posterior.

En Saturno no hay muestreo horario. El Hermite diario da 151408 km y 1981-08-26 04:08:02 TDB, con rapidez 16.99 km/s y energía saturnocéntrica negativa (−106.19 km²/s²). Este resultado **no es una reconstrucción física fiable del periapsis** ni evidencia de captura. NASA NSSDCA da 161081.1 km y 03:24:57 UTC: el interpolante difiere unos −9674 km y +42 min tras conversión temporal aproximada. La tabla de comparación conserva el resultado, claramente marcado como insuficiente. No se puede cerrar una validación fina histórica de Saturno con estos TXT diarios.

N predice mínimo de 4.3412 millones de km, el 21 de agosto, y sí entra a SOI de Saturno. P corregida lo reproduce con pocos km de diferencia en distancia mínima para f=1. K con f=1 falla el encuentro: 70.811 millones de km, fuera de SOI. Las tablas de errores frente a Voyager son descriptivas contra un interpolante aproximado, no incertidumbres científicas certificadas.

## Qué causa eta

En la entrada SOI J, eta solar≈0.454, eta de Saturno≈6.54e−5 y eta del residuo de aceleración del origen≈0.853. En periapsis, ≈7.84e−7, 3.86e−10 y 1.03e−4, respectivamente. Las normas no se suman escalarmente; la eta total usa la suma vectorial. Por tanto no se debe llamar «marea solar» a toda la diferencia P0−N ni inferir solo lunas a partir del residuo.

## Diseño validado en su nivel de hipótesis

Hohmann desde LEO de 500 km: impulso 6.288 km/s, vuelo 997.726 d, exceso joviano 5.6433 km/s. Se adopta rp_J=714920 km (10 radios), superior al mínimo declarado 142984 km. El giro exterior de 115.917° conserva la rapidez asintótica; Δε_sol=105.886 km²/s², con residuo de la identidad del producto escalar ≈2.84e−14.

Se optimiza la fase de Saturno circular a 9.58 UA. La primera solución solar aproxima 500000 km y una fase de 24.8967° al encuentro con J. Al incluir Sol+Saturno y el término indirecto, se reajusta la fase a 24.9203° y se obtiene rp_S=500000.0001 km en unos 1163.123 d después del impulso asintótico J. La búsqueda detiene los ensayos que alcanzan el radio físico de Saturno; la solución aceptada es exterior. La fase inicial no ajustada falla por unos 279 millones de km en el modelo solar.

Esto responde a la existencia de un diseño **idealizado de cónicas asintóticas**, no a una misión continua de precisión con efemérides reales. No se afirma que el impulso instantáneo en la posición nominal de J sea una integración atravesando el planeta: es la idealización asintótica del documento original. La animación ajustada realiza otro experimento con propagación continua y Lambert como semilla. Para validación de navegación se necesitarían estados y maniobras de misión y una geometría de salida reconstruida.

## Pendientes explícitos

* Efemérides horarias de Voyager y Saturno alrededor del 26 de agosto de 1981 para una validación histórica fina; su procedencia sigue siendo trayectoria aproximada de diseño.
* SPICE/solución de navegación e historia de maniobras para separar error gravitacional de correcciones reales. Las fuentes narrativas no bastan para reconstruir impulsos.
* Comparación con baricentros de sistemas y satélites para atribuir el residuo del origen a cuerpos concretos. No se inventan datos faltantes.
* Contrastar el diseño ajustado Tierra–Júpiter–Saturno con un plano B y una época de navegación reconstruidos. La animación usa un blanco geométrico exterior de diseño.
* Completar las interacciones previas y la revisión independiente del equipo en la [bitácora 4D](BITACORA_4D.md), y ensayar con la [guía de defensa](DEFENSA.md).
