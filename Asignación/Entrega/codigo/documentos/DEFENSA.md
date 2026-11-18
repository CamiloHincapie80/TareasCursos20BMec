# Guía para la defensa oral

Las diapositivas siguen el orden de esta guía. Para preparar el equipo de dos o tres personas, registren quién asumió modelado, simulación y validación y cuándo intercambiaron esas funciones. Si son dos personas, compartan la validación. Cada integrante debe poder explicar los tres aspectos y resolver preguntas sobre el código.

1. **Plantear el reto.** Explicar cómo una fuerza central puede girar un vector velocidad y, al cambiar de marco, alterar la energía solar de la sonda.
2. **Dibujar el mecanismo.** En una figura heliocéntrica, marcar la Tierra, los dos puntos tangentes y Júpiter. En el marco de Júpiter, indicar periapsis y las dos velocidades asintóticas. En el de Saturno, señalar el mínimo y la rama saliente.
3. **Comparar con un mismo estado.** Mostrar el cuadro que mide el cambio en el mínimo saturniano, definir el umbral de 1000 km y explicar por qué se escoge una entrada distinta para cada factor de frontera.
4. **Defender la evidencia.** Mostrar la diferencia P–N como control de integración, el efecto de descartar fuerzas en K y la limitación de Horizons diario cerca de Saturno.

## Recorrido por el código

| Archivo | Qué debe poder explicar cada integrante |
| --- | --- |
| [horizons.py](../src/horizons.py) | Cómo se validan unidades, escala temporal y marcos. Por qué Hermite usa los vectores de velocidad. |
| [dinamica.py](../src/dinamica.py) | Término indirecto solar, continuidad de estado, detección de esferas y reinicios en nudos. |
| [encuentros.py](../src/encuentros.py) | Cómo se pasa de muestras a un mínimo refinado y por qué se distingue un mínimo en el borde. |
| [transferencia_leo.py](../src/transferencia_leo.py) | Cómo altura e impulso fijan el exceso terrestre, y cómo el giro joviano fija la velocidad de salida. |
| [optimizacion.py](../src/optimizacion.py) | Qué parámetro se ajusta para alcanzar Saturno circular y cómo se distingue blanco de periapsis. |
| [ajuste_saturno.py](../src/ajuste_saturno.py) | Cómo Lambert genera la semilla y qué residuos se minimizan en las integraciones con objetivos históricos. |
| [validacion.py](../src/validacion.py) | Qué se modifica en cada prueba de convergencia y cómo se diagnostican las fuerzas omitidas. |
| [main.py](../main.py) y [reporte.py](../herramientas/reporte.py) | De qué datos salen las cifras finales y dónde se guardan las figuras y tablas. |

## Preguntas para ensayar

- **¿Qué se conserva?** En un tramo central, la energía específica y el momento angular vectorial relativos al atractor. En una dispersión ideal aislada, también coincide el módulo de ambas velocidades asintóticas. La energía heliocéntrica puede variar porque se suma la velocidad planetaria a vectores relativos de distinta dirección.
- **¿Por qué conviene P?** Cambiar únicamente de coordenadas deja la fuerza física de N intacta. El residuo P–N ayuda a medir el efecto numérico antes de interpretar cambios de K como físicos.
- **¿Por qué K nominal histórico puede fallar aunque K ajustado llegue?** El primero hereda el estado de N y luego descarta fuerzas. El ajustado optimiza su propia velocidad de salida en presencia de esas omisiones. Ambos experimentos responden preguntas distintas.
- **¿Cómo se observa Saturno?** En el diseño circular el atractor está en la ecuación posterior a Júpiter. Se grafica antes y después del periapsis. En K, Saturno gobierna dentro de su esfera y el Sol vuelve a gobernar al salir.
- **¿Qué mediría un error histórico de alta precisión?** Estados de navegación con mayor resolución y un historial de maniobras, además de una dinámica planetaria más completa.
