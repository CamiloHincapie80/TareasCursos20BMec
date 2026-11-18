# Hoja de hipótesis y decisiones de modelado

La comparación principal parte de un mismo estado de la sonda al entrar en la región de Júpiter. El diseño, en cambio, busca condiciones iniciales que permitan llegar a Saturno. La distancia de encuentro se mide desde el centro de cada planeta.

| Decisión | Modelo elegido | Alternativa | Justificación | Impacto esperado |
| --- | --- | --- | --- | --- |
| Sistema físico | Sonda de prueba y planetas puntuales | Integrar un sistema planetario completo | Aislar las fuerzas esenciales | Se omiten lunas, achatamiento, radiación y maniobras. |
| Fuerzas de K | Sol fuera y planeta dentro de cada esfera | Fuerzas simultáneas de N | Evaluar cónicas empalmadas | El error se acumula hasta Saturno. |
| Frontera de K | Radio de Laplace instantáneo por un factor | Radio fijo o esfera de Hill | Medir sensibilidad a la frontera | Cambian la entrada y la duración de cada tramo. |
| Marco de P | Traslación con la aceleración completa del origen | Fuerza planetaria y marea solar en P0 | Controlar el efecto del cambio de coordenadas | P y N deben converger. P0 permite medir el efecto de las omisiones. |
| Datos de Voyager 2 | Estado de noviembre de 1978. Datos horarios en Júpiter y diarios en Saturno | Integrar desde lanzamiento con gravedad terrestre y resolver ambos encuentros por hora | Aprovechar los datos disponibles | Saturno histórico tiene menor resolución y requiere contraste externo. |
| Diseño circular | Periapsis joviano de 10 radios y fase libre de Saturno | Ajustar época y geometría con efemérides | Exponer el mecanismo de dispersión | El resultado circular depende de la fase elegida. |
| Criterio de equivalencia | Cambio menor de 1000 km en el mínimo saturniano | Error de posición a fecha fija | Medir el efecto sobre el objetivo del reto | El umbral depende del objetivo adoptado. |

**Predicciones físicas.** La rapidez jovicéntrica asintótica se conserva en el sobrevuelo ideal. El giro seleccionado aumenta la energía solar. Al acercarse a Júpiter, la razón de marea solar a fuerza joviana disminuye aproximadamente como el cubo de la distancia. Si solo cambia el marco, las trayectorias deben converger al reducir el paso y la tolerancia.

## Experimentos y condiciones de contorno

| Experimento | Estado inicial y frontera | Observables |
| --- | --- | --- |
| Histórico, `main.py todo` | Estado de Horizons del 1 de noviembre de 1978 TDB para N. Para cada $f$, transformar al marco de Júpiter el estado de N en su primer cruce entrante. P, P0 y K empiezan allí. Integración principal hasta el 1 de octubre de 1981, con una extensión separada hasta diciembre para registrar la salida de Saturno. | Estado a fecha fija, mínimos de distancia y épocas, cambios de marco e invariantes. |
| Circular, `src/transferencia_leo.py` y `src/optimizacion.py` | LEO nominal de 500 km, impulso prógrado tangencial. Órbitas heliocéntricas de radios $1.49598\times10^8$ km, $7.78479\times10^8$ km y 9,58 UA. Periapsis joviano 714 920 km. La salida asintótica rotada se suma a la velocidad orbital de Júpiter. Se resuelve la fase saturniana para un mínimo de 500 000 km. | Energía ganada, giro, fecha del mínimo saturniano y curvatura en el marco de Saturno. |
| Animación ajustada a efemérides, `src/ajuste_saturno.py` | Inicio del 20 de agosto de 1977 a las 14:29 UTC. Se aproxima la posición planetaria hasta el comienzo del archivo del 5 de septiembre por propagación de dos cuerpos hacia atrás. Dos arcos Lambert generan una semilla prógrada. K ajusta la velocidad inicial. El continuo ajusta geometría y época y vuelve a propagar desde la Tierra. El blanco está desplazado un millón de km en la normal eclíptica respecto al centro de Saturno. | Residuo al blanco, mínimo al centro, pasos por las esferas. Voyager se muestra desde la primera muestra disponible. |

El blanco de la animación y el mínimo geométrico son observables distintos. Hohmann es el arco circular coplanar del diseño y usa una fase joviana prescrita. La variante ajustada usa efemérides y solución Lambert como semilla. El informe prioriza el experimento histórico controlado y el diseño circular.

## Ecuaciones y discretización

El estado que integra `System.run` es $y=(\mathbf r,\mathbf v)$ con $y'=(\mathbf v,\mathbf a)$ en el marco seleccionado. DOP853 construye cada paso por etapas de Runge–Kutta explícitas y de tamaño adaptativo. Se controla el error local con escalas de la forma `atol + rtol * abs(y)`. La propagación histórica usa `rtol=1e-12`, `atol=1e-8` y un paso máximo de 12 h. La trayectoria ajustada a efemérides usa `1e-11`, `1e-7` y 3 días de paso máximo. Los límites de paso y tolerancias se prueban aparte. Las configuraciones se reflejan en `src/dinamica.py`, `src/ajuste_saturno.py` y `src/validacion.py`.

La frontera es el cero de la distancia relativa menos el radio de Laplace. La dirección del evento distingue entrada y salida. En el cruce se resta o suma el estado del centro planetario antes de continuar. Los nudos de Hermite son límites de segmento para P en marco planetario debido al salto en la segunda derivada del interpolante. El mínimo se aísla con una malla y se refina por minimización escalar usando la solución densa. Un mínimo en el borde se marca explícitamente.

En el diseño circular, Brent resuelve la fase objetivo y la propagación Sol–Saturno termina las trayectorias de prueba que alcanzan el radio del planeta. En la animación ajustada, mínimos cuadrados reduce residuos de posición y restricciones jovianas, con diferencias centrales en los refinamientos sensibles. La malla dibujada por la interfaz es un muestreo de la solución, mientras que tablas de mínimos y eventos usan la solución densa.

Las derivaciones del cambio de marco y de los invariantes están en [FISICA.md](FISICA.md). Las pruebas de resolución y el contraste histórico están en [VALIDACION.md](VALIDACION.md).
