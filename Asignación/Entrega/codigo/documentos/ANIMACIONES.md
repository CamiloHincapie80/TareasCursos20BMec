# Transferencia de Hohmann y comparación orbital

```bash
./abrir_animaciones.sh
```

El lanzador usa el entorno local o el intérprete indicado en la variable
PYTHON. Desde Entrega/codigo, también puede utilizar el entorno del proyecto
si ambas carpetas siguen juntas. La primera preparación sin caché puede
tardar varios minutos.

Para abrir la ventana 3D, ejecuta el lanzador desde una terminal de la
sesión gráfica. `--comprobar --sin-pantalla` comprueba la escena y sus
controles con renderizado EGL y se cierra al terminar.

## Transferencia de Hohmann

La primera pestaña muestra una órbita circular terrestre, el impulso tangencial
en el perigeo de salida, la cónica terrestre y la transferencia heliocéntrica
en el modelo circular de las diapositivas. El impulso se fecha el **20 de agosto
de 1977 a las 14:29 UTC**, al igual que el lanzamiento de Voyager 2.

La altura y el impulso son independientes. **Usar Δv de Hohmann** calcula el
valor tangente para la altura elegida. Con 500 km es **6,287952 km/s** y el
vuelo heliocéntrico nominal dura **997,73 días**. Un valor distinto puede
mantener la sonda ligada a la Tierra, dejar el afelio por debajo de la
órbita de Júpiter o cambiar la hora y el ángulo de cruce. La fase de Júpiter
en la escena circular se mantiene en su valor nominal para mostrar el efecto
de modificar el impulso. El control de reproducción recorre la geometría por
etapas.

En el caso nominal, la reproducción pasa a una hipérbola jovicéntrica de
dos cuerpos con periapsis de 10 radios de Júpiter. El giro asintótico es
115,92° y conserva los 5,643 km/s de rapidez relativa. La velocidad de
salida se suma al vector de velocidad orbital de Júpiter y se propaga
con Sol y Saturno circular. La fase de Saturno corresponde al diseño
idealizado que deja 500.000 km de distancia mínima a su centro. La vista
saturnocéntrica muestra los cinco días anteriores y posteriores al mínimo. La
curva solar continúa después del encuentro con la misma integración. Alterar
el impulso conserva la fase nominal de Júpiter y muestra el desfase al
cruce. El sobrevuelo dibujado requiere la geometría nominal.

## Comparación hasta Saturno

Una escena 3D y un reloj compartido superponen tres trayectorias:

- **Cian: gravedad continua ajustada**, transferencia asistida con dispersión gravitatoria continua en Júpiter, dirigida al Saturno histórico.
- **Coral: esferas de influencia**, una sola trayectoria de cónicas empalmadas (K).
- **Violeta: Voyager 2**, datos originales de Horizons en su marco heliocéntrico J2000.

Los antiguos nombres `hohmann`, `esferas` y `ambas` se aceptan como alias de la misma vista. El factor de las esferas puede configurarse con `./abrir_animaciones.sh --factor 0.5` (0.1, 0.25, 0.5 o 1).

## Salida desde la Tierra y comparación

Ambos modelos parten el **20 de agosto de 1977** desde la posición terrestre estimada para el lanzamiento. Para cubrir el intervalo entre ese día y el comienzo de los archivos planetarios, el 5 de septiembre, se prolongan hacia atrás los estados de Tierra, Júpiter y Saturno mediante dinámica solar de dos cuerpos. Sus velocidades iniciales se ajustan por separado para buscar un encuentro con Júpiter y después Saturno. No se integra el escape terrestre ni se modela la gravedad terrestre: la posición de la Tierra representa la salida en el modelo heliocéntrico. K no copia ningún prefijo de N.

El cálculo inicial resuelve dos transferencias Lambert prógradas de cero revoluciones y ajusta las fechas del encuentro para igualar las velocidades hiperbólicas jovianas de entrada y salida, con un periapsis preliminar de 10 radios jovianos. Ese diseño sirve de semilla, no de resultado final certificado.

- **Esferas (K):** se ajusta mediante disparo la velocidad de salida terrestre de su integración continua de cónicas empalmadas. Se usa el Sol fuera de las SOI y Júpiter o Saturno dentro, conservando posición y velocidad al cambiar de marco. Se apunta a un punto separado un millón de km del centro de Saturno en la fecha objetivo. La aceleración durante la aproximación determina el mínimo real, que se mide después con la solución densa.
- **Gravedad continua ajustada:** se ajustan la orientación del estado en periapsis, el exceso hiperbólico y las fechas joviana y saturniana mediante integración en ambos sentidos desde el periapsis. Se impone un radio de 10 radios jovianos y el objetivo exterior de Saturno. Con la velocidad terrestre convergida se propaga una única solución heliocéntrica Sol+Júpiter+Saturno, sin impulsos ni cambios de fuerza en la SOI. El sobrevuelo es una hipérbola perturbada por el Sol; sus elementos jovicéntricos son osculantes. El cálculo circular clásico original permanece intacto en `src/hohmann.py`.

Se comprueba que ambos modelos tengan un mínimo interior fuera del radio de Saturno y dentro de su SOI. También se verifica que no intersecten Júpiter.

La escena dura **10 años (3652.5 días)**, del 20 de agosto de 1977 al 20 de agosto de 1987 aproximadamente. Las trayectorias se siguen integrando después de Saturno; no se prolongan con posiciones inmóviles. La duración también garantiza al menos dos años después de la fecha objetivo. Se añaden muestras horarias alrededor de los encuentros y los instantes exactos de cambio de marco. En el sobrevuelo cian se añaden muestras cada cinco minutos durante cuatro días alrededor del periapsis.

Los planetas de la segunda pestaña y el objetivo de los dos modelos usan las mismas efemérides, con la aproximación inicial indicada arriba. Voyager conserva coordenadas y épocas originales, sin giro ni desplazamiento. Su primer estado disponible ya está a unos 13.34 millones de km de la Tierra: no se inventa un lanzamiento. La curva de Voyager se muestra desde el comienzo de la cobertura local el 5 de septiembre.

## Controles

| Control | Función |
|---|---|
| Reproducir / Pausar / Espacio | Controla el reloj común |
| Barra de tiempo / Inicio | Busca un instante o vuelve al comienzo |
| Días/s / Repetir | Velocidad temporal y repetición |
| Botones por color | Muestran u ocultan cada una de las tres trayectorias |
| Marco solar / Tierra / Júpiter / Saturno | Cambia el origen común de la escena |
| Cámara inicial / R | Restablece el encuadre |
| Arrastrar / rueda | Gira / acerca la cámara 3D |
| Recorridos completos | Muestra u oculta la trayectoria futura tenue |
| Órbitas planetarias | Muestra u oculta las curvas planetarias, conservando sus marcadores |
| Esferas de influencia | Muestra u oculta las fronteras de K |
| Ir a un evento… | Salta al lanzamiento, al cruce terrestre o al mínimo en Saturno |

Las coordenadas visuales están en UA. Marcadores y textos están ampliados para facilitar su lectura; las SOI conservan su escala. Cambiar el marco resta la posición planetaria de cada instante en float64 antes de enviar los datos a OpenGL.

## Tres tarjetas de datos

Cada trayectoria tiene una tarjeta con su color y estos valores:

- Velocidad heliocéntrica actual, en km/s.
- Energía orbital específica solar actual, ε = v²/2 − μ☉/r, en km²/s². No incluye el potencial planetario ni necesita la masa de la sonda.
- Distancia mínima al **centro** de Júpiter durante el intervalo completo, en km.
- Ángulo de dispersión δ, en grados.
- Distancia mínima al centro de Saturno y fecha de ese mínimo, para comprobar la llegada.

Los mínimos y δ son resultados del encuentro completo, no valores acumulados hasta el instante mostrado. Los mínimos de los modelos se refinan con sus soluciones densas y las efemérides históricas. Voyager usa los interpolantes horarios de la sonda y Júpiter del 1 al 20 de julio de 1979; su mínimo saturniano se estima con la serie diaria y no se certifica como periapsis de navegación.

Para los tres modelos, δ = 2 asin(1/e) se calcula con la hipérbola osculante jovicéntrica en el mínimo: es una estimación de dos cuerpos, no el cambio de dirección heliocéntrico. Para la trayectoria de gravedad continua también se registran a = −μ_J/v∞², rp = a(1−e) y b = |a|√(e²−1). El periapsis y la velocidad determinan estos parámetros; el ángulo ya no es un giro prescrito. Los estados y fechas ajustados quedan registrados en los metadatos de la caché, que se invalida al cambiar el código físico.

Las posiciones y velocidades se interpolan visualmente entre muestras. Voyager conserva la serie horaria en Júpiter y la diaria fuera. Su tramo anterior a 1989 está identificado por Horizons como trayectoria aproximada; el muestreo diario de Saturno no certifica periapsis ni rapidez históricos.

## Preparación y comprobaciones

```bash
./abrir_animaciones.sh --preparar
./abrir_animaciones.sh --recalcular --preparar
.venv/bin/python -m unittest discover -s tests -v
./abrir_animaciones.sh --comprobar --sin-pantalla
```

Las cachés se guardan en `resultados/animaciones` con hashes de datos, código, factor y condiciones de cálculo. La comprobación sin pantalla requiere EGL y Qt offscreen: verifica una sola escena con las tres trayectorias, las tarjetas, visibilidad, marcos y el reloj común. No genera capturas. Abrir la ventana interactiva requiere una terminal con acceso al escritorio y OpenGL.

Las funciones de comparación histórica N/P/P0/K anteriores siguen disponibles para los cálculos científicos y sus pruebas; la interfaz muestra únicamente K como modelo de esferas de influencia.

### Salida en el marco de la Tierra

Selecciona **Marco de la Tierra** para centrar la vista en el planeta con ejes
paralelos a los eclípticos J2000. La escena muestra la esfera de influencia
terrestre de Laplace a su escala nominal, independiente del factor aplicado a
Júpiter y Saturno. Usa **Ir a un evento…** para saltar al cruce saliente de cada
modelo o al inicio de los datos. Reduce el avance para seguir las primeras horas.

Las dos trayectorias ajustadas parten de la posición heliocéntrica de la
Tierra con velocidades de salida prescritas. La fecha mostrada corresponde
al cruce geométrico de la esfera de influencia, pues la simulación comparativa
omite la órbita de estacionamiento y la gravedad terrestre. Los datos locales de
Voyager 2 empiezan el 5 de septiembre de 1977, cuando la sonda ya se encuentra
fuera de la esfera, por lo que aparece «sin datos» el día del lanzamiento. Una vez disponible su trayectoria, se indica que el cruce fue anterior a los datos.

Los tres botones con los nombres de las trayectorias empiezan activados. Pulsa
cualquiera para ocultar o volver a mostrar su recorrido, marcador y etiqueta.
