# Bitácora de interacción con IA

**Equipo:** Santiago Delgado Cobos (2241400) y Angel Camilo Hincapie Quintero (2241945).

**Contacto:** hincapiec837@gmail.com.

**Herramienta:** asistente de programación Codex.
**Alcance de este registro:** interacciones visibles en la sesión de trabajo. Los espacios pendientes corresponden a interacciones previas y verificaciones humanas sin registro disponible en esta sesión.

Metodología de referencia. [Material 4D indicado en el enunciado](https://github.com/nunezluis/MisCursos/blob/main/MisMateriales/Presentaciones/MC_S01_InteraccionHumanoIA.pdf).

## 1. Delegación

El equipo pidió asistencia para documentar los modelos, identificar las variables de ajuste, programar vistas comparativas y auditar el flujo de entrega. El asistente implementó los cambios y ejecutó comprobaciones numéricas y gráficas. Las decisiones previas de autoría y la aceptación final de los supuestos quedan a cargo del equipo.

## 2. Descripción

| Petición visible del equipo (fragmento) | Tarea concreta |
| --- | --- |
| «cómo están los modelos, como se ajustan los parámetros» | Distinguir el ajuste a un blanco saturniano de la comparación con Voyager 2. |
| «añade el marco de referencia de la tierra» | Mostrar el primer cruce saliente de la esfera de influencia y registrar cuándo faltan datos de Voyager. |
| «que esté la orbita circular y se haga la transferencia tal como explican las diapositivas» | Construir el arco circular de estacionamiento y el impulso tangencial, con altura e impulso editables. |
| «necesito que el modelo de transferencia de hohmann y scattering stén bien» | Verificar tangencia y giro hiperbólico, y usar la velocidad de salida en la propagación posterior. |
| «por qué en saturno no se continua graficando la física?» | Distinguir la integración del recorte de la gráfica y mostrar la salida. |
| «auditar, limpiar y estructurar el directorio actual» | Elaborar el reporte UIS y verificar una compilación con fuentes y datos únicamente. |

## 3. Discernimiento

- La premisa de error histórico casi nulo no se sostuvo. Se compararon distancias al centro, épocas y posiciones a una misma fecha. Objetivos de ajuste similares explican cercanía entre algunos mínimos de los modelos, mientras que ambos pueden alejarse del valor histórico.
- La trayectoria ajustada por Lambert y disparo continuo se había presentado junto a Hohmann. A partir de la corrección del equipo se implementó la construcción nominal desde LEO, diferenciando los dos experimentos por sus ecuaciones y condiciones iniciales.
- La gráfica nominal de Saturno terminaba en periapsis aunque la integración seguía. La observación del equipo llevó a renderizar ambos lados del mínimo en el marco planetario, usando la solución ya integrada.
- Se trató el mínimo histórico obtenido con datos diarios como estimación de interpolación y se contrastó con el valor publicado por NASA. La discrepancia P–N se interpretó como residuo numérico al comprobar la identidad de fuerzas entre marcos.

## 4. Diligencia

| Evidencia | Acción verificable | Responsable registrado |
| --- | --- | --- |
| Geometría nominal | Pruebas de radios extremos, tangencias y tiempo de vuelo. | Asistente |
| Asistencia y salida de Saturno | Comprobar módulo asintótico, ganancia de energía y continuidad de trayectorias. | Asistente |
| Reproducibilidad | Verificar las seis sumas de entrada, regenerar figuras y PDF y ejecutar pruebas. | Asistente |
| Control de calidad del modelo | Señalar que faltaba el arco de Hohmann y que la animación se cortaba en Saturno. | Equipo, vía mensajes de la sesión |
| Revisión manual independiente | [Agregar cálculo, datos o líneas revisadas y evidencia, si se realizó.] | [Integrante y fecha] |

### Registro para sesiones anteriores o nuevas

| Fecha y persona | Delegación y criterio de éxito | Descripción o prompt | Discernimiento y decisión | Diligencia, corrección manual y evidencia |
| --- | --- | --- | --- | --- |
| [completar] | [completar] | [cita o enlace] | [aceptado, modificado o descartado y por qué] | [cálculo, prueba, archivo o commit] |

## Revisión de acuerdo con el enunciado

- **Delegación.** El equipo encargó a Codex contrastar el trabajo con el apartado general, el reto III.3 y los cinco entregables. El criterio fue que el reporte se pudiera leer en español y que cada conclusión tuviera evidencia.
- **Descripción.** «Los documentos no se entienden». La revisión se centró en explicar qué fuerza usa cada modelo, qué magnitud se conserva y cómo se mide el error de las cónicas.
- **Discernimiento.** Los mínimos de K y N ajustados por separado son cercanos, pero sus fechas difieren casi 16 días. Esa diferencia obligó a precisar en qué experimento se mide el error de la aproximación. Para esa medición se conservó el protocolo con un estado común en cada frontera joviana.
- **Diligencia.** Codex regeneró las tablas, verificó las seis sumas SHA-256 y ejecutó 23 pruebas. La comprobación independiente de las fórmulas de dispersión se puede repetir con `make comprobar`. La referencia histórica se contrastó con la tabla de sobrevuelo de NASA. El equipo puede registrar su propia revisión y la rotación de roles en la fila anterior.
