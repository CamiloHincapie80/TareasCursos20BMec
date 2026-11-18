# Formulación corregida

## Sistema y escalas

Sonda sin masa dinámica: seis coordenadas de estado (tres grados de libertad), cartesianas eclípticas J2000, km y s, tiempo TDB. No hay vínculos holónomos adicionales. N prescribe los centros de Júpiter y Saturno; no conserva el momento ni la energía de un sistema aislado porque no integra sus reacciones. La Tierra se reserva para el diseño analítico. Sol y planetas son masas puntuales; no se incluyen satélites, achatamiento, presión de radiación ni maniobras. Las distancias de periapsis se refieren al centro, no a las nubes.

## Newton y origen acelerado

Si x es una posición inercial, r=x−x_sol y R_i=x_i−x_sol, la aceleración solar debida a los planetas incluidos es Σ μ_i R_i/R_i³. Su resta produce

`a_N(t,r) = −μ_sol r/r³ + Σ_i μ_i[−(r−R_i)/|r−R_i|³ − R_i/R_i³]`.

El primer sumando del corchete es directo; el segundo es la aceleración opuesta del origen. El signo indirecto negativo del código original es correcto. Con efemérides reales la suma truncada no es la aceleración completa del Sol: N es una aproximación gravitacional heliocéntrica declarada, no una referencia baricéntrica exacta. Sin efeméride solar baricéntrica no se añade una corrección supuestamente exacta.

Para cualquier planeta prescrito k, la identidad cinemática obliga a

`rho=r−R_k; u=v−R'_k; rho''=a_N(t,R_k+rho)−R''_k`.

Esta es P corregida, idéntica físicamente a N. R'_k y R''_k son las derivadas del mismo Hermite de posición. Su segunda derivada tiene saltos en los nudos: se detiene y reinicia la integración P en cada nudo mientras el marco es planetocéntrico, sin cambiar estado ni fuerzas físicas.

Una segunda derivación separa las fuerzas:

`rho'' = −μ_k rho/rho³ + [g_sol(R_k+rho)−g_sol(R_k)]`

`+ Σ_{i≠k}[g_i(R_k+rho−R_i)−g_i(R_k−R_i)] + c_k(t)`

`c_k = g_{sol+k}(R_k) + Σ_{i≠k}[g_i(R_k−R_i)+g_i(R_i)] − R''_k`,

con `g_i(z)=−μ_i z/|z|³`. El término c_k mide la diferencia entre la aceleración del planeta prescrita y la del sistema simplificado, incluyendo la reacción del Sol y cuerpos omitidos. No puede identificarse únicamente con las lunas sin un experimento adicional. P0 conserva solo la atracción central y el corchete solar. Incluso en dos cuerpos la aceleración heliocéntrica del planeta contiene μ_sol+μ_k, no solamente μ_sol.

## Lagrange y Hamilton

Por unidad de masa, sea

`U(r,t)=μ_sol/r + Σ_i μ_i[1/|r−R_i|−r·R_i/R_i³]`.

Entonces `L=|v|²/2+U`, `p=v`, `H=|p|²/2−U`. Euler–Lagrange da `r''=∇U`, y Hamilton da `r'=∂H/∂p`, `p'=−∂H/∂r`. En general `dH/dt=−∂U/∂t≠0`.

Al trasladar r=R_k+rho, el término cruzado R'_k·rho' se elimina restando la derivada total de R'_k·rho. Descartando funciones solo temporales,

`L_k=|u|²/2+U(R_k+rho,t)−rho·R''_k`;

`H_k=|p_k|²/2−U(R_k+rho,t)+rho·R''_k`.

Su gradiente reproduce P exactamente. Esta construcción demuestra qué términos faltan en el lagrangiano planetocéntrico original.

## Conservaciones y transferencia de energía

En dos cuerpos, `ε=|v|²/2−μ/r` y **el vector** `h=r×v` son constantes; también lo es el vector de excentricidad. En K se conservan por región, no a través de todas las fronteras. En P/N, `dε_J/dt=u·a_pert` y `dh_J/dt=rho×a_pert`. La energía kepleriana heliocéntrica satisface `dε_sol/dt=v·(a_N−g_sol)`; no es H, que incluye potenciales planetarios e indirectos.

En dispersión aislada hiperbólica, `ε_J=v_inf²/2>0`, `e=1+rp v_inf²/μ_J`, `b²=rp²+2 μ_J rp/v_inf²`, `δ=2 asin(1/e)`. A distancia finita la rapidez es `sqrt(v_inf²+2μ_J/rho)`, no v_inf. El ángulo entre velocidades finitas no es automáticamente δ. En una trayectoria perturbada δ calculada en periapsis es un diagnóstico **osculante**, no la deflexión acumulada real.

En el encuentro instantáneo con posición y V_J iguales a ambos lados,

`Δε_sol = V_J·(u_out−u_in) + (|u_out|²−|u_in|²)/2`.

El último término desaparece solo si se conserva la rapidez asintótica. Por ello una sonda puede ganar energía solar sin ganar energía asintótica planetaria. No se aplica la identidad instantánea a dos fechas separadas por 132 días con posiciones y V_J distintas.

## Frontera y predicciones anteriores al cálculo

Laplace usa `a(μ_p/μ_sol)^(2/5)`, asociado a comparación de perturbaciones en órbitas casi keplerianas, no igualdad puntual de fuerzas ni una discontinuidad física. Se conserva como variante del original el radio **instantáneo** R(t) en lugar del semieje; está rotulado así en las tablas. El factor f cambia las fuerzas descartadas de K y el tiempo de acumulación; no hay razón para que el error sea monótono. P no debería depender físicamente de f.

Antes de integrar se predice: Hohmann llega con exceso retrógrado relativo a Júpiter; una rotación hacia componente radial exterior puede aumentar energía solar y alcanzar Saturno; la fase es necesaria además de la energía. K debe conservar ε y h regionales. P debe converger a N. La marea solar disminuye aproximadamente como rho³ respecto a la atracción joviana al acercarse al planeta, pero el residuo c_k tiene otra dependencia. A mayor b o rp, δ tiende a cero; en el límite matemático frontal δ tiende a π, aunque el planeta impone un límite físico antes de rho=0.

## Separación de errores

* P−N: integración, interpolación densa y representación numérica del cambio de marco; misma aceleración física y mismas CI.
* P0−P: omisiones de fuerzas y movimiento del origen, además del residuo numérico.
* K−N: fuerzas eliminadas dentro **y fuera** de las SOI, frontera elegida y propagación posterior. Las tablas hasta Saturno no aíslan únicamente la marea dentro de Júpiter.
* Datos horarios excluidos frente a Hermite diario: interpolación de Voyager; no error de N.
* Cambio de paso y método: integración a efeméride fija. Planetas cada dos días: sensibilidad a resolución de efemérides, no convergencia del integrador.
* N sin Saturno: sensibilidad a cuerpo omitido. N contra Voyager: mezcla de física incompleta, estado inicial, trayectoria de diseño aproximada y datos; no medición pura del error K.

Se declara equivalencia **para la distancia mínima a Saturno** si `|rp_modelo−rp_N|<1000 km` en el intervalo 1981-01-01 a 1981-10-01. No implica equivalencia del estado completo. No es una precisión histórica certificada. La elección permite resolver diferencias grandes de K/P0 manteniendo un margen sobre el residuo numérico P−N.
