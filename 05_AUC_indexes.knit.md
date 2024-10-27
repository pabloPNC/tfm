# (PART) Evaluación de índices {-}

En esta sección del trabajo se estudiarán los distintos índices aplicados a los marcadores de la sección anterior. En primer lugar, se describirá el fundamento teórico y limitaciones de cada uno de ellos. Tras ello serán aplicados al conjunto de datos utilizando las distintas variables respuesta. Finalmente se evaluará la capacidad para seleccionar los clasificadores de mayor rendimiento, además de contrastar sus resultados con marcadores de utilidad aplicados actualmente en clínica.

# Evaluación de índices {#index-evaluation}



Para esta sección necesitaremos el conjunto de datos generado en la anterior sesión. El conjunto de datos contiene los niveles de diversos genes expresados en células luminales de próstata de pacientes con cancer.



## AUCs

### Introducción

En la literatura, 2 de las aproximaciones estadísticas más usadas para evaluar y comparar el rendimiento de una función o clasificador binario han sido las curvas $ROC$ (Característica Operativa del Receptor) y el área bajo dichas curvas ($AUC$).

El principal objetivo de un clasificador binario es discrimar instancias con una condición de interés de aquellas que no lo posean. Su precisión viene determianada por 2 métricas: su sensibilidad y su especificidad o la tasa verdaderos positivos ($TPR = Sensitivity$) y la tasa de falsos positivos ($FPR = 1 - Specificity$).

Cuando se trabaja con marcadores continuos, estos pueden ser dicotomizados en instancias positivas y negativas seleccionando una de las puntuaciones del marcador como un punto de corte, o umbral de decisión. Las curva $ROC$ representa los valores de $TPR$ y $FPR$ para cada uno de los umbrales de decisión de un marcador. Dado que estos puntos se calculan para todo un rango de umbrales, el $AUC$ actúa como una métrica global, mostrando los intercambios entre ambos parámetros de la curva, la ganancia o perdida de una métrica frente a otra.

### Cálculo de AUCs para diagnóstico

Comprendido el fundamento teórico es posible realizar los cálculos. Para ello se utilizará la función `pROC::auc`, la cual calculará el $AUC$ utilizando una variable para hacer la predicción (`predictor`) y una variable que contenga la respuesta a dicha predicción, es decir una variable respuesta o gold standard (`response`). La función requiere que la variable respuesta sea de tipo `factor` y utilice solamente 2 clases, donde la primera se corresponda con los casos sin y la segunda con los casos con la condición de estudio. 

En el conjunto de datos con el que se trabaja, la variable respuesta toma valores de `0` y `1` para indicar los casos con la condición de estudio y sin ella. Sin embargo, si se comprueba la ordenación de esta clases, se encuentran invertidas.


```r
levels(roc_data$disease)
```

```
## [1] "1" "0"
```

Por este motivo antes de realizar el cálculo de los $AUC$s se regorganizarán los valores de `disease`. Para esta tarea se utilizará la función `forcats::fct_relevel`, la cual tomará como entrada un `factor` y un vector indicando el nuevo orden de las clases de este.


```r
roc_data$disease <- fct_relevel(roc_data$disease, c("0", "1"))
levels(roc_data$disease)
```

```
## [1] "0" "1"
```

Ordenadas las clases en la variable respuesta, es posible aplicar la función `pROC::auc`. Dado que se buscará aplicar la métrica a todos los genes, se obtendrá un conjunto de datos temporal que contendrá solo los genes (`diagnostic_dataset`) y se utilizará la función `purrr::map_dbl` sobre este. La función iterará sobre las columnas (marcadores) y aplicará la función `pROC::auc` sobre cada una de ellas utilizando la variable `disease` como respuesta.


```r
diagnostic_dataset <- roc_data %>% select(
    -c("sample", "gleason_score", "malignancy", "disease", "prognostic")
)

auc_diagnostic <- map_dbl(
    diagnostic_dataset,
    auc,
    response = roc_data$disease,
    direction = "<",
    quiet = TRUE,
    .progress = TRUE
)
```

```
## ■■■■■■■ 20% | ETA: 4s ■■■■■■■■■ 26% | ETA: 4s ■■■■■■■■■■ 31% | ETA: 3s
## ■■■■■■■■■■■■ 38% | ETA: 3s ■■■■■■■■■■■■■■ 44% | ETA: 2s ■■■■■■■■■■■■■■■■■ 54% |
## ETA: 2s ■■■■■■■■■■■■■■■■■■■ 61% | ETA: 1s ■■■■■■■■■■■■■■■■■■■■■ 66% | ETA: 1s
## ■■■■■■■■■■■■■■■■■■■■■■ 71% | ETA: 1s ■■■■■■■■■■■■■■■■■■■■■■■■ 76% | ETA: 1s
## ■■■■■■■■■■■■■■■■■■■■■■■■■ 81% | ETA: 1s ■■■■■■■■■■■■■■■■■■■■■■■■■■■ 87% | ETA:
## 1s ■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 91% | ETA: 0s ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
## 97% | ETA: 0s
```

### Cálculo de AUCs para pronóstico

De forma similar a la sección anterior, es posible calcular las $AUC$s pero esta vez en un contexto de pronóstico para la enfermedad. Sin embargo, esta sección requerirá de un filtrado previo: hasta ahora se estaba trabajando tanto con tejidos sanos como enfermos, pero para este contexto solamente se necesarios casos con la enfermedad.

Para realizar este primer filtrado se generará un conjunto de datos temporal con los pacientes enfermos (`disease == 1`). A continuación, se reorganizarán las clases de la variable (`prognostic`) y eliminarán las que no son utilizadas (`Normal`, aquellas correspondientes a pacientes sanos).


```r
disease <- roc_data %>% filter(disease == 1)
disease[["prognostic"]] <- fct_drop(disease[["prognostic"]])
disease[["prognostic"]] <- fct_relevel(disease[["prognostic"]], c("0", "1"))
```

Realizado el filtrado, se realizará un procedimiento similar para calcular las $AUC$s en este nuevo contexto.


```r
prognostic_dataset <- disease %>% select(
    -c("sample", "gleason_score", "malignancy", "disease", "prognostic")
)

auc_prognostic <- map_dbl(
    prognostic_dataset,
    auc,
    response = disease$prognostic,
    direction = "<",
    quiet = TRUE,
    .progress = TRUE
)
```

```
## ■■■■■■■■■■■■■■ 43% | ETA: 1s ■■■■■■■■■■■■■■■■ 48% | ETA: 1s
## ■■■■■■■■■■■■■■■■■■■■■■ 71% | ETA: 1s ■■■■■■■■■■■■■■■■■■■■■■■■ 77% | ETA: 1s
## ■■■■■■■■■■■■■■■■■■■■■■■■■■ 82% | ETA: 1s ■■■■■■■■■■■■■■■■■■■■■■■■■■■ 88% | ETA:
## 0s ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 94% | ETA: 0s ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
## 99% | ETA: 0s
```

## pAUC

### Introducción

En ciertos contextos el estudio de la curva $ROC$ en su totalidad, no es de interés. En aplicaciones médicas y en otros campos, suelen seleccionarse ciertos rangos de la curva debido al interés por unos valores concretos de $TPR$ o $FPR$.

Por ejemplo en un contexto de diagnóstico, un marcador con un bajo $FPR$ (alta especifidad) es relevante a la hora de confirmar la presencia de una condición ("rule in"), dado que es muy improbable que se de un falso positivo. Mientras que un marcador con un alto $TPR$ (alta sensibilidad) es relevante para descartar la presencia de una condición ("rule out"), dado que es muy improbable que se de un falso negativo.

En este tipo de situaciiones el área parcial de la curva $ROC$ ($pAUC$) atrae mucha más atención como métrica, ya que sintetiza la información en la región de interés en lugar de en toda la curva.

### Cálculo de pAUCs en regiones de alta especificidad

Un diagnóstico temprano de ciertas enfermedades, como es el caso del cáncer de prostata, es gran importancia dado que puede mejorar el pronóstico de un paciente en gran medida [@early-detection]. Sin embargo, el tratamiento de estas enfermedades puede ser agresivo y por tanto debe considerarse en gran medida. En este contexto, es conveniente reducir casos en los que se trate a pacientes sanos, es decir el número de falsos positivos, por lo que se requieren condiciones de alta especificidad.

Con el fin de trabajar en estas condiciones se seleccionarán arbitrariamente rangos de $TPR$ contenidos en un intervalo de $(0.9, 1)$. Para realizar estos cálculos se usará de nuevo la función `pROC::auc` indicando que se trabajará sobre el rango indicado (`partial.auc = c(0.9, 1)`) y sobre la especificidad (`partial.auc.focus = "spec"`).


```r
pauc_sp_diagnostic <- map_dbl(
    diagnostic_dataset,
    auc,
    partial.auc = c(0.9, 1),
    partial.auc.focus = "spec",
    response = roc_data$disease,
    direction = "<",
    quiet = T,
    .progress = T
)
```

De forma similar, determinar el pronóstico de un paciente es relevante para seleccionar un tratamiento más o menos agresivo. Cometer falsos positivos en este contexto implica consecuencias graves para el paciente, por este motivo, también se considerarán condiciones de alta especificidad para determinar un pronóstico.


```r
pauc_sp_prognostic <- map_dbl(
    prognostic_dataset,
    auc,
    partial.auc = c(0.9, 1),
    partial.auc.focus = "spec",
    response = disease$prognostic,
    direction = "<",
    quiet = T,
    .progress = T
)
```

```
##  ■■■■■■■■■■■■■                     40% |  ETA:  2s
```

### Cálculo de pAUCs en regiones de alta sensibilidad

De forma similar, es posible traladar el caso a condiciones de alta sensibilidad. Una alta sensibilidad se requiere en enfermedades que son fatales si no son tratadas, es decir, cometer un falso negativo implica la muerte del paciente. El mal pronóstico de la enfermedad implica graves consecuencias, por lo que al realizarlo, se debe considerar una alta sensibilidad.

De este modo para calcular los índices, simplemente se realizará el procedimiento anterior pero esta vez indicando que se utilizará la sensibilidad (`partial.auc.focus = "sens"`).


```r
pauc_sn_prognostic <- map_dbl(
    prognostic_dataset,
    auc,
    partial.auc = c(0.9, 1),
    partial.auc.focus = "sens",
    response = disease$prognostic,
    direction = "<",
    quiet = T,
    .progress = T
)
```

Finalmente y aunque quizá menos interesante desde el punto de vista clínico, también es relevante un contexto de alta sensibilidad para el diagnóstico de la enfermedad. Como se ha explicado antes, un diagnóstico temprano puede ser relevante para evitar la muerte del paciente. Por este motivo, también debería considerarse esta métrica para realizar los cálculos.


```r
pauc_sn_diagnostic <- map_dbl(
    diagnostic_dataset,
    auc,
    partial.auc = c(0.9, 1),
    partial.auc.focus = "sens",
    response = roc_data$disease,
    direction = "<",
    quiet = T,
    .progress = T
)
```

```
##  ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1s
```

## Índices derivados de pAUC

### Problemas del pAUC

Aunque el $pAUC$ puede ser útil para trabajar en regiones de interés en la curva $ROC$, no está libre de fallas. Este índice ha sido cuestionado por su falta de interpretabilidad, ya que aunque algunos marcadores presentan valores altos de $AUC$ y rinden bien en las zonas seleccionadas, presentan valores de $pAUC$ muy cercanos a 0. Esto puede comprobarse con los 10 marcadores de mayor $AUC$ en condiciones de alta especificidad:


```r
enframe(auc_diagnostic) %>%
    arrange(desc(value)) %>%
    left_join(
        enframe(pauc_sp_diagnostic),
        join_by(name == name)
    ) %>%
    dplyr::rename(
        auc = value.x,
        pauc = value.y
    ) %>%
    .[0:10,]
```

```
## # A tibble: 10 × 3
##    name              auc   pauc
##    <chr>           <dbl>  <dbl>
##  1 ENSG00000105707 0.944 0.0746
##  2 ENSG00000181378 0.925 0.0749
##  3 ENSG00000124839 0.923 0.0680
##  4 ENSG00000006625 0.921 0.0553
##  5 ENSG00000242110 0.916 0.0672
##  6 ENSG00000159263 0.910 0.0528
##  7 ENSG00000138028 0.909 0.0505
##  8 ENSG00000173638 0.907 0.0513
##  9 ENSG00000245750 0.906 0.0549
## 10 ENSG00000176058 0.904 0.0605
```

Como se observa, los $pAUC$ no superan el $0.5$; valor límite en el $AUC$ con el que se consideraría que un clasificador rinde igual que uno completamente aleatorio.  Además, en ciertos escenarios, el $pAUC$ no puede distinguir que clasificador presentaría un mejor rendimiento. Por ejemplo, cuando se comparan clasificadores con $pAUC$s iguales y curvas $ROC$ que se cruzan [@pauc-comparison].

Para lidiar con estos defectos, distintas transformaciones del $pAUC$ se han desarrollado para trabajar tanto en condiciones de alta sensibilidad como de alta especificidad, las cuales serán tratadas en las siguientes secciones.

### Índices para la especificidad

#### SpAUC

Una de las primeras transformaciones descritas es el **área parcial bajo la curva estandarizada** ($SpAUC$), propuesto por @spauc. Esta transformación trabaja sobre un rango concreto de especificidad, $(FPR_1, FPR_2)$, donde los valores que puede tomar el $pAUC$ se escalan a valores dentro del intervalo $(0.5, 1)$ haciendo que la métrica sea más interpretable.

<!-- TODO -->

Esta transformación trabaja con sobre un rango concreto de especificidad $(FPR_1, FPR_2)$, donde los límites superior e inferior para el $pAUC$ se escalan a los valores en el intervalo $(0.5, 1)$, haciendo que la métrica sea más interpretable. A pesar de ello, la métrica todavía presenta algunas desventajas.

El límite superior para el área que engloba el $pAUC$, esta definida por el rectángulo formado por los puntos: $(FPR_1, 0), (FPR_2, 0), (FPR_1, 1)$ y $(FPR_2, 1)$. Mientras que el límite inferior está definido por el trapecio con puntos: $(FPR_1, 0), (FPR_2, 0), (FPR_1, FPR_1)$ y $(FPR_2, FPR_2)$. Este será el trapecio que coincide con la diagonal de la curva. Si calculamos las áreas de los límites, el área $(A)$ que puede tomar el $pAUC$ está definida de la siguiente forma:

$$ 
\dfrac{1}{2}(FPR_2 - FPR_1)(FPR_2 + FPR_1) \le A \le FPR_2 - FPR_1
$$ 

Realizando la siguiente tranformación:

$$
\dfrac{1}{2}\left[1 + \dfrac{A - \text{min}}{\text{max} - \text{min}} \right]
$$ El índice toma un valor de $1$ cuando $A = \text{max}$ y un valor de $0.5$ cuando $A = \text{min}$, permitiendo ver el $pAUC$ de forma escalada al área total.

El problema que presenta la métrica, es que asume el mismo límite inferior para todas las curvas. En la práctica es común trabajar con curvas que atraviesen la diagonal, por lo que el límite inferior no está bien definido. Además, aunque resuelve el problema de facilitar la interpretación, no solventa el problema de 2 curvas con igual $pAUC$ que se cruzan.

Conociendo sus ventajas y desventajas, calcularemos el índice para contrastarlo con otros métodos. Para ello utilizaremos también la función `pROC::auc` indicando el argumento `partial.auc.correct = TRUE` para que aplique la corrección descrita por McClish. Además indicaremos `allow.invalid.partial.auc.correct = FALSE`, de esta forma la función no calculará el índice (aunque sea posible) si la curva está descrita por debajo del area mínima, puesto que el límite inferior no está bien definido, en su lugar devolverá un `NA`.


```r
spauc_sp_diagnostic <- map_dbl(
    diagnostic_dataset,
    auc,
    partial.auc = c(0.9, 1),
    partial.auc.focus = "spec",
    partial.auc.correct = TRUE,
    allow.invalid.partial.auc.correct = FALSE,
    response = roc_data$disease,
    direction = "<",
    quiet = TRUE,
    .progress = TRUE
)
```


```r
spauc_sp_prognostic <- map_dbl(
    prognostic_dataset,
    auc,
    partial.auc = c(0.9, 1),
    partial.auc.focus = "spec",
    partial.auc.correct = TRUE,
    allow.invalid.partial.auc.correct = FALSE,
    response = disease$prognostic,
    direction = "<",
    quiet = TRUE,
    .progress = TRUE
)
```

```
##  ■■■■■■■■■■■■                      37% |  ETA:  2s
```

#### TpAUC

Para solventar los problemas que el $SpAUC$ todavía presenta, [Vivo et al. han propuesto](https://link.springer.com/article/10.1007/s11634-017-0295-9) un índice alternativo para rangos especificidades concretos, el "tighter partial area index" $TpAUC$. Este nuevo índice, encuentra nuevos límites todavía más ajustados y teniendo en cuenta la forma que puede adoptar la curva $ROC$, es decir, siendo compatible con curvas que pasan por debajo de la diagonal o son inadecuadas.

Los límites iniciales propuestos por este índice parten de los establecidos por el $SpAUC$:

$$
\dfrac{1}{2}(FPR_2 - FPR_1)(FPR_2 + FPR_1) \le A \le FPR_2 - FPR_1
$$

Sin embargo, Vivo et al. proponen que estos límites se pueden ajustar a las regiones del area parcial donde esté definida la curva. Es decir, ajustandose a los puntos de la curva presentes en la región de interes, donde $TPR_i = TPR(FPR_i)$ para $i = 1,2$. De esta forma los límites se ajustan a:

$$
TPR_1(FPR_2 - FPR_1) \le A \le TPR_2(FPR_2-FPR_1)
$$

Este índice inferior se ajusta incluso cuando las curvas son inadecuadas en el límite inferior. De esta forma, cuando las curvas son adecuadas y no presentan ganchos, la siguiente expresión puede utilizarse:

$$
\text{max}\left\{TPR_1(FPR_2 - FPR_1), \dfrac{1}{2}(FPR_2 - FPR_1)(FPR_2 + FPR_1)\right\} \le A \le TPR_2(FPR_2-FPR_1)
$$

Finalmente, se debe considerar que las curvas $ROC$ normalmente son concavas, salvando la presencia de irregularidades (e.g. ganchos). Una curva $ROC$ concava, presenta un ratio de probabilidad positivo ($PLR(t) = ROC(t)/t$) decreciente. De esta forma el $PLR$ alcaza valores mínimos en el extremo superior de la curva.

Trasladado a este caso, si se puede asegurar que la curva se mantiene concava en la región de interés, es decir $PLR$ alcanza su mínimo valor en $FPR_2$ un límite inferior más ajustado se puede determinar, correspondiéndose al trapecio con puntos: $(FPR_1,0), (FPR_2, 0), (FPR_1, TPR_2)$ y $(FPR_2, TPR_2)$. Desarrollando esta área obtenemos:

$$
\dfrac{1}{2}(TPR_1+TPR_2)(FPR_2-FPR_1) \le A \le TPR_2(FPR_2-FPR_1)
$$

Conociendo los límites inferior y superior, y usándolos según las condiciones descritas, aplicando la transformación descrita por McClish para encapsular los valores entre $(0.5,1)$:

$$
\dfrac{1}{2}\left[1 + \dfrac{A - \text{min}}{\text{max} - \text{min}} \right]
$$

Se obtiene finalmente el $TpAUC$.

De forma similar, procederemos a calcularlo para contrastarlo con las demás métricas. El algoritmo para el cálculo de este índice está implementado en el paquete `ROCpAI` de [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/ROCpAI.html) con la función `tpAUC`. El paquete además incluye funciones para el cálculo de sus estabilidades, aunque para este caso nos bastará con la función principal.

La función utiliza como entrada un dataset donde la primera columna se corresponda con la respuesta y las siguientes los predictores.


```r
tpauc_sp_diagnostic_assay <- tpAUC(
    dataset = diagnostic_dataset %>% mutate(
        disease = roc_data$disease,
        .before = 1
    ),
    low.value = 0,
    up.value = 0.1,
    plot = FALSE
)
```


```r
tpauc_sp_prognostic_assay <- tpAUC(
    dataset = prognostic_dataset %>% mutate(
        prognostic = disease$prognostic,
        .before = 1
    ),
    low.value = 0,
    up.value = 0.1,
    plot = FALSE
)
```

Dado que la salida de la función es un objeto con diversas métricas, utilizaremos la función `assay` para extraer los $SpAUC$ almacenados en en `St_pAUC`. Finalmente nombraremos cada uno de los valores con su respectivo identificador.


















































































































