# (PART) Evaluación de índices {-}

# Evaluación de índices {#index-evaluation}



En esta sección del trabajo se estudiarán los distintos índices aplicados a los marcadores de la sección anterior. En primer lugar, se describirá el fundamento teórico y limitaciones de cada uno de ellos. Tras ello, estos índices serán aplicados al conjunto de datos utilizando las distintas variables respuesta. Finalmente, se evaluará su capacidad para seleccionar los clasificadores de mayor rendimiento, además de contrastar sus resultados con marcadores de utilidad aplicados actualmente en la práctica.

Para esta sección necesitaremos el conjunto de datos generado en la anterior sesión. El conjunto de datos contiene los niveles de diversos genes expresados en células luminales de próstata de pacientes con cáncer.



## AUCs

### Introducción

En la literatura, 2 de las aproximaciones estadísticas más usadas para evaluar y comparar el rendimiento de una función o clasificador binario han sido las curvas $ROC$ (Característica Operativa del Receptor) y el área bajo dichas curvas ($AUC$).

El principal objetivo de un clasificador binario es discriminar instancias con una condición de interés. Su precisión viene determinada por 2 métricas: su sensibilidad y su especificidad o, de forma equivalente la tasa verdaderos positivos ($TPR = Sensitivity$) y de falsos positivos ($FPR = 1 - Specificity$).

Cuando se trabaja con marcadores continuos, estos pueden ser dicotomizados en instancias positivas y negativas seleccionando una de las puntuaciones del marcador como un punto de corte, o umbral de decisión. Una curva $ROC$ representa los valores de $TPR$ y $FPR$ para cada uno de los umbrales de decisión de un marcador. Dado que estos puntos se calculan para todo un rango de umbrales, el $AUC$ actúa como una métrica global mostrando los intercambios entre ambos parámetros de la curva, es decir, mostrando la ganancia o perdida de una métrica frente a otra.

### Cálculo de AUCs para diagnóstico

Comprendido el fundamento teórico es posible realizar los cálculos. Para ello, se utilizará la función `pROC::auc`, la cual calculará el $AUC$ utilizando una variable para hacer una predicción (`predictor`) y una variable que contenga la respuesta a dicha predicción, es decir la condición de interés, variable respuesta o gold standard (`response`). La función requiere que la variable respuesta sea de tipo `factor` y utilice solamente 2 valores, donde el primero se corresponda con los casos sin la condición de estudio y el segundo con ella. 

En el conjunto de datos con el que se trabaja, la variable respuesta toma valores de `0` y `1` para indicar los casos con la condición de estudio y sin ella. Sin embargo, si se comprueba el orden de estos valores, se observa que está invertido.


```r
levels(roc_data$disease)
```

```
## [1] "1" "0"
```

Por este motivo, antes de realizar el cálculo de los $AUC$s, se reorganizarán los valores de `disease`. Para esta tarea se utilizará la función `forcats::fct_relevel`, la cual tomará como entrada un `factor` y un vector indicando el nuevo orden para sus clases.


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

### Cálculo de AUCs para pronóstico

De forma similar a la sección anterior, es posible calcular las $AUC$s pero esta vez en un contexto de pronóstico para la enfermedad. Sin embargo, esta sección requerirá de un filtrado previo: hasta ahora se estaba trabajando tanto con tejidos sanos como con enfermos, sin embargo, para este contexto solamente serán necesarios casos con la enfermedad.

Para realizar este primer filtrado, se generará un conjunto de datos temporal con los pacientes enfermos (`disease == 1`). A continuación, se reorganizarán las clases de la variable (`prognostic`) y eliminarán las que no son utilizadas (`Normal`, aquellas correspondientes a pacientes sanos).


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

## pAUC

### Introducción

En ciertos contextos, el estudio de la curva $ROC$ en su totalidad no es de interés. En aplicaciones médicas y en otros campos, suelen seleccionarse ciertos rangos de la curva debido al interés por unos valores concretos de $TPR$ o $FPR$.

Por ejemplo, en un contexto de diagnóstico, un marcador con un bajo $FPR$ (alta especificidad) es relevante a la hora de confirmar la presencia de una condición ("rule in"), dado que es muy improbable que se de un falso positivo. Mientras que un marcador con un alto $TPR$ (alta sensibilidad) es relevante para descartar la presencia de una condición ("rule out"), dado que es muy improbable que se de un falso negativo.

En este tipo de situaciones, el área parcial de la curva $ROC$ ($pAUC$) atrae mucha más atención como métrica, ya que sintetiza la información en la región de interés en lugar de en toda la curva.

### Cálculo de pAUCs en regiones de alta especificidad

Un diagnóstico temprano de ciertas enfermedades, como es el caso del cáncer de próstata, es de gran importancia dado que puede mejorar el pronóstico de un paciente en gran medida [@early-detection]. Sin embargo, el tratamiento de estas enfermedades puede ser agresivo y por tanto debe considerarse cuidadosamente. En este contexto, es conveniente reducir los casos en los que se trata a pacientes sanos, es decir el número de falsos positivos, por ello se requieren condiciones de alta especificidad.

Con el fin de trabajar en estas condiciones se seleccionarán arbitrariamente rangos de $1 - FPR$ contenidos en un intervalo de $(0.9, 1)$. Para realizar estos cálculos se usará de nuevo la función `pROC::auc` indicando que se trabajará sobre el rango indicado (`partial.auc = c(0.9, 1)`) y sobre la especificidad (`partial.auc.focus = "spec"`).


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

### Cálculo de pAUCs en regiones de alta sensibilidad

De forma similar, es posible trasladar el caso a condiciones de alta sensibilidad. Una alta sensibilidad se requiere en enfermedades que son fatales si no son tratadas, es decir, cometer un falso negativo implica la muerte del paciente. De esta forma, cuando un mal pronóstico para la enfermedad implica graves consecuencias, se debe considerar una alta sensibilidad.

De este modo, para calcular los índices, simplemente se realizará el procedimiento anterior pero esta vez indicando que se utilizará la sensibilidad (`partial.auc.focus = "sens"`).


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

Finalmente, desde el punto de vista clínico también es relevante un contexto de alta sensibilidad para el diagnóstico de la enfermedad. Como se ha explicado antes, un diagnóstico temprano puede ser relevante para evitar la muerte del paciente. Por este motivo, también debería considerarse esta métrica para realizar los cálculos.


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

## Índices derivados de pAUC

### Problemas del pAUC

Aunque el $pAUC$ puede ser útil para trabajar en regiones de interés en la curva $ROC$, no está libre de fallas. Este índice ha sido cuestionado por su falta de interpretabilidad, ya que aunque algunos marcadores presentan valores altos de $AUC$ y rinden bien en las zonas seleccionadas, presentan valores de $pAUC$ muy cercanos a 0. Esto puede comprobarse, por ejemplo, con los 10 marcadores de mayor $AUC$ en condiciones de alta especificidad:


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
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> name </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000105707 </td>
   <td style="text-align:right;"> 0.9436868 </td>
   <td style="text-align:right;"> 0.0746093 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000181378 </td>
   <td style="text-align:right;"> 0.9252605 </td>
   <td style="text-align:right;"> 0.0749387 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000124839 </td>
   <td style="text-align:right;"> 0.9229620 </td>
   <td style="text-align:right;"> 0.0680432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000006625 </td>
   <td style="text-align:right;"> 0.9207784 </td>
   <td style="text-align:right;"> 0.0553402 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000242110 </td>
   <td style="text-align:right;"> 0.9157217 </td>
   <td style="text-align:right;"> 0.0672311 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000159263 </td>
   <td style="text-align:right;"> 0.9103203 </td>
   <td style="text-align:right;"> 0.0528042 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000138028 </td>
   <td style="text-align:right;"> 0.9091710 </td>
   <td style="text-align:right;"> 0.0505440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000173638 </td>
   <td style="text-align:right;"> 0.9071407 </td>
   <td style="text-align:right;"> 0.0513178 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000245750 </td>
   <td style="text-align:right;"> 0.9061063 </td>
   <td style="text-align:right;"> 0.0549264 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000176058 </td>
   <td style="text-align:right;"> 0.9044591 </td>
   <td style="text-align:right;"> 0.0604658 </td>
  </tr>
</tbody>
</table>


Como se observa, los $pAUC$ no superan el $0.5$; valor límite en el $AUC$ con el que se considera que un clasificador rinde igual que uno completamente aleatorio. Además, en ciertos escenarios, el $pAUC$ no puede distinguir que clasificador presentaría un mejor rendimiento. Por ejemplo, cuando se comparan clasificadores con $pAUC$s iguales y curvas $ROC$ que se cruzan [@pauc-comparison].

Para lidiar con estas limitaciones, distintas transformaciones del $pAUC$ se han desarrollado para trabajar en condiciones de alta sensibilidad y especificidad.

### Índices para la especificidad

#### SpAUC {#spauc}

Una de las primeras transformaciones descritas es el **área parcial bajo la curva estandarizada** ($SpAUC$), propuesto por @spauc. Esta transformación trabaja sobre un rango concreto de especificidad, $(FPR_1, FPR_2)$, donde los valores que puede tomar el $pAUC$ se escalan a valores dentro del intervalo $(0.5, 1)$ haciendo que la métrica sea más interpretable.

El límite superior para el área del $pAUC$ está definido por el rectángulo formado por los puntos: $(FPR_1, 0), (FPR_2, 0), (FPR_1, 1)$ y $(FPR_2, 1)$. Mientras que, el límite inferior está definido por el trapecio con puntos: $(FPR_1, 0), (FPR_2, 0), (FPR_1, FPR_1)$ y $(FPR_2, FPR_2)$, el cual coincide con la diagonal de la curva. 

Si se desarrollan las áreas de los límites superior e inferior, el área que puede tomar $pAUC$ ($A$) está definida de la siguiente forma:

\begin{equation}
    \dfrac{1}{2}(FPR_2 - FPR_1)(FPR_2 + FPR_1) \le A \le FPR_2 - FPR_1
    (\#eq:spauc-bounds)
\end{equation}

Finalmente, aplicando la siguiente transformación \@ref(eq:mcclish-transformation), se logra que el índice tome un valor de $1$ cuando $A = \text{max}$ y un valor de $0.5$ cuando $A = \text{min}$, escalando el $pAUC$ al área total.

\begin{equation}
    \dfrac{1}{2}\left[1 + \dfrac{A - \text{min}}{\text{max} - \text{min}} \right]
    (\#eq:mcclish-transformation)
\end{equation}

A pesar de la mejora en su interpretabilidad, la métrica todavía presenta algunas limitaciones. En la práctica, las curvas $ROC$ suelen atravesar la diagonal, por lo que, el límite inferior definido en \@ref(eq:spauc-bounds) no estará bien definido en estos casos. Además, el índice no solventa el problema de comparar 2 clasificadores con igual $pAUC$ y curvas $ROC$ que se cruzan.

Conociendo las ventajas y limitaciones de la métrica, esta se aplicará al conjunto de datos para contrastarla con otras más adelante. Para ello se utilizará la función `pROC::auc` indicando el argumento `partial.auc.correct = TRUE`, que aplicará la transformación descrita por McClish. Además, se indicará el parámetro `allow.invalid.partial.auc.correct = FALSE`, de esta forma, la función no calculará el índice cuando la curva se encuentra descrita por debajo del área mínima, puesto que el límite inferior no estará bien definido, devolviendo en su lugar un `NA`.


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

#### TpAUC {#tpauc}

Para solventar los problemas que presenta el $SpAUC$, @tpauc han propuesto una transformación alternativa para regiones de alta especificidad, el denominado *tighter partial area index* ($TpAUC$). Este nuevo índice define límites más ajustados que el anterior y considera la forma que adopta una curva $ROC$, es decir, es compatible con curvas definidas por debajo de la diagonal o inadecuadas.

Los límites iniciales propuestos por el $TpAUC$ parten de los establecidos por el $SpAUC$, \@ref(eq:spauc-bounds). Sin embargo, Vivo et al. proponen que estos límites pueden ajustarse a las regiones donde la curva se encuentre definida. Es decir, estos límites se ajustan a los puntos de la curva presentes en la región de interés, donde $TPR_i = TPR(FPR_i)$ para $i = 1,2$. Estos límites ajustados se definen de la siguiente forma:

\begin{equation}
    TPR_1(FPR_2 - FPR_1) \le A \le TPR_2(FPR_2-FPR_1)
    (\#eq:tpauc-bounds)
\end{equation}

El nuevo límite inferior es más consistente, ajustándose a la forma de la curva incluso cuando estas no son adecuadas, es decir, están definidas por debajo de la línea del azar. Sin embargo, cuando las curvas $ROC$ son adecuadas, se pueden describir límites más ajustados. Por definición, estas curvas son siempre no decrecientes, por lo que los siguientes límites podrían ser descritos:

\begin{equation}
    \text{max}\left\{TPR_1(FPR_2 - FPR_1), \dfrac{1}{2}(FPR_2 - FPR_1)(FPR_2 + FPR_1)\right\} \le \\
    A \le TPR_2(FPR_2-FPR_1)
    (\#eq:tpauc-proper-bounds)
\end{equation}

Finalmente, se pueden encontrar límites más ajustados cuando las curvas son cóncavas. Una curva $ROC$ cóncava presenta un ratio de probabilidad positivo ($PLR(t) = ROC(t)/t$) decreciente, es decir, este alcanza valores mínimos en el extremo superior de la curva. Así, si esta condición se cumple en la región de interés, se puede describir un nuevo límite inferior. Este límite más ajustado se corresponde con el trapecio de puntos: $(FPR_1,0), (FPR_2, 0), (FPR_1, TPR_2)$ y $(FPR_2, TPR_2)$, que considerando su área se obtiene la siguiente expresión:

\begin{equation}
    \dfrac{1}{2}(TPR_1+TPR_2)(FPR_2-FPR_1) \le A \le TPR_2(FPR_2-FPR_1)
    (\#eq:tpauc-concave-bounds)
\end{equation}

Descritos los límites y condiciones donde se aplican, es posible utilizar la transformación descrita por McClish \@ref(eq:mcclish-transformation) para encapsular los valores en el rango $(0.5, 1)$. Tras aplicar la transformación el resultado obtenido se corresponde con el $TpAUC$.

De forma similar al $SpAUC$, el $TpAUC$ se aplicará en el conjunto de datos para más tarde contrastarlo con las demás métricas. El algoritmo para calcular el $TpAUC$ se encuentra implementado en el paquete `ROCpAI` de [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/ROCpAI.html), concretamente, en la función `ROCpAI::tpAUC`. 

La función `tpAUC` requiere de indicar la región de la curva en la que aplicar los cálculos (`low.value` y `up.value`) y el conjunto de datos. Este último requerirá que la primera columna se corresponda con la variable respuesta, mientras que el resto con los demás predictores.


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

La salida de la función devolverá un objeto con diversas métricas. Utilizando la función `assay` se podrán extraer los $TpAUC$ almacenados en `St_pAUC`. Finalmente, los identificadores de cada clasificador se asignarán a cada uno de los resultados.


```r
tpauc_sp_diagnostic <- as.numeric(assay(tpauc_sp_diagnostic_assay)$St_pAUC)
tpauc_sp_prognostic <- as.numeric(assay(tpauc_sp_prognostic_assay)$St_pAUC)
names(tpauc_sp_diagnostic) <- colnames(diagnostic_dataset)
names(tpauc_sp_prognostic) <- colnames(prognostic_dataset)
```

### Índices para la sensibilidad

#### NpAUC {#npauc}

El índice *normalized partial area* ($NpAUC$), propuesto por @npauc, fue uno de los primeros enfocados en regiones de alta sensibilidad. Este índice es conceptualmente similar al $SpAUC$, salvo que en este caso son $TPR_1 = TPR_0$ y $TPR_2 = 1$ los que definen sus límites superior e inferior. Concretamente, estos límites están definidos como:

+ El rectángulo de puntos $(1,1), (1, TPR_0), (0, 1)$ y $(0, TPR_0)$, es decir, de longitud $1$ y altura dado por la banda $1 - TPR_0$.
+ El triángulo de puntos $(TPR_0, TPR_0), (1, TPR_0)$ y $(1,1)$.

Estos límites, de forma desarrollada, quedan expresados de la siguiente forma:

\begin{equation}
    \dfrac{1}{2}(1 - TPR_0)^2 \le A \le 1 - TPR_0
    (\#eq:npauc-bounds)
\end{equation}

Para realizar la normalización, cada término se divide entre el límite superior ($1 - TPR_0$) generando la expresión: 

\begin{equation}
    \dfrac{1}{2}(1-TPR_0) \le \dfrac{A}{1 - TPR_0} \le 1
    (\#eq:npauc-normalized-bounds)
\end{equation}

Así, el $NpAUC$ se corresponde con el término central, quedando definido como:

\begin{equation}
        NpAUC = \dfrac{A}{1 - TPR_0}
    (\#eq:npauc-index)
\end{equation}

Debido a esta transformación, el $NpAUC$ se puede interpretar como el valor promedio de la especificidad para los valores de sensibilidad por encima de $TPR_0$. Sin embargo, a pesar de ser más interpretable el índice sigue presentando desventajas. En primer lugar, el $NpAUC$ puede tomar valores inferiores a $0.5$, dificultando su interpretación. En segundo lugar, no es capaz de distinguir clasificadores con igual $pAUC$ y curvas $ROC$ que se cruzan, independientemente de su forma.

Igual que en los casos anteriores, el $NpAUC$ se aplicará al conjunto de datos para contrastarlo más adelante. El índice no está implementado en ningún paquete conocido, de modo que para aplicarlo de forma sencilla se ha incluido en un paquete propio, `ROCnGO`. Este paquete cuenta con la función `ROCnGO::np_auc` que permite calcular el índice con un procedimiento similar a los anteriores.


```r
npauc_sn_diagnostic <- map_dbl(
    diagnostic_dataset,
    \(gene) {
        np_auc(
            predictor = gene,
            response = roc_data$disease,
            lower_tpr = 0.9
        )
    },
    .progress = TRUE
)
```


```r
npauc_sn_prognostic <- map_dbl(
    prognostic_dataset,
    \(gene) {
        np_auc(
            predictor = gene,
            response = disease$prognostic,
            lower_tpr = 0.9
        )
    },
    .progress = TRUE
)
```

#### FpAUC {#fpauc}

De forma similar al $TpAUC$ con respecto al $SpAUC$, @fpauc propusieron el *fitted partial area index* ($FpAUC$) para solventar las desventajas del $NpAUC$ en rangos de alta sensibilidad. Este nuevo índice buscaba:

1. Ser interpretable como métrica de diagnóstico. Esto incluye ser igual a $AUC$ cuando $TPR_0 = 0$.
2. Aplicarse a cualquier forma de la curva $ROC$, incluyendo las inadecuadas.
3. Distinguir dos clasificadores con $pAUC$ iguales y que sus curvas $ROC$ se crucen.

Para realizar el cálculo del $FpAUC$ se requiere definir nuevos límites para el $pAUC$ en condiciones de alta sensibilidad, para ello, se partirá de los límites iniciales definidios en $NpAUC$ \@ref(eq:npauc-bounds). En primer lugar, un nuevo límite superior puede ser definido. Este nuevo límite puede ajustarse a la zona donde se encuentre definida la curva, es decir, al rectángulo de puntos $(FPR_0, TPR_0), (FPR_0, 1), (1, TPR_0)$ y $(1,1)$ donde $FPR_0 = FPR(TPR_0)$. Desarrollando el área de este rectángulo se obtiene la siguiente expresión \@ref(eq:tpauc-upper-bound).

\begin{equation}
    \dfrac{1}{2}(1 - TPR_0)^2 \le A \le (1-FPR_0)(1 - TPR_0)
    (\#eq:tpauc-upper-bound)
\end{equation}

Como se ha explicado anteriormente en la sección \@ref(tpauc), las curvas $ROC$ son funciones monótonas no decrecientes. Esto no asegura su concavidad, ya que en ocasiones pueden pasar por debajo de la línea del azar y mostrar rendimientos peores al de un predictor aleatorio. La concavidad de una curva $ROC$ puede determinarse por métricas como el ratio de probabilidad negativa ($NLR = (1 - ROC(t))/(1 - t)$). El $NLR$ se corresponde con la pendiente de una línea recta que pasa por un punto de la curva hasta el punto $(1,1)$, por tanto, un $NLR$ decreciente se corresponderá con curvas $ROC$ cóncavas. 

De esta forma, si se asegura la concavidad en la región de alta sensibilidad, o el $NLR$ es decreciente en dicha zona, un nuevo límite inferior puede ser definido. Este nuevo límite más ajustado se define por el triángulo de puntos $(FPR_0, TPR_0), (1, TPR_0)$ y $(1,1)$, que desarrollado se puede expresar de la siguiente forma:

\begin{equation}
    \dfrac{1}{2}(1 - FPR_0)(1 - TPR_0) \le A \le (1-FPR_0)(1 - TPR_0)
    (\#eq:tpauc-concave-bound)
\end{equation}

Finalmente, las curvas $ROC$ pueden presentar ganchos (regiones donde no la curva $ROC$ no es cóncava) en el extremo superior de la región de interés, en estas situaciones el límite inferior definido en \@ref(eq:tpauc-upper-bound) tampoco se ajusta a la curva. De esta forma, si existe un punto en la región de interés que cumpla $NLR(t) > \text{max}\ \{1, NLR_0\}$, no se podrá encontrar un límite inferior que se ajuste y por tanto:

\begin{equation}
    0 \le A \le (1-FPR_0)(1 - TPR_0)
    (\#eq:tpauc-upper-hook-bound)
\end{equation}

Conociendo todos los casos posibles, y aplicando la transformación descrita por McClish \@ref(eq:mcclish-transformation) se obtiene el $FpAUC$.

De forma similar a los casos anteriores, el $FpAUC$ se calculará para contrastarlo con el resto. Al igual que el $NpAUC$, este índice no se encuentra implementado en ningún paquete conocido, por lo que para facilitar su aplicación se ha implementado en el paquete `ROCnGO`. Para realizar los cálculos sobre el conjunto de datos, se utilizará la función `ROCnGO::fp_auc`:


```r
fpauc_sn_diagnostic <- map_dbl(
    diagnostic_dataset,
    \(gene) {
        fp_auc(
            predictor = gene,
            response = roc_data$disease,
            lower_tpr = 0.9
        )
    },
    .progress = TRUE
)
```


```r
fpauc_sn_prognostic <- map_dbl(
    prognostic_dataset,
    \(gene) {
        fp_auc(
            predictor = gene,
            response = disease$prognostic,
            lower_tpr = 0.9
        )
    },
    .progress = TRUE
)
```

## Comparación de índices

Hasta ahora, para identificar posibles marcadores de diagnóstico y pronóstico se han utilizado distintos índices en condiciones de alta sensibilidad y especificidad. Sin embargo, como se ha explicado previamente, cada índice presenta una serie de ventajas y desventajas. Por consiguiente, para evaluar cual de todos los índices presenta un mejor rendimiento para seleccionar clasificadores, se compararán cada uno de sus resultados.

### Agrupar por condiciones

Para facilitar su acceso, los resultados se agruparán en función de las condiciones que se han medido (alta sensibilidad o especificidad) y su finalidad (diagnóstico o pronóstico).


```r
indexes_sn_diagnostic <- tibble(
    identifier = names(auc_diagnostic),
    auc = auc_diagnostic,
    pauc = pauc_sn_diagnostic,
    npauc = npauc_sn_diagnostic,
    fpauc = fpauc_sn_diagnostic
)
```


```r
indexes_sn_prognostic <- tibble(
    identifier = names(auc_prognostic),
    auc = auc_prognostic,
    pauc = pauc_sn_prognostic,
    npauc = npauc_sn_prognostic,
    fpauc = fpauc_sn_prognostic
)
```


```r
indexes_sp_diagnostic <- tibble(
    identifier = names(auc_diagnostic),
    auc = auc_diagnostic,
    pauc = pauc_sp_diagnostic,
    spauc = spauc_sp_diagnostic,
    tpauc = tpauc_sp_diagnostic
)
```


```r
indexes_sp_prognostic <- tibble(
    identifier = names(auc_prognostic),
    auc = auc_prognostic,
    pauc = pauc_sp_prognostic,
    spauc = spauc_sp_prognostic,
    tpauc = tpauc_sp_prognostic
)
```

### Marcadores para diagnóstico

#### Métricas generales

Antes de analizar en detalle los distintos índices, es conveniente describir algunas métricas generales. Para ello es posible utilizar la función `ROCnGO::summarize_dataset`.


```r
metrics_diagnostic_tpr <- summarize_dataset(
    diagnostic_dataset,
    predictors = NULL,
    response = roc_data$disease,
    ratio = "tpr",
    threshold = 0.9
)
```


```r
metrics_diagnostic_fpr <- summarize_dataset(
    diagnostic_dataset,
    predictors = NULL,
    response = roc_data$disease,
    ratio = "fpr",
    threshold = 0.1
)
```

El conjunto de datos utilizado para diagnóstico contiene la expresión de 2651 genes provenientes de 554 muestras extraídas de pacientes con cáncer de próstata, de las cuales, 52 se corresponden con tejido sano y 502 con tejido canceroso. Los genes de estudio se corresponden con los expresados diferencialmente en experimentos de scRNA-seq. 

De los 2651 genes seleccionados, 1722 (64.93%) presentan un $AUC \ge 0.5$. En regiones de alta sensibilidad ($TPR \ge 0.9$), de este número 184 presentan curvas $ROC$ cóncavas (10.68% y 6.94% del total), 148 presentan algún tipo de gancho sin atravesar la línea del azar (8.59% y 5.58% del total) y 1390 atraviesan en algún punto dicha línea (80.72% y 52.43% del total).

Para regiones de alta especificidad ($FPR \le 0.1$), de los predictores con un $AUC \ge 0.5$, 405 presentan curvas cóncavas (23.51%, un 15.28% del total), 1032 curvas parcialmente adecuadas (59.93%, un 38.92% del total) y 285 contienen algún punto bajo la línea de azar (16.55%, un 10.75% del total).

#### Índices de sensbilidad: $NpAUC$ y $FpAUC$

Tras describir algunas métricas generales de los índices, se explorarán algunas medidas de tendencia central para realizar análisis más detallados. Con estas medidas será posible comprobar, entre otras cosas, el rango de valores que toman los índices.


```r
summary(indexes_sn_diagnostic)
```

```
##   identifier             auc               pauc              npauc        
##  Length:2651        Min.   :0.07099   Min.   :0.000000   Min.   :0.00000  
##  Class :character   1st Qu.:0.44546   1st Qu.:0.003134   1st Qu.:0.03134  
##  Mode  :character   Median :0.56491   Median :0.009447   Median :0.09447  
##                     Mean   :0.55931   Mean   :0.010910   Mean   :0.10910  
##                     3rd Qu.:0.68095   3rd Qu.:0.016338   3rd Qu.:0.16338  
##                     Max.   :0.94369   Max.   :0.060244   Max.   :0.60244  
##                                                                           
##      fpauc       
##  Min.   :0.5000  
##  1st Qu.:0.6952  
##  Median :0.7659  
##  Mean   :0.7517  
##  3rd Qu.:0.8164  
##  Max.   :0.9957  
##  NA's   :145
```

Como puede observarse, el $pAUC$ toma valores inferiores a $0.5$, el $NpAUC$ toma valores tanto valores inferiores como superiores a dicho valor y finalmente el $FpAUC$ solo presenta valores superiores. A simple vista, se puede observar como el $pAUC$ se convierte en un métrica de difícil interpretación, al contrario que el $NpAUC$ y $FpAUC$ que tienen límites definidos. Por otra parte, aunque tanto el $NpAUC$ como el $FpAUC$ están comprendidos por debajo de $1$ por su límite superior, no presentan el mismo límite inferior (secciones \@ref(npauc) y \@ref(fpauc)).

Para ilustrar esta situación se utilizará el gen `ENSG00000169347`, el clasificador de mayor $FpAUC$ en la región de interés. Sus métricas pueden ser rápidamente visualizadas utilizando la función `ROCnGO::summarize_predictor`.


```r
summarize_predictor(roc_data, ENSG00000169347, disease, "tpr", 0.9)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.4502184 </td>
   <td style="text-align:right;"> 0.0134041 </td>
   <td style="text-align:right;"> 0.1340408 </td>
   <td style="text-align:right;"> 0.9957314 </td>
   <td style="text-align:left;"> Concave </td>
  </tr>
</tbody>
</table>

El clasificador toma aproximadamente el valor máximo para el $FpAUC$, sin embargo, presenta valores reducidos de $AUC$, $pAUC$ e incluso $NpAUC$. Para explicar estos resultados, la curva $ROC$ del clasificador será representada con la función `ROCnGO::plot_roc_curve`. Además, para facilitar la explicación, también se representará el límite inferior del $NpAUC$.


```r
plot_roc_curve(roc_data, response = disease, predictor = ENSG00000169347) +
    add_chance_line() +
    add_tpr_threshold_line(0.9) +
    add_npauc_normalized_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000169347,
        threshold = 0.9
    )
```

<img src="05_AUC_indexes_files/figure-html/Predictor ENSG00000169347 con límite inferior de NpAUC-1.png" width="672" />

Como se puede observar, el límite inferior del $NpAUC$ se corresponde con la mitad del área de la banda $(TPR_0, 1)$. Dado que la curva $ROC$ queda descrita por debajo de esta área, el $NpAUC$ toma valores inferiores a $0.5$ y en consecuencia dificultando su interpretación. 

Este no es el caso del límite inferior definido por el $FpAUC$, el cual se ajusta a la forma de la curva. Para visualizarlo de forma detalla, se graficará la curva $ROC$ en la región de alta sensibilidad utilizando la función `ROCnGO::plot_partial_roc_curve`.


```r
plot_partial_roc_curve(
    roc_data, 
    response = disease,
    predictor = ENSG00000169347,
    ratio = "tpr",
    threshold = 0.9
    ) +
    add_chance_line() +
    add_tpr_threshold_line(0.9) +
    add_fpauc_concave_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000169347,
        threshold = 0.9
    )
```

<img src="05_AUC_indexes_files/figure-html/Predictor con límite inferior de FpAUC-1.png" width="672" />

Ampliando la figura en la región de interés, se observa como el límite inferior del $FpAUC$ se ajusta a la región definida por la curva $ROC$, en lugar de a toda la banda $(TPR_0, 1)$. Así, el $FpAUC$ toma valores mayores a $0.5$, siendo una métrica más interpretable. Este punto puede contrastarse con `ENSG00000105707`, el predictor de mayor valor para $NpAUC$.


```r
summarize_predictor(roc_data, ENSG00000105707, disease, "tpr", 0.9)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.9436868 </td>
   <td style="text-align:right;"> 0.0602436 </td>
   <td style="text-align:right;"> 0.6024364 </td>
   <td style="text-align:right;"> 0.8642639 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
</tbody>
</table>



```r
plot_roc_curve(roc_data, response = disease, predictor = ENSG00000105707) +
    add_chance_line() +
    add_tpr_threshold_line(0.9) +
    add_npauc_normalized_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000105707,
        threshold = 0.9
    )
```

<img src="05_AUC_indexes_files/figure-html/Predictor ENSG00000105707 con límite inferior de NpAUC-1.png" width="672" />


```r
plot_partial_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000105707,
        ratio = "tpr",
        threshold = 0.9
    ) +
    add_chance_line() +
    add_tpr_threshold_line(0.9)
```

<img src="05_AUC_indexes_files/figure-html/Predictor ENSG00000105707 con límite inferior de FpAUC-1.png" width="672" />

En este caso, la mayor parte de la curva $ROC$ está contenida por encima del límite definido por $NpAUC$. Esto confirma el razonamiento anterior, ya que se obtienen valores superiores a $0.5$. Por otra parte, el $FpAUC$ no puede definir ningún límite inferior puesto que la curva presenta un gancho bajo la línea de azar en el extremo superior. A pesar de ello, el $FpAUC$ considera esta posibilidad por lo que sigue devolviendo valores interpretables.

Otra característica de interés de los índices es su capacidad para distinguir predictores con el mismo $pAUC$ y curvas $ROC$ que se cruzan. De todo el conjunto de datos, 633 predictores presentan un $pAUC$ idéntico a otro, formando 246 grupos de $pAUC$ iguales. Para evaluar el rendimiento de los índices en este contexto se utilizarán los genes `ENSG00000065361`y `ENSG00000130653`, aquellos con el mayor valor de $pAUC$ entre los grupos formados.


```r
bind_rows(
    "ENSG00000065361" = summarize_predictor(
        roc_data, ENSG00000065361, disease, "tpr", 0.9
    ),
    "ENSG00000130653" = summarize_predictor(
        roc_data, ENSG00000130653, disease, "tpr", 0.9
    ),
    .id = "Identifier"
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000065361 </td>
   <td style="text-align:right;"> 0.7715676 </td>
   <td style="text-align:right;"> 0.0274287 </td>
   <td style="text-align:right;"> 0.2742875 </td>
   <td style="text-align:right;"> 0.8565737 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000130653 </td>
   <td style="text-align:right;"> 0.7384692 </td>
   <td style="text-align:right;"> 0.0274287 </td>
   <td style="text-align:right;"> 0.2742875 </td>
   <td style="text-align:right;"> 0.8351422 </td>
   <td style="text-align:left;"> Partially proper </td>
  </tr>
</tbody>
</table>

```r
plot_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000065361
) +
    add_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000130653
    ) +
    add_chance_line() +
    add_threshold_line(0.9, "tpr")
```

<img src="05_AUC_indexes_files/figure-html/Comparación de marcadores ENSG00000065361 y ENSG00000130653-1.png" width="672" />

Ambos predictores presentan $AUC$ relativamente altos y sus curvas se cruzan en diversos puntos dentro de la región de alta sensibilidad. Dado que el $NpAUC$ simplemente realiza una normalización utilizando la banda de alta sensibilidad \@ref(eq:npauc-normalized-bounds), la métrica no es capaz de seleccionar un predictor mejor, dando el mismo valor para ambas ($NpAUC = 0.2742875$). Por otra parte, el $FpAUC$ adapta sus límites a la forma de la curva $ROC$ por ello es capaz de seleccionar al de mejor rendimiento, `ENSG00000065361`. Siguiendo este ejemplo, se deduce que el $FpAUC$ presenta una mayor capacidad de discriminación entre predictores con igual $pAUC$, lo cual puede contrastarse con todo el conjunto de datos.


```r
metrics_diagnostic_tpr$data %>%
    filter(pauc != 0) %>%
    group_by(pauc) %>%
    filter(n() > 1) %>%
    mutate(
        npauc_equal = all(
            map_lgl(np_auc, \(x) near(x, np_auc[[1]]))
        ),
        fpauc_equal = all(
            map_lgl(fp_auc, \(x) near(x, fp_auc[[1]]))
        )
    ) %>%
    summarize(
        npauc_equal = all(npauc_equal),
        fpauc_equal = all(fpauc_equal)
    ) %>%
    summarize(
        npauc_equal = sum(npauc_equal),
        fpauc_equal = sum(fpauc_equal)
    ) %>%
    mutate(
        npauc_distintct = 246 - npauc_equal,
        fpauc_distintct = 246 - fpauc_equal,
    ) %>%
    relocate(
        sort(names(.))
    )
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> fpauc_distintct </th>
   <th style="text-align:right;"> fpauc_equal </th>
   <th style="text-align:right;"> npauc_distintct </th>
   <th style="text-align:right;"> npauc_equal </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 199 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 246 </td>
  </tr>
</tbody>
</table>

De los 246 grupos formados, el $NpAUC$ no es capaz de desempatar en ningún grupo, mientras que el $FpAUC$ es capaz de seleccionar un mejor clasificador en 47 de ellos. Estos cálculos sobre todo el conjunto de datos refuerzan la afirmación de que el $FpAUC$ tiene una mejor capacidad de discriminación que el $NpAUC$.

#### Índices de especificidad: $SpAUC$ y $TpAUC$

Al igual que con los índices de sensibilidad, en primer lugar, se explorarán algunas métricas de tendencia central.


```r
summary(indexes_sp_diagnostic)
```

```
##   identifier             auc               pauc              spauc       
##  Length:2651        Min.   :0.07099   Min.   :0.000000   Min.   :0.5000  
##  Class :character   1st Qu.:0.44546   1st Qu.:0.004635   1st Qu.:0.5278  
##  Mode  :character   Median :0.56491   Median :0.011952   Median :0.5645  
##                     Mean   :0.55931   Mean   :0.015553   Mean   :0.5813  
##                     3rd Qu.:0.68095   3rd Qu.:0.023314   3rd Qu.:0.6213  
##                     Max.   :0.94369   Max.   :0.074939   Max.   :0.8681  
##                                                          NA's   :711     
##      tpauc       
##  Min.   :0.0000  
##  1st Qu.:0.7347  
##  Median :0.7841  
##  Mean   :0.7763  
##  3rd Qu.:0.8312  
##  Max.   :1.0000  
## 
```

Observando las métricas el $pAUC$ solo toma valores inferiores a $0.5$ lo cual dificulta su interpretabilidad, a diferencia del $SpAUC$ y $TpAUC$ que pueden tomar valores superiores. Cabe destacar que el $TpAUC$ toma valores iguales a 0 cuando teóricamente debería estar comprendido en $(0.5, 1)$. 


```r
indexes_sp_diagnostic %>%
    filter(tpauc == 0)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> spauc </th>
   <th style="text-align:right;"> tpauc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000169116 </td>
   <td style="text-align:right;"> 0.1625038 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000135378 </td>
   <td style="text-align:right;"> 0.1314358 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000148985 </td>
   <td style="text-align:right;"> 0.5433650 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000171867 </td>
   <td style="text-align:right;"> 0.1755287 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000085662 </td>
   <td style="text-align:right;"> 0.1333895 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000167711 </td>
   <td style="text-align:right;"> 0.2537542 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000101856 </td>
   <td style="text-align:right;"> 0.0867300 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000185432 </td>
   <td style="text-align:right;"> 0.1069951 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000121691 </td>
   <td style="text-align:right;"> 0.2978854 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000148344 </td>
   <td style="text-align:right;"> 0.3502145 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000091138 </td>
   <td style="text-align:right;"> 0.2779651 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000186526 </td>
   <td style="text-align:right;"> 0.4682616 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000156642 </td>
   <td style="text-align:right;"> 0.1963684 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000170954 </td>
   <td style="text-align:right;"> 0.2587726 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

Una exploración de los datos muestra como estos valores del $TpAUC$ se corresponden con valores de $pAUC$ iguales a $0$, lo cual se puede comprobar de forma gráfica representando algunos de ellos.


```r
plot_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000169116
) +
    add_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000135378
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Representación de predictores con pAUC igual a 0-1.png" width="672" />

Dado que las curvas de ambos predictores toman valores de $TPR = 0$, no hay un area definida debajo de la curva. Ampliando la imagen esto puede observarse con mayor claridad.


```r
plot_partial_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000169116,
    ratio = "fpr",
    threshold = 0.1
) + add_partial_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000135378,
        ratio = "fpr",
        threshold = 0.1
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Ampliación de predictores con pAUC igual a 0-1.png" width="672" />

Tras comprobar el origen de los valores atípicos, es posible recalcular las métricas filtrando los valores de $pAUC = 0$. De esta forma, $TpAUC$ y $SpAUC$ se encuentran limitados por encima de $0.5$ convirtiéndose en métricas fácilmente interpretables.


```r
summary(indexes_sp_diagnostic %>% filter(pauc != 0))
```

```
##   identifier             auc               pauc               spauc       
##  Length:2637        Min.   :0.07099   Min.   :7.660e-06   Min.   :0.5000  
##  Class :character   1st Qu.:0.44622   1st Qu.:4.766e-03   1st Qu.:0.5278  
##  Mode  :character   Median :0.56624   Median :1.198e-02   Median :0.5645  
##                     Mean   :0.56098   Mean   :1.564e-02   Mean   :0.5813  
##                     3rd Qu.:0.68162   3rd Qu.:2.337e-02   3rd Qu.:0.6213  
##                     Max.   :0.94369   Max.   :7.494e-02   Max.   :0.8681  
##                                                           NA's   :697     
##      tpauc       
##  Min.   :0.5000  
##  1st Qu.:0.7354  
##  Median :0.7846  
##  Mean   :0.7804  
##  3rd Qu.:0.8314  
##  Max.   :1.0000  
## 
```

A pesar de tener buena interpretabilidad, el $SpAUC$ presenta valores de `NA`. Como se ha explicado anteriormente, el $SpAUC$ tiene su límite definido con el trapecio de puntos $(FPR_1, 0), (FPR_2, 0), (FPR_1, FPR_1)$ y $(FPR_2, FPR_2)$, o en otras palabras, el que coincide con la línea del azar. Este límite no se ajusta a curvas definidas por debajo de la línea (sección \@ref(spauc)), por esta razón la función devuelve estos valores. Esto es posible comprobarlo gráficamente con el predictor `ENSG00000181754`, un predictor que aun presentan una forma cóncava en la región de interés pasa por debajo de la línea del azar.


```r
summarize_predictor(
    roc_data,
    ENSG00000181754,
    disease,
    "fpr",
    0.1
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.4311408 </td>
   <td style="text-align:right;"> 0.0049425 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.8215741 </td>
   <td style="text-align:left;"> Concave </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000181754
) +
    add_spauc_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000181754,
        lower_threshold = 0,
        upper_threshold = 0.1
    ) +
    add_tpauc_concave_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000181754,
        lower_threshold = 0,
        upper_threshold = 0.1
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Plot predictor ENSG00000181754-1.png" width="672" />

Como se observa en la figura, el límite inferior del $SpAUC$ se encuentra definido por encima de la curva $ROC$, por otra parte el límite del $TpAUC$ se adecua a la región y forma de la curva. Para contrastar este ajuste es posible seleccionar otro predictor con una forma de curva distinta, por ejemplo, `ENSG00000111667`.


```r
summarize_predictor(
    roc_data,
    ENSG00000111667,
    disease,
    "fpr",
    0.1
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.413768 </td>
   <td style="text-align:right;"> 0.0049954 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.7985348 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000111667
) +
    add_spauc_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000111667,
        lower_threshold = 0,
        upper_threshold = 0.1
    ) +
    add_tpauc_under_chance_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000111667,
        lower_threshold = 0,
        upper_threshold = 0.1
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Plot predictor ENSG00000111667-1.png" width="672" />

En este caso, la curva $ROC$ queda definida por debajo de la línea de azar y el $TpAUC$ no puede establecer un área mínima, por lo que establece el límite inferior a $0$. Por otra parte, el límite inferior del $SpAUC$ no se adecua y no puede ser calculado. Con estas observaciones, se deduce que tanto el $SpAUC$ y $TpAUC$ son métricas interpretables cuando se trabaja con curvas $ROC$ adecuadas al quedar comprendidas entre $0.5$ y $1$. Sin embargo, el $TpAUC$ se convierte en una métrica más adecuada en condiciones generales, ya que sus límites se ajustan a la forma de la curva $ROC$ a diferencia del $SpAUC$, imposibilitando su uso en curvas comprendidas debajo de la línea del azar.

Una vez considerada la interpretabilidad de las métricas, se puede contrastar su capacidad para distinguir entre predictores con curvas que se cruzan e igual $pAUC$. Del conjunto de datos, 374 predictores presentan un $pAUC$ igual a otro, formando 174 grupos de $pAUC$ iguales. Para evaluar el rendimiento de los índices utilizaremos los predictores `ENSG00000175336` y `ENSG00000169733`, aquellos que presentan el mayor valor de $pAUC$ en este subgrupo de predictores.


```r
bind_rows(
    "ENSG00000175336 " = summarize_predictor(
        roc_data,
        ENSG00000175336,
        disease,
        "fpr",
        0.1
    ),
    "ENSG00000169733" = summarize_predictor(
        roc_data,
        ENSG00000169733,
        disease,
        "fpr",
        0.1
    ),
    .id = "Identifier"
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000175336 </td>
   <td style="text-align:right;"> 0.8286086 </td>
   <td style="text-align:right;"> 0.0475329 </td>
   <td style="text-align:right;"> 0.7238576 </td>
   <td style="text-align:right;"> 0.8695316 </td>
   <td style="text-align:left;"> Partially Proper </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000169733 </td>
   <td style="text-align:right;"> 0.8126724 </td>
   <td style="text-align:right;"> 0.0475329 </td>
   <td style="text-align:right;"> 0.7238576 </td>
   <td style="text-align:right;"> 0.9123511 </td>
   <td style="text-align:left;"> Partially Proper </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000175336
) +
    add_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000169733
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Comparación de marcadores ENSG00000175336 y ENSG00000169733-1.png" width="672" />

Comprobando sus $AUC$, ambos predictores presentan un buen rendimiento global. El $SpAUC$ se mantiene igual para ambos predictores ($SpAUC = 0.7238576$) mientras que el $TpAUC$ toma valores distintos, lo que permite identificar el predictor de mejor rendimiento. Esta diferencia puede deberse a que el $TpAUC$ adapta sus límites a los puntos de la curva (sección \@ref(tpauc)), a diferencia del $SpAUC$ que los mantienen iguales.

Esta capacidad de discriminación entre predictores se puede contrastar con el resto del conjunto de datos:


```r
metrics_diagnostic_fpr$data %>%
    filter(pauc != 0) %>%
    group_by(pauc) %>%
    filter(n() > 1) %>%
    mutate(
        spauc_equal = all(
            map_lgl(sp_auc, \(x) near(x, sp_auc[[1]]))
        ),
        tpauc_equal = all(
            map_lgl(tp_auc, \(x) near(x, tp_auc[[1]]))
        )
    ) %>%
    summarize(
        spauc_equal = all(spauc_equal),
        tpauc_equal = all(tpauc_equal)
    ) %>%
    dplyr::count(spauc_equal, tpauc_equal)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> spauc_equal </th>
   <th style="text-align:left;"> tpauc_equal </th>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> 34 </td>
  </tr>
</tbody>
</table>

De los 374 predictores con un $pAUC > 0$ y al menos igual a otro predictor (174 grupos), se observa como el $SpAUC$ no es capaz de diferenciar el predictor de mayor rendimiento en ningún caso. Además, el índice solo toma valores para 95 de los grupos, no estando definido para el resto. Por otro lado, el $TpAUC$ es capaz de distinguir un mejor predictor en 138 de los grupos, siendo incapaz de hacerlo solamente en 36 de ellos.

Los resultados sobre el conjunto de datos refuerzan que ambos índices son fácilmente interpretables. Sin embargo, muestran que el $SpAUC$ no es aplicable cuando la curva $ROC$ pasa por debajo de la línea del azar. Además, el $SpAUC$ no es capaz de discriminar entre ciertos predictores. Por otra parte, el $TpAUC$ resulta más útil en la práctica al poder aplicarse a cualquier tipo de curva y ser capaz de discriminar entre predictores con igual $pAUC$.

### Marcadores para pronóstico

#### Métricas generales

En esta sección se repetirán los procedimientos realizados previamente, pero en esta ocasión serán aplicados a marcadores de pronóstico.


```r
metric_prognostic_tpr <- summarize_dataset(
    prognostic_dataset,
    predictors = NULL,
    response = disease$prognostic,
    ratio = "tpr",
    threshold = 0.9
)
```


```r
metric_prognostic_fpr <- summarize_dataset(
    prognostic_dataset,
    predictors = NULL,
    response = disease$prognostic,
    ratio = "fpr",
    threshold = 0.1
)
```

El conjunto de datos utilizado para pronóstico contiene la expresión de 2652 genes provenientes de 502 muestras extraídas de pacientes con la enfermedad. Los genes de estudio se corresponde con los diferencialmente expresados en experimentos de scRNA-seq.

De los 2651 genes seleccionados, 1015 (38.27%) presentan un $AUC \ge 0.5$. En regiones de alta sensibilidad ($TPR \ge 0.9$), de este número 173 presentan curvas $ROC$ cóncavas (17.04% y 6.52% del total), 296 presentan algún tipo de gancho sin llegar a atravesar la línea del azar (29.16% y 11.16% del total) y 546 atraviesan en algún punto dicha línea (52.43% y 20.59% del total).

Para regiones de alta especificidad ($FPR \le 0.1$), de los predictores con un $AUC \ge 0.5$, 127 presentan curvas cóncavas (12.51%, un 4.79% del total), 376 curvas parcialmente adecuadas (37.04%, un 14.18% del total) y 511  contienen algún punto bajo la línea de azar (50.34%, un 19.27% del total).

#### Índices de sensibilidad: $NpAUC$ y $FpAUC$

De forma similar al procedimiento en condiciones de diagnóstico, se explorarán algunas medidas de tendencia central.


```r
summary(metric_prognostic_tpr$data)
```

```
##   identifier             auc              pauc              np_auc       
##  Length:2651        Min.   :0.2485   Min.   :0.000000   Min.   :0.00000  
##  Class :character   1st Qu.:0.4289   1st Qu.:0.003098   1st Qu.:0.03098  
##  Mode  :character   Median :0.4767   Median :0.004921   Median :0.04921  
##                     Mean   :0.4778   Mean   :0.005477   Mean   :0.05477  
##                     3rd Qu.:0.5249   3rd Qu.:0.007376   3rd Qu.:0.07376  
##                     Max.   :0.7377   Max.   :0.024389   Max.   :0.24389  
##                                                                          
##      fp_auc       curve_shape       
##  Min.   :0.5000   Length:2651       
##  1st Qu.:0.6950   Class :character  
##  Median :0.7434   Mode  :character  
##  Mean   :0.7354                     
##  3rd Qu.:0.7848                     
##  Max.   :0.9758                     
##  NA's   :1
```

En este contexto, tanto el $pAUC$ como el $NpAUC$ toman valores inferiores a $0.5$ mientras que el $FpAUC$ toma superiores. Estos resultados probablemente se deban a los límites definidos por el $NpAUC$ y $FpAUC$. El $NpAUC$ define su límite inferior como la mitad del área de la banda $(TPR_0, 1)$ \@ref(eq:npauc-normalized-bounds), en consecuencia, una curva $ROC$ contenida por este límite resultará en valores de $NpAUC$ inferiores a $0.5$. Para visualizarlo, se escogerá el predictor `ENSG00000185615` que posee el mayor valor de $NpAUC$ en los datos. 


```r
summarize_predictor(
    prognostic_dataset %>% mutate(prognostic = disease$prognostic),
    ENSG00000185615,
    prognostic,
    "tpr",
    0.9
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.7143536 </td>
   <td style="text-align:right;"> 0.0243888 </td>
   <td style="text-align:right;"> 0.2438877 </td>
   <td style="text-align:right;"> 0.7743255 </td>
   <td style="text-align:left;"> Partially proper </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    prognostic_dataset %>% mutate(prognostic = disease$prognostic),
    response = prognostic,
    predictor = ENSG00000185615
    ) +
    add_chance_line() +
    add_tpr_threshold_line(0.9) +
    add_fpauc_partially_proper_lower_bound(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000185615,
        threshold = 0.9
    ) +
    add_npauc_normalized_lower_bound(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000185615,
        threshold = 0.9
    ) 
```

<img src="05_AUC_indexes_files/figure-html/Predictor ENSG00000185615 con límite inferior de NpAUC-1.png" width="672" />

Como se observa en la figura, la curva $ROC$ queda definida por debajo del límite definido por $NpAUC$ pero por encima del definido por $FpAUC$. De esta forma, se confirma que el $FpAUC$ presenta mejor interpretabilidad, además de ser aplicable independientemente de la forma de la curva.

Además de su interpretabilidad, es interesante evaluar la capacidad para distinguir predictores con igual $pAUC$ y cuyas curvas se cruzan en algún punto en la región de interés. De todo el conjunto de datos, 241 predictores presentan un $pAUC$ idéntico a al menos otro, formando 114 grupos en función de esta métrica. Para evaluar el rendimiento de los índices se utilizarán los predictores `ENSG00000272333` y `ENSG00000096070`, aquellos de mayor $pAUC$ en los grupos formados.


```r
bind_rows(
    "ENSG00000272333" = summarize_predictor(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        ENSG00000272333,
        prognostic,
        "tpr",
        0.9
    ),
    "ENSG00000096070" = summarize_predictor(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        ENSG00000096070,
        prognostic,
        "tpr",
        0.9
    ),
    .id = "Identifier"
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000272333 </td>
   <td style="text-align:right;"> 0.5947269 </td>
   <td style="text-align:right;"> 0.0131597 </td>
   <td style="text-align:right;"> 0.1315975 </td>
   <td style="text-align:right;"> 0.7782804 </td>
   <td style="text-align:left;"> Partially proper </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000096070 </td>
   <td style="text-align:right;"> 0.5004503 </td>
   <td style="text-align:right;"> 0.0131597 </td>
   <td style="text-align:right;"> 0.1315975 </td>
   <td style="text-align:right;"> 0.8087762 </td>
   <td style="text-align:left;"> Concave </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    prognostic_dataset %>% mutate(prognostic = disease$prognostic),
    response = prognostic,
    predictor = ENSG00000272333
) +
    add_roc_curve(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000096070
    ) +
    add_chance_line() +
    add_threshold_line(0.9, "tpr")
```

<img src="05_AUC_indexes_files/figure-html/Comparación de marcadores ENSG00000272333 y ENSG00000096070-1.png" width="672" />

En este caso, los predictores no presentan un buen rendimiento global, a diferencia de en la región de interés. En esta región ambas curvas se cruzan en varios puntos y mientras que el $NpAUC$ no es capaz de distinguir el predictor de mayor rendimiento ($NpAUC = 0.13159748$), el $FpAUC$ identifica que `ENSG00000096070` ($FpAUC = 0.8087762$) rinde mejor que `ENSG00000272333` ($FpAUC = 0.7782804$). De forma similar a las condiciones de diagnóstico, dado que el $FpAUC$ adapta sus límites a la forma y  puntos de la curva, puede distinguir entre predictores con igual $pAUC$ cuyas curvas se cruzan. Esta capacidad de discriminación puede contrastarse con el resto del conjunto de datos.


```r
metric_prognostic_tpr$data %>%
    filter(pauc != 0) %>%
    group_by(pauc) %>%
    filter(n() > 1) %>%
    mutate(
        npauc_equal = all(
            map_lgl(np_auc, \(x) near(x, np_auc[[1]]))
        ),
        fpauc_equal = all(
            map_lgl(fp_auc, \(x) near(x, fp_auc[[1]]))
        )
    ) %>%
    summarize(
        npauc_equal = all(npauc_equal),
        fpauc_equal = all(fpauc_equal)
    ) %>%
    dplyr::count(npauc_equal, fpauc_equal)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> npauc_equal </th>
   <th style="text-align:left;"> fpauc_equal </th>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> 76 </td>
  </tr>
</tbody>
</table>

El $FpAUC$ es capaz de desempatar entre los predictores en 38 de los 114 grupos formados, mientras que el $NpAUC$ no es capaz de hacerlo en ninguno de ellos.

#### Índices de especificidad: $SpAUC$ y $TpAUC$

Igual que en los casos anteriores, para realizar la evaluación de los índices se explorarán algunas métricas generales en primer lugar.


```r
summary(metric_prognostic_fpr$data)
```

```
##   identifier             auc              pauc               sp_auc      
##  Length:2651        Min.   :0.2485   Min.   :0.0002866   Min.   :0.5000  
##  Class :character   1st Qu.:0.4289   1st Qu.:0.0036109   1st Qu.:0.5065  
##  Mode  :character   Median :0.4767   Median :0.0054696   Median :0.5148  
##                     Mean   :0.4778   Mean   :0.0063296   Mean   :0.5196  
##                     3rd Qu.:0.5249   3rd Qu.:0.0082330   3rd Qu.:0.5276  
##                     Max.   :0.7377   Max.   :0.0277573   Max.   :0.6198  
##                                                          NA's   :1167    
##      tp_auc       curve_shape       
##  Min.   :0.5000   Length:2651       
##  1st Qu.:0.7088   Class :character  
##  Median :0.7492   Mode  :character  
##  Mean   :0.7443                     
##  3rd Qu.:0.7855                     
##  Max.   :0.9322                     
## 
```

De nuevo, el $pAUC$ solo toma valores inferiores a $0.5$ a diferencia del $SpAUC$ y $TpAUC$ que toman valores superiores. Además, el $SpAUC$ toma valores de `NA`, los cuales se corresponden para casos donde el límite inferior \@ref(eq:spauc-bounds) no se cumple, como se ha visto previamente. A diferencia del escenario anterior, el $TpAUC$ no toma valores atípicos por debajo de $0.5$, esto se debe a que no hay predictores con $pAUC = 0$.

Para contrastar la interpretabilidad de ambos índices se utilizará arbitrariamente el predictor `ENSG00000110013`, que no presenta valor para el $SpAUC$ y uno relativamente alto para el $TpAUC$.


```r
summarize_predictor(
    prognostic_dataset %>% mutate(prognostic = disease$prognostic),
    ENSG00000110013,
    prognostic,
    "fpr",
    0.1
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.3827888 </td>
   <td style="text-align:right;"> 0.0008352 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.9322034 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    prognostic_dataset %>% mutate(prognostic = disease$prognostic),
    response = prognostic,
    predictor = ENSG00000110013
) +
    add_spauc_lower_bound(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000110013,
        lower_threshold = 0,
        upper_threshold = 0.1
    ) +
    add_tpauc_under_chance_lower_bound(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000110013,
        lower_threshold = 0,
        upper_threshold = 0.1
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Plot predictor ENSG00000196363-1.png" width="672" />

Como se observa en la figura, la curva $ROC$ se encuentra contenida bajo el límite inferior de $SpAUC$ y por tanto obteniendo un valor de `NA`. Por otra parte, el límite de $TpAUC$ ($0$, dado que $TPR = 0$ para toda la región) se ajusta a la forma de la curva y por tanto puede calcularse.

Al igual que en los casos anteriores, de forma general $SpAUC$ y $TpAUC$ son métricas interpretables dado que se comprenden entre valores de $0.5$ y $1$. Sin embargo, el $SpAUC$ define límites que no se ajustan a ciertas curvas y por tanto imposibilitando su aplicación, mientras que el $TpAUC$ se ajusta a la forma de la curva, pudiéndose aplicar a cualquiera y por tanto haciéndola una métrica más idónea.

Una vez considerada la interpretabilidad de ambas métricas, se contrastará su capacidad para discriminar predictores con curvas que se cruzan y de igual $pAUC$. Del conjunto de datos, 917 predictores presentan un $pAUC$ igual al menos a otro, formándose 407 grupos en función de esta métrica. Para estudiar la capacidad de discriminación, se utilizarán los predictores `ENSG00000177943` y `ENSG00000136783`, aquellos de mayor $pAUC$ entre los grupos formados.


```r
bind_rows(
    "ENSG00000177943" = summarize_predictor(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        ENSG00000177943,
        prognostic,
        "fpr",
        0.1
    ),
    "ENSG00000136783" = summarize_predictor(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        ENSG00000136783,
        prognostic,
        "fpr",
        0.1
    ),
    .id = "Identifier"
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000177943 </td>
   <td style="text-align:right;"> 0.6587407 </td>
   <td style="text-align:right;"> 0.0172194 </td>
   <td style="text-align:right;"> 0.5643124 </td>
   <td style="text-align:right;"> 0.8111201 </td>
   <td style="text-align:left;"> Partially Proper </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000136783 </td>
   <td style="text-align:right;"> 0.6306722 </td>
   <td style="text-align:right;"> 0.0172194 </td>
   <td style="text-align:right;"> 0.5643124 </td>
   <td style="text-align:right;"> 0.7711047 </td>
   <td style="text-align:left;"> Partially Proper </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
         prognostic_dataset %>% mutate(prognostic = disease$prognostic),
         response = prognostic,
         predictor = ENSG00000177943
) + 
    add_roc_curve(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000136783
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Comparación de marcadores ENSG00000177943 y ENSG00000136783-1.png" width="672" />

Ambos predictores presentan un rendimiento ligeramente mejor que un clasificador aleatorio a nivel global, pero en la región de interés presentan un rendimiento mucho mejor. Los predictores presentan el mismo $SpAUC$ ($0.5643124$), a diferencia del $TpAUC$ que toma valores distintos e indica que `ENSG00000177943` presenta un mejor rendimiento.


```r
plot_partial_roc_curve(
         prognostic_dataset %>% mutate(prognostic = disease$prognostic),
         response = prognostic,
         predictor = ENSG00000177943,
         ratio = "fpr",
         threshold = 0.1
) + 
    add_partial_roc_curve(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000136783,
        ratio = "fpr",
        threshold = 0.1
    ) +
    add_tpauc_partially_proper_lower_bound(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000177943,
        upper_threshold = 0.1,
        lower_threshold = 0
    ) +
    add_tpauc_partially_proper_lower_bound(
        prognostic_dataset %>% mutate(prognostic = disease$prognostic),
        response = prognostic,
        predictor = ENSG00000136783,
        upper_threshold = 0.1,
        lower_threshold = 0
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Ampliación de comparación de marcadores ENSG00000177943 y ENSG00000136783-1.png" width="672" />

Para este caso, ambos predictores poseen el mismo límite inferior para $TpAUC$, sin embargo sus límites superiores será distintos dado que poseen distinto $TPR_2$, \@ref(eq:tpauc-proper-bounds). Por este motivo, el índice toma valores distintos estableciendo que el predictor `ENSG00000177943` ($TpAUC = 0.8111201$) rinde mejor que `ENSG00000136783` ($TpAUC = 0.7711047$). Esta capacidad de discriminación se puede contrastar en todo el conjunto de datos.


```r
metric_prognostic_fpr$data %>%
    filter(pauc != 0) %>%
    group_by(pauc) %>%
    filter(n() > 1) %>%
    mutate(
        spauc_equal = all(
            map_lgl(sp_auc, \(x) near(x, sp_auc[[1]]))
        ),
        tpauc_equal = all(
            map_lgl(tp_auc, \(x) near(x, tp_auc[[1]]))
        )
    ) %>%
    summarize(
        spauc_equal = all(spauc_equal),
        tpauc_equal = all(tpauc_equal)
    ) %>%
    dplyr::count(spauc_equal, tpauc_equal)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> spauc_equal </th>
   <th style="text-align:left;"> tpauc_equal </th>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 184 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 138 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> 56 </td>
  </tr>
</tbody>
</table>

De los 407 grupos de grupos formados, el $SpAUC$ no es capaz de diferenciar un marcador que rinda mejor en ninguno de ellos. Además, no es capaz de aplicarse debido a la forma de la curva en 194 de ellos. Por otra parte el $TpAUC$ es capaz de diferenciar un mejor predictor en 323 grupos de los 407 totales. Este hecho refuerza que el $TpAUC$ es una mejor métrica a la hora de evaluar el rendimiento de un predictor, al ser aplicable a cualquier tipo de curva y discriminar entre predictores con igual $pAUC$, a diferencia del $SpAUC$.

## Resultados

### Rendimiento de métricas

Tras haber contrastado distintas transformaciones del $pAUC$ en condiciones de alta sensibilidad y especificidad, se observa como el $TpAUC$ y $FpAUC$ presentan un mejor rendimiento a la hora de seleccionar predictores relevantes. Estos índices solventan limitaciones descritas por el propio $pAUC$ y otras transformaciones derivadas, como el $SpAUC$ y el $NpAUC$, lo que hace que se conviertan en métricas más idóneas en la práctica. En primer lugar, ambas métricas presentan límites inferiores y superiores que se ajustan a las curvas $ROC$ de un predictor, independientemente de su forma. En segundo lugar, son capaces de discriminar entre predictores que presentan un mismo $pAUC$ y cuyas curvas $ROC$ se cruzan en la región de interés, algo que no había sido posible hasta ahora. Finalmente, pero no menos importante, comprenden sus valores entre $0.5$ y $1$, lo cual facilita su interpretación a la hora de seleccionar predictores.

Como se ha comprobado en las secciones anteriores, estas limitaciones son comunes en la práctica por lo que el $TpAUC$ y $FpAUC$ presentan un valor añadido como métricas de selección de predictores.

### Relevancia clínica

Aunque teóricamente el $TpAUC$ y $FpAUC$ solventan las limitaciones de otros índices, en la práctica, estas también deben de ser capaces de encontrar marcadores clínicamente relevantes. Esta capacidad puede comprobarse con el conjunto de datos generado y predictores ya conocidos para el diagnóstico y estratificación de pacientes en cáncer de próstata.

Desde la primera vez que fue descrito, el antígeno específico de próstata (PSA) ha sido utilizado en los planes de diagnóstico y tratamiento de pacientes con cáncer de próstata. Desde su aplicación, la mortalidad y el número de diagnósticos en estados avanzados de la enfermedad se ha reducido, sin embargo, este biomarcador presenta ciertas limitaciones [@tumor-heterogeneity]. Aunque el PSA es específico del tejido, este puede presentarse en lesiones benignas, es decir, generando falsos positivos. Además, a la hora de realizar una estratificación de riesgo de un paciente requiere de utilizar otras técnicas y métricas adicionales.

El nombre de PSA es utilizado en el ámbito clínico, pero en otros ámbitos también es conocido como kallikrein related peptidase 3 (*KLK3*), una proteasa expresada en las células del tejido. Este biomarcador se encuentra entre los descritos en el conjunto de datos con el identificador `ENSG00000142515`.


```r
inner_join(
    metrics_diagnostic_tpr$data %>%
        filter(identifier == "ENSG00000142515") %>% select(-curve_shape),
    metrics_diagnostic_fpr$data %>%
        filter(identifier == "ENSG00000142515") %>% select(-curve_shape),
    join_by(identifier, auc),
    suffix = c("_sens", "_spec")
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc_sens </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:right;"> pauc_spec </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000142515 </td>
   <td style="text-align:right;"> 0.649479 </td>
   <td style="text-align:right;"> 0.0191005 </td>
   <td style="text-align:right;"> 0.1910052 </td>
   <td style="text-align:right;"> 0.7640208 </td>
   <td style="text-align:right;"> 0.0135228 </td>
   <td style="text-align:right;"> 0.544857 </td>
   <td style="text-align:right;"> 0.7856116 </td>
  </tr>
</tbody>
</table>

Además de *KLK3*, entre los datos se pueden encontrar otros marcadores de diagnóstico conocidos, como *AMACR* (`ENSG00000242110`) y *FOLH1* (`ENSG00000086205`). El rendimiento de estos marcadores puede comparase con los diversos índices, para así, contrastarlos en función relevancia clínica.


```r
known_markers <- c(
    "ENSG00000242110",
    "ENSG00000086205",
    "ENSG00000142515"
)

metrics_diagnostic_fpr$data %>%
    filter(identifier %in% known_markers)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000086205 </td>
   <td style="text-align:right;"> 0.8114848 </td>
   <td style="text-align:right;"> 0.0414879 </td>
   <td style="text-align:right;"> 0.6920416 </td>
   <td style="text-align:right;"> 0.8457328 </td>
   <td style="text-align:left;"> Partially Proper </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000142515 </td>
   <td style="text-align:right;"> 0.6494790 </td>
   <td style="text-align:right;"> 0.0135228 </td>
   <td style="text-align:right;"> 0.5448570 </td>
   <td style="text-align:right;"> 0.7856116 </td>
   <td style="text-align:left;"> Partially Proper </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000242110 </td>
   <td style="text-align:right;"> 0.9157217 </td>
   <td style="text-align:right;"> 0.0672311 </td>
   <td style="text-align:right;"> 0.8275320 </td>
   <td style="text-align:right;"> 0.8231707 </td>
   <td style="text-align:left;"> Concave </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000142515
) + 
    add_roc_curve(roc_data, response = disease, predictor = ENSG00000242110) + 
    add_roc_curve(roc_data, response = disease, predictor = ENSG00000086205) + 
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Marcadores de diagnóstico-1.png" width="672" />

Tras representar las curvas $ROC$ de los clasificadores se puede observar como *KLK3* no presenta un buen rendimiento a nivel global. Sin embargo, este predictor presenta mejor rendimiento en la región de interés, lo cual puede explicar su uso hasta la fecha. A pesar de ello, como se ha explicado previamente, *KLK3* presenta un gran número de falsos positivos, lo cual, explica que otros predictores como *FOLH1* y *AMACR* presentan mejores puntuaciones en la región de interés.


```r
plot_partial_roc_curve(
    roc_data,
    response = disease,
    predictor = ENSG00000142515,
    ratio = "fpr",
    threshold = 0.1
) + 
    add_partial_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000242110,
        ratio = "fpr",
        threshold = 0.1
    ) + 
    add_partial_roc_curve(
        roc_data,
        response = disease,
        predictor = ENSG00000086205,
        ratio = "fpr",
        threshold = 0.1
    ) + 
    add_tpauc_concave_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000242110,
        upper_threshold = 0.1,
        lower_threshold = 0
    ) +
    add_tpauc_partially_proper_lower_bound(
        roc_data,
        response = disease,
        predictor = ENSG00000086205,
        upper_threshold = 0.1,
        lower_threshold = 0
    ) +
    add_chance_line() +
    add_threshold_line(0.1, "fpr")
```

<img src="05_AUC_indexes_files/figure-html/Marcadores de diagnóstico en alta especificidad-1.png" width="672" />

Si nos centramos en esta zona, el $SpAUC$ califica al *AMACR* como el predictor de mejor rendimiento ($SpAUC = 0.828$), pero además, puntúa relativamente bajo al resto, es decir a *KLK3* y *FOLH1* ($SpAUC = 0.545$ y $SpAUC = 0.692$). Por otra parte, el $TpAUC$ califica a *FOLH1* como el predictor de mayor rendimiento ($TpAUC = 0.846$) y al *AMACR* como valores muy similares ($TpAUC = 0.823$).

Contrastando estos resultados con otros estudios, tanto el *AMACR* como el *FOLH1* se encuentran sobrexpresados en pacientes con la enfermedad [@amacr-marker; @folh1-marker]. Sin embargo, *FOLH1* presenta un mayor interés, por además de su uso a nivel de diagnóstico, como diana terapéutica [@folh1-therapeutic]. Esto refuerza la interpretación dada por $TpAUC$, donde tanto *AMACR* y *FOLH1* presentan potencial como clasificadores, y en concreto este último, un rendimiento ligermante superior.

Si se selecciona alguno los primeros marcadores identificados por el $TpAUC$, puede observarse que se encuentran relacionados con la enfermedad, por ejemplo *ESRP2* (`ENSG00000103067`), al cual estudios recientes lo describen por tener un papel relevante en el desarrollo de la enfermedad [@esrp2-marker].


```r
metrics_diagnostic_fpr$data %>%
    filter(identifier == "ENSG00000103067") %>%
    select(-curve_shape)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000103067 </td>
   <td style="text-align:right;"> 0.8339335 </td>
   <td style="text-align:right;"> 0.0499464 </td>
   <td style="text-align:right;"> 0.7365598 </td>
   <td style="text-align:right;"> 0.9497378 </td>
  </tr>
</tbody>
</table>

Por otra parte, *KLK3* presenta un rendimiento peor que un clasificador aleatorio para distinguir entre pacientes con una enfermedad más o menos agresiva.


```r
inner_join(
    metric_prognostic_fpr$data %>%
        filter(identifier == "ENSG00000142515") %>%
        select(-curve_shape),
    metric_prognostic_tpr$data %>%
        filter(identifier == "ENSG00000142515") %>%
        select(-curve_shape),
    join_by(identifier, auc),
    suffix = c("_spec", "_sens")
)
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc_spec </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:right;"> pauc_sens </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000142515 </td>
   <td style="text-align:right;"> 0.3177106 </td>
   <td style="text-align:right;"> 0.0008106 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.779661 </td>
   <td style="text-align:right;"> 0.0011463 </td>
   <td style="text-align:right;"> 0.0114632 </td>
   <td style="text-align:right;"> 0.6690821 </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    prognostic_dataset %>% mutate(prognostic = disease$prognostic),
    response = prognostic,
    predictor = ENSG00000142515
) +
    add_chance_line()
```

<img src="05_AUC_indexes_files/figure-html/ROC KLK3 estratificación de riesgo-1.png" width="672" />

La curva $ROC$ del predictor queda definida por debajo de la línea de azar, confirmando un rendimiento peor que un clasificador completamente aleatorio. Cabe destacar, que el cálculo de la curva se ha realizado considerando como condición de interés padecer una forma de la enfermedad más agresiva, y por tanto es posible que el *KLK3* rinda mejor identificando formas menos agresivas de la enfermedad. Para comprobarlo, es posible invertir la forma de la curva, lo cual resulta equivalente a invertir la condición de interés. 


```r
inverse_prognostic_dataset <- prognostic_dataset %>%
    mutate(
        prognostic = disease$prognostic
    )

inverse_prognostic_dataset$prognostic <- fct_recode(
    inverse_prognostic_dataset$prognostic,
    "1" = "0",
    "0" = "1"
)

inverse_prognostic_dataset$prognostic <- fct_relevel(
    inverse_prognostic_dataset$prognostic,
    "0"
)

tpr_metrics <- summarize_predictor(
    inverse_prognostic_dataset,
    predictor = ENSG00000142515,
    response = prognostic,
    ratio = "tpr",
    threshold = 0.9
)
fpr_metrics <- summarize_predictor(
    inverse_prognostic_dataset,
    predictor = ENSG00000142515,
    response = prognostic,
    ratio = "fpr",
    threshold = 0.1
)


inner_join(
    fpr_metrics,
    tpr_metrics,
    join_by(auc),
    suffix = c("_spec", "_sens")
) %>%
    mutate(
        identifier = "ENSG00000142515",
        .before = 1
    ) %>%
    select(
        !starts_with("curve_shape")
    )
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc_spec </th>
   <th style="text-align:right;"> sp_auc </th>
   <th style="text-align:right;"> tp_auc </th>
   <th style="text-align:right;"> pauc_sens </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000142515 </td>
   <td style="text-align:right;"> 0.6822894 </td>
   <td style="text-align:right;"> 0.0151036 </td>
   <td style="text-align:right;"> 0.5531767 </td>
   <td style="text-align:right;"> 0.7319498 </td>
   <td style="text-align:right;"> 0.0152133 </td>
   <td style="text-align:right;"> 0.152133 </td>
   <td style="text-align:right;"> 0.7315611 </td>
  </tr>
</tbody>
</table>


```r
plot_roc_curve(
    inverse_prognostic_dataset,
    response = prognostic,
    predictor = ENSG00000142515
) +
    add_chance_line()
```

<img src="05_AUC_indexes_files/figure-html/ROC inversa de KLK3-1.png" width="672" />

Al realizar la inversión, las diversas métricas presentan ciertas mejoras a excepción del $TpAUC$ que se ve ligeramente reducido. A pesar de ello, las métricas obtenidas distan de ser ideales para el ámbito clínico.

Para realizar un mejor contraste, estas métricas se compararán con otros predictores usados en tests de estratificación de riesgo [@risk-tests], como por ejemplo: Decipher [@decipher] y Oncotype [@oncotype].


```r
decipher <- tibble::tribble(
    ~name, ~ensemble,
    "LASP1", "ENSG00000002834",
    "IQGAP3", "ENSG00000183856",
    "NFIB", "ENSG00000147862",
    "S1PR4", "ENSG00000125910",
    "THBS2", "ENSG00000186340",
    "ANO7", "ENSG00000146205",
    "PCDH7", "ENSG00000169851",
    "MYBPC1", "ENSG00000196091",
    "EPPK1", "ENSG00000261150",
    "PBX1", "ENSG00000185630",
    "NUSAP1", "ENSG00000137804",
    "ZWILCH", "ENSG00000174442",
    "UBE2C", "ENSG00000175063",
    "CAMK2N1", "ENSG00000162545",
    "RABGAP1", "ENSG00000011454",
    "TNFRSF19", "ENSG00000127863"
)

oncotype <- tibble::tribble(
    ~name, ~ensemble,
    "BGN", "ENSG00000182492",
    "COL1A1", "ENSG00000108821",
    "SFRP4", "ENSG00000106483",
    "FLNC", "ENSG00000128591",
    "GSN", "ENSG00000148180",
    "TPM2", "ENSG00000198467",
    "GSTM2", "ENSG00000213366",
    "TPX2", "ENSG00000088325",
    "LAMB3", "ENSG00000196878",
    "FAM13C", "ENSG00000148541",
    "KLK2", "ENSG00000167751",
    "AZGP1", "ENSG00000160862",
    "SRD5A2", "ENSG00000277893",
    "DUSP1", "ENSG00000120129",
    "FOS", "ENSG00000170345"
)

metric_prognostic_tpr$data %>%
    arrange(desc(fp_auc)) %>%
    filter(
        identifier %in% oncotype$ensemble |
        identifier %in% decipher$ensemble
    )
```
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> identifier </th>
   <th style="text-align:right;"> auc </th>
   <th style="text-align:right;"> pauc </th>
   <th style="text-align:right;"> np_auc </th>
   <th style="text-align:right;"> fp_auc </th>
   <th style="text-align:left;"> curve_shape </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000196091 </td>
   <td style="text-align:right;"> 0.3068206 </td>
   <td style="text-align:right;"> 0.0008597 </td>
   <td style="text-align:right;"> 0.0085974 </td>
   <td style="text-align:right;"> 0.7536232 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000127863 </td>
   <td style="text-align:right;"> 0.4276427 </td>
   <td style="text-align:right;"> 0.0015230 </td>
   <td style="text-align:right;"> 0.0152297 </td>
   <td style="text-align:right;"> 0.7246377 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000146205 </td>
   <td style="text-align:right;"> 0.3072628 </td>
   <td style="text-align:right;"> 0.0005371 </td>
   <td style="text-align:right;"> 0.0053713 </td>
   <td style="text-align:right;"> 0.6980676 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000160862 </td>
   <td style="text-align:right;"> 0.3320233 </td>
   <td style="text-align:right;"> 0.0002522 </td>
   <td style="text-align:right;"> 0.0025219 </td>
   <td style="text-align:right;"> 0.6859903 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000167751 </td>
   <td style="text-align:right;"> 0.3659052 </td>
   <td style="text-align:right;"> 0.0018865 </td>
   <td style="text-align:right;"> 0.0188651 </td>
   <td style="text-align:right;"> 0.6739130 </td>
   <td style="text-align:left;"> Hook under chance </td>
  </tr>
</tbody>
</table>

Entre los datos, es posible encontrar marcadores utilizados en estos tests, los cuales presentan puntuaciones relativamente bajas, a excepción del $FpAUC$ que las puntúa en mejor medida. Contrastando estos resultados, se puede observar como el $FpAUC$ es más representativo de los marcadores actualmente utilizados en tests de estratificación de riesgos. Es más, entre los predictores con mayor $FpAUC$ se encuentra `ENSG00000112210` (*RAB23*), el cual se encuentra relacionado con la progresión de la enfermedad [@rab23-marker].

## Discusión y conclusiones

El desarrollo de nuevas tecnologías permite el estudio simultáneo de un gran número de marcadores, donde los análisis basados en $AUC$ son habitualmente utilizados. En el ámbito clínico, regiones concretas de las curvas $ROC$ son más relevantes, principalmente para tratar con los costes y consecuencias de diagnósticos erróneos. De esta forma, las métricas que evalúan estas regiones, como el $pAUC$, están ganando mayor interés. A pesar de ello, esta métrica presenta limitaciones para su uso debido a su interpretabilidad, por lo que transformaciones alternativas para trabajar en estas regiones resultan de mayor utilidad.

En este trabajo se ha analizado el rendimiento de métricas derivadas del $pAUC$ en regiones de interés, para ello se ha utilizado un conjunto de datos generado con distintas muestras de cáncer de próstata. Las diferentes métricas han sido evaluadas en función de su capacidad para encontrar marcadores de diagnóstico y estratificación de riesgo, escenarios donde condiciones de alta especificidad y sensibilidad son requeridas.

Entre las métricas descritas, el $SpAUC$ fue la primera transformación del $pAUC$ descrita, la cual permite trabajar en regiones de alta especificidad. Esta métrica es fácilmente interpretable al comprender sus valores entre $0.5$ y $1$, sin embargo, no es aplicable a regiones de la curva $ROC$ descritas por debajo de la línea de azar. Estas regiones son comunes en la práctica, lo cual limita su uso en gran medida. Este problema ha sido resuelto por el $TpAUC$ la cual, además de tomar valores entre $0.5$ y $1$, puede aplicarse a regiones que pasen por debajo de la línea del azar.

De forma similar, el $NpAUC$ fue propuesto para trabajar en regiones de alta sensibilidad. Esta métrica puede tomar valores máximos de $1$, pero al no estar limitada por debajo, también puede tomar valores inferiores a $0.5$, lo que dificulta su interpretación y aplicación a ciertos predictores. Estos problemas han sido resueltos por el $FpAUC$ el cual, además de contener sus valores entre $0.5$ y $1$, se puede aplicar independientemente de la forma de la curva.

Además de presentar una mejor interpretabilidad, el $TpAUC$ y $FpAUC$ han sido capaces de distinguir entre predictores que presentan el mismo $pAUC$ y cuyas curvas $ROC$ se cruzan, algo que hasta ahora no había sido posible con las anteriores métricas. De esta forma, tanto el $TpAUC$ como el $FpAUC$ poseen una mayor capacidad de discriminación, permitiendo la selección de otros marcadores que convencionalmente se habrían descartado.

Finalmente, para considerar la relevancia clínica de las métricas, los resultados obtenidos se han contrastado con predictores relevantes en el ámbito clínico. Con respecto a marcadores de diagnóstico, el $TpAUC$ otorga puntuaciones similares a los 2 marcadores de mayor interés actualmente, a diferencia del $SpAUC$ que penaliza a uno de ellos. En condiciones de estratificación de riesgo, el $NpAUC$ presenta valores muy bajos para predictores que actualmente se utilizan en tests de estratificación de riesgo, a diferencia del $FpAUC$ que les otorga puntuaciones mayores. 

Estos resultados revelan la utilidad práctica del $TpAUC$ y $FpAUC$ frente al resto. Su facilidad para ser interpretadas y capacidad de discriminar entre marcadores altamente específicos o sensibles incentiva su uso en el ámbito clínico, permitiendo la selección de marcadores que puedan ayudar en la toma de decisiones.
