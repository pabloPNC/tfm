# (PART) Generación del conjunto de datos {-}

# Generación de la matriz de conteo {#count-matrix-generation}

En esta sección se trabajará con datos de experimentos de RNA-seq, los cuales han sido extraídos de TCGA, para generar una matriz de conteo con genes expresados en la enfermedad. Esta matriz será procesada en pasos posteriores para finalmente ser utilizada para evaluar el rendimiento de los diversos índices.



## Descripción de los datos

A la hora de trabajar con ellos, los datos descargados se encontrarán en el directorio `./tcga_data/`. En este se podrán encontrar distintos subdirectorios y ficheros que se utilizarán durante el procesamiento de los datos:

+ **clinical**\. Directorio que contiene archivos de metadatos sobre los donantes de tejido. En su interior se encuentra el fichero `clinical.tsv` con variables de interés para el trabajo.
+ **gdc_download**\. Directorio que contiene, a su vez, distintos subdirectorios, cada uno de ellos con un archivo de lecturas para una muestra de tejido.
+ **gdc_sample_sheet**\. Fichero `*.tsv` que contiene información sobre los archivos descargados. En su mayoría, contiene identificadores de las muestras y archivos descargados, tipo de experimento, tipo de tumor, etc. Mayoritariamente es de utilidad para relacionar los archivos de metadatos con su respectiva muestra.

### Preprocesamiento del fichero clinical y gdc_sample_sheet {#clinical-preprocessing}

Dado que se trabaja con un gran número de ficheros de distinto tipo (datos clínicos, metadatos y datos de experimentos), es posible simplificar la tarea uniendo los ficheros de metadatos (`gdc_sample_sheet.tsv`) y datos clínicos (`clinical.tsv`) en un solo.

Para esta tarea, se ha creado un script en python (`tcga_metadata_processing.py`)[^2]. En concreto, los pasos realizados por el script son los siguientes:

1. Lectura de los archivos
1. Renombrado de las columnas para usar un estilo común
    + Cambio de espacios por `_`
    + Uso de minúsculas
1. Selección de variables de interés y muestras que las describan
    + Puntuación de Gleason primaria
    + Puntuación de Gleason secundaria
1. Fusión de archivos utilizando identificadores del paciente
    + `case_id`
    + `case_submitter_id`
1. Transformación de variables
    + Transformación a tipo numérico
    + Cálculo de nuevas variables
        + Puntuación de Gleason[^3]

Realizados estos pasos, el script genera el archivo `file_clinical_sheet.csv`, el cual contiene la información necesaria de cada archivo y las variables clínicas de cada paciente del que proviene la muestra. Finalmente, para trabajar con estos archivos, se guardaran las rutas en las que se encuentran.


```r
tcga_data_folder <- "./tcga_data/gdc_download_20240611_200434.803317"
file_clinical_path <- "./file_clinical_sheet.csv"
```

## Creación de matriz de conteo

Como se ha descrito en el apartado anterior, se está trabajando con un gran número de archivos de lecturas. De forma similar a metadatos y datos clínicos, los archivos de lecturas pueden combinarse para facilitar su acceso a ellos, pero para realizar este proceso se ha de entender la estructura que presentan.

Los archivos de lecturas con los que se está trabajando destinan sus 6 primeras líneas a una serie de métricas y nombres de las variables. Tras ellas, el resto de filas representan un gen con los valores de dichas variables.

```
## gene-model: GENCODE v36
gene_id	gene_name	gene_type	unstranded	stranded_first	stranded_second	tpm_unstranded	fpkm_unstranded	fpkm_uq_unstranded
N_unmapped			1242436	1242436	1242436			
N_multimapping			2589737	2589737	2589737			
N_noFeature			1331220	16814349	16833920			
N_ambiguous			3211136	753902	752634			
ENSG00000000003.15	TSPAN6	protein_coding	2743	1379	1364	72.3539	20.8836	20.1004
ENSG00000000005.6	TNMD	protein_coding	12	5	7	0.9728	0.2808	0.2702
ENSG00000000419.13	DPM1	protein_coding	745	398	347	73.8513	21.3158	20.5164
```
### Importar archivos de lecturas {#import-read-files}

Conociendo la estructura de los archivos es posible planificar una forma de automatizar el proceso para todos ellos. En primer lugar, se deberán de leer y para ello se define la función `to_tsv`. La función simplemente lee el archivo en la ruta especificada (`file_path`) omitiendo las 6 primeras lineas y utilizando un tabulador como elemento de separación. Tras ello, re-asigna el nombre de las columnas y devuelve la matriz generada.



```r
to_tsv <- function(file_path) {
    column_names <- c("gene_id", "gene_name", "gene_type", "unstranded",
                      "stranded_first", "stranded_second", "tpm_unstranded",
                      "fpkm_unstranded","fpkm_uq_unstranded")
    counts <- read_tsv(file_path, skip = 6, col_names = FALSE, show_col_types = FALSE)
    colnames(counts) <- column_names
    return(counts)
}
```

### Selección de variables y renombrado {#process-read-files}

Una vez es posible leer los conjuntos de datos es posible decidir que variables utilizar. Como se ha observado antes, los archivos de lecturas presentan las siguientes variables: 

+ **gene_id** y **gene_name**\. Indican el acrónimo para el gen y su identificador.
+ **gene_type**\. Pequeña descripción del tipo de gen con el que se trata.
+ **unstranded**\. Número de lecturas que se alinean con la región de un gen, sin tener en cuenta la hebra de ADN en la que se encuentra.
+ **stranded_first** y **stranded_second**\. Número de lecturas que se alinean a la región de una hebra de ADN o su complementaria.
+ **tpm_unstranded**, **fpkm_unstranded** y **fpkm_uq_unstranded**\. Diversas normalizaciones de las lecturas en crudo. Respectivamente: transcritos por millón, fragmentos por kilobase por millón de lecturas y cuartil superior de FPKM.

Al realizar análisis resulta más conveniente trabajar con identificadores frente a nombres, por lo que estos serán seleccionados. Con respecto a las métricas, se plantea la cuestión de si utilizar lecturas `unstranded` o `paired` (`stranded_first` y `stranded_second`). Habitualmente las segundas se suelen utilizar para experimentos donde es relevante conocer la hebra de la que proviene un transcrito, por ejemplo para el descubrimiento *de novo* de transcritos o estudio de isoformas [@best-practices], por este motivo se escogerán las primeras.

La segunda cuestión que se plantea es la normalización de los datos que se debe usar. Las lecturas normalizadas son de utilidad cuando se compara la expresión de distintos genes, ya que pueden normalizar parámetros que influyen en esta expresión, e.g., longitud del transcrito, tamaño de librería, etc. Puesto que se trabaja con muestras de distintos pacientes, resulta adecuado utilizar la métrica de transcritos por millón (TPM). Esta métrica normaliza las diferencias en la composición de transcritos en el denominador en lugar de simplemente dividir por el número de lecturas en la librería [@best-practices], como se observa en \@ref(eq:normalization-metrics).

$$
C_g = \text{Número de lecturas alineadas al gen } g \\
L_g = \text{Longitud del gen } g \\
N = \text{Número total de genes}
$$

\begin{equation}
    FPKM = \dfrac{C_g*10^9}{(\sum^N_{i=1} C_i)L_g} \\
    TPM = \dfrac{(C_g*10^3/L_g)*10^6}{\sum^N_{i=1}(C_i*10^3/L_i)}
    (\#eq:normalization-metrics)
\end{equation}

Identificadas las variables a utilizar, es posible definir una función para automatizar el proceso. La función `rename_col` utilizará como entrada una matriz de conteo (`data`) y un nombre que se le provea (`sample_id`), con ellos renombrará la columna `tpm_unstranded` al nombre especificado, más adelante se utilizará para indicar el identificador de la muestra. Finalmente, seleccionará dicha columna junto con la que contiene los identificadores de los genes.


```r
rename_col <- function(data, sample_id) {
    count <- data %>%
        dplyr::rename("{ sample_id }" := tpm_unstranded) %>%
        dplyr::select(gene_id, all_of(sample_id))
    return(count)
}
```

### Unión de archivos

Estando automatizada la lectura de archivos y selección de variables de interés, solamente falta automatizar el proceso para unir estos archivos. Para esta tarea, se define `join_count_matrices`. 

La función tomará como entrada 2 matrices de conteo (`count_matrix_a` y `count_matrix_b`) y utilizará los identificadores de los genes (`gene_id`) para realizar su unión. De esta forma, si hay genes identificados en la segunda matriz que no están en la primera, se añadirá una nueva fila para dicho gen. En caso de que un gen exista en la primera matriz, las cuentas se colocarán en la fila adecuada. De esta forma la función devolverá una matriz donde cada fila se corresponderá a un gen concreto, mientras que cada columna a una muestra.


```r
join_count_matrices <- function(count_matrix_a, count_matrix_b) {
    full_join(
        count_matrix_a,
        count_matrix_b,
        by = c("gene_id"), keep = FALSE
    )
}
```

### Construcción de la matriz de conteo

Creadas las funciones para el flujo de trabajo, simplemente basta con poder aplicarlas a cada uno de los archivos de lecturas. Para que la tarea sea más sencilla, se realizará en primer lugar el preprocesado de los archivos y tras ello su unión.

En primer lugar se define la función `transform_counts`, la cual simplemente aplica las funciones `to_tsv` y `rename_col`, definidas en las secciones \@ref(import-read-files) y \@ref(process-read-files).


```r
transform_counts <- function(file_path, sample_id) {
    count_matrix <- to_tsv(file_path) %>%
        rename_col(sample_id)
    return(count_matrix)
}
```

La función toma como entrada la ruta de un archivo de lecturas y un identificador de muestra, por lo que para aplicar esta función a todos los archivos, se requerirá iterar sobre cada una de las rutas y sus identificadores. Esta información, está disponible en el archivo `file_clinical_sheet.csv` previamente generado (sección \@ref(clinical-preprocessing)), por lo que se importará en primer lugar. 

Tras ello utilizando la función `stringr::str_c` es posible concatenar la ruta del directorio (`tcga_data_folder`) con el nombre de cada uno de los archivos (`file_id` y `file_name`) obteniendo así la ruta completa de cada uno de ellos. Estas rutas, se guardarán en una nueva columna (`full_path`) creada utilizando la función `dplyr::mutate`.

Teniendo las rutas completas e identificadores, se usará la función `purrr::map2` para iterar simultáneamente sobre 2 inputs y pasarlos en una función. En este caso, se aplicará sobre las rutas (`file_clinical$full_path`) y los identificadores de las muestras (`file_clinical$sample_id`) para pasarlos a la función `transform_counts`. De esta forma se obtendrá una lista de archivos de lecturas procesados.


```r
file_clinical <- as.data.frame(read_csv(file_clinical_path, col_names = TRUE))
file_clinical <- file_clinical %>%
    mutate(
        full_path = str_c(tcga_data_folder, file_id, file_name, sep = "/")
    )

count_matrices <- map2(
    file_clinical$full_path,
    file_clinical$sample_id,
    transform_counts,
    .progress = T
)
```

Preprocesados los archivos de lecturas, el siguiente pasó será realizar la unión de todos archivos. Para ello se utilizará la función `purrr:reduce` la cual combina los elementos de un vector en un valor único. De forma detallada, `reduce` actúa de la siguiente forma:

1. Toma una función como input que es aplicada a los 2 primeros elementos de un vector.
1. La función se aplica al nuevo resultado y al siguiente elemento del vector.
1. El paso anterior se repite hasta que no queden más elementos en el vector.

De esta forma, cuando se aplica `reduce` con la función `join_count_matrices`, se realizarán uniones de los ficheros de forma consecutiva hasta obtener una única matriz de conteo.


```r
merged_matrices <- as.data.frame(purrr::reduce(count_matrices, join_count_matrices))
rownames(merged_matrices) <- merged_matrices$gene_id
```

## Procesado de la matriz de conteo

Aunque en el paso anterior se ha generado una matriz de conteo, se debe tener en cuenta los posteriores análisis a realizar. En ellos, se deberán realizar una serie de cálculos sobre la matriz, por ello, será necesario que esta cumpla con algunas condiciones:

1. Los genes o biomarcadores deberán encontrarse en las columnas de la matriz, dado que serán las variables de estudio.
1. La matriz debe de contener una variable binaria que represente la variable a predecir e.g. si la muestra proviene de un paciente enfermo o no.

### Pivotado de la matriz

Como se ha descrito al comienzo de la sección, la estructura de los datos es relevante para realizar los análisis aguas abajo. Por este motivo, en esta sección se realizará una reestructuración de los mismos utilizando las funciones `tidyr::pivot_longer` y `tidyr::pivot_wider`. 

La primera, compactará todas la columnas especificadas (`cols = !gene_id`, es decir todas menos `gene_id`) en una nueva columna (`names_to = "sample"`) y su valores en otra especificada (`values_to = "counts"`). Compactada la información de las muestras, podemos utilizar la segunda para crear nuevas columnas con el nombre de los genes (`names_from = gene_id`) que usen los valores previamente guardados (`values_from = counts`).


```r
tcga_data <- merged_matrices %>%
    pivot_longer(
        cols = !gene_id,
        names_to = "sample",
        values_to = "counts"
    ) %>%
    pivot_wider(
        names_from = gene_id,
        values_from = counts
    )
```

De esta forma, se logra que los genes se correspondan con las columnas y cada entrada con la muestra de un paciente.


|sample           | ENSG00000000003.15| ENSG00000000005.6| ENSG00000000419.13|
|:----------------|------------------:|-----------------:|------------------:|
|TCGA-KK-A8IG-01A |            29.0736|            0.2262|            76.9529|
|TCGA-G9-A9S7-01A |            50.6646|            0.0459|           106.8449|
|TCGA-EJ-7789-11A |           111.8758|            1.2765|           114.1114|
|TCGA-YL-A8SA-01A |            51.1902|            0.0597|            90.8608|
|TCGA-XK-AAJ3-01A |            36.0334|            0.0000|            63.9224|

### Cálculo de Gold Standards

Una vez se ha reestructurado la matriz es posible incluir las variables a predecir, que se usarán en posteriores análisis. Estas variables, representarán patrones que suceden en la realidad o también conocidos como "Gold Standards". En este caso los patrones que se incluirán en los datos serán variables binarias que, posteriormente, permitirán evaluar los genes seleccionados como predictores.

Estas variables se establecerán a partir de la puntuación de Gleason, por lo que, el primer paso será incluirla en la matriz. Dado que el orden de las muestras es el mismo en la matriz que en la tabla de metadatos, es posible incluir la columna sin mayor dificultad.


```r
tcga_data[["gleason_score"]] <- file_clinical[["gleason_score"]]
```

Incluidas las puntuaciones de Gleason es posible generar las demás variables. En primer lugar se agruparán las muestras en función de la malignidad del tumor, es decir en función de su puntuación de Gleason:

+ **Puntuación de 5 o menor**\. Se podría considerar un cáncer de bajo grado. En el conjunto de datos utilizado, todos los tejido en este rango serán tejidos sanos (e.g. hiperplásias)
+ **Puntuación de 6 o 7**\. Cáncer de grado medio. En el conjunto de datos serán referidos como casos de malignidad baja.
+ **Puntuación de 8 o mayor**\. Cáncer de grado alto. En el conjunto de datos serán referidos como casos de malignidad alta.

Con la función `dplyr::mutate` se creará la variable `malignancy` que ,utilizando la función `dplyr::case_when`, clasificará las distintas muestras en `Normal`, `Low Malignancy` y `High Malignancy` según los criterios establecidos. 


```r
tcga_data <- tcga_data %>%
    mutate(
        malignancy = case_when(
            gleason_score >= 8 ~ "High Malignancy",
            gleason_score == 6 | gleason_score == 7 ~ "Low Malignancy",
            gleason_score <= 5 ~ "Normal"
        )
    )
```

Es importante considerar que en pasos posteriores se deberá trabajar con variables binarias, por ello que se requerirá realizar transformaciones adicionales. En primer lugar, se buscará trabajar con biomarcadores de diagnóstico. Los biomarcadores de este tipo deberían ser capaces de diferenciar o clasificar entre pacientes sanos y enfermos, por ello, se utilizará la información de la variable `malignancy` para generar una variable que represente estas condiciones. Esta variable será `disease` que, utilizando la función `forcats::fct_collapse`, compactará los valores de `malignancy` de `High Malignancy` y `Low Malignancy` en 1 (enfermos), y valores de `Normal` en 0 (sanos). 


```r
tcga_data <- tcga_data %>%
    mutate(
        disease = fct_collapse(
            malignancy,
            "0" = "Normal", 
            "1" = c("High Malignancy", "Low Malignancy")
        )
    )
```

En segundo lugar, se tratará de buscar biomarcadores de pronóstico o estratificación de riesgo. Los biomarcadores de este tipo deberían ser capaces de distinguir entre pacientes donde la enfermedad es menos probable que progrese frente a pacientes donde es más probable que lo haga, es decir un cáncer de crecimiento lento y menos agresivo frente a uno de crecimiento rápido y mayor agresividad. En este caso, se realizará un procedimiento similar creando la variable `prognostic`, donde 0 se corresponderá con valores de `Low Malignancy` y de 1 con `High Malignancy`[^4].


```r
tcga_data <- tcga_data %>%
    mutate(
        prognostic = fct_collapse(
            malignancy,
            "0" = "Low Malignancy",
            "1" = "High Malignancy"
        )
    )
```

Una vez generada la matriz con los genes y variables a predecir, el conjunto de datos está preparado para realizar los análisis pertinentes. Sin embargo, dado que anteriormente se han seleccionado biomarcadores relevantes de la enfermedad (sección \@ref(cancer-cells-selection)), será posible utilizarlos para filtrar el conjunto de datos, de esta forma se logra:

+ **Aumento de precisión**\. Los experimentos de RNA-seq muestran perfiles transcripcionales promedios de todas las células de la muestra. Esto incluye células epiteliales, fibroblastos, células musculares, etc. Por este motivo, realizar un filtrado con los biomarcadores seleccionados en experimentos de scRNA-seq, permite seleccionar solamente aquellas células de interés para la enfermedad.
+ **Reducción del tamaño del conjunto de datos**\. Los análisis aguas abajo son exigentes a nivel computacional, por lo que una reducción del conjunto de datos acelera en gran medida el proceso.

Así, en pasos posteriores, se realizará este filtrado usando los datos provenientes de experimentos de scRNA-seq, para finalmente realizar el análisis de los biomarcadores seleccionados de forma más óptima.



[^2]: El uso de un script de Python para esta tarea ha sido completamente arbitrario. Para algunos de los pasos este lenguaje utiliza una sintaxis más sencilla, sin embargo, estos también podrían haberse realizado con R.

[^3]: La puntuación de Gleason se obtiene a través de la suma de la puntuación primaria y la secundaria:
    $$\text{P. Gleason} = \text{P. Gleason Primaria} + \text{P. Gleason Secundaria}$$

[^4]: En la creación de la variable `prognostic` solamente se deben considerar pacientes enfermos. Sin embargo, debido a que los datos provienen de la variable `malignancy` es posible encontrar encontrar casos de tejido sano, de esta forma la variable podrá tomar los valores: `Normal`, `0` y `1`). Para mantener los datos en una matriz común, la variable `prognostic`se mantendrá en esta forma, pero a la hora de realizar los cálculos las entradas con valor `Normal` serán filtradas.
