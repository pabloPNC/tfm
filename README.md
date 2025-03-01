# Descripción

Este repositorio contiene los scripts utilizados para realizar el análisis de diversos índices de áreas parciales, los cuales son utilizados como métricas del rendimiento de un clasificador en condiciones de alta sensibilidad o especificidad.

# Contenidos

En el repositorio se encuentran disponibles una serie de archivos que han sido utilizados para llevar a cabo este análisis.

En primer lugar, se pueden encontrar una serie de archivos `*.Rmd`, los cuales contienen descritas las distintas secciones del análisis completo además del código fuente en R utilizado durante el trabajo.

Para facilitar su lectura, estos documentos, así como los resultados del código utilizado, han sido compilados en archivos `*.html` los cuales se pueden encontrar bajo el directorio `_book/`. 

Finalmente, en el repositorio también se encuentra disponible un script en python (`tcga_metadata_processing.py`). Este pequeño script contiene el código utilizado para procesar las variables clínicas de muestras extraídas de *The Cancer Genome Atlas* (TCGA) durante el trabajo.

# Estructura del análisis

Este análisis tiene como finalidad estudiar distintos índices de áreas parciales como métricas del rendimiento de un clasificador, además de identificar cual de ellos presenta un mejor desempeño identificando biomarcadores. Así, los distintos pasos seguidos en el análisis se encuentran divididos en 3 secciones principales:

1. **Búsqueda de marcadores**. Utilizando datos de experimentos de scRNA-seq, en esta sección se identifican células de tejido responsable de la enfermedad de estudio, cáncer de próstata, y se seleccionan marcadores genéticos que puedan estar relacionadas con ella.
1. **Generación de una matriz de cuentas**. En esta sección se genera el conjunto de datos sobre el que analizar los distintos índices de áreas parciales. Para ello, se extraen perfiles de expresión de muestras del TCGA filtrados con los marcadores identificados en la sección anterior.
1. **Evaluación de índices**. En esta última sección se describe el fundamento matemático de los distintos índices a evaluar. Finalmente, estos son aplicados al conjunto de datos y sus resultados son contrastados con marcadores ya conocidos en la práctica.

