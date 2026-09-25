#' Carga de librerías
## ELIMINEN las verificacionees personales a medida que corren el codigo, no es parte del codigo en realidad (view, dim, length...)

library(SummarizedExperiment) #ESTA LIBRERIA LA USARON?
library(TCGAbiolinks) # LA USARON?
library(dplyr) # LA USARON??
library(edgeR)
library(limma)

expr <- read.delim(
  "7TK4C4-expression-matrix.tsv",
  header = TRUE,
  check.names = FALSE
)

#' Copia de la matriz
matrix_exp <- expr

#' Creación de matriz de conteos
counts_matrix <- matrix_exp[, grepl("_count$", colnames(matrix_exp))]
rownames(counts_matrix) <- matrix_exp$gene_id

## VISUALIZAR
View(counts_matrix)
dim(counts_matrix)
head(counts_matrix)

## POR QUE LA DOBLE COMILLA? ESO TE LO VUELVE COMENTARIO Y NO SE PUEDE EJECTUAR
## CREO QUE NO ES NECESARIA LA CPM, BORRAR ESTA SECCION SI NO LA USARON
## EN REALIDAD NO LA USARON ASI QUE BORRENLA Y LA CREAN EN LA SECCION QUE SEA UTIL PERO PARA 
## DGE NO SE USO, SOLO GENERARA CONFUSION, CPM SE USA EN HIPATHIA COLOCARLA CUANDO SE USE
# O BIEN DESPUES DE dge <- calcNormFactors(dge) Y COLOCAR UN COMENTARIO DE QUE SE VA A USAR LUEGO PARA HIPATHIA
'
"# <- ESTAS
#' Creación de matriz de cpm
cpm_matrix <- matrix_exp[, grepl("_cpm$", colnames(matrix_exp))]
rownames(cpm_matrix) <- matrix_exp$gene_id
View(cpm_matrix)
dim(cpm_matrix)
head(cpm_matrix)"  # <- ESTAS


# CPM MATRIX DE HIPATHIA
cpm_matrix <- cpm(
  dge,
  log = FALSE
)
'


#' Verificación de que counts sea igual al orden de las filas
all(colnames(counts_matrix) == 
      colnames(matrix_exp)[grepl("_count$", colnames(matrix_exp))])

# SE PRECISA EN LOS COMMENTS ESTE YA ES CODIGO OFICIAL NO COMO EN CLASE QUE DETALLAMOS MUCHO
# Guardar matriz de conteos 
saveRDS(counts_matrix, "matriz_conteos.rds")


#' Asignación de nombre de la variable a las réplicas trabajadas
colnames(counts_matrix)
## CAMBIO: antes de asignar las condiciones, comprobar cuántas muestras
## existen realmente en la matriz.
ncol(counts_matrix)
colnames(counts_matrix)


## ASIGNAR CONDICIONES 
## CUANDO NECESITAMOS VARIAS CONDICIONES O LISTAS SE VE MAL LISTAR TODO Y POCO REPRODUCIBLE, LO MEJOR
## ES HACER UN LOOP. EN ESTE CASO VEO QUE LA MATRIZ NO TIENE LOS NOMBRES DE LAS CONDICIONES ENTONCES PODEMOS HACER 
# MINI TABLA DE METADATA QUE SEA MAS AGRADABLE AL LECTOR
sample_info <- data.frame(
  sample = colnames(counts_matrix),
  condition = c(
    rep("vehicle_6", 3),
    rep("etoh_6", 6),
    "lipof_6",
    "control_6",
    rep("lipof_6", 4),
    rep("control_24", 5),
    "control_6",
    "control_24",
    rep("vehicle_24", 6),
    rep("etoh_24", 3),
    "control_6",
    rep("etoh_24", 3),
    rep("lipof_24", 6),
    rep("control_6", 3),
    rep("vehicle_6", 3)
  )
)

#### DEL CODIGO ANTERIOR ES MEJOR AUN SI LOS METADATOS ESTAN EN OTRO SHEET Y LO SUBIMOS A GITHUB ASI PODEMOS HACER:
sample_info <- read.csv("sample_metadata.csv")
sample_info$condition <- with(
  sample_info,
  paste(treatment, time, sep = "_")
)
## PERO UDS DECIDEN COMO PREFIEREN. SI QUIEREN DEJAR SU LISTA AL FINAL FUNCIONA SOLO NO ES BUENA PRACTICA.


## CAMBIO: en lugar de escribir una segunda lista de condiciones,
## se comprueba directamente que exista una condición por muestra.
nrow(groups)
ncol(counts_matrix)

stopifnot(nrow(groups) == ncol(counts_matrix))

## hay muchísimo código que se puede eliminar



# Confimación de número de muestras
## Verificar metadatos y matriz de conteos

# CAMBIO: se reemplazaron varias comprobaciones redundantes por una
# verificación directa de que el número y orden de las muestras coinciden.
stopifnot(
  nrow(groups) == ncol(counts_matrix),
  groups$sample == colnames(counts_matrix)
)

# Verificar valores faltantes y negativos
anyNA(counts_matrix)
any(counts_matrix < 0)




## Verificar datos y condiciones

is.numeric(counts_matrix)
table(groups$condition)
head(rownames(counts_matrix))


## Definir condiciones experimentales
groups$condition <- factor(groups$condition)
class(groups$condition) ## es interno en el codigo original lo podriamos quitar, opcional
levels(groups$condition) ## same

## Construir MATRIZ DE DISE˜NO no diseño de matriz
design <- model.matrix(~0 + condition, data = groups)
colnames(design) <- levels(groups$condition)



# Creación de las comparaciones a evaluar en edgeR
contrast <-limma::makeContrasts(
  etoh_6_vs_control_6 = etoh_6 - control_6,
  etoh_24_vs_control_24 = etoh_24 - control_24,
  lipof_6_vs_vehicle_6 = lipof_6 - vehicle_6,
  lipof_24_vs_vehicle_24 = lipof_24 - vehicle_24,
  vehicle_6_vs_control_6 = vehicle_6 - control_6,
  vehicle_24_vs_control_24 = vehicle_24 - control_24,
  etoh_6_vs_etoh_24 = etoh_6 - etoh_24,
  lipof_6_vs_lipof_24 = lipof_6 - lipof_24,
  
  levels = colnames(design)
  )


# Creación del objeto DGEList
dge <- DGEList(
  counts = counts_matrix,
  group = groups$condition
)


## Filtrar genes de baja expresión
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]


# Normalización con TMM
dge <- calcNormFactors(dge)
dge$samples


# Estimación de la dispersión
dge <- estimateDisp(dge, design)
dge$common.dispersion
dge$tagwise.dispersion[1:10]
plotBCV(dge)


# Ajustar el modelo estadístico
ajuste <- glmQLFit(dge, design)

# Verificación
dim(dge)
dge$samples$norm.factors


# Ejecución de comparaciones ## VALE LA PENA HACER UN LOOP NO ES BUENA PRACTICA EJECUTAR EL MISMO CODIGO MIL VECES
## Ejecutar los contrastes
qlf_results <- lapply(
  colnames(contrast),
  function(x) {
    glmQLFTest(
      ajuste,
      contrast = contrast[, x]
    )
  }
)

names(qlf_results) <- colnames(contrast)

# Obtención de los resultados de las comparaciones
## Obtener tablas de resultados
results <- lapply(
  qlf_results,
  function(x) topTags(x, n = Inf)$table
)

## Añadir anotación génica
annotation <- matrix_exp[, c("gene_id", "gene_name", "gene_biotype")]

## HACER UN LOOP! SE VE MEJOR, LIMITA ERRORES Y LO VUELVE REPRODUCIBLE
# EJEMPLO: YA TENIAN UN ERROR  NO COINCIDE ETOH24 CON ETOH6
# EN  
''# etoh 24
  result_etoh24$gene_name <- matrix_exp$gene_name[
  match(rownames(result_etoh6), matrix_exp$gene_id) ## OJO! ERROR AQUI NO COINCIDE ETOH24 CON ETOH6
]
result_etoh6$gene_name <- matrix_exp$gene_name[
  match(rownames(result_etoh6), matrix_exp$gene_id)
]
''

results <- lapply(
  results,
  function(result) {
    result$gene_name <- annotation$gene_name[
      match(rownames(result), annotation$gene_id)
    ]
    
    result$gene_biotype <- annotation$gene_biotype[
      match(rownames(result), annotation$gene_id)
    ]
    
    result
  }
)

# Selección de genes diferencialmente expresados
## SAME, NO CORRER MUCHAS VECES EL MISMO CODIGO
## Seleccionar genes diferencialmente expresados
DEGs <- lapply(
  results,
  function(result) {
    result[
      result$FDR < 0.05 &
      abs(result$logFC) >= 1,
    ]
  }
)


DEGs[["etoh_6_vs_control_6"]]
DEGs[["lipof_24_vs_vehicle_24"]]

dim(DEG_lipof) # 340 X 7


## NOTAS EXTRA
## ESTA bien el approach con EDGER solol recordar y conocer bien el metodo que están usando, la razon y las diferencias.
## El análisis se realiza con edgeR utilizando los conteos crudos,
## normalización TMM y el modelo quasi-likelihood (glmQLFit).
## limma se utiliza únicamente para definir los contrastes mediante
## makeContrasts(); por tanto, no es necesario aplicar voom().
