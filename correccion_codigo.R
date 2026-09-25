#' Carga de librerías
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
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
" # <- ESTAS
#' Creación de matriz de cpm
cpm_matrix <- matrix_exp[, grepl("_cpm$", colnames(matrix_exp))]
rownames(cpm_matrix) <- matrix_exp$gene_id
View(cpm_matrix)
dim(cpm_matrix)
head(cpm_matrix)"  # <- ESTAS


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
groups <- data.frame(
  sample = colnames(counts_matrix),
  condition = c(
    "vehicle_6",
    "vehicle_6",
    "vehicle_6",
    "etoh_6",
    "etoh_6",
    "etoh_6",
    "etoh_6",
    "etoh_6",
    "etoh_6",
    "lipof_6",
    "control_6",
    "lipof_6",
    "lipof_6",
    "lipof_6",
    "lipof_6",
    "control_24",
    "control_24",
    "control_24",
    "control_24",
    "control_24",
    "control_6",
    "control_24",
    "vehicle_24",
    "vehicle_24",
    "vehicle_24",
    "vehicle_24",
    "vehicle_24",
    "vehicle_24",
    "etoh_24",
    "etoh_24",
    "etoh_24",
    "control_6",
    "etoh_24",
    "etoh_24",
    "etoh_24",
    "lipof_24",
    "lipof_24",
    "lipof_24",
    "lipof_24",
    "lipof_24",
    "lipof_24",
    "control_6",
    "control_6",
    "control_6",
    "vehicle_6",
    "vehicle_6",
    "vehicle_6"
  )
)

## CAMBIO: en lugar de escribir una segunda lista de condiciones,
## se comprueba directamente que exista una condición por muestra.
nrow(groups)
ncol(counts_matrix)

stopifnot(nrow(groups) == ncol(counts_matrix))

# Confimación de número de muestras
length(colnames(counts_matrix))
length(c("Vehicle_6",
         "Vehicle_6",
         "Vehicle_6",
         "etoh_6",
         "etoh_6",
         "etoh_6",
         "etoh_6",
         "etoh_6",
         "etoh_6",
         "lipof_6",
         "control_6",
         "lipof_6",
         "lipof_6",
         "lipof_6",
         "lipof_6",
         "control_24",
         "control_24",
         "control_24",
         "control_24",
         "control_24",
         "control_6",
         "control_24",
         "vehicle_24",
         "vehicle_24",
         "vehicle_24",
         "vehicle_24",
         "vehicle_24",
         "vehicle_24",
         "etoh_24",
         "etoh_24",
         "etoh_24",
         "control_6",
         "etoh_24",
         "etoh_24",
         "etoh_24",
         "lipof_24",
         "lipof_24",
         "lipof_24",
         "lipof_24",
         "lipof_24",
         "lipof_24",
         "control_6",
         "control_6",
         "control_6",
         "vehicle_6",
         "vehicle_6",
         "vehicle_6"))



# Comprobación de nombres asignados correctamente
data.frame(
  counts_matrix = colnames(counts_matrix),
  groups = groups$sample,
  condition = groups$condition
  ) ## GROUP ES IGUAL A SAMPLE? SAMPLE DEBERIA SER GROUP$SAMPLE NO?

##' Comprobación de nombres asignados correctamente
data.frame(
  sample = groups$sample,
  condition = groups$condition
)

all(groups$sample == colnames(counts_matrix))


# Verificación de valores negativos ## FALTANTES (NAs) Y NEGATIVOS
any(is.na(counts_matrix))
any(counts_matrix < 0)

# Verificación de que sean números ## VERIFICACI´ON DE QUE DATOS SEAN NUM´ERICOS
is.numeric(counts_matrix[,1])

# Verificar cuántas muestras hay en cada condición
table(groups$condition)

# Verificar que los gene_ID sean rows ## FILAS O ESPA˜NOL O INGLES, NO ES CODIGO PERSONAL DEBE TENER FORMATO
head(rownames(counts_matrix))


# Cambio de objetos a factor
groups$condition <- factor(groups$condition)
class(groups$condition)
levels(groups$condition)

# Diseño de matriz de comparación de grupos
design <- model.matrix(~0 + condition, data = groups)
colnames(design) <- levels(groups$condition)
design
colnames(design)

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

# Evaluación de la creación de las comparaciones
View(contrast)
colnames(contrast)


# Creación del objeto DGEList
dge <- DGEList(
  counts = counts_matrix,
  group = groups$condition
)
dge


# Eliminación de genes con baja expresión
keep <- filterByExpr(dge, design)
summary(keep)
nrow(dge)

dge <- dge[keep, , keep.lib.sizes = FALSE]
nrow(dge)
dim(dge)


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
View(ajuste)
ajuste


# Verificación
dim(dge)
dge$samples$norm.factors


# Ejecución de comparaciones
qlf_etoh6 <- glmQLFTest(ajuste, 
                        contrast = contrast[, "etoh_6_vs_control_6"])
qlf_etoh24 <- glmQLFTest(ajuste, 
                        contrast = contrast[, "etoh_24_vs_control_24"])
qlf_lipof6 <- glmQLFTest(ajuste, 
                        contrast = contrast[, "lipof_6_vs_vehicle_6"])
qlf_lipof24 <- glmQLFTest(ajuste, 
                        contrast = contrast[, "lipof_24_vs_vehicle_24"])
qlf_vehicle6 <- glmQLFTest(ajuste, 
                        contrast = contrast[, "vehicle_6_vs_control_6"])
qlf_vehicle24 <- glmQLFTest(ajuste, 
                        contrast = contrast[, "vehicle_24_vs_control_24"])
qlf_etoh <- glmQLFTest(ajuste, 
                        contrast = contrast[, "etoh_6_vs_etoh_24"])
qlf_lipof <- glmQLFTest(ajuste, 
                        contrast = contrast[, "lipof_6_vs_lipof_24"])

# Obtención de los resultados de las comparaciones
result_etoh6 <- topTags(qlf_etoh6, n = Inf)$table
result_etoh24 <- topTags(qlf_etoh24, n = Inf)$table
result_lipof6 <- topTags(qlf_lipof6, n = Inf)$table
result_lipof24 <- topTags(qlf_lipof24, n = Inf)$table
result_vehicle6 <- topTags(qlf_vehicle6, n = Inf)$table
result_vehicle24 <- topTags(qlf_vehicle24, n = Inf)$table
result_etoh <- topTags(qlf_etoh, n = Inf)$table
result_lipof <- topTags(qlf_lipof, n = Inf)$table

# Anotación génica en base a la matriz de expresión

## HACER UN LOOP! SE VE MEJOR, LIMITA ERRORES Y LO VUELVE REPRODUCIBLE
# etoh 6
result_etoh6$gene_name <- matrix_exp$gene_name[
  match(rownames(result_etoh6), matrix_exp$gene_id)
]
result_etoh6$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_etoh6), matrix_exp$gene_id)
]

# etoh 24
result_etoh24$gene_name <- matrix_exp$gene_name[
  match(rownames(result_etoh6), matrix_exp$gene_id) ## OJO! ERROR AQUI NO COINCIDE ETOH24 CON ETOH6
]
result_etoh24$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_etoh24), matrix_exp$gene_id)
]

# lipofilico 6
result_lipof6$gene_name <- matrix_exp$gene_name[
  match(rownames(result_lipof6), matrix_exp$gene_id)
]
result_lipof6$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_lipof6), matrix_exp$gene_id)
]

# lipofilico 24
result_lipof24$gene_name <- matrix_exp$gene_name[
  match(rownames(result_lipof24), matrix_exp$gene_id)
]
result_lipof24$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_lipof24), matrix_exp$gene_id)
]

# Vehicle 6
result_vehicle6$gene_name <- matrix_exp$gene_name[
  match(rownames(result_vehicle6), matrix_exp$gene_id)
]
result_vehicle6$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_vehicle6), matrix_exp$gene_id)
]

# Vehicle 24
result_vehicle24$gene_name <- matrix_exp$gene_name[
  match(rownames(result_vehicle24), matrix_exp$gene_id)
]
result_vehicle24$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_vehicle24), matrix_exp$gene_id)
]

# etoh
result_etoh$gene_name <- matrix_exp$gene_name[
  match(rownames(result_etoh), matrix_exp$gene_id)
]
result_etoh$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_etoh), matrix_exp$gene_id)
]

# lipofilic
result_lipof$gene_name <- matrix_exp$gene_name[
  match(rownames(result_lipof), matrix_exp$gene_id)
]
result_lipof$gene_biotype <- matrix_exp$gene_biotype[
  match(rownames(result_lipof), matrix_exp$gene_id)
]

# Selección de genes diferencialmente expresados
DEG_etoh6 <- result_etoh6[
  result_etoh6$FDR < 0.05 &
    abs(result_etoh6$logFC) >= 1,
    ]
dim(DEG_etoh6) # 1451 X 7

DEG_etoh24 <- result_etoh24[
  result_etoh24$FDR < 0.05 &
    abs(result_etoh24$logFC) >= 1,
]
dim(DEG_etoh24) # 0 X 7

DEG_lipof6 <- result_lipof6[
  result_lipof6$FDR < 0.05 &
    abs(result_lipof6$logFC) >= 1,
]
dim(DEG_lipof6) # 1032 X 7

DEG_lipof24 <- result_lipof24[
  result_lipof24$FDR < 0.05 &
    abs(result_lipof24$logFC) >= 1,
]
dim(DEG_lipof24) # 221 X 7

DEG_vehicle6 <- result_vehicle6[
  result_vehicle6$FDR < 0.05 &
    abs(result_vehicle6$logFC) >= 1,
]
dim(DEG_vehicle6) # 3 X 5

DEG_vehicle24 <- result_vehicle24[
  result_vehicle24$FDR < 0.05 &
    abs(result_vehicle24$logFC) >= 1,
]
dim(DEG_vehicle24) # 31 X 7

DEG_etoh <- result_etoh[
  result_etoh$FDR < 0.05 &
    abs(result_etoh$logFC) >= 1,
]
dim(DEG_etoh) # 183 X 7

DEG_lipof <- result_lipof[
  result_lipof$FDR < 0.05 &
    abs(result_lipof$logFC) >= 1,
]
dim(DEG_lipof) # 340 X 7
