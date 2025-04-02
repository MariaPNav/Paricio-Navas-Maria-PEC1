# PEC1. Anàlisi de dades òmiques

## Maria Paricio Navas

library(SummarizedExperiment) # Es carreguen diferents llibreries
library(tidyverse)
library(readr)
library(pheatmap)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(cluster)
library(RColorBrewer)
library(matrixStats)
library(reshape2)
library(mixOmics)
library(dplyr)

# 1. Obtenció i estructura de les dades
cachexia_df <- read_csv("human_cachexia.csv") # Es carrega el fitxer csv amb les dades [1]

glimpse(cachexia_df) # Es visualitzen les primeres files

## Identificació de dades d'expressió i metadades

data_matrix <- as.matrix(cachexia_df[, 3:ncol(cachexia_df)]) # Es crea la matriu d'expressió
data_matrix <- t(data_matrix)  # Es transposa per a que els metabòlits estiguin en files i les mostres en columnes (per a poder crear l'objecte de SummarizedExperiment)
head(data_matrix) # S'observen les primeres files de la matriu

rownames(data_matrix) <- colnames(cachexia_df)[3:ncol(cachexia_df)]  # S'asignen els noms dels metabòlits a les files
colnames(data_matrix) <- cachexia_df[[1]]  # S'asignen els noms de les mostres (Patient ID) a les columnes

metab_metadata <- data.frame(Metabolito = rownames(data_matrix)) # Es creen les metadades de les files (metabòlits)
rownames(metab_metadata) <- rownames(data_matrix) # S'asignen els noms de les files a les metadades

pacients_metadata <- data.frame( # Es creen les metadades de les columnes (dels pacients) 
  patient_id = colnames(data_matrix), # S'asigna el nom a les columnes (Patient ID)
  condition = rep(c("control", "cachexia"), length.out = ncol(data_matrix)) # S'afegeix la condició
  )
rownames(pacients_metadata) <- pacients_metadata$patient_id

## Creació de l'objecte de classe Summarized Experiment:

se <- SummarizedExperiment( # Es crea un objecte de clase SummarizedExperiment[2]
  assays = list(counts = data_matrix),
  rowData = metab_metadata,
  colData = pacients_metadata)
dim(se) # Es visualitza la dimensió de l'objecte
se # Es visualitza el propi SummarizedExperiment per observar que tot sigui correcte (nom de les columnes, files, etc)

save(se, file = "summarized_experiment_cachexia.Rda") # Es guarda el SummarizedExperiment (se) com a Rda

# 2. Preprocessament de les dades

summary(t(assay(se))) # Es fa un resum estadístic del SummarizedExperiment amb la seva transposada, perquè s'ha canviat l'ordre previament

boxplot(t(assay(se))) # Es fa un diagrama de caixes del SummarizedExperiment per a observar les variançes en les dades

png("grafic1.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
boxplot(t(assay(se)), las = 2, main = "Boxplot sense normalitzar")
dev.off()

dades_log <- log1p(assay(se))  # S'aplica una transformació logarítmica (log(x+1)) per a estabilitzar la variança [3]

distancies <- rowSums((dades_log) - matrixStats::colMeans2(dades_log)^2) # Es calcula la mitjana dels quadrats de cada columna de la matriu x, es a dir, per a cada variable/metabòlit estima la seva potència quadràtica mitja [4] 
llindar <- quantile(distancies, 0.99) # Es defineix el llindar
outliers <- which(distancies > llindar) # Es calculen els outliers en els que les distàncies superin el llindar
outliers # Es mostren els outliers

se <- se[, -outliers] # S'eliminen els outliers de l'objecte SummarizedExperiment (se)
dades_log <- dades_log[, -outliers]  # S'eliminen també els outliers de les dades transformades

normalize_quantils <- function(x) { # Es crea una funció per a normalitzar per quantils de forma manual amb rank i sort [5]
  ranked <- apply(x, 2, rank, ties.method = "min")
  sorted <- apply(x, 2, sort)
  mean_values <- apply(sorted, 1, mean)
  norm_matrix <- apply(ranked, 2, function(r) mean_values[r])
  return(norm_matrix)}

dades_lognorm <- normalize_quantils(dades_log) # Es normalitzen les dades transformades

head(dades_lognorm) # Es mostren les primeres files de les dades transformades i normalitzades

rownames(dades_lognorm) <- rownames(assay(se)) # S'asignen els noms de les files per a que siguin iguals a les del SummarizedExperiment (se)
colnames(dades_lognorm) <- colnames(assay(se)) # S'assignen els noms de les columnes per a que siguin iguals a les del SummarizedExperiment (se)

png("grafic2.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
boxplot(dades_lognorm, las = 2, main = "Boxplot Transformació Log + Normalització per quantils", col = "lightblue") # Es visualitza la transformació i normalització amb un diagrama de caixes
dev.off()

# 3. Control de qualitat

png("grafic3.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
hist(as.vector(dades_lognorm), # Es fa un histograma com a control de qualitat per a veure una distribució global dels valors d'expressió després del preprocessament
     main = "Histograma de valors d'expressió després de preprocessar les dades",
     xlab = "Expressió (log)",
     col = "lightblue", breaks = 30)
dev.off()

png("grafic4.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
qqnorm(as.vector(as.matrix(dades_lognorm)),  # Es grafica un qqplot per a comprovar i les dades segueixen una distribució normal [6]
       main = "QQ-Plot de les dades preprocessades")
qqline(as.vector(as.matrix(dades_lognorm)), col = "orange")
dev.off()

# 4. Anàlisi exploratori

## Estadística descriptiva

mitjanes <- rowMeans(dades_lognorm) # Es calculen les mitjanes de cada metabòlit
variançes <- apply(dades_lognorm, 1, var) # Es calculen les variançes de cada metabòlit

taula_estadistica <- data.frame( # Es crea un dataframe amb el nom del metabòlit i les dues estadístiques calculades
  Metabòlit = rownames(dades_lognorm),
  Mitjana = round(mitjanes, 3),
  Variança = round(variançes, 5))

taula_estadistica <- taula_estadistica[order(-taula_estadistica$Variança), ] # Es mostren els metabòlits més variables primer (variança descendent)

head(taula_estadistica, 10) # S'observen els metabòlits més variables de tots

## PCA (Anàlisi de Components Principals)

pca_res <- prcomp(t(dades_lognorm), scale. = TRUE) # S'aplica la PCA a la matriu trasposada perquè es necessiten els metabòlits a les columnes i les mostres a les files [7]
summary(pca_res)
var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)
cum_var <- cumsum(var_exp)

png("grafic5.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
plot(cum_var, type = "b", xlab = "Nº de components principals", ylab = "Variança acumulada",
     main = "Anàlisi PCA") # Es mostra un gràfic amb els diferents PCA marcant límits al 95% de la variança en taronja i el 99% de la variança en verd
abline(h = 0.95, col = "orange", lty = 2)
abline(h = 0.99, col = "green", lty = 2)
dev.off()

which(cum_var >= 0.95)[1]  # És calcula el component mínim que explica al menys el 95% de la variança. S'observa que es PC38. Es redueix molt la dimensionalitat, ja que es passa de tenir 63 variables a tenir-ne 38.

loadings <- pca_res$rotation[, 1:2]  # Per a veure els metabòlits que més participen en la variança total es miren les carregues en els primers components [7]

var_metabolitos_pc1 <- sort(abs(loadings[, 1]), decreasing = TRUE) # Es veuen els metabòlits que més contribueixen a la variança total
head(var_metabolitos_pc1, 10)

## Heatmap dels metabólits més variables

vars <- apply(dades_lognorm, 1, var) # Es calcula la variança de cada metabòlit

var_metabs <- order(vars, decreasing = TRUE)[1:10] # Es seleccionen els 10 metabòlits més variables

mat_var <- dades_lognorm[var_metabs, ] # Es crea un subconjunt de la matriu de dades amb els més variables

png("grafic6.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
pheatmap(mat_var,  # Es grafica un heatmap [8]
         annotation_col = as.data.frame(colData(se)),  # Anotació de condicions
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         scale = "row",
         show_rownames = TRUE,
         main = "Heatmap dels 10 metabólits més variables")
dev.off()

## Clustering jeràrquic de les mostres

dist_mat <- dist(dades_lognorm, method = "euclidean") # Es calculen les distàncies euclidianes [9]
hc <- hclust(dist_mat, method = "average") # Es calcula l'agrupament amb el mètode d'average linkage [9]

png("grafic7.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
plot(hc, main = "Dendrograma")
dev.off()

## Boxplots d'expressió per cada metabòlit en cada grup

vars <- apply(dades_lognorm, 1, var) # Es calcula la variança per metabòlit
var_metabs <- names(sort(vars, decreasing = TRUE))[1:10] # Es seleccionen els 10 més variables ordenant amb la funció sort()

expr_var <- dades_lognorm[var_metabs, ] # Es busca l'expressió dels metabòlits més variables
expr_var_p <- melt(expr_var) # Es preparen les dades (més amples)
colnames(expr_var_p) <- c("Metabolit", "Mostra", "Expressió")
expr_var_p$Condicio <- colData(se)$condition[match(expr_var_p$Mostra, colnames(dades_lognorm))] # S'afegeix la condició

png("grafic8.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
ggplot(expr_var_p, aes(x = Condicio, y = Expressió, fill = Condicio)) + # Es dissenya el gràfic de caixes per cada metabòlit diferenciant l'expressió en cada grup
  geom_boxplot(outlier.shape = 21) +
  facet_wrap(~ Metabolit, scales = "free_y") +
  theme_minimal() +
  labs(title = "Expressió per grups dels metabòlits més variables")
dev.off()

# 5. Anàlisi estadístic

grups <- colData(se)$condition # Es crea un vector per cada grup

resultats_ttest <- data.frame( # Es crea una taula per guardar els resultats
  Metabolit = rownames(dades_lognorm), # Es mostrarà el metabòlit corresponent
  p_value = NA, # Es mostrarà el p-valor
  mitjn_cachexia = NA, # Es mostrarà la mitjana d'expressió del metabòlit dels pacients de caquèxia
  mitjn_control = NA) # Es mostrarà la mitjana d'expressió del metabòlit dels individus control

for (i in 1:nrow(dades_lognorm)) { # S'aplica el t-test per cada metabòlit (iterant a la fila) [10]
  expr <- dades_lognorm[i, ]
  grup_cachexia <- expr[grups == "cachexia"]
  grup_control <- expr[grups == "control"]
  
  test <- t.test(grup_cachexia, grup_control)
  
  resultats_ttest$p_value[i] <- test$p.value
  resultats_ttest$mitjn_cachexia[i] <- mean(grup_cachexia)
  resultats_ttest$mitjn_control[i] <- mean(grup_control)}

write.csv(resultats_ttest, "resultats_ttest.csv", row.names = FALSE) # Es crea un csv amb la taula per a mostrar-la en l'informe
resultats_ttest <- resultats_ttest[order(resultats_ttest$p_value), ] # S'ordenen els resultats per p-valor

head(resultats_ttest, 10) # S'observen els 10 primers, amb el p-valor més baix

# 6. Anàlisi multivariant

X <- t(dades_lognorm) # Es trasposen les dades per tenir les mostres en files i els metabòlits en columnes

Y <- factor(colData(se)$condition) # S'afegeixen els grups per condició (controls o caquèxics)

splsda_model <- splsda(X, Y, ncomp = 2, keepX = c(10, 10)) # S'aplica el model sPLS-DA amb dos components i es seleccionen 10 variables per component [11]

png("grafic9.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe

plotIndiv(splsda_model, comp = c(1,2),  # Es grafica la separació dels grups (controls i caquèxics)
          group = Y, 
          legend = TRUE, 
          title = "sPLS-DA: Separació per grup (controls i pacients de caquèxia)")
dev.off()

png("grafic10.png", width = 800, height = 600) # Es guarda en format .png el gràfic per a mostrar-ho a l'informe
plotVar(splsda_model, comp = 1:2, # Es fa un Biplot on es mostren els metabòlits i les mostres
        cex = 0.8, 
        title = "Metabòlits més rellevants")
dev.off()

loadings <- selectVar(splsda_model, comp = 1) # Es defineixen els loadings [11]

metab_import <- data.frame( # Es crea un dataframe amb els metabòlits i les seves càrregues
  Metabolit = rownames(loadings$value),  # Nom del metabòlit
  load = loadings$value[, 1]             # Valor de la càrrega (load)
)

metab_import_ordre<- metab_import[order(abs(metab_import$load), decreasing = TRUE), ] # S'ordenen per a que apareguin de més a menys rellevant
write.csv(metab_import_ordre, "metab_import_ordre.csv", row.names = FALSE) # Es crea un csv amb la taula per a mostrar-la en l'informe
head(metab_import_ordre, 10) # Es mostren els 10 primers

# 7. Referències

## 1.	GitHub. nutrimetabolomics/metaboData: https://github.com/nutrimetabolomics/metaboData 

## 2.	Bioconductor. SummarizedExperiment: https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html 

## 3.	RPubs. Vega,M. (2023). Análisis multivariante en R: https://rpubs.com/Matiu9714/1035673 

## 4.	R Core Team. matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors): https://cran.r-project.org/web/packages/matrixStats/matrixStats.pdf 

## 5.	Sanderson, S. P. (2024). Exploring data with R: https://www.spsanderson.com/steveondata/posts/2024-03-28/index.html 

## 6.	RDocumentation. qqnorm: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/qqnorm 

## 7.	Gil, C. PCA en R: https://rpubs.com/cristina_gil/pca 

## 8.	RDocumentation. pheatmap: https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap 

## 9.	Arquez, M. PCA Example in R: https://rpubs.com/arquez9512/597881 

## 10.	R-Coder. Prueba t en R: https://r-coder.com/t-test-en-r/ 
  
## 11.	RDocumentation. mixOmics - sPLS-DA: https://www.rdocumentation.org/packages/mixOmics/versions/6.3.2/topics/splsda 

## 12.	Repositori de GitHub on es troba aquest anàlisi: https://github.com/MariaPNav/Paricio-Navas-Maria-PEC1 