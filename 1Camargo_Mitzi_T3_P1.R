---
title: "Camargo_Mitzi_T3_P1.R"
author: "Mitzi Naomi Camargo Arellano"
date: "22/3/2022"
output: html_document
---

## WGCNA

*(a) Instala la biblioteca **WGCNA***

Primero instalamos algunos otros paquetes para que nos permita poder instalar "WCNA"

```{r}
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
BiocManager::install(c("GO.db", "preprocessCore", "impute"));

install.packages("WGCNA")
library("WGCNA")
```

Sigue el tutorial , es decir, en un Rmarkdown ve corriendo lo que marca el tutoral, 1. Networks analysis of liver expression data . . . para construir una red de co expresión, a partir de los datos del tutorial responde las siguientes preguntas:

1\. Preliminares y entrada de datos

```{r}
options(stringsAsFactors = FALSE)
```

1.1 Cargamos la base de datos de higado femenino

```{r}
femData = read.csv ( "LiverFemale3600.csv" )
```

1.2 Observamos el conjunto de datos

```{r}
dim(femData)
nombres(femData)
cabeza(femData)
```

1.3 Para mantener los datos que contienen la expresión del gen y los nombres de los genes y queden como índice en el marco de datos

```{r}
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH;
rownames(datExpr0) = names(femData)[-c(1:8)];
```

1.4 Comprobamos si teniamos genes que tuvieran valores faltantes

```{r}
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
```

1.5 Agrupamos la matriz transpuesta con el fin de identificar valores atípicos en la matriz

```{r}
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,2,1,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 15, col = "red"); # A partir de esta última línea pudimos trazar una línea de corte
```

1.6 Identificamos el valor atípico

```{r}
clust = cutreeStatic ( sampleTree , cutHeight =  15 , minSize =  10 ) 
table ( clust )
```

1.7 Eliminamos el valor atípico y construimos un marco de datos principal

```{r}
keepSamples = (clust==1) #Conservamos el cluster 1, porque aquí están las muestras que nos interesan
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
```

1.8 Se introducen los datos clínicos, se preparan y se depuran

```{r}
traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)
#Eliminamos las columnas con la información inecesaria
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)
#Se forma un marco de datos análogo a los datos de expresión que contendrán los rasgos clínicos. 
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Mice);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
```

1.9 Repetimos la agrupación de muestras junto con un mapa de calor para los datos fenotípicos

```{r}
#Se reagruparon las muestras
sampleTree2 = hclust(dist(datExpr), method = "average")
#Se convierten los rasgos a una representación con calor: El "blanco" nos indica bajo, el "rojo" alto y el gris "entrada faltante"
traitColors = numbers2colors(datTraits, signed = FALSE);#Se trazó el dendograma y los colores se plasman abajo
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Dendograma de muestra y mapa de calor de rasgos")
```

1.10 El análisis se guarda en un archivo RData

```{r}
( datExpr , datTraits , file = "  FemaleLiver-01-dataInput.RData" )
```

2.  Construcción automática de redes y detección de módulos

    ```{r}
    # Se permiten subprocesos dentro de WGCNA. Para acelerar cálculos
    ## Línea de código necesaria!
    allowWGCNAThreads ()
    #Se cargan los datoas guardados en la primera parte
    # En lnames se guardan los nombres de las variables cargadas
    lnames = load(file = "FemaleLiver-01-dataInput.RData");
    lnames
    ```

2.1 **Forma más conveniente y automática de detectar módulos y construir una red con WGCNA.**

Mediante este método se identifica una potencia a la que se eleva la matriz de correlación para calcular la matriz de adyacencia de la red en base al criterio de aproximación libre de escala.

```{r}
# Se elijió un conjunto de potencias soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# LLamamos a la función de análisis de topología de la red
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Ploteamos los resultados
par(mfrow = c(1,2));
# Índice de ajuste de topología sin escala 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#Corte en la R^2 de h
abline(h=0.90,col="red")
#Conectividad media en la función
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

2.2 **Contribución de la red**

Se eligieron 6 como las potencias más bajas que constituyen una topología sin escala. Se indica la función que genere módulos de tamaño 30 y que fusione módulos más del 25% similares y guarden la matrzi de superposición topológica en un objeto.

```{r}
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
```

2.3 Cada modulo se irá con su tamaño

```{r}
tabla ( net $ colores )
```

2.4 Dendograma resultante de la construcción de los módulos con la agrupación de genes

```{r}
# Se convierten etiquetas en colores para plotear
mergedColors = labels2colors(net$colors)
# Se grafica el dendrograma y los colores el módulo se muestra abajo
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

2.5 Se guardan los resultados como archivo RData

```{r}
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")
```

3.  **Relacionar módulos con información externa e identificar genes importantes**

    ```{r}
    lnames = load(file = "FemaleLiver-01-dataInput.RData");
    #La variable lnames contiene los nombres de las variables cargadas.
    lnames #Para cargar la red que se guardo de la segunda parte
    lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
    lnames
    ```

    3.3 **Cuantificación de asociaciones módulo-rasgo**

    Se identifican los módulos significativamente asociados con los rasgos clínicos medidos. Se buscan asociaciones más significativas.

    ```{r}
    #Se define el número de genes y muestras
    nGenes = ncol(datExpr);
    nSamples = nrow(datExpr);
    # Se recalcular los ME con etiquetas de color 
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, datTraits, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    ```

    3.4 **Visualizar la asociación módulo-rasgo.** Se realiza una interpretación más entendible a partir del código de colores. Se representa cada gen propio del módulo y su coeficiente de correlación.

    ```{r}
    # Se muestran las correlaciones y sus valores p 
    textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                               signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8, 1, 1));
    # Se muestran los valores de correlación dentro de un diagrama de mapa de calor 
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = names(datTraits),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    ```

    ```{r}
    names(datExpr)[moduleColors=="magenta"]
    ```

    3.5 **Archivo de anotación de sonda proporcionado por el fabricante para facilitar la anotación funcional.**

    ```{r}
    annot = read.csv(file = "GeneAnnotation.csv");
    dim(annot)
    names(annot)
    probes = names(datExpr)
    probes2annot = match(probes, annot$substanceBXH)
    # Número de sondas sin anotación:
    sum(is.na(probes2annot))
    # El resultado arrojaría 0.
    ```

    3.6 **Se recopilo toda la información de genes significativos relacionados con el peso corporal.**

    ```{r}
    # Se crea un formato de datos inicial
    geneInfo0 = data.frame(substanceBXH = probes,
                          geneSymbol = annot$gene_symbol[probes2annot],
                          LocusLinkID = annot$LocusLinkID[probes2annot],
                          moduleColor = moduleColors,
                          geneTraitSignificance,
                          GSPvalue)
    # Se ordenan módulos por su importancia para el peso
    modOrder = order(-abs(cor(MEs, weight, use = "p")));
    # Se agrega información de membresía del módulo en el orden elegido 
    for (mod in 1:ncol(geneModuleMembership))
    {
      oldNames = names(geneInfo0)
      geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                             MMPvalue[, modOrder[mod]]);
      names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                           paste("p.MM.", modNames[modOrder[mod]], sep=""))
    }
    # Se ordenan los genes en la variable geneInfo primero por color de módulo, luego por geneTraitSignificance
    geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
    geneInfo = geneInfo0[geneOrder, ]
    ```

    3.7 Los resultados se guardan en un archivo de salida para su posterior análisis.

    ```{r}
    geneInfo , archivo =  "geneInfo.csv" )
    ```

    **4. Se interconecta el análisis de redes con otros datos, como la anotación funcional y la ontología de genes**

```{r}
# Cargamos los datos de expresión y rasgo guardados en la primera parte 
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#La variable lnames tiene nombres de las variables cargadas. 
lnames
# Se carga datos de red guardados en la segunda parte
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames
```

```{r}
# Leer la anotación de la sonda
annot = read.csv(file = "GeneAnnotation.csv");
# Se hace match entre las sondas en el conjunto de datos con los ID de sonda en el archivo de anotación
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Se obtienen los ID de Locuis Link correspondientes
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Se elijen los módulos interesantes 
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
  # Seleccionamos las sondas de módulo 
  modGenes = (moduleColors==module)
  # Obtenemos los códigos de ID de entrada
  modLLIDs = allLLIDs[modGenes];
  # Se escriben en un archivo
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# Como antecedentes en el análisis de enriquecimiento, se usan todas las sondas en el análisis. 
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)
```

```{r}
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);
```

```{r}
#anRichment(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);  # Does not work yet.
```

```{r}
tab = GOenr$bestPTerms[[4]]$enrichment
```

```{r}
names(tab)
```

```{r}
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
```

```{r}
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Se redondean las columnas numéricas a 2 decimales:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Se trunca el nombre del término a un máximo de 40 caracteres 
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Se acortan los nombres de las columnas: 
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Se establece el ancho de la salida de R. 
options(width=95)
# Se muestra la tabla de enriquecimiento: 
screenTab
```

**5. Exportación de redes a software externo**

```{r}
# Se cargan los datos de expresión y rasgo guardados en la primera parte 
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#La variable lnames contiene los nombres de las variables cargadas. 
lnames
# También cargamos los datos de red guardados en la segunda parte..
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames
```

```{r}
# Se calcula otra vez la superposición topológica si es necesario
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Se lee en el archivo de anotación 
annot = read.csv(file = "GeneAnnotation.csv");
# Se seleccionan módulos 
modules = c("brown", "red");
# Se seleccionan sondas de módulo 
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Se selecciona el mod de superposición topológica correspondiente
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Exportamos la red a archivos de lista de nodos y borde Cytoscape puede leer 
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.5,
  nodeNames = modProbes,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[inModule]);
```

***Custionario***

i\. En dónde está el objeto matriz de expresión y de qué dimensión es?

ii\. ¿Por qué se eliminan datos que son demasiado distintos ? (Vean la gráfica Sample

clustering to detect outliers)

iii\. ¿Qué criterio utilizan para generar la red, es decir explica el significado de la

variable softpower?

iv\. ¿Por qué crees que genes que pertenecen al mismo clúster sin relevantes.

v\. Discute algunos de los resultados que te parezcan interesantes de los clústers y su

relación con los meta- datos ( datos de loas hembras ratones).
