---
title: "Camargo_Mitzi_T3_P3.R"
author: "Mitzi Naomi Camargo Arellano"
date: "22/3/2022"
output: html_document
---

## Red de co-gustos

*A partir de las dos versiones de la red de co-gustos: la dicotómica (D) y la pesada (P).*

*Elaboramos un programa en R que calculara:*

1.  Dibuja la red con al menos tres layouts (D y P).

Nota. Tiene que cargarse la librería "igraph". En el momento se le indicara como.

1.  1 Cargamos la base de datos con "read.csv"

    ```{r}
    red_cogustos <- read.csv("Red de co-gustos 2022 - Hoja 1 (1).csv")
    red_cogustos #Imprimimos
    ```

    1.2 Usamos row.names para cohesionar y que no se hagan duplicados

```{r}
row.names(red_cogustos) <- red_cogustos[ ,1]
red_cogustos
```

1.3 Quitamos las lineas que no utilizamos.

```{r}
red_cogustos <- red_cogustos[ ,-1] #Quitamos la primera columna 
red_cogustos #Visualizamos
mtzred_cogustos1.0 <- cor(t(red_cogustos)) #Hacemos una correlación y transpoción
View(mtzred_cogustos1.0) #Visualizamos
```

1.  4 Quitamos los 1 de la diagonal, ya que compara a la persona consigo mismo

    ```{r}
    diag(mtzred_cogustos1.0) <- 0
    View(mtzred_cogustos1.0) #Visualizamos ****
    ```

    1.5 Hacemos una matriz de correlación para la **red pesada**

    ```{r}
    matriz_adyacencia_rpesada <- (1+mtzred_cogustos1.0)/2 #Se ajustan los valores conforme los datos de correlación
    View(matriz_adyacencia_rpesada) #Visualizo mis datos para observar que se realizó el ajuste


    diag(matriz_adyacencia_rpesada) <-0 #Cambio los 1 de la diagonal. Ya que compararían a la persona consigo mismo**
    View(matriz_adyacencia_rpesada) #Observo la matriz y verifico que los 1 hayan sido sustituidos.
    ```

    Cargamos la librería "igraph".

    ```{r}
    library(igraph)
    ```

    1.  6 Dibujamos 3 layout para la **red pesada**

```{r}
red_pesada <- graph_from_adjacency_matrix(matriz_adyacencia_rpesada, mode = "undirected",weighted = TRUE) 
E(red_pesada)$color <-"black" #A partir del color negro podrmos ver que tan conectados se encuentran los nodos. Entre más intensidad mayor conexión. 
plot(red_pesada, edge.width=E(red_pesada)$weight)

# Primer estilo
plot(red_pesada, layout = layout.davidson.harel)
# Segundo estilo
plot(red_pesada, layout = layout.circle)
# Tercer estilo
plot(red_pesada, layout = layout.kamada.kawai)
```

1.  7 Hacemos una matriz de correlación para la **red dicotomica**

    ```{r}
    matriz_adyacencia_rdicotomica <- (1+mtzred_cogustos1.0)/2 # Realizamos la matriz de adyacencia
    View(matriz_adyacencia_rdicotomica) #Visualizo mis datos para observar que se realizó el ajuste


    diag(matriz_adyacencia_rdicotomica) <-0 #Asigno 0 a los valores de la diagonal
    View(matriz_adyacencia_rdicotomica) #Observo la matriz y verifico que los 1 hayan sido sustituidos.


    matriz_adyacencia_rdicotomica1.0<-ifelse(matriz_adyacencia_rdicotomica>0.5,1,0) 
    View(matriz_adyacencia_rdicotomica1.0)
    ```

    1.8 Cargamos la librería "igraph".

    ```{r}
    library(igraph)
    ```

1.9 Dibujamos 3 layout para la **red dicotomica**

```{r}
red_dicotomica <- graph_from_adjacency_matrix (matriz_adyacencia, mode = "undirected")

# Primer estilo
plot(red_dicotomica, layout = layout.fruchterman.reingold.grid)
# Segundo estilo
plot(red_dicotomica, layout = layout.graphopt)
# Tercer estilo
plot(red_dicotomica, layout = layout.gem)
```

2.  *Distribucion de conectividades (D)*

```{r}
distribucion_conectividades <- degree(red_dicotomica) 
distribucion_conectividades #Imprimimos
```

3.  *Nodos más conectados (D)*

    ```{r}
    nodos_mas_conectados <- degree(red_dicotomica) #Me muestra el número de conexiones
    sort(nodos_mas_conectados, decreasing = TRUE) #Uso sort para orden de mayor a menor el número de conexiones 
    ```

<!-- -->

4.  *Los nodos más importantes con al menos tres medidas de centralidad (D)*

```{r}
i1 <- sort(eccentricity(red_dicotomica))[1:3]
i1
i2 <- sort(degree(red_dicotomica), decreasing = TRUE)[1:3]
i2 
i3 <- sort(closeness(red_dicotomica), decreasing= TRUE)[1:3]
i3 
```

5.  *Los clústers obtenidos con al menos tres métodos de clusterización D y P*

    Clústers para la **red pesada \*\*\*\*\***

```{r}
View(matriz_adyacencia_rpesada)

# Primer método
pesada_cluster1<- cluster_edge_betweenness(red_pesada)
table(membership(pesada_cluster1))

plot(pesada_cluster1,red_pesada)

# Segundo método

pesada_cluster2 <- cluster_fast_greedy(red_pesada) 
table(membership(pesada_cluster2))
plot(pesada_cluster2,red_pesada)

# Tercer método

pesada_cluster3 <- cluster_label_prop(red_pesada) 
table(membership(pesada_cluster2))
plot(pesada_cluster3,red_pesada)
```

Clústers para la **red dicotomica \*\*\*\*\***

```{r}
View(mat_ady_1)

# Primer método
dicotomica_cluster1.0 <- cluster_edge_betweenness(red_dicotomica,directed = FALSE)
table(membership(dicotomica_cluster1.0))
plot(dicotomica_cluster1.0,red_dicotomica)

# Segundo método

dicotomica_cluster2.0 <- cluster_fast_greedy(red_dicotomica)
table(membership(dicotomica_cluster2.0))
plot(dicotomica_cluster2.0,red_dicotomica)

# Tercer método

dicotomica_cluster3.0 <-cluster_label_prop(red_dicotomica)
table(membership(dicotomica_cluster3.0))
plot(dicotomica_cluster3.0,red_dicotomica)
```

6.  *Discute si las redes D y P son dirigidas o no.*

    No son dirigidas porque las interacciones no necesariamente son devueltas, si hay conexiones bidireccionales pero no en la misma proporción. La conexión entre ij no implica la conexión entre ji.

7.  *¿Cómo podrías encontrar clicas, si las hay?*

    Las clicas representan subgrupos que se encuentran más conectados entre sí que no pueden formar un cluster nuevo ya que no se encuentran lo suficientemente conectadas. Podríamos observarlo con los nodos más agrupados que se encuentran en los clusters.
