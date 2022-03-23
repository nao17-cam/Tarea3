---
title: "Camargo_Mitzi_T3_P2.R"
author: "Mitzi Naomi Camargo Arellano"
date: "21/3/2022"
output: html_document
---

## Red booleana

*Cargamos la librería "BoolNet"*

```{r}
library(BoolNet)
```

#### *Generamos una red de regulación transripcional (con 5 nodos y, al menos, 7 interacciones regulatorias, con al menos 3 inhibiciones)*

![](imagen1.jpg)

1.  *Hacemos las tablas de verdad de las reglas*

    Tabla 1.

    | B   | D   | B AND D |
    |-----|-----|---------|
    | 0   | 0   | 0       |
    | 0   | 1   | 0       |
    | 1   | 0   | 0       |
    | 1   | 1   | 1       |

    Tabla 2.

    | B   | C   | D   | B AND D AND NOT C |
    |-----|-----|-----|-------------------|
    | 0   | 0   | 0   | 0                 |
    | 0   | 0   | 1   | 0                 |
    | 0   | 1   | 0   | 0                 |
    | 0   | 1   | 1   | 0                 |
    | 1   | 0   | 0   | 0                 |
    | 1   | 0   | 1   | 0                 |
    | 1   | 1   | 0   | 1                 |
    | 1   | 1   | 1   | 0                 |

    Tabla 3.

    | A   | D   | A AND NOT D |
    |-----|-----|-------------|
    | 0   | 0   | 0           |
    | 0   | 1   | 0           |
    | 1   | 0   | 1           |
    | 1   | 1   | 0           |

    Tabla 4.

    | A   | E   | A AND E |
    |-----|-----|---------|
    | 0   | 0   | 0       |
    | 0   | 1   | 0       |
    | 1   | 0   | 0       |
    | 1   | 1   | 1       |

    Tabla 5.

    | A   | B   | A AND NOT B |
    |-----|-----|-------------|
    | 0   | 0   | 0           |
    | 0   | 1   | 0           |
    | 1   | 0   | 1           |
    | 1   | 1   | 0           |

2.  *Escribimos las reglas para generar una red de tipo Boolnet*. *Lo pasamos a un archivo de texto plano.*

    *Uso **loadNetwork** para convertir de texto a red booleana*

    ```{r}
    red_boo1.0 <- loadNetwork("archivitito.txt")
    red_boo1.0 #Me imprime
    ```

    *Otra forma para ver las reglas booleanas*

    ```{r}
    reg_boo <- truthTableToSymbolic(red_boo1.0)
    reg_boo # Para visualizar mis reglas
    ```

3.  *Encuentramos todos los atractores de la red. Usamos la función "getAttractors".*

    ```{r}
    todos_atractores <- getAttractors(red_boo1.0)
    todos_atractores #Imprimo
    ```

4.  *¿Cuál sería el estado final más probable?*

    Todos los genes se apaguen. Solo tenemos un atractor con 32 estados.

5.  *Dibuja todos los estados y sus atractores*

    ```{r}
    plotAttractors(todos_atractores)
    ```
