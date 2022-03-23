---
title: "Camargo_Mitzi_T3_P4.R"
author: "Mitzi Naomi Camargo Arellano"
date: "22/3/2022"
output: html_document
---

## Red de señalización

*A partir de la red de señalización mostrada en la figura 1*

![](imagen2.jpg)

*(a) Elabora una representación verbal de la red*

Empezando por "Depolar"

Depolar es activado por KEV, Ca2+, AnionEM y es reprimido por H+ATPase y KOUT. CalM es reprimido por Depolar. KEV es activado por Ca2+. H+ATPase es reprimido por Ca2+. KOUT es activado por Depolar. KAP es activado por Depolar. Closure es activado por KOUT, KAP y AnionEM. AnionEM es activado por Ca2+. Ca2+ es activado por CalM, CIS y es reprimido por Ca2+ATPase. Ca2+ATPase es activado por Ca2+. NOS es activado por Ca2+. NO es activado por NOS. GC es activado por NO. cGMP es activado por GC. ADPRc es activado por NO. cADPR es activado por ADPRc. CIS es activado por cGMP y cADPR. PLC es activado por Ca2+. InsP3 es activado por PLC.

(*b) Elabora una tabla de posibles reglas de la red*

+----------------+--------------------+----------------+------------------------------------------------+
| Nodo           | Activadores        | Represores     | Regla                                          |
+================+====================+================+================================================+
| Depolar (A)    | KEV, Ca2+, AnionEM | H+ATPase, KOUT | KEV & Ca2+ OR AnionEM ADN NOT H+ATPase OR KOUT |
+----------------+--------------------+----------------+------------------------------------------------+
| CalM (B)       | \-                 | Depolar        | Depolar                                        |
+----------------+--------------------+----------------+------------------------------------------------+
| KEV (C)        | Ca2+               | \-             | Ca2+                                           |
+----------------+--------------------+----------------+------------------------------------------------+
| H+ATPase (D)   | \-                 | Ca2+           | Ca2+                                           |
+----------------+--------------------+----------------+------------------------------------------------+
| KOUT (E)       | Depolar            | \-             | Depolar                                        |
+----------------+--------------------+----------------+------------------------------------------------+
| KAP (F)        | Depolar            | \-             | Depolar                                        |
+----------------+--------------------+----------------+------------------------------------------------+
| Closure (G)    | KOUT, KAP, AnionEM | \-             | KOUT AND KAP OR NOT AnionEM                    |
+----------------+--------------------+----------------+------------------------------------------------+
| AnionEM (H)    | Ca2+               | \-             | Ca2+                                           |
+----------------+--------------------+----------------+------------------------------------------------+
| Ca2+ (I)       | CalM, CIS          | Ca2+ATPase     | CalM AND CIS AND NOT Ca2+ATPase                |
+----------------+--------------------+----------------+------------------------------------------------+
| Ca2+ATPase (J) | Ca2+               | \-             | Ca2+                                           |
+----------------+--------------------+----------------+------------------------------------------------+
| NOS (K)        | Ca2+               | \-             | Ca2+                                           |
+----------------+--------------------+----------------+------------------------------------------------+
| NO (L)         | NOS                | \-             | NOS                                            |
+----------------+--------------------+----------------+------------------------------------------------+
| GC (M)         | NO                 | \-             | NO                                             |
+----------------+--------------------+----------------+------------------------------------------------+
| cGMP (N)       | GC                 | \-             | GC                                             |
+----------------+--------------------+----------------+------------------------------------------------+
| ADPRc (O)      | NO                 | \-             | NO                                             |
+----------------+--------------------+----------------+------------------------------------------------+
| cADPR (Q)      | ADPRc              | \-             | ADPRc                                          |
+----------------+--------------------+----------------+------------------------------------------------+
| CIS (R)        | cGMP, cADPR        | \-             | cGMP AND cADPR                                 |
+----------------+--------------------+----------------+------------------------------------------------+
| PLC (S)        | Ca2+               | \-             | Ca2+                                           |
+----------------+--------------------+----------------+------------------------------------------------+
| InsP3 (T)      | PLC                | \-             | PLC                                            |
+----------------+--------------------+----------------+------------------------------------------------+

*(c) Encuentra y discute biológicamente el significado de los atractores (Usa BoolNet)*

```{r}
#Cargamos la librería
library("BoolNet")
# Subimos el archivo plano
regulacion <- loadNetwork("archivo2022.txt")
regulacion
```

```{r}
#Generamos los atractores 
atractores2022 <- getAttractors(regulacion)
```

La red es una representación simplificada de la que se muestra en los siguientes papers Li S, Assmann SM, Albert R. Predicting essential components of signal transduction networks: a dynamic model of guard cell abscisic acid signaling. PLoS Biol 2006;4:e312. Saadatpour A, Albert I, Albert R. Attractor analysis of asynchronous Boolean models of signal transduction networks. J Theor Biol 2010;4:641--56.
