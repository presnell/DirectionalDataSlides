Introduction
============

Examples
--------

Wish to analyze data in which response is a "direction":

-   2d directional data are called *circular* data
-   3d directional data are called *spherical* data
-   not all "directional" data are directions in the usual sense
-   "directional" data may also arise in higher dimensions

Wind Directions
---------------

-   Recorded at Col de la Roa, Italian Alps
-   n = 310 (first 40 listed below)
-   Radians, clockwise from north
-   Source: Agostinelli (CSDA 2007); also R package `circular`

### Data

    ascii(matrix(wind[1:40], ncol=5, byrow=TRUE), digits=2,
          include.rownames=FALSE, include.colnames=FALSE)

### Plot

``` {#windDataPlot}
require("circular")
par(mar=c(0,0,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plot(windc, cex=1.5, axes=FALSE,
     bin=360, stack=TRUE, sep=0.035, shrink=1.3)
axis.circular(at=circular(seq(0, (7/4)*pi, pi/4),
                  template="geographics"),
              labels=c("N","NE","E","SE","S","SW","W","NW"),
              cex=1.4)
ticks.circular(circular(seq(0, (15/8)*pi, pi/8)),
               zero=pi/2, rotation="clock",
               tcl=0.075)
```

![](file:Plots/wind.png)

Arrival Times at an ICU
-----------------------

-   24-hour clock times (format `hrs.mins`)
-   n = 254 (first 32 listed below)
-   Source: Cox & Lewis (1966); also Fisher (1993) and R package
    `circular`

### Data

    ascii(matrix(fisherB1[1:32], ncol=4, byrow=TRUE), digits=2,
          include.rownames=FALSE, include.colnames=FALSE)

### Plot

``` {#icuDataPlot}
## Note that pch=17 does not work properly here.
par(mar=c(0,0,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plot(fisherB1c, cex=1.5, axes=TRUE,
     bin=360, stack=TRUE, sep=0.035, shrink=1.3)
```

![](file:Plots/icu.png)

Primate Vertebrae
-----------------

-   Orientation of left superior facet of last lumbar vertebra in
    humans, gorillas, and chimpanzees
-   Source: Keifer (2005 UF Anthropology MA Thesis)

![as superior articulate process.](file:Pictures/Gray93.png)

Plot of Human Data
------------------

![chimpanzees (red), 16 gorillas (green) and 19 humans
(blue).](file:Pictures/vertebraeOnSphere.png)

Butterfly Migrations
--------------------

-   Direction of travel observed for 2649 migrating butterflies in
    Florida
-   Source: Thomas J Walker, University of Florida, Dept of Entomology
    and Nematology
-   Other variables:
    -   site: 23 locations in Florida
    -   observer: Thomas Walker (tw) or James J. Whitesell (jw)
    -   species: cloudless sulphur (cs), gulf fritillary (gf),
        long-tailed skipper (lt)
    -   distance to coast (km)
    -   date and time of observation
    -   percentage of sky free of clouds
    -   quality of sunlight: (b)right, (h)aze, (o)bstructed, (p)artly
        obstructed
    -   presence/absence and direction (N, NE, E, SE, S, SW, W, NW) of
        wind
    -   temperature

Why is the Analysis of Directional Data Different?
--------------------------------------------------

-   First three observations from the wind directions data:
    `paste(round(wind[1:3], 2), collapse=", ")`{.r .rundoc-block
    rundoc-language="R"}
-   The mean of these three numbers is `round(mean(wind[1:3]), 2)`{.r
    .rundoc-block rundoc-language="R"} {{{results(`2.47`)}}}
-   What do you think?

``` {#meanAnglePlot}
par(mar=c(0,0,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plot(windc[1:3], cex=2, lwd=1.5, axes=TRUE, ticks=TRUE, tcl=0.05)
points(circular(mean(wind[1:3]), units="radians", template="geographics"),
       pch=8, cex=4) 
```

![](file:Plots/meanAngle.png)

Graphical Display of Directional Data
=====================================

Graphical Display of Circular Data (in R)
-----------------------------------------

-   Have already seen simple dot plots for circular data, e.g., for the
    wind data:

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code" rundoc-eval="no"}
<<windConvert>>
<<windDataPlot>>
```

Graphical Display of Circular Data (in R) (ctd)
-----------------------------------------------

-   and for the ICU data:

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code" rundoc-eval="no"}
<<icuDataPlot>>
```

-   and one more ...

Graphical Display of Circular Data (in R) (ctd)
-----------------------------------------------

``` {#antsDataPlot}
par(mar=c(0,0,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plot(fisherB10c$set1, units="degrees", zero=pi/2,
     rotation="clock", pch=16, cex=1.5)
ticks.circular(circular(seq(0, (11/6)*pi, pi/6)),
               zero=pi/2, rotation="clock", tcl=0.075)
points(fisherB10c$set2, zero=pi/2,
       rotation="clock", pch=16, col="darkgrey",
       next.points=-0.1, cex=1.5)
points(fisherB10c$set3, zero=pi/2,
       rotation="clock", pch=1,
       next.points=0.1, cex=1.5)
```

![three different experimental conditions:](file:Plots/ants.png)

Graphical Display of Circular Data (in R) (ctd)
-----------------------------------------------

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code" rundoc-eval="no"}
<<antsDataPlot>>
```

Circular Histograms
-------------------

-   [Circular
    histograms](https://www.google.com/search?q=R+circular+histogram)
    exist (see Fisher and Mardia and Jupp) but is there a ready-made
    function in R?

Rose Diagrams
-------------

-   Invented by [Florence
    Nightingale](https://en.wikipedia.org/wiki/Florence_Nightingale)
    (elected first female member of the Royals Statistical Society in
    1859; honorary member of ASA)
-   [Nightingale's rose in
    R](https://github.com/jennybc/r-graph-catalog/tree/master/figures/fig05-14_nightingale-data)
    (see also [this
    post](http://www.r-bloggers.com/going-beyond-florence-nightingales-data-diagram-did-flo-blow-it-with-wedges/)
    and the [R graph
    catalog](http://shiny.stat.ubc.ca/r-graph-catalog/))
-   Note that radii of segments are proportional to *square root* of the
    frequencies (counts), so that areas are proportional to frequencies.
    Is this the right thing to do?
-   Rose diagrams suffer from the same problems as histograms. The
    impression conveyed may depend strongly on:
    -   the binwidth of the cells
    -   the choice of starting point for the bins

Adding a Rose Diagram to the Plot of Wind Directions
----------------------------------------------------

``` {#windRosePart .r .rundoc-block rundoc-language="R" rundoc-exports="code" rundoc-eval="no"}
rose.diag(windc, bins=16, col="darkgrey",
          cex=1.5, prop=1.35, add=TRUE)
```

Adding a Rose Diagram to the Plot of Wind Directions
----------------------------------------------------

``` {#windRose}
<<windDataPlot>>
<<windRosePart>>
```

![(segment radii are proportional to square roots of
counts).](file:Plots/windRose.png)

Changing the Binwidth
---------------------

### Fewer/Wider Bins

``` {#windRoseWideBins}
<<windDataPlot>>
<<windRoseWideBinsPart>>
```

![](file:Plots/windRoseWide.png)

### Narrow Bins

``` {#windRoseNarrowBins}
<<windDataPlot>>
<<windRoseNarrowBinsPart>>
```

![](file:Plots/windRoseNarrow.png)

Changing the Radii
------------------

-   I think that the default "radii proportional to counts" is generally
    best, but this is not always obvious. The scale certainly makes a
    big difference however.

``` {#windRoseLinearPart .r .rundoc-block rundoc-language="R" rundoc-exports="code" rundoc-eval="no"}
rose.diag(windc, bins=16, col="darkgrey",
          radii.scale="linear",
          cex=1.5, prop=2.4, add=TRUE)
```

Changing the Radii
------------------

``` {#windRoseLinear}
<<windDataPlot>>
<<windRoseLinearPart>>
```

![(segment radii proportional to
counts).](file:Plots/windRoseLinear.png)

Kernel Density Estimates
------------------------

``` {#windKdensPart .r .rundoc-block rundoc-language="R" rundoc-exports="code" rundoc-eval="no"}
lines(density.circular(windc, bw=40), lwd=2, lty=1)
```

Kernel Density Estimates
------------------------

``` {#windKdens}
par(mar=c(0,0,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plot(windc, cex=1.5, axes=FALSE,
     bin=360, stack=TRUE, sep=0.035, shrink=1.7)
axis.circular(at=circular(seq(0, (7/4)*pi, pi/4),
                  template="geographics"),
              labels=c("N","NE","E","SE","S","SW","W","NW"),
              cex=1.4)
ticks.circular(circular(seq(0, (15/8)*pi, pi/8)),
               ## zero=pi/2, rotation="clock",
               tcl=0.075)
<<windRosePart>>
<<windKdensPart>>
```

![and kernel density estimate.](file:Plots/windKdens.png)

Spherical Data
--------------

-   Are there any canned routines for plotting spherical data in R?

Basic Summary Statistics
========================

Mean Direction and Mean Resultant Length
----------------------------------------

-   First three observations from the wind directions data:

<!-- -->

    theta <- wind[1:3]
    x <- sin(theta)
    y <- cos(theta)
    ascii(cbind(theta, x, y), digits=2,
          include.rownames=FALSE, include.colnames=TRUE)

-   resultant (sum of direction vectors): (`round(xsum, 3)`{.r
    .rundoc-block rundoc-language="R"}, `round(ysum, 3)`{.r
    .rundoc-block rundoc-language="R"})

-   mean vector: $(\bar{x}, \bar{y}) = $ (`round(xbar, 3)`{.r
    .rundoc-block rundoc-language="R"}, `round(ybar, 3)`{.r
    .rundoc-block rundoc-language="R"})

-   resultant length (Euclidean norm of resultant): R =
    `round(resultantLength, 3)`{.r .rundoc-block rundoc-language="R"}

-   mean resultant length: $\bar{R} = $
    `round(meanResultantLength, 3)`{.r .rundoc-block
    rundoc-language="R"}

-   mean direction: $(\bar{x}, \bar{y})/\bar{R} = $
    (`round(meanDirection[1], 3)`{.r .rundoc-block rundoc-language="R"},
    `round(meanDirection[2], 3)`{.r .rundoc-block rundoc-language="R"})

-   $\tilde{\theta} = $ `round(meanDirectionRadians, 3)`{.r
    .rundoc-block rundoc-language="R"}

Plot
----

``` {#meanDirection}
par(mar=c(0,0,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plot(windc[1:3], cex=2, lwd=1.5, axes=TRUE, ticks=TRUE, tcl=0.05)
points(circular(meanDirectionRadians, units="radians", template="geographics"),
       pch=8, cex=4) 
```

![and their sample mean direction.](file:Plots/meanDirection.png)

Aside: Generating from the Uniform Distribution on the Sphere
=============================================================

Generating Random Points on the Sphere
--------------------------------------

-   Wish to generate a random "direction" in d-dimensions; i.e., an
    observation from the uniform distribution in the $d-1$ sphere.
-   Usual way: let X  ∼ N~d~(0, I) and return U = X/||X||.
-   An alternative rejection sampler:
    -   Repeat until ||X|| &lt;= 1
        -   Let X be uniformly distributed on the cube \[-1,1\]^d^
    -   Return U = X/||X||
-   What is the acceptance rate for the rejection sampler:
    -   Volume of the $d - 1$ sphere is $\pi^{d/2}/\Gamma(d/2 + 1)$
    -   Volume of \[-1,1\]^d^ is 2^d^
    -   Acceptance rate is $(\pi^{1/2}/2)^d/\Gamma(d/2 + 1)$
    -   Curse of dimensionality

<!-- -->

    accRate <- function(d) ((sqrt(pi)/2)^d)/gamma(d/2 + 1)
    d <- 2:10
    ## ar <- matrix(accRate(d), nrow=1,
    ##              dimnames=list("accept rate", "d"=d))
    ar <- rbind("dimension"=d, "accept rate (%)"= 100*accRate(d))
    ascii(ar, digits=0, include.rownames=TRUE, include.colnames=FALSE)

Code for Timing Results
-----------------------

``` {#runifSphereR .r .rundoc-block rundoc-language="R" rundoc-exports="code"}
runifSphere <- function(n, dimension, method=c("norm", "cube", "slownorm")) {
    method <- match.arg(method)
    if (method=="norm") {
        u <- matrix(rnorm(n*dimension), ncol=dimension)
        u <- sweep(u, 1, sqrt(apply(u*u, 1, sum)), "/")
    } else if (method=="slownorm") {
        u <- matrix(nrow=n, ncol=dimension)
        for (i in 1:n) {
            x <- rnorm(dimension)
            xnorm <- sqrt(sum(x^2))
            u[i,] <- x/xnorm
        }
    } else {
        u <- matrix(nrow=n, ncol=dimension)
        for (i in 1:n) {
            x <- runif(dimension, -1, 1)
            xnorm <- sqrt(sum(x^2))
            while (xnorm > 1) {
                x <- runif(dimension, -1, 1)
                xnorm <- sqrt(sum(x^2))
            }
            u[i,] <- x/xnorm
        }
    }
    u
}
```

Easy fix for Borel's paradox in 3-d
-----------------------------------

Take longitude $\phi \sim U(0,2\pi)$ independent of latitude
$\theta = \arcsin(2U-1)$, $U \sim U(0,1)$.

Rotationally Symmetric Distributions
====================================

Comparison of Projected Normal and Langevin Distributions
---------------------------------------------------------

One way that we might compare the $\nlangevin(\mu, \kappa)$ and
$\npn(\gamma\mu, I)$ distributions by choosing *κ*and *γ*to give the
same mean resultant lengths and comparing the densities of the cosine of
the angle *θ*between $U$ and $\mu$.

Of course matching mean resultant lengths is not necessarily the best
way to compare these families of distributions.

$d = 2$
-------

``` {#PNvLvMF2}
par(mar=c(2,2,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plotPNvLvMF(2)
```

![](file:Plots/PNvLvMF2.png)

$d = 3$
-------

``` {#PNvLvMF3}
par(mar=c(2,2,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plotPNvLvMF(3)
```

![](file:Plots/PNvLvMF3.png)

$d = 4$
-------

``` {#PNvLvMF4}
par(mar=c(2,2,0,0)+0.1, oma=c(0,0,0,0)+0.1)
plotPNvLvMF(4)
```

![](file:Plots/PNvLvMF4.png)

Regression
==========

Gould's Model
-------------

A.k.a., the [barber
pole](https://commons.wikimedia.org/wiki/File:Barber-pole-01.gif#)
model.

Gould's Model: Likelihood
-------------------------

Calculate the (profile) log-likelihood for Gould (1969 Biometrics) model
for simple (single predictor) regression with an intercept. For fixed
"slope" *β*, this function "profiles out" (maximizes over) the
"intercept" term and optionally the concentration parameter *κ*.

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
loglklhd.gould <- function(beta, theta, x, do.kappa=FALSE) {
    res <- sapply(beta,
                  function(b, th, x) {
                      sqrt(sum(cos(th - b*x))^2
                           + sum(sin(th - b*x))^2)
                  },
                  th=theta, x=x)
    if (do.kappa) {
        n <- length(theta)
        kappa <- sapply(res/n, imrlLvMF, dimen=2)
        res <- n*log(constLvMF(kappa, dimen=2)) + kappa*res
    }
    res
}
```

Gould's Model with Equally Spaced X
-----------------------------------

``` {#gouldLatticeXPlot1}
<<gouldLatticeXData>>
<<gouldPlot>>
```

``` {#gouldLatticeXPlot2}
<<gouldPlot>>
```

Gould's Model with Equally-Spaced X: Kappa Not Profiled Out
-----------------------------------------------------------

![*κ*not profiled out.](file:Plots/gouldLatticeX1.png)

Gould's Model with Equally-Spaced X: Kappa Profiled Out
-------------------------------------------------------

![*κ*profiled out.](file:Plots/gouldLatticeX2.png)

Gould's Model with Random X: Data Generation
--------------------------------------------

``` {.r}
alpha <- 0
beta <- 1
kappa = 2.5
x <- rnorm(10)
mu <- as.circular((alpha + beta*x) %% (2*pi))
theta <- as.circular(mu + rvonmises(length(mu), mu=0, kappa=kappa))
```

``` {#gouldRandomXPlot1}
<<gouldPlot>>
```

``` {#gouldRandomXPlot2}
<<gouldPlot>>
```

Gould's Model with Random X: Kappa Not Profiled Out
---------------------------------------------------

![*κ*not profiled out.](file:Plots/gouldRandomX1.png)

Gould's Model with Random X: Kappa Profiled Out
-----------------------------------------------

![*κ*profiled out.](file:Plots/gouldRandomX2.png)
