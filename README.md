This package provides various indices, like Crude Migration Rate, different Gini indices or the Coefficient of Variation among others, to show the (un)equality of migration.

# Installation

Most recent stable version can be installed directly from CRAN:

```
install.packages('migration.indices')
```

Or please use `devtools` for the development version: [![Build Status](https://travis-ci.org/daroczig/migration.indices.png?branch=master)](https://travis-ci.org/daroczig/migration.indices)

```
library(devtools)
install_github('migration.indices', 'daroczig')
```

# Examples:

The below migration matrix is a demo dataset based on

 * David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262

```r
library(migration.indices)
data(migration.hyp)
res <- migration.gini(migration.hyp)
res
```

Resulting in:

```
Total Flows Gini Index:              0.2222222 
Rows Gini Index:                     0 
Standardized Rows Gini Index:        0 
Columns Gini Index:                  0.05555556 
Standardized Columns Gini Index:     25 
Exchange Gini Index:                 0.05555556 
Standardized Exchange Gini Index:    25 
In-migration Field Gini Index:       vector 
Weighted In-migration Gini Index:    0.1222222 
Out-migration Field Gini Index:      vector 
Weighted Out-migration Gini Index:   0 
Migration-weighted Mean Gini Index:  0.06111111 
```

Where *In* and *Out-migration Field Gini Index* are vectors:

```rout
> res$migration.gini.in
[1] 0.2000000 0.5000000 0.3333333
> res$migration.gini.out
[1] 0 0 0
```

## World dataset

Using the [Global Bilateral Migration Database](http://data.worldbank.org/data-catalog/global-bilateral-migration-database) (2000) of World Bank:

```
data(migration.world)
migration.gini(migration.world)
```

Results in:

```
Total Flows Gini Index:              0.9855174 
Rows Gini Index:                     0.004272225 
Standardized Rows Gini Index:        0.4335008 
Columns Gini Index:                  0.004067787 
Standardized Columns Gini Index:     0.4127565 
Exchange Gini Index:                 1.382263e-05 
Standardized Exchange Gini Index:    0.001402575 
In-migration Field Gini Index:       vector 
Weighted In-migration Gini Index:    0.004241788 
Out-migration Field Gini Index:      vector 
Weighted Out-migration Gini Index:   0.004171575 
Migration-weighted Mean Gini Index:  0.004206681 
```
