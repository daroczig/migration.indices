This package provides different indices to show (un)equality of migration.

# Installation

Please use `devtools`:

```
library(devtools)
install_github('migration.indeces', 'daroczig')
```

# Examples:

The below migration matrix is a demo dataset based on

 * David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262

```
data(migration.hyp)
res <- migration.gini(m)
res
```

Resulting in:

```
Total Flows Gini Index:	                 0.2222222
Rows Gini Index:                         0
Standardized Rows Gini Index:            0
Columns Gini Index:                      0.05555556
Standardized Columns Gini Index:         25
Exchange Gini Index:                     0.05555556
Standardized Exchange Gini Index:        25
In-migration Field Gini Index:           vector
Out-migration Field Gini Index:          vector
```

Where *In* and *Out-migration Field Gini Index* are vectors:

```
> res$migration.gini.in
[1] 0.1000000 0.2500000 0.1666667
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
Total Flows Gini Index:            0.9855174
Rows Gini Index:                   0.004272225
Standardized Rows Gini Index:      0.4335008
Columns Gini Index:                0.004067787
Standardized Columns Gini Index:   0.4127565
Exchange Gini Index:               1.382263e-05
Standardized Exchange Gini Index:  0.001402575
In-migration Field Gini Index:     vector
Out-migration Field Gini Index:    vector
```