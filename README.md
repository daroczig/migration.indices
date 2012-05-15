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
m   <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
res <- migration.gini(m)
res
```

Resulting in:

```
Total Flows Gini Index:			 0.2222222
Rows Gini Index:			 0
Standardized Rows Gini Index:		 0
Columns Gini Index:			 0.05555556
Standardized Columns Gini Index:	 25
Exchange Gini Index:			 0.05555556
Standardized Exchange Gini Index:	 25
In-migration Field Gini Index:		 vector
Out-migration Field Gini Index:		 vector
```

Where *In* and *Out-migration Field Gini Index* are vectors:

```
> res$migration.gini.in
[1] 0.1000000 0.2500000 0.1666667
> res$migration.gini.out
[1] 0 0 0
```