# Global Equation Solver

## Definition

Global Equation Solver (GES) is a global complex zeros and poles finding algorithm analyzing logarithm of the function. For a complex function $f = f(z)$, algorithm approximates complex points $z^*$, such that $f(z^\*) = 0$ or $[f(z^\*)]^{-1}=0$. The function $f$ is assumed to be meromorphic in the inspected domain, however, output of the algorithm can be beneficial in other cases as well.

The theoretical basis of the algorithm and some applications are described in following articles:
1. Viktor A. Frantsuzov, Anton V. Artemyev, "A global argument-based algorithm for finding complex zeros and poles to investigate plasma kinetic instabilities", Journal of Computational and Applied Mathematics, [link](http://dx.doi.org/10.1016/j.cam.2024.116217)
2. **WIP**

Refer to these publications if the algorithm is used in a scientific work.

## Representation

The solver is represented as a class in MATLAB (see examples).

$f$ can be set as an anonymous function:
```
f = @(z) ...;
```

The inspected domain has to be represented as a matrix:
```
domain = ...
    [
        x_min y_min
        x_max y_max
    ];
```

The class can be initialized as follows:
```
sol = GES( ...
          f, ...
          domain ...
        );
```

Output contains following general properties:
```
sol.Domain %analyzed domain
sol.DomainNorm %domain normalization
sol.Tol %approximation relative tolerance level
sol.Func %analyzed function
sol.DT %Delaunay triangulation
sol.FuncEval %function evaluations
sol.PointNum %triangulation points number
sol.CandPoint %candidate points
sol.CandRegion %candidate regions
```

The algorithm is affected by 4 optional parameters: approximation relative tolerance level (limits time consumption), maximum number of triangulation points (limits time consumption), minimum number of triangulation points (prevents early halting) and maximum value of the absolute value flow to consider (prevents early halting).
```
tol = ...; %positive scalar with value less than 1  (default 1 / 2^10)
point_num_max = ...; %positive integer (default 0, no limit)
point_num_min = ...; %positive integer (default 25)
prop_max = ...; %positive scalar (default Inf)
```

The initial mesh can be refined by adding points to it:
```
add_point = ...; %vector (default empty)
```

Additionally, batching can be performed to limit memory consumption:
```
batch_size = ...; %nonnegative integer (default 0, no batching)
```

Displaying of the progress can be turned on optionally:
```
display = ...; %either "off" or "on" (default "off")
```

Complete algorithm initialization can be represented as:
```
sol = GES( ...
          f, ...
          domain, ...
          Tol=tol, ...
          PointNumMax=point_num_max, ...
          PointNumMin=point_num_min, ...
          PropMax=prop_max, ...
          AddPoint=add_point, ...
          BatchSize=batch_size, ...
          Display=display ...
        );
```

Some class properties can be altered after the initialization:
```
sol.Tol = ...;
sol.PointNumMax = ...;
sol.PointNumMin = ...;
sol.PropMax = ...;
sol.BatchSize = ...;
sol.Display = ...;

sol.fitTriang; %refit the triangulation
```

Fast visualization of the solution can be done through the class method:
```
sol.visTriang;
```

The type of region to visualize can be set through an optional argument:
```
region_type = ... %subset of [0 -1 1] (unknown, poles and zeros, respectively, default [-1 1])
sol.visTriang(region_type);
```


