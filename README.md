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
domain =...
    [
        x_min y_min
        x_max y_max
    ];
```

The approximation error is controlled through a relative tolerance level:

```
tol = ...
```

Total time consumption can be limited by maximum number of triangulation points:

```
point_num_max = ...;
```

The class can be initialized as follows:

```
sol = GES(...
          f,...
          domain,...
          tol,...
          point_num_max...
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

The algorithm is affected by 2 optional parameters: minimum number of triangulation points for halting (before reaching maximum, default 25) and maximum value of the absolute value flow to consider (default 1).

```
point_num_min = ...;
prop_max = ...;
```

Additionally, batching can be performed to limit memory consumption:

```
batch_size = ...;
```

Displaying of the progress can be turned on optionally:

```
display = ...;
```

Complete algorithm initialization can be represented as:

```
sol = GES(...
          f,...
          domain,...
          tol,...
          point_num_max,...
          PointNumMin=point_num_min,...
          PropMax=prop_max,...
          BatchSize=batch_size,...
          Display=display...
        );
```

Fast visualization of the solution can be done with class method:

```
sol.visTriang;
```

Maximum number of triangulation points can be altered after the initialization:

```
sol.PointNumMax = ...;
sol.fitTriang;
```
