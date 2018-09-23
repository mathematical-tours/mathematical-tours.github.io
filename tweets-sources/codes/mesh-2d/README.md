## `MESH2D: Delaunay-based mesh generation in MATLAB`

`MESH2D` is a `MATLAB` / `OCTAVE`-based unstructured mesh-generator for two-dimensional polygonal geometries, providing a range of relatively simple, yet effective two-dimensional meshing algorithms. `MESH2D` includes variations on the "classical" Delaunay refinement technique, a new "Frontal"-Delaunay refinement scheme, a non-linear mesh optimisation method, and auxiliary mesh and geometry pre- and post-processing facilities. 

<p align="center">
  <img src = "../master/poly-data/lake-1.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "../master/poly-data/lake-2.png">
</p>

Algorithms implemented in `MESH2D` are "provably-good" - ensuring convergence, geometrical and topological correctness, and providing guarantees on algorithm termination and worst-case element quality bounds. Support for user-defined "mesh-spacing" functions and "multi-part" geometry definitions is also provided, allowing `MESH2D` to handle a wide range of complex domain types and user-defined constraints. `MESH2D` typically generates very high-quality output, appropriate for a variety of finite-volume/element type applications.

`MESH2D` is a simplified version of my <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> mesh-generation algorithm (a `C++` code). `MESH2D` aims to provide a straightforward `MATLAB` / `OCTAVE` implementation of these Delaunay-based triangulation and mesh optimisation techniques. 

## `Code Structure`

`MESH2D` is a pure `MATLAB` / `OCATVE` package, consisting of a core library + associated utilities:

    MESH2D::
    ├── MAIN-DIR. -- core MESH2D library functions. See REFINE2, SMOOTH2, TRIDEMO, etc.
    ├── aabb-tree -- support for fast spatial indexing, via tree-based data-structures.
    ├── geom-util -- geometry processing, repair, etc.
    ├── hfun-util -- mesh-spacing definitions, limiters, etc.
    ├── hjac-util -- solver for Hamilton-Jacobi eqn's.
    ├── mesh-ball -- circumscribing balls, orthogonal balls etc.
    ├── mesh-cost -- mesh cost/quality functions, etc.
    ├── mesh-file -- mesh i/o via ASCII serialisation.
    ├── mesh-util -- meshing/triangulation utility functions.
    ├── poly-data -- polygon definitions for demo problems, etc.
    └── poly-test -- fast inclusion test for polygons.

## `Starting Out`

After downloading and unzipping the current <a href="https://github.com/dengwirda/mesh2d/archive/master.zip">repository</a>, navigate to the installation directory within <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> and run the set of examples contained in `tridemo.m`:
```
tridemo( 0); % a very simple example to get everything started.
tridemo( 1); % investigate the impact of the "radius-edge" threshold.
tridemo( 2); % Frontal-Delaunay vs. Delaunay-refinement refinement.
tridemo( 3); % explore impact of user-defined mesh-size constraints.
tridemo( 4); % explore impact of "hill-climbing" mesh optimisations.
tridemo( 5); % assemble triangulations for "multi-part" geometries.
tridemo( 6); % assemble triangulations with "internal" constraints.
tridemo( 7); % investigate the use of "quadtree"-type refinement.
tridemo( 8); % explore use of custom, user-defined mesh-size functions.
tridemo( 9); % larger-scale problem, mesh refinement + optimisation. 
tridemo(10); % medium-scale problem, mesh refinement + optimisation. 
```

Note that a call to `libpath()` should be included in `MESH2D` scripts to ensure the `MATLAB` / `OCTAVE` path variable includes all `MESH2D` sub-directories.

For <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> users, performance can be improved by compiling elements of the `MESH2D` library. See `mkoctfile` for additional information (note: requires a `-dev` installation of <a href="https://www.gnu.org/software/octave">`OCTAVE`</a>).

## `References!`

If you make use of `MESH2D` please include a reference to the following! `MESH2D` is designed to provide a simple and easy-to-understand implementation of Delaunay-based mesh-generation techniques. For a much more advanced, and fully three-dimensional mesh-generation library, see the <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> package. `MESH2D` makes use of the <a href="https://github.com/dengwirda/aabb-tree">`AABBTREE`</a> and <a href="https://github.com/dengwirda/find-tria">`FINDTRIA`</a> packages to compute efficient spatial queries and intersection tests. 

`[1]` - Darren Engwirda, <a href="http://hdl.handle.net/2123/13148">Locally-optimal Delaunay-refinement and optimisation-based mesh generation</a>, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, September 2014.

`[2]` - Darren Engwirda, Unstructured mesh methods for the Navier-Stokes equations, Honours Thesis, School of Aerospace, Mechanical and Mechatronic Engineering, The University of Sydney, November 2005.


