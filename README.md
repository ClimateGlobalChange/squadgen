```
\file  README.md
\author  Paul Ullrich
\version May 8, 2024
```

Copyright 2024 Paul Ullrich

This file is distributed as part of the Tempest source code package.
Permission is granted to use, copy, modify and distribute this
source code and its documentation under the terms of the GNU General
Public License.  This software is provided "as is" without express
or implied warranty.

# Spherical Quadrilateral Mesh Generator (SQuadGen)

## REQUIRED EXTERNAL LIBRARIES

libnetcdf

## BUILDING

Modify `src/Makefile` to specify NetCDF paths

```sh
make all
```

## INSTALLING FROM CONDA PACKAGE

An alternative to building SQuadGen from source is to install a pre-built Linux
verison from a conda package.

If you don't already have a conda base environment set up, you can install
[miniforge3](https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge3).

To add SQuadGen to an existing conda environment, activate the environment and
run:
```sh
conda install squadgen
```

To create a new conda environment with SQuadGen in it (and any other packages
you like, run:
```sh
conda create -n <env_name> squadgen <package> <package>
```
Then, activate the environment with:
```sh
conda activate <env_name>
```
The `SQuadGen` executable will be in your path.

## SYNOPSIS
```sh
  ./SQuadGen <Parameter List>
```

## PARAMETERS

```
  --grid_type <string> ["CS"] (Options: CS | ICO | OCT1 | OCT2)
  --refine_type <string> ["LOWCONN"] (Options: LOWCONN | CUBIT | LOWCONNOLD)
  --refine_level <integer> [2]
  --resolution <integer> [10]
  --refine_file <string> [""]
  --refine_rect <string> [""]
  --output <string> [""]
  --loadcsrefinementmap <bool> [false]
  --smooth_type <string> ["NONE"] (Options: NONE | SPRING | PRESSURE)
  --smooth_dist <integer> [1] (Smooth distance, -1 = smooth entire mesh)
  --smooth_iter <integer> [10]
  --lon_base <double> [-180.000000]
  --lat_base <double> [0.000000]
  --x_rotate <double> [0.000000]
  --y_rotate <double> [0.000000]
  --lon_ref <double> [0.000000]
  --lat_ref <double> [0.000000]
  --orient_ref <double> [0.000000]
  --tessellate <integer> [0]
  --subcellres <integer> [0]
  --invert <bool> [false]
  --block_refine <bool> [false]
```

## DESCRIPTION

  `SQuadGen` is a mesh generation utility that uses a cubed-sphere base mesh to
  generate quadrilateral meshes with user-specified enhancements.  In order to
  determine where enhancement is desired, the user provides a PNG file which
  corresponds to a latitude-longitude grid.  Raster values with higher brightness
  (whiter values) are tagged for refinement.  The algorithm uses a basic paving
  technique and supports two paving stencil types: Low-connectivity (`LOWCONN`)
  and `CUBIT`-type transition regions.

## OPTIONS

`--grid_type (CS) | ICO`
<dd>

  Type of base mesh to use for refinement: cubed-sphere (`CS`) or icosahedral
  (`ICO`).  If `ICO` is specified, generate a basic icosahedral flag mesh.  All
  refinement criteria will be ignored in this case.
</dd>

`--refine_type (LOWCONN) | CUBIT | LOWCONNOLD`
<dd>

  Type of paving stencil to use `LOWCONN`, `CUBIT` or `LOWCONNOLD`.  `LOWCONN` provides
  a few enhancements over the previous `LOWCONN` old template by removing
  spurious elements with high aspect ratios (typically near interior corners).
</dd>

`--refine_level <integer>`
<dd>

  Number of levels of refinement to use in grid generation. Each refinement
  level corresponds to a 2x refinement of the mesh.
</dd>

`--resolution <integer>`
<dd>

  Base resolution of the mesh. For cubed-sphere meshes this corresponds to the
  number of faces along each edge of the cubed-sphere.
</dd>

`--refine_file <filename>`
<dd>

  Filename of the PNG file to use for specifying refinement regions.
</dd>

`--refine_rect <string>`
<dd>

  An argument specifying rectangular regions on the cubed-sphere where
  refinement should be applied. Each refine_rect should be separated by a
  semi-colon and consists of five arguments:

  ```
  <lon1>,<lat1>,<lon2>,<lat2>,<level>
  ```

  where `<lon1>,<lat1>` are the longitude-latitude coordinates of the lower-right
  corner of the refinement region, `<lon2>,<lat2>` are the longitude-latitude
  coordinates of the upper-right corner of the refinement region, and `<level>`
  is the number of levels of refinement.
</dd>

`--output <filename>`
<dd>

  Filename for the output EXODUS-format grid file.
</dd>

`--loadcsrefinementmap`
<dd>

  (ADVANCED) If specified, the refinement map will be reloaded from the
  previously generated refine_map.dat file.  This option allows for manual
  editing of the cubed-sphere refine map.
</dd>

`--smooth_type (NONE) | SPRING`
<dd>

  If `SPRING` is specified, use spring-dynamics type mesh smoothing once the
  refined mesh has been generated.
</dd>

`--smooth_dist <integer>`
<dd>

  By default (smooth dist 1) only nodes which are part of the transition region
  will be smoothed.  Larger values of smooth_dist includes nodes of distance
  smooth_dist away from the transition region (by edge connectivity).  If
  smooth_dist is set to -1 then the entire mesh will be subjected to smoothing.
</dd>

`--smooth_iter <integer>`
<dd>

  Number of smoothing iterations to perform.  A larger value will typically
  result in a smoother mesh.
</dd>

`--lon_base <double>`
<dd>

  Longitude coordinate (in degrees) of the left-most edge of the PNG image.
</dd>

`--lon_shift <double>`
<dd>

  Amount of rotation (in degrees) to apply to the base mesh.
</dd>

`--lon_ref <double>`
<dd>

  The central longitude of the first panel of the cubed-sphere (in degrees).
</dd>

`--lat_ref <double>`
<dd>

  The central latitude of the first panel of the cubed-sphere (in degrees).
</dd>

`--orient_ref <double>`
<dd>

  The orientation of the cubed sphere at `<lon_ref>,<lat_ref>` in degrees.
  A value of 0 means the cubed-sphere is aligned perfectly zonal at the
  reference point. A value of 45 means the cubed-sphere is inclined by 45
  degrees at this point.
</dd>

`--invert`
<dd>

  Invert the base image (so black identifies regions of refinement).
</dd>

`--block_refine`
<dd>

  Turn on block-wise tagging of elements for refinement.
</dd>
