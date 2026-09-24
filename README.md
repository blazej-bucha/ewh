# Equivalent water heights at the Earth's surface

This repository offers a Python code to compute equivalent water heights from
GRACE and GRACE-FO monthly gravitational models.  The equivalent water heights
are computed at the Earth's surface to avoid errors of spherical and
ellipsoidal approximations, thereby allowing precise mass loss estimates for
the Greenland and West Antarctic ice sheets.

This repository calculates equivalent water heights from a single monthly GRACE 
solution.  To actually study the Earth's mass transport from GRACE and 
GRACE-FO, you would additionally have to

* consider glacial isostatic adjustment,

* apply corrections to degree-1 spherical-harmonic coefficients,

* replace the GRACE and GRACE-FO coefficients V2,0 and V3,0 by the same
  coefficients from satellite laser ranging,

* process all GRACE and GRACE-FO monthly gravitational models to get monthly 
  snapshots and

* compare the monthly snapshots with respect to a long-term mean model.

This repository does not implement those steps.


# Repository structure

* `./src/ewh.py` -- Python module to compute the surface density and equivalent
  water heights at the Earth's surface using the method described in
  [Citation](#citation).  You should build your mass transport estimators on
  routines from this file.

* `./src/example_run.py` -- A simple Python example showing how to actually
  calculate equivalent water heights at the Earth's surface using the routines
  from `./src/ewh.py` and data from `./data/input`.

* `./data/input/earth-shape.tbl` -- Spherical-harmonic coefficients of the
  Earth's shape up to degree 150.  Synthesis of these coefficients yields the
  full spherical radius of surface points, that is, the distance from the
  origin of the coordinate system (geocenter) to the Earth's surface.  These
  coefficients were used in the paper from the [Citation](#citation) section.

* `./data/input/GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gfc` -- An example 
  monthly GRACE gravitational model developed at the Center for Space Research 
  (CRS).

* `./data/output/ewh-reference.png` -- A reference plot of equivalent water
  heights created by `./src/example_run.py` with the default settings.

* `./data/output/ewh-reference.tbl` -- Reference spherical-harmonic
  coefficients of equivalent water heights created by `./src/example_run.py`
  with the default settings.


# Dependencies

The requirements are

* [numpy](https://numpy.org) to work with arrays,

* [matplotlib](https://matplotlib.org) to plot the outputs and

* [pyharm](https://charmlib.org) to perform spherical-harmonic transforms.

Install the requirements by:

`pip install numpy matplotlib pyharm`


# Example run

Enter the folder with the source files:

```bash
cd src
```

Run the test example:

```bash
python3 example_run.py
```

On success, you should find two files in `data/output`:

* `ewh.png` -- global plot of equivalent water heights for a single month,

* `ewh.tbl` -- spherical-harmonic coefficients of equivalent water heights.

Compare these files with the reference ones:

* `./data/output/ewh-reference.png`,

* `./data/output/ewh-reference.tbl`.


# Citation

If you use the modelling method or its implementation from this repository,
please consider to provide the following reference in your work:

* Bucha, B., **under review**.  Determination of Mass Transport at the Earth's 
  Surface From GRACE and GRACE-FO.

Full reference will be provided as soon as possible.

