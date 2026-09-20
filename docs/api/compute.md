# Compute

Trajectory and structure analyses. Import with `from molpy.compute import...`.

Numerical kernels live in the high-performance backend. The public types are identity-style for a stable Python import path — there is no second science
implementation in molpy. Compose **raw Computes** with **Fits** (and an optional
SI scale) the same way the Rust API does.

Like freud’s [API modules](https://freud.readthedocs.io/en/stable/), each
`molpy.compute` module has its own page under [Compute](../compute/index.md)
with an overview table and full signatures. This page is the **index** plus the
shared contract / result types.

!!! note "Analysis units (LAMMPS *real*)"
    Length **Å**, charge **e**, **time fs**, volume Å³, temperature K.
    Vibrational spectra take `dt_fs` in femtoseconds and report cm⁻¹.
    GROMACS trajectories are nm-native — scale lengths ×10 before analysis.
    MSD / Einstein routes need **unwrapped** coordinates.

## Architecture: raw Compute → Fit → scale

| Layer | Role | Examples |
|-------|------|----------|
| Raw Compute | Correlation / MSD / ACF curve only | `EinsteinConductivity`, `GreenKuboConductivity`, `DebyeRelaxation`, `MSD` |
| Fit | Integrate or slope-fit the curve | `CumulativeTrapezoid`, `LinearFit`, `DebyeFit`, `EinsteinHelfandSpectrum`, `GreenKuboSpectrum` |
| Scale | MD → SI prefactor in your script | $1/(6 V k_B T)$, $1/(3 V k_B T)$, $1/d$ |

Two all-in-one recipe classes do exist today — `IonicConductivity` and
`DielectricSusceptibility`, both in `molpy.compute.dielectric` — and both
**violate** the three-layer rule above: a single call unwraps the trajectory,
collects the dipole series, runs the raw Compute, fits the slope *and* applies
the SI prefactor, so the fit window and the unit conversion that decide the
published number are buried in a library default instead of being visible in
your script. They are routed to a `/mol:refactor` that will split them into
primitives and delete the façade. **Do not use them in new code**: compose
`EinsteinConductivity` → `LinearFit` → your own prefactor, as the
[PMSD](../compute/pmsd.md) and [Dielectric](../compute/dielectric.md) pages
show.

Self-diffusion uses `MSD` (Einstein) and `Acf` / `signal.acf_fft` (Green–Kubo);
see the [MSD](../compute/msd.md) and [VACF](../compute/vacf.md) guides.

## Module index

| Module | Primary exports | Guide |
|--------|-----------------|-------|
| `neighborlist` | `NeighborList` | [NeighborList](../compute/neighborlist.md) |
| `rdf` | `RDF` | [RDF](../compute/rdf.md) |
| `density` | `LocalDensity`, `GaussianDensity` | [Density](../compute/density.md) |
| `diffraction` | `StaticStructureFactorDebye` | [Diffraction](../compute/diffraction.md) |
| `pmft` | `PMFTXY` | [PMFT](../compute/pmft.md) |
| `distribution` | distance / angle / dihedral / combined DF | [Distribution](../compute/distribution.md) |
| `spatial` | `SpatialDistribution` | [Spatial](../compute/spatial.md) |
| `order` | Steinhardt family | [Order](../compute/order.md) |
| `environment` | `BondOrder` | [Environment](../compute/environment.md) |
| `shape` | COM, gyration, inertia, $R_g$ | [Shape](../compute/shape.md) |
| `cluster` | `Cluster`, `ClusterCenters`, `ClusterProperties` | [Cluster](../compute/cluster.md) |
| `decomposition` | `DescriptorRow`, `Pca`, `KMeans` | [Decomposition](../compute/decomposition.md) |
| `hbond` | `HBonds`, `HBondCriterion` | [HBond](../compute/hbond.md) |
| `voronoi` | radical Voronoi tessellation | [Voronoi](../compute/voronoi.md) |
| `msd` | `MSD` | [MSD](../compute/msd.md) |
| `pmsd` | `EinsteinConductivity` | [PMSD](../compute/pmsd.md) |
| `jacf` | `GreenKuboConductivity` | [JACF](../compute/jacf.md) |
| `onsager` | `Onsager` | [Onsager](../compute/onsager.md) |
| `persist` | `Persist` | [Persist](../compute/persist.md) |
| `van_hove` | `VanHove` | [Van Hove](../compute/van_hove.md) |
| `reorientation` | `LegendreReorientation` | [Reorientation](../compute/reorientation.md) |
| `dielectric` | dielectric raw/fit helpers | [Dielectric](../compute/dielectric.md) |
| `spectra` | VDOS / IR / Raman / VCD / ROA | [Spectra](../compute/spectra.md) |
| `signal` | `acf_fft`, windows, frequency grid | [Signal](../compute/signal.md) |

---

## Shared types

### The `Compute` contract

`Compute` is not a base class to inherit from. It is a structural
`typing.Protocol` owned by the molrs backend and re-exported here
(`molpy.compute.Compute is molrs.compute.Compute`): any class that defines a
`compute(...)` method satisfies it, with no subclassing and no registration.
Writing one is covered in
[Adding a Compute Operation](../developer/extending-compute.md).

::: molpy.compute.Compute

### Result types

::: molpy.compute.result
