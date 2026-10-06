# AxioNyx results

## $z=5$ CDM density and ULDM wavefunction cutouts

**Export $z=5$ density cutouts for the MW-like halo candidates**
We check the z=5 snapshot to find MW-like halo candidates.
The source snapshot is `run/axionyx_sim_plt00568` in the parent case directory.


### File structure
There are 36 H5 files in the `f=0.01_m=1e-24_z5` folder which means we find 36 MW-like halo candidates (Halo selection creterion is below.)

#### H5 file context

Each cutout contains four aligned $51^3$ comoving-grid fields spanning
$9.9609375\,\mathrm{cMpc}$. The original density and coordinate datasets retain
their values and types when the wavefunction is added.

For each $z=0$ candidate, we trace the CDM particles within $3R_{200c}$ back
to $z\simeq5$ and centre the cutout on the grid cell containing their periodic
centroid. The fixed $51^3$ box includes the surrounding environment, so its
boundaries do not represent a $z=5$ halo boundary.


For example, `f=0.01_m=1e-24_z5/halo_0271_z5_cutout.h5` contains:

```text
halo_0271_z5_cutout.h5
├── rho_cdm        (51, 51, 51), float32
├── rho_uldm       (51, 51, 51), float32
├── psi_uldm_real  (51, 51, 51), float64
├── psi_uldm_imag  (51, 51, 51), float64
└── offset_cmpc    (51,),        float64
```

These are datasets inside each HDF5 file

| Dataset | Meaning | Unit |
| --- | --- | --- |
| `rho_cdm` | Physical CDM density within the grid | physical $M_\odot\,\mathrm{kpc}^{-3}$ |
| `rho_uldm` | AxioNyx `AxDens`, converted to physical density | physical $M_\odot\,\mathrm{kpc}^{-3}$ |
| `psi_uldm_real` | `AxRe`, rescaled to match the exported density normalization | $(M_\odot\,\mathrm{kpc}^{-3})^{1/2}$ |
| `psi_uldm_imag` | `AxIm`, rescaled to match the exported density normalization | $(M_\odot\,\mathrm{kpc}^{-3})^{1/2}$ |
| `offset_cmpc` | Cell-centre offsets, shared by the three coordinate axes | comoving Mpc |


HDF5 attributes

| Attribute | Meaning | Unit |
| --- | --- | --- |
| `grid_center_cmpc_z5` | Centre of cell `[25,25,25]` | cMpc |
| `target_particle_count` | Number of tracked CDM particles selected within $3R_{200c}$ at $z=0$ | Count |
The comoving cell width is

$$
\Delta x=\frac{50\,\mathrm{cMpc}}{256}=0.1953125\,\mathrm{cMpc}.
$$

For indices $i,j,k=0,\ldots,50$,

$$
\mathrm{offset\_cmpc}[i]=(i-25)\Delta x,
$$

All four three-dimensional arrays use `(x, y, z)` index order. The values
`rho_cdm[i,j,k]`, `rho_uldm[i,j,k]`, `psi_uldm_real[i,j,k]`, and
`psi_uldm_imag[i,j,k]` refer to the same cell, with comoving coordinates

$$
(x_i,y_j,z_k)=\bigl(x_c+\mathrm{offset}[i],\,
y_c+\mathrm{offset}[j],\,z_c+\mathrm{offset}[k]\bigr),
$$

wrapped periodically into the $50\,\mathrm{cMpc}$ simulation box. Here
$(x_c,y_c,z_c)$ is the HDF5 attribute `grid_center_cmpc_z5`. Thus
`offset_cmpc` is a coordinate displacement in cMpc, not a density correction.
It is one shared coordinate vector: use index `i` for the x offset, `j` for
the y offset, and `k` for the z offset. `offset_cmpc[25] = 0`, so `[25,25,25]`
is the central cell. The first and last cell centres are at offsets
$-4.8828125$ and $+4.8828125\,\mathrm{cMpc}$; the full side length includes
the outer half-cell on each end.



### Wavefunction normalization and use

The fields obey

$$
\rho=\mathrm{AxDens}=\mathrm{AxRe}^2+\mathrm{AxIm}^2,
$$

where the comoving density is in $M_\odot\,\mathrm{cMpc}^{-3}$. 

Run this example from the cutouts repository root:

```python
import h5py
import numpy as np

with h5py.File("f=0.01_m=1e-24_z5/halo_0271_z5_cutout.h5", "r") as f:
    rho_cdm = f["rho_cdm"][:]
    rho_uldm = f["rho_uldm"][:]
    psi = f["psi_uldm_real"][:] + 1j * f["psi_uldm_imag"][:]
    offset = f["offset_cmpc"][:]
    center = f.attrs["grid_center_cmpc_z5"]
    box = np.full(3, 50.0)  # Parent box width, cMpc.
    a = 0.166666852930862  # Source snapshot scale factor.

np.testing.assert_allclose(np.abs(psi)**2, rho_uldm, rtol=2e-7, atol=0)

i, j, k = 25, 26, 25
position_cmpc = (center + np.array([offset[i], offset[j], offset[k]])) % box
cdm_density_at_cell = rho_cdm[i, j, k]
uldm_field_at_cell = psi[i, j, k]
physical_offset_kpc = offset * a * 1000
```


### Halo selection

Here's the halo selection creterion:

1. $0.8\leq M_{200c}/(10^{12}M_\odot)\leq1.2$;
2. host halo (`Type = 0`);
3. no halo with $M_{\rm vir}\geq M_{\rm vir,host}/3$ inside $3R_{200c}$.

The exported IDs are 13, 14, 40, 63, 77, 97, **107**, 186, 198, **199**,
203, 223, 224, 226, 228, 229, 233, 239, 241, 242, 247, 248, 249, 251,
252, 253, 254, 256, 257, 259, 268, **271**, 272, **273**, 277, and 286.

IDs **107, 199, 271, and 273** are the higher-priority relaxed candidates:
they additionally satisfy $X_{\rm off}/R_{\rm vir}<0.07$ and
$|\eta-1|<0.35$, where $\eta\equiv2T/|U|$.

## Results of level4-AMR simulation of MW-like halo 271

$m_{\rm ULDM}=10^{-24}\,\mathrm{eV}$ and $f=0.01$.
Snapshots: $z\simeq5,3,1$.
This is the snapshot of a AMR zoom in simulation. The zoom in region contains a target MW-like halo.

The boundary of the snapshot is defined by the smallest cube that contains particles that within $3R_{200}$ of the MW-like halo at $z=0$.

### h5 file structure

Directory: `f=0.01_m=1e-24_amr_zoomin_halo271/`.
Files contain consecutive rows, with up to 1,000,000 cells per part and less than 100 MB per file.

| Folder | Files | Total cells |
| --- | ---: | ---: |
| `z5/` | 105 | 104,096,403 |
| `z3/` | 87 | 86,174,195 |
| `z1/` | 42 | 41,421,736 |


| Dataset | Shape | Meaning / unit |
| --- | --- | --- |
| `cells/center_cmpc` | $(N,3)$ | comoving cell-center coordinates, Mpc |
| `cells/width_cmpc` | $(N,)$ | comoving cubic cell side length,  Mpc |
| `cells/level` | $(N,)$ | Actual AMR level |
| `cells/rho_cdm` | $(N,)$ | physical CDM density,  $M_\odot\,\mathrm{kpc}^{-3}$ |
| `cells/rho_uldm` | $(N,)$ | physical ULDM density,  $M_\odot\,\mathrm{kpc}^{-3}$ |
| `cells/psi_uldm_real` | $(N,)$ | Wavefunction real part, $(M_\odot\,\mathrm{kpc}^{-3})^{1/2}$ |
| `cells/psi_uldm_imag` | $(N,)$ | Wavefunction imaginary part, same unit |



### Reading

```python
from pathlib import Path
import h5py

folder = Path("f=0.01_m=1e-24_amr_zoomin_halo271/z5")
for path in sorted(folder.glob("*_part*.h5")):
    with h5py.File(path, "r") as f:
        cells = f["cells"]
        for start in range(0, len(cells["level"]), 65536):
            rows = slice(start, start + 65536)
            center = cells["center_cmpc"][rows]
            width = cells["width_cmpc"][rows]
            level = cells["level"][rows]
            rho_cdm = cells["rho_cdm"][rows]
            rho_uldm = cells["rho_uldm"][rows]
            psi = cells["psi_uldm_real"][rows] + 1j * cells["psi_uldm_imag"][rows]
```
