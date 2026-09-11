# $z=5$ CDM density and ULDM wavefunction cutouts

The notebook cell **Export $z=5$ density cutouts for the MW-like candidates**
calls `../export_z5_cutouts.py` to create one HDF5 file per halo under `z5/`.
The source snapshot is `run/axionyx_sim_plt00568` in the parent case directory.
Each cutout contains four aligned $51^3$ native-grid fields spanning
$9.9609375\,\mathrm{cMpc}$. The original density and coordinate datasets retain
their values and types when the wavefunction is added.

For each $z=0$ candidate, we trace the CDM particles within $3R_{200c}$ back
to $z\simeq5$ and centre the cutout on the grid cell containing their periodic
centroid. The fixed $51^3$ box includes the surrounding environment; its
boundaries do not represent a $z=5$ halo boundary.

## File structure and coordinates

For example, `z5/halo_0271_z5_cutout.h5` contains:

```text
halo_0271_z5_cutout.h5
├── rho_cdm        (51, 51, 51), float32
├── rho_uldm       (51, 51, 51), float32
├── psi_uldm_real  (51, 51, 51), float64
├── psi_uldm_imag  (51, 51, 51), float64
└── offset_cmpc    (51,),        float64
```

These are datasets inside one HDF5 file, not separate files.

| Dataset | Meaning | Unit |
| --- | --- | --- |
| `rho_cdm` | CDM particle mass deposited with periodic, cell-centred CIC | physical $M_\odot\,\mathrm{kpc}^{-3}$ |
| `rho_uldm` | Native AxioNyx `AxDens`, converted to physical density | physical $M_\odot\,\mathrm{kpc}^{-3}$ |
| `psi_uldm_real` | Native `AxRe`, rescaled to match the exported density normalization | $(M_\odot\,\mathrm{kpc}^{-3})^{1/2}$ |
| `psi_uldm_imag` | Native `AxIm`, rescaled to match the exported density normalization | $(M_\odot\,\mathrm{kpc}^{-3})^{1/2}$ |
| `offset_cmpc` | Cell-centre offsets, shared by the three coordinate axes | comoving Mpc |

The native comoving cell width is

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

Both density fields use physical $M_\odot\,\mathrm{kpc}^{-3}$ at $z=5$;
positions and box sizes are comoving. File attributes record the continuous
particle centroid, the actual centre cell, comoving cell width
`cell_width_cmpc` in cMpc, redshift, target count, and the corresponding $z=0$
halo properties.

## Wavefunction normalization and use

The native fields obey

$$
\rho_{\rm native}=\mathrm{AxDens}=\mathrm{AxRe}^2+\mathrm{AxIm}^2,
$$

where the native density is in $M_\odot\,\mathrm{cMpc}^{-3}$. For the scale
factor $a$ recorded in each file, the exported fields are

$$
\rho_{\rm ULDM}=\mathrm{AxDens}\,\frac{a^{-3}}{10^9},\qquad
\psi=(\mathrm{AxRe}+i\,\mathrm{AxIm})\sqrt{\frac{a^{-3}}{10^9}}.
$$

Thus $\lvert\psi\rvert^2=\rho_{\rm ULDM}$, up to the float32 rounding of the
stored density. The export preserves the native phase and uses the same
periodic cell selection as the density, without interpolation, smoothing, or
phase rotation. The complex field is the nonrelativistic envelope, not the
rapidly oscillating real relativistic scalar field. Dataset attributes identify
the native source fields; file attributes record the conversion factor,
normalization, phase convention, source snapshot, axis order, and parent box
width.

Run this example from the cutouts repository root:

```python
import h5py
import numpy as np

with h5py.File("z5/halo_0271_z5_cutout.h5", "r") as f:
    rho_cdm = f["rho_cdm"][:]
    rho_uldm = f["rho_uldm"][:]
    psi = f["psi_uldm_real"][:] + 1j * f["psi_uldm_imag"][:]
    offset = f["offset_cmpc"][:]
    center = f.attrs["grid_center_cmpc_z5"]
    box = f.attrs["parent_box_width_cmpc"]
    a = f.attrs["scale_factor"]

np.testing.assert_allclose(np.abs(psi)**2, rho_uldm, rtol=2e-7, atol=0)

i, j, k = 25, 26, 25
position_cmpc = (center + np.array([offset[i], offset[j], offset[k]])) % box
cdm_density_at_cell = rho_cdm[i, j, k]
uldm_field_at_cell = psi[i, j, k]
physical_offset_kpc = offset * a * 1000
```



The snapshot has $z\simeq4.9999933$, with physical cell width
$\simeq32.5521\,\mathrm{kpc}$. Sub-cell correlations are not resolved. These cutouts
come from the Level-0 parent run. Different
candidate windows may overlap and should not automatically be treated as
independent realizations.


## Halo selection

All 36 objects passed the hard $z=0$ cuts:

1. $0.8\leq M_{200c}/(10^{12}M_\odot)\leq1.2$;
2. host halo (`Type = 0`);
3. no halo with $M_{\rm vir}\geq M_{\rm vir,host}/3$ inside $3R_{200c}$.

The exported IDs are 13, 14, 40, 63, 77, 97, **107**, 186, 198, **199**,
203, 223, 224, 226, 228, 229, 233, 239, 241, 242, 247, 248, 249, 251,
252, 253, 254, 256, 257, 259, 268, **271**, 272, **273**, 277, and 286.

IDs **107, 199, 271, and 273** are the higher-priority relaxed candidates:
they additionally satisfy $X_{\rm off}/R_{\rm vir}<0.07$ and
$|\eta-1|<0.35$, where $\eta\equiv2T/|U|$.
