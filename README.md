# $z=5$ CDM and ULDM density cutouts

The notebook cell **Export $z=5$ density cutouts for the MW-like candidates**
creates one HDF5 file per halo under `z5/`. Each file contains two aligned
$51^3$ native-grid fields spanning $9.9609375\,\mathrm{cMpc}$:

- `rho_cdm`: CDM particle mass deposited with periodic, cell-centred CIC;
- `rho_uldm`: the AxioNyx `AxDens` field without additional interpolation;
- `offset_cmpc`: a length-51 coordinate array giving the cell-centre offsets
  along each spatial axis.

The native comoving cell width is

$$
\Delta x=\frac{50\,\mathrm{cMpc}}{256}=0.1953125\,\mathrm{cMpc}.
$$

For indices $i,j,k=0,\ldots,50$,

$$
\mathrm{offset\_cmpc}[i]=(i-25)\Delta x,
$$

and the cell `rho_cdm[i,j,k]` or `rho_uldm[i,j,k]` has comoving coordinates

$$
(x_i,y_j,z_k)=\bigl(x_c+\mathrm{offset}[i],\,
y_c+\mathrm{offset}[j],\,z_c+\mathrm{offset}[k]\bigr),
$$

wrapped periodically into the $50\,\mathrm{cMpc}$ simulation box. Here
$(x_c,y_c,z_c)$ is the HDF5 attribute `grid_center_cmpc_z5`. Thus
`offset_cmpc` is a coordinate displacement in cMpc, not a density correction.

Both density fields use physical $M_\odot\,\mathrm{kpc}^{-3}$ at $z=5$;
positions and box sizes are comoving. File attributes record the continuous
particle centroid, the actual centre cell, comoving cell width
`cell_width_cmpc` in cMpc, redshift, target count, and the corresponding $z=0$
halo properties.

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
$|\eta-1|<0.35$, where $\eta\equiv2T/|U|$. The other 32 are retained for halo-to-halo comparisons,
but at least one of these two relaxation diagnostics lies outside the adopted
reference range. This does not prove that they are unsuitable: the parent
$195\,\mathrm{kpc}$ grid does not resolve their centres well, and a single
$z=0$ catalog cannot by itself rule out a recent major merger. Merger histories
are therefore required before treating those objects as clean relaxed halos.
