# Scientific exchange contract · version 1.0

The runtime accepts UTF-8 JSON records up to 4,000,000 bytes. Non-finite numbers,
unknown schema versions, malformed curves/tables and incompatible input maps
are rejected. All numerical quantities use the explicit units stated in their
field definitions; display-unit conversions occur inside the owning model.

## Material: `polymer-material/1.0`

Required: `schema`, `sample_id`, `name`, `chemistry` and `conditions`.
`parameters` and `history` are lists and default to empty lists.
Chemistry entries are text, at most 4,000 characters each. Sample identifiers
and names contain 1–200 characters. Temperature, if present, is positive kelvin.

```json
{
  "schema": "polymer-material/1.0",
  "sample_id": "example-01",
  "name": "Illustrative polymer sample",
  "chemistry": {
    "repeat_unit": "unspecified",
    "architecture": "linear",
    "tacticity": "unspecified",
    "sequence": "unspecified",
    "end_groups": "unspecified",
    "notes": "Demonstration metadata; no calibrated chemistry is implied."
  },
  "conditions": {
    "temperature_K": 298.15,
    "state": "unspecified",
    "solvent": "unspecified"
  },
  "parameters": [],
  "history": []
}
```

An optional `molar_mass_distribution` contains equal-length lists
`molar_mass_g_mol` and `number_fraction`, positive masses, normalized nonnegative
number fractions and `basis: "number_fraction"`. Do not relabel mass fractions
as number fractions. Chain Lab converts explicitly selected source conventions.

Derived parameter records carry `key`, `value`, `unit`, `quantity`, `origin`
and `model_scope`. Model calculations append derived records and lineage; they
cannot overwrite sample identity. Repeated parameters can represent different
conditions or models. They are not an unqualified global source of coefficients.

## Calculation: `polymer-result/1.0`

The record contains `run_id`, `created_utc`, `application` (id/name/version),
`material`, exact `inputs`, `results`, `references`, an interpretation statement
and `checksum_sha256`. Results include metrics, plots, tables, notes, warnings,
material updates and explicit handoffs. Plot X/Y arrays must match and contain
at most 1,200 finite points per series. Table rows must match their columns.

The checksum is SHA-256 of UTF-8 JSON with sorted keys, compact separators,
Unicode preserved, non-finite numbers forbidden, and the checksum field omitted.
It is an integrity check, not authentication. Changing a calculated value requires
a new calculation or a separate documented measurement record; it must not be
presented as the original run.

## Workspace: `polymer-session/1.0`

Contains the current `material` and up to 50 complete calculation `runs`.
Each run is independently checked. The same 4 MB document limit applies.
Saving an excessive workspace reports an error; individual runs can still be
exported. Imported records remain inspectable even if their application is not
installed. Unknown application code is never dynamically imported from a file.

## Transfers

`POST /api/exchange/transfer` accepts a completed result and target application
ID. It applies only an explicit handoff emitted by the originating scientific
module. Target field names must exist in the receiving module's contract.
The receiving model validates values and units under its own input definitions.
Without an explicit map, only the material record and lineage transfer.

Current physical mapping: Network → Transport. `thickness_um` is the **full
equilibrium film thickness**, `diffusivity_m2_s` is the user-supplied effective
solute diffusivity, with `mode="release"`, `temperature_K`, `duration_s` and
`n_points`. The geometry is fixed and fully equilibrated before diffusion.

There is no tracer-D → Fickian-D map and no MWD → effective network density map.
These omissions preserve scientific meaning rather than substitute heuristics.

## Local API

- `GET /api/catalog`: installed app metadata and model field definitions.
- `GET /api/health`: runtime version and installed app IDs.
- `POST /api/calculate/{id}`: `{"inputs": {...}, "material": {...}}`.
- `POST /api/exchange/validate`: a material/result/workspace JSON record.
- `POST /api/exchange/transfer`: `{"document": {...}, "target": "transport"}`.

Unknown inputs and invalid scientific states return readable 422 responses.
Missing applications return 404. Foreign origins/hosts are rejected; the CLI
serves loopback only. The runtime has no remote job submission or credential
storage. This API is intended for local research workflows.

## Version 0.2 molecular-view extension

The backward-compatible `results.molecular_views` list stores numerical 2D drawing coordinates and explicit chemical graph connectivity. Each view has `title`, `caption`, optional `smiles`, `atoms` and `bonds`. Atoms carry unique integer `index`, `element`, finite `x`/`y`, optional chemical `label`, formal charge, isotope, chirality and integer `unit_index`. Bonds reference atom indices with `begin`/`end`, numeric order 1, 1.5, 2 or 3, and optional `stereo`/`direction` metadata. The runtime rejects duplicate indices, invalid endpoints, duplicate bonds, non-finite coordinates and invalid orders before rendering or accepting imported result records.

The GUI constructs SVG elements from this data; arbitrary SVG/HTML supplied by a calculation or import is not inserted. A molecular drawing does not imply equilibrium geometry, force-field assignment or a stereochemical assignment beyond the supplied template. Existing 1.0 result files without molecular views remain readable.

Material `chemistry` metadata now includes optional `repeat_smiles`, `stereochemistry_source`, `branch_topology`, `functional_groups`, `synthesis_route` and `conversion_basis` strings. Existing records lacking them remain valid. These fields preserve evidence; text entries are not interpreted as physical model coefficients.
