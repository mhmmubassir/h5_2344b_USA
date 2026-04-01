# MD phase 1 setup (organized, scripts unchanged)

This folder organizes the first phase of your Amber MD run **without changing any script contents**.

## Folder layout

- `scripts/`
  - `rosetta_to_amber_glynamingv8_finalfor_rosettaoriname.py`
  - `create_inputs.sh`
  - `multi_run_create_inputs.sh`
  - `submit.CPU.sapelo2.sh`
  - `submit.GPU.sapelo2.sh`
- `z.template/`
  - `leapTop.in`
  - `leapBot.in`
  - `1.min.in` … `10.produ.in`
- `y.input_files/`
  - input `.pdb` and `.bond` files referenced by `systemList.txt`
- `systemList.txt`

## Confirmed tleap outputs

From `leapTop.in` + `leapBot.in`, tleap writes:

- `CPLX.pdb`
- `CPLX_Neut_Sol.parm7`
- `CPLX_Neut_Sol.rst7`

It also writes intermediate unsolvated and neutralized files:

- `CPLX.parm7`, `CPLX.rst7`, `CPLX.pdb`
- `CPLX_Neut.parm7`, `CPLX_Neut.rst7`, `CPLX_Neut.pdb`
- `CPLX_Neut_Sol.parm7`, `CPLX_Neut_Sol.rst7`, `CPLX_Neut_Sol.pdb`

## How these original scripts expect to be run

The original scripts assume this relative layout:

- run from the phase-1 root directory
- `systemList.txt` is in the root
- `y.input_files/<name>.pdb` and `y.input_files/<name>.bond` exist
- `create_inputs.sh` is run from a directory where `../z.template/` is visible from each generated run folder

Because of that, the easiest way to use the existing scripts unchanged is:

1. Keep `systemList.txt` in the root.
2. Keep the input structures in `y.input_files/`.
3. Keep `z.template/` in the root.
4. Keep `create_inputs.sh` and `multi_run_create_inputs.sh` in the root **or** call them from there.

## Example logical workflow

1. Prepare a Rosetta-converted Amber/GLYCAM PDB.
2. Place `name.pdb` and `name.bond` in `y.input_files/`.
3. Add `name` to `systemList.txt`.
4. Run `multi_run_create_inputs.sh <replicates> systemList.txt`.
5. Each generated folder runs tleap, creates the Amber topology/restart, then submits CPU equilibration.
6. CPU equilibration submits GPU production.

## Example files included here

This organized package includes example inputs for:

- `9dip_23.amber_renamegly.pdb`
- `9dip_23.amber_renamegly.bond`

matching the included `systemList.txt`.
