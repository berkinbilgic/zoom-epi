# Zoom EPI

MATLAB/[Pulseq](https://pulseq.github.io/) code that generates a 1 mm isotropic
3D-EPI fMRI sequence with a reduced (zoomed) field of view. The sequence also
emits digital trigger outputs for synchronizing transcranial magnetic
stimulation (TMS) hardware with the acquisition.

## Repository layout

```
sequence/
├── script_main_3d_epi_1mm_zoom_R3x2_tms_v150_switch.m   main script (run this)
├── writeEPI_bb_v10_acs_ice_tms_switch.m                  builds the EPI and GRE-ACS sequences
├── writeEPI_bb_v5_acs_ice.m                              single-TR scan for trajectory loading
├── utils/gradSpectrum_bb.m                               gradient spectrum vs. forbidden frequencies
├── requirements.txt                                      Python packages for the CAIPI pattern
└── external/                                             pinned third-party submodules
```

## Requirements

- MATLAB with the Signal Processing Toolbox and the Curve Fitting Toolbox.
  Tested with R2024b and R2026a.
- Python 3, available as `python3` on the `PATH` that MATLAB's `system()`
  sees, with numpy, scipy and matplotlib. Verified with Python 3.9.25,
  numpy 1.23.5, scipy 1.9.3 and matplotlib 3.4.3.

## Getting started

```bash
git clone --recurse-submodules https://github.com/berkinbilgic/zoom-epi.git
cd zoom-epi/sequence
python3 -m pip install -r requirements.txt   # only if the packages are missing
```

In MATLAB, change the current folder to `zoom-epi/sequence` and run the first
section of `script_main_3d_epi_1mm_zoom_R3x2_tms_v150_switch.m`. All paths in
the script are relative to `sequence/`. The `.seq` files are written to
`sequence/<dd-Mmm-yyyy>/`, which git ignores.

If you cloned without `--recurse-submodules`, run
`git submodule update --init` first.

## Default protocol

These are the values set in the main script. The protocol is configured by
editing the script.

| Parameter | Value |
|---|---|
| Resolution | 1 mm isotropic |
| FOV / matrix | 220 × 66 × 64 mm / 220 × 66 × 64 |
| Encoding | 3D EPI; `Ry = 1`, `Rz = 1`, no partial Fourier (TE/TR presets for `Ry` = 2, 3 are in the script) |
| TE / volume TR | 47 ms / 7.04 s (110 ms per partition) |
| Flip angle | Ernst angle for T1 = 1.25 s (23.7°) |
| Excitation | Sinc, 8 ms, time-bandwidth product 30, 3 kHz RF frequency offset; fat saturation on |
| Gradient limits | 28 mT/m, 120 T/m/s (Siemens Skyra) |
| Readout | 3.8 µs ADC dwell (rounded to the 0.1 µs Siemens raster), 1.1 ms echo spacing |

**Slab selection.** `use_gz_excitation = 0` (default) plays the RF pulse
without the scanner's slab-select gradient: the external local coil (ARES)
provides the spatial selection, and a 5.25 ms delay takes the place of the Gz
rephaser. Set `use_gz_excitation = 1` for a conventional Gz slab select with
rephaser.

**Triggers.** `use_trig = 1` adds a 1 ms digital output pulse on `ext1` before
every excitation, followed by a 5.25 ms delay. With the 0.1 ms RF dead time,
the rising edge leads the RF pulse by 6.35 ms. `trig_delay_switch = 1` adds one
more trigger and a 500 ms delay after each EPI volume, for driving the TMS
switch.

**Two configurations.** The script is set up for the EPI-only time series:
one frame, no dummies, no ghost reference and no GRE ACS. Repeat it with the
interpreter's repetitions on the scanner. A commented block in the script
switches to a calibration configuration: one dummy frame, then two frames, the
first of which is the ghost-correction reference (`doRefscan = true`). That
block leaves `use_gre_acs = 0`; set it to 1 to also include the GRE ACS scan.

## Outputs

| File | Contents |
|---|---|
| `3d-epi_RyxRz_<Ry>x<Rz>del<shift>_pf<pf>[_use_gre_acs]_rfFreq_<Hz>_Hz_scaleGrad_acs_<s>.seq` | EPI sequence, with the GRE ACS scan prepended when `use_gre_acs = 1` |
| `3d-epi_..._only_gre_acs_ice_rfFreq_<Hz>_Hz_scaleGrad_acs_<s>.seq` | Standalone GRE ACS scan with ICE labels for Siemens online reconstruction; its ky/kz order is flipped relative to the ACS inside the EPI file. Empty when `use_gre_acs = 0`. |

### CAIPI sampling pattern

`writeEPI_*.m` builds the ky-kz sampling pattern by running

```
python3 external/3DEPI/caipi/skippedcaipi_sampling.py <ny> <nz> <Ry> <Rz> <caipiShiftZ> 1
```

This writes `caipi.mat` to the current folder, which the function then loads.
If the Python call fails, the function prints the command and waits for Enter,
then loads whatever `caipi.mat` is in the current folder. Make sure that file
was generated for the current parameters before you continue.

## Diagnostic sections

The later sections of the main script are optional checks.

- **seqeyes** needs the seqeyes viewer, which is not included. Edit its
  `addpath` line to point to your copy.
- **ICE labels** is only meaningful with `use_gre_acs = 1`.
- **Gradient spectrum** (`utils/gradSpectrum_bb.m`) and **PNS**
  (`external/safe_pns_prediction`) need the Siemens `MP_GPA_*.asc`
  gradient-system file for your scanner. These files are vendor-confidential,
  so they are not included and `*.asc` is git-ignored. Edit the `.asc` paths in
  `gradSpectrum_bb.m` and in the PNS section to point to your local copies.
- **Trajectory** uses `writeEPI_bb_v5_acs_ice.m` to write a single-TR scan
  with a ghost-reference frame for loading the k-space trajectory in
  reconstruction. It also saves the sequence object as a `.mat` file.

## Known messages

- Pulseq warns `Unknown or undefined RF pulse parameter 'use'` for the
  fat-saturation pulse, once per partition.
- The timing check reports `dwell time is not an integer multiple of
  sys.adcRasterTime` for every ADC block. `sys` sets `adcRasterTime = 2 µs`,
  while the dwell is chosen on the 0.1 µs Siemens raster. The `.seq` file is
  still written.

Both messages are also produced by the original working copy of this code.

## Reproducibility

In a given MATLAB environment, the code in this repository writes `.seq`
files that are byte-identical to those of the working-directory code it was
taken from. Across MATLAB versions or sessions, a few samples of the
fat-saturation RF shape can differ in the last printed digit (about 1e-9 of the
peak), so the `.seq` signature hash can differ even though the sequences are
functionally identical.

## Third-party code

Each toolbox is a git submodule pinned to the commit whose files match the
copies this sequence was developed with.

| Submodule (`sequence/external/`) | Upstream | Version | Used for |
|---|---|---|---|
| `pulseq` | [pulseq/pulseq](https://github.com/pulseq/pulseq) | `v1.5.0` | Pulseq MATLAB toolbox |
| `toppe` | [toppeMRI/toppe](https://github.com/toppeMRI/toppe) | `74b0389` | SLR design of the fat-sat pulse, `systemspecs`, `vararg_pair` |
| `PulCeq` | [HarmonizedMRI/PulCeq](https://github.com/HarmonizedMRI/PulCeq) | `ab76fe8` | `trap4ge` |
| `SMS-EPI` | [HarmonizedMRI/SMS-EPI](https://github.com/HarmonizedMRI/SMS-EPI) | `23f7854` | `rf2pulseq` |
| `3DEPI` | [HarmonizedMRI/3DEPI](https://github.com/HarmonizedMRI/3DEPI) | `08e8322` | Skipped-CAIPI sampling pattern (Python) |
| `safe_pns_prediction` | [filip-szczepankiewicz/safe_pns_prediction](https://github.com/filip-szczepankiewicz/safe_pns_prediction) | `c1c46c7` | PNS prediction (SAFE model) |

The EPI generator is adapted from `writeEPI.m` in the day-2 SMS-EPI tutorial
of the ISMRM virtual meeting series *Vendor-Agnostic Pulse Sequence Programming
with Pulseq* (November 15–17, 2023;
[pulseq/ISMRM-Virtual-Meeting--November-15-17-2023](https://github.com/pulseq/ISMRM-Virtual-Meeting--November-15-17-2023)). The CAIPI pattern
uses the skipped-CAIPI code of R. Stirnberg; see Stirnberg & Stöcker, *Magn
Reson Med* 2021;85:1540–1551, https://doi.org/10.1002/mrm.28486.

## License

The code outside `sequence/external/` is released under the MIT License (see
[LICENSE](LICENSE)). Each submodule is covered by its own upstream license.
