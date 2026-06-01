[![License](https://img.shields.io/github/license/UPSY-group/UPSY-models)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-Published-blue.svg)](https://doi.org/10.5194/gmd-18-3635-2025)
[![Paper](https://img.shields.io/badge/Paper-Published-blue.svg)](https://doi.org/10.5194/egusphere-2026-930)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/UPSY-group/UPSY-models)
[![Latest Commit](https://img.shields.io/github/last-commit/UPSY-group/UPSY-models)](https://github.com/UPSY-group/UPSY-models/commits)
[![Workflow Status](https://github.com/UPSY-group/UPSY-models/actions/workflows/UFE_test_suite.yml/badge.svg)](https://github.com/UPSY-group/UPSY-models/actions)
[![Workflow Status](https://github.com/UPSY-group/UPSY-models/actions/workflows/UPSY_test_suite.yml/badge.svg)](https://github.com/UPSY-group/UPSY-models/actions)

![UPSY_logo](UPSY_logo.png)

Welcome to the Utrecht Polar SYstem (UPSY) models repo! Here you will find the UPSY modelling toolkit, the
Utrecht Finite Volume Ice-Sheet Model (UFEMISM), and the One-Layer Antarctic Model for Dynamical
Downscaling of Ice–Ocean Exchanges (LADDIE).

See https://github.com/UPSY-group/UPSY-models/wiki/Getting-started for how to set up your model.

### Python tools
Some tools are available to plot model output on its native mesh. To use these, you need
to install [Miniforge3](https://conda-forge.org/download/) and run:
```
conda env create -f environment.yml
conda install -n base -c conda-forge conda-libmamba-solver
conda config --set solver libmamba
conda activate upsy
python -m pip install -e . --no-deps --no-build-isolation
```
To the package can then be loaded by `import upsy`,
You can also try out these commands in the terminal:
```
upsy-diagnose-run rundir
upsy-plot-2dfigure rundir
upsy-plot-3dfigure rundir
```
For additional help, try `upsy-plot-2dfigure -h`
