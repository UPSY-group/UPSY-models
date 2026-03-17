![UPSY_logo](UPSY_logo.png)

Welcome to the Utrecht Polar SYstem (UPSY) models repo! Here you will find the UPSY modelling toolkit, the
Utrecht Finite Volume Ice-Sheet Model (UFEMISM), and the One-Layer Antarctic Model for Dynamical
Downscaling of Ice–Ocean Exchanges (LADDIE).

See https://github.com/UPSY-group/UPSY-models/wiki/Getting-started for how to set up your model.

### Python tools
Some tools are available to plot model output on its native mesh. To use these, you need
to install [Miniforge3](https://conda-forge.org/download/) and run:
```
conda env create -y -n upsy -f environment.yml
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
