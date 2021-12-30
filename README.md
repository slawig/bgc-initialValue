# bgc-initialValue

Marine ecosystem models are important to identify the processes that affects for example the global carbon cycle.

For many marine ecosystem models, the existence and uniqueness of periodic solutions has not yet been analytically proven. Thus, we implemented a generator for different initial concentrations and started simulations with these different initial concentrations.



## Installation

To clone this project with **git** run:
>git clone https://github.com/slawig/bgc-initialValue.git



## Usage

The project consists of the Python packages initialValue and Python scripts in the directory InitialValueSimulation to run a simulation using different initial concentrations.


The Python packages util is available in the repository https://github.com/slawig/bgc-util.git.



### Python package initialValue

This package summarizes the functions to generate different initial concentrations for a marine ecosystem model and to start the spin-up calculations with different (generated) initial concentration.



### Python scripts

Python scripts exist for the application to start the simulations and evaluate them.

The scripts to use different initial concentrations are available in the directory `InitialValueSimulation`. The script `InitialValue.py` starts the spin-up calculations using different initial concentrations which can be identified with parameters. The script `InitialValue_Plot.py` can be used to visualize the results.



## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
