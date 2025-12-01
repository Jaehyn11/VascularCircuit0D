# README

## User customize in source code (IMPORTANT!)

ctrl+f "setme" and change the working directory.
(both in `solver_serial.cpp` and `plot.py`)

## Input file - RLC serial chain circuit. User input parameters.

`VascularCircuit0D/input/input.inp`

## Compilation/Build:

`g++ -O2 solver_serial.cpp -o solver_serial.exe`

## Running 0D Simulation:

`./solver_serial.exe`

### Output

`VascularCircuit0D/output.csv`

## Post-processing:

### Run:

`python plot.py`

### Output:

1. $Q_i:$ `VascularCircuit0D/flow.png`

2. $P_i:$ `VascularCircuit0D/pressure.png`
