# Feynman Decision Diagram: Classical Data Structure for Quantum Circuits

## Build

```commandline
git clone url-of-this-repo
git submodule update --init --recursive
make build test
```

## Execution

```commandline
./build/src/cudd_circuit_bdd -h
./build/src/cudd_circuit_bdd -g gate_sets/default_2.json -s 1 -t 0 -f circuits/ghz.qasm
./build/src/cudd_circuit_bdd -g gate_sets/default_2.json -s 1 -t 2 -f circuits/ghz.qasm
./build/src/cudd_circuit_bdd -g gate_sets/default_2.json -s 1 -t 4 -f circuits/ghz.qasm -f circuits/ccx.qasm
```
