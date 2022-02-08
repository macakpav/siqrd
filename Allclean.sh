#!/bin/bash
cd "${0%/*}" || exit                                # Run from this directory
rm -f *.aux *.log plot_simulation.pdf plot_estimation.pdf
cd cpp
make clean
cd ..

