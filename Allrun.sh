#!/bin/bash
cd "${0%/*}" || exit                                # Run from this directory
cd cpp
sed 's/optimize=.*/optimize=true/' Makefile
make allrun
cd ..
./plot_simulation.sh
./plot_estimation.sh
