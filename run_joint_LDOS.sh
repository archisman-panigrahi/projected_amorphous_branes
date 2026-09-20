#!/usr/bin/env bash

set -euo pipefail

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_directory"

m_0_values=(-1.0 1.0)
L=100
x_values=(0.05 0.10 0.15 0.20 0.25)
periodic_values=(0 1)

for m_0 in "${m_0_values[@]}"; do
    for x in "${x_values[@]}"; do
        for periodic in "${periodic_values[@]}"; do
            echo "Running m_0=$m_0, proportion_sites=$x, x_periodic=$periodic, y_periodic=$periodic"

            julia --startup-file=no PTB_chern_BHZ_random_sites_joint_LDOS.jl \
                m_0="$m_0" \
                L="$L" \
                proportion_sites="$x" \
                x_periodic="$periodic" \
                y_periodic="$periodic"
        done
    done
done
