#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

m_0=-1.75
L=100
x_values=(0.03 0.04 0.05 0.15)
periodic_values=(0)
disorder_strength_values=(0.005 0.01 0.5)
site_configs=(1 2 3 4 5)
disorder_configs=(1 2 3 4 5)

total_runs=$((${#x_values[@]} * ${#periodic_values[@]} * ${#site_configs[@]} * ${#disorder_configs[@]} * ${#disorder_strength_values[@]}))
current_run=0

for x in "${x_values[@]}"; do
    for periodic in "${periodic_values[@]}"; do
        for site_config in "${site_configs[@]}"; do
            site_seed=$((1000 + site_config))
            for disorder_config in "${disorder_configs[@]}"; do
                disorder_seed=$((2000 + 100 * site_config + disorder_config))
                # Reuse each configuration across strengths; only its amplitude changes.
                for disorder_strength in "${disorder_strength_values[@]}"; do
                    ((current_run += 1))
                    echo "[$current_run/$total_runs] Running proportion_sites=$x, x_periodic=$periodic, y_periodic=$periodic, site_seed=$site_seed, disorder_seed=$disorder_seed, disorder_strength=$disorder_strength"

                    julia_args=(
                        m_0="$m_0"
                        L="$L"
                        proportion_sites="$x"
                        x_periodic="$periodic"
                        y_periodic="$periodic"
                        disorder_strength="$disorder_strength"
                        site_seed="$site_seed"
                        disorder_seed="$disorder_seed"
                    )
                    if [[ -n "${RUN_ID:-}" ]]; then
                        julia_args+=(run_id="$RUN_ID")
                    fi

                    julia --startup-file=no PTB_chern_BHZ_random_joint_LDOS_disorder_chern.jl "${julia_args[@]}"
                done
            done
        done
    done
done

echo "Completed all $total_runs runs."
