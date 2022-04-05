rule basic_simulation:
    output:
        clusterAllocations="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/{seed}/simple/clusterAllocations.csv",
        clusterAllocationsNovar="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/{seed}/novar/clusterAllocations.csv",
        clusterAllocationsLatents="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/{seed}/latents/clusterAllocations.csv",
        pca="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/pca.png",
        summary="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/results.csv",
        rds="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/results.rds",
    script:
        "scripts/basic_snakemake.R"
