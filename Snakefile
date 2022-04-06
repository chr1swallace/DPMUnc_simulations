localrules: test

wildcard_constraints:
    seed="\d+",
    n="\d+",
    d="\d+",
    k="\d+",

rule simplex_simulations:
    input:
        expand("simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/results.csv",
               seed=range(20),
               d=2,
               noise_factor=[5, 7, 10, 13],
               var_latents=[4, 9, 12, 16]),
        expand("simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/results.csv",
               seed=range(20),
               d=[3, 5, 10, 20],
               noise_factor=[5, 10],
               var_latents=[4, 12]),

rule test:
    input:
        expand("results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/results.csv",
               seed=range(20),
               n=50,
               d=4,
               k=[2, 4, 8],
               var_latents=["1e-1", "5e-1", 1, 2, 3],
               var_means=[1, 2, 4, 6, 9],
               noise_factor=["1e-1", "5e-1", 1, 2, 4])

rule simplex_simulation:
    output:
        clusterAllocations="simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/{seed}/simple/clusterAllocations.csv",
        clusterAllocationsNovar="simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/{seed}/novar/clusterAllocations.csv",
        clusterAllocationsLatents="simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/{seed}/latents/clusterAllocations.csv",
        summary="simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/results.csv",
        rda="simplex_results/{d}_{var_latents}_{noise_factor}/{seed}/results.rda",
    script:
        "scripts/simplex_simulation.R"

rule basic_simulation:
    output:
        clusterAllocations="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/{seed}/simple/clusterAllocations.csv",
        clusterAllocationsNovar="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/{seed}/novar/clusterAllocations.csv",
        clusterAllocationsLatents="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/{seed}/latents/clusterAllocations.csv",
        summary="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/results.csv",
        rda="results/{n}_{d}_{k}_{var_latents}_{var_means}_{noise_factor}_{seed}/results.rda",
    script:
        "scripts/basic_snakemake.R"

#snakemake --snakefile DPMUnc.rules -k -j 1000 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -c {cluster.cpus-per-task} -t {cluster.time} ut {cluster.error} -J {cluster.job} "
