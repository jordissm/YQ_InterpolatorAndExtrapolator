# Yield interpolator/extrapolator in $Y_Q$

This is a program to perform interpolation and extrapolation predictions for yield (and yield ratios) of charged particles in heavy-ion collisions using arbitrary-order polynomials.

## How to build and install the project
### Prerequisites

| Software   | Required version |
|------------|-----------------|
| cmake      |                 |
| yaml-cpp   |                 |
| OpenMP     |                 |

### Compilation and installation
On the project directory, run
```shell
bash build.sh -DUSE_OMP=OFF
```
where the `-DUSE_OMP` can be set to `ON` if OpenMP is available on your system.

## Running the program

### Input

The input is divided between a single YAML and one or more data files. Data files need to be placed on the `data` subdirectory, where example files are provided.

### Execution (local)
```shell
./YQ_InterpolatorAndExtrapolator
```

### Execution (SLURM Cluster)
```shell
sbatch submit_job.sh
```

### Output

## License
