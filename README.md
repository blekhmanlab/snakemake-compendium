# HMC Processing Pipeline

## What is this?
This Snakemake pipeline defines the steps taken to process 16S rRNA gene amplicon sequencing data from a single BioProject, for inclusion in the Human Microbiome Compendium. It is retrieved automatically by the [Compendium Manager](https://github.com/blekhmanlab/compendium) software as part of the compilation process, but this can also be used to process individual projects for which you would like the output to match HMC file formats.

## Setup and configuration

**Installation.** When running this pipeline, either from the command line or as part of a batch job, the dependencies in the `requirements.txt` file must be installed, using commands such as these:

```sh
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

**Config file.** The file at `config/config.yaml` must be updated with two absolute file paths:

* **`resources_path`** describes a directory storing Singularity container images for use in the pipeline (see "Containers," below)
* **`cmanager_path`** describes the directory where Compendium Manager resources are stored. If you are using the Compendium Manager software, this is the root directory of that repository. This must be specified even if you are not using the manager; it also indicates where files should be downloaded for processing.

**Slurm job submission.** The `run_snakemake.example.slurm` file is a suggestion for how one could start this pipeline by submitting it as a batch job on a Slurm system. This is not necessary for the pipeline to run, but if you are going to use a `run_snakemake.slurm` file, it should be modified to fit your computing environment.

**Snakefile.** The Snakefile should also be modified to reflect the specifics of your computing environment. The current implementation of this pipeline assumes it has been launched by the Compendium Manager software on a computing cluster managed by the Slurm workload scheduler, but these are not necessary for the pipeline to work. If you are not using the Compendium Manager, the three conditions at the top of the Snakefile (`onstart`, `onsuccess`, `onerror`) should be removed. These submit Compendium Manager jobs and create files used for tracking progress but are unneeded for running a single instance of the pipeline.

**Accession list.** The pipeline expects a file called `SraAccList.txt` to be included at its root. Each line in this text file should list a single SRA run accession (e.g. `SRR12345`); the file should contain all samples to be processed from a single project.

### Containers

To avoid compatability issues, all steps in the pipeline are currently executed within Singularity containers, but the included scripts should work within either the Singularity or Docker runtimes.

**Singularity.** The containers used here are available as Docker containers, both of which should both be [converted to Singularity images](https://docs.sylabs.io/guides/3.0/user-guide/build_a_container.html) by running commands such as these:

```sh
singularity build sra-tools_3.0.1.sif docker://ncbi/sra-tools:3.0.1
singularity build dada2_1260a.sif docker://blekhmanlab/dada2:1.26.0a
```

The pipeline will look for these SIF files in the directory specified by the `resources_path` configuration value in `config/config.yaml`.

**Docker.** If your system supports Docker rather than Singularity, the included scripts should work as intended. The containers are available directly from Docker Hub:

* [`ncbi/sra-tools`](https://hub.docker.com/r/ncbi/sra-tools), version 3.0.1
* [`blekhmanlab/dada2`](https://hub.docker.com/r/blekhmanlab/dada2), version 1.26.0a

The only required changes will be to the steps invoking the `singularity exec` command, which should be modified to use the equivalent Docker commands instead.

## Executing

Snakemake itself can be invoked using the [conventional approach](https://snakemake.readthedocs.io/en/stable/executing/cli.html) as defined in the documentation:

```sh
snakemake --jobs 45 --slurm --default-resources slurm_account=richard slurm_partition=hello
```

Note that the Snakemake `--use-singularity` flag would enable containers to be defined in each of the rules, rather than invoking Singularity commands directly in the rules themselves. Inconsistencies in how this has run led us to side-step this feature for now.
