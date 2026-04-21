# SCRATCH-TCR

## Introduction
SCRATCH-TCR is a comprehensive and scalable **Nextflow DSL2** pipeline for end-to-end **single-cell T-cell receptor (TCR) analysis**, unifying VDJ quality control, TCR-aware integration with gene expression, multi-algorithm clonotype clustering, and repertoire profiling within a reproducible framework that enables biologically meaningful insights into T-cell phenotypes, clonal expansion, diversity, and immune dynamics.


## Prerequisites
Before running the subworkflow, ensure you have the following installed:
- [Nextflow](https://www.nextflow.io/) (version 21.04.0 or higher)
- [Java](https://www.oracle.com/java/technologies/javase-downloads.html) (version 8 or higher)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) for containerized execution
- [Git](https://git-scm.com/)

## Installation
Clone the repository to your local machine:
```bash
git clone https://github.com/WangLab-ComputationalBiology/SCRATCH-TCR.git
cd SCRATCH-TCR
```


## Overview

### SCRATCH-TCR starts from:

- **Cell Ranger VDJ outputs** (`outs/` directories)
- an **annotated single-cell object** (for example a Seurat object)
- a **sample sheet**

## The pipeline then performs:

1. **VDJ QC**
2. **T-cell integration**
3. Optional downstream analyses:
   - **TCRi**
   - **CoNGA**
   - **GLIPH2**
   - **TCRdist3**
   - **GIANA**
   - **Consensus clustering**
   - **Repertoire analysis**
4. **Master summary report**

The workflow is orchestrated from `main.nf`, where VDJ QC runs first, T-cell integration creates the common baseline, downstream modules run conditionally, consensus clustering aggregates selected clustering outputs, and master summary collects outputs across modules. :contentReference[oaicite:3]{index=3} :contentReference[oaicite:4]{index=4}

---

## Key capabilities

- Modular **DSL2** architecture
- Containerized execution with **Docker** or **Singularity**
- Support for **LSF** profile-based execution
- TCR-centered downstream analysis from annotated GEX + VDJ data
- Configurable metadata mapping for:
  - labels
  - sample IDs
  - patient IDs
  - condition
  - timepoint
  - batch
- Optional execution of major downstream methods through parameter toggles
- Publication/report-style module outputs via Quarto-based reporting paths configured in `nextflow.config` :contentReference[oaicite:5]{index=5}

---

## Workflow structure

### Main stages

#### 1. VDJ QC
Performs QC and filtering of VDJ contigs and generates curated summary tables and figures for downstream use and final reporting. :contentReference[oaicite:6]{index=6}

#### 2. T-cell integration
Builds the common integrated T-cell baseline using QC-filtered VDJ data and the annotated object. All major downstream modules operate from this shared output. :contentReference[oaicite:7]{index=7}

#### 3. Downstream modules
The following modules can be enabled or disabled independently:

- **TCRi**
- **CoNGA**
- **GLIPH2**
- **TCRdist3**
- **GIANA**
- **Repertoire**

These run after T-cell integration and produce module-specific reports and exports. :contentReference[oaicite:8]{index=8} :contentReference[oaicite:9]{index=9}

#### 4. Consensus clustering
Consensus clustering combines outputs from GLIPH2, TCRdist3, and GIANA when enabled. :contentReference[oaicite:10]{index=10}

#### 5. Master summary
Aggregates selected outputs across the workflow into a final summary report. :contentReference[oaicite:11]{index=11}

---

## Parameters
### Mandatory Parameters

The root workflow checks for the following required inputs:

- `--input_annotated_object`: Path to annotated object, such as an .RDS
- `--input_vdj_contigs`: Path or glob to VDJ output directories
- `--sample_sheet`: Path to sample sheet CSV

If these are not supplied, the workflow exits with an error.

## Basic usage

`Run with Singularity'

nextflow run main.nf -profile singularity \
  --input_vdj_contigs "VDJ/*/outs" \
  --input_annotated_object "project_annotated_object.RDS" \
  --sample_sheet "samplesheet.csv" \
  --project_name "SCRATCH_TCR_run" \
  --outdir "results"
  
`Run with Docker`

```bash
nextflow run main.nf -profile docker \
  --input_vdj_contigs "VDJ/*/outs" \
  --input_annotated_object "project_annotated_object.RDS" \
  --sample_sheet "samplesheet.csv" \
  --project_name "SCRATCH_TCR_run" \
  --outdir "results"
```

`Run selected modules only`

```bash
nextflow run main.nf -profile singularity \
  --input_vdj_contigs "VDJ/*/outs" \
  --input_annotated_object "project_annotated_object.RDS" \
  --sample_sheet "samplesheet.csv" \
  --run_tcri false \
  --run_conga false \
  --run_gliph2 true \
  --run_tcrdist3 true \
  --run_giana true \
  --run_consensus true \
  --run_repertoire false \
  --run_master_summary true
```

Module toggles are exposed in nextflow.config

### Important parameters
#### Metadata mapping

The workflow supports configurable metadata fields:

`--label_col`
`--sample_col`
`--patient_col`
`--condition_col`
`--timepoint_col`
`--batch_col`
`--Global plotting / embedding`


### Global resource:

max_cpus
max_memory
max_time

These are applied through a helper function that constrains module requests.


### Module-specific configuration

The config defines detailed parameter blocks for:

VDJ QC
T-cell integration
TCRi
CoNGA
GLIPH2
TCRdist3
GIANA
Consensus clustering
Repertoire
Master summary

This makes the pipeline flexible while preserving a consistent workflow entry point.


## Outputs

On successful completion, the workflow reports result locations for:

VDJ_QC
TCell_Integration
TCRi
CoNGA
GLIPH2
TCRdist3
GIANA
Repertoire
Consensus_Clustering
Master_Summary

## Configuration
The subworkflow can be configured using the `nextflow.config` file. Modify this file to set default parameters, profiles, and other settings. An institution profile should be created whenever running the pipeline in an HPC environment, please refer to [Step-by-step guide to writing an institutional profile](https://nf-co.re/docs/tutorials/use_nf-core_pipelines/config_institutional_profile)

## Output
Upon successful completion, the results will be available in a directory named after your project (`<project_name>`). You can open the report in your browser:
```plaintext
Done! Open the following report in your browser -> <path/to/launchDir>/<project_name>/report/index.html
```

## Documentation
For more detailed documentation and advanced usage, refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) and the comments within the subworkflow script (`main.nf`).

## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

## License
This project is available under the GNU General Public License v3.0. See the LICENSE file for more details.

## Contact
For questions or issues, please contact:
- sazaidi@mdanderson.org
- lwang22@mdanderson.org
 
