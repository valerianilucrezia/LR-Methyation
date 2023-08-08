# LR-Methyation
Analysis scripts for comparison between EPIC and Nanopore methylation data


## Requirements
 `R`

Required libraries:
*  `ggplot2 `
*  `optparse `
*  `dplyr `
*  `patchwork `
*  `ggpubr `
*  `ggplogo `
*  `stringr`
*  `tidyr`
*  `gridExtra`

## How to run
### Analysis
```
Rscript ./scripts/main.R  \ 
        --nanopore $np_file \
        --epic $epic_file \
        --output $output_directory \
        --sample $sample_name 
```

### Run Full analysis with a **job array** in a cluster with SLURM scheduler
Modify `run_analysis.sbatch` based on number of samples and your resources, then:

```
sbatch run_analysis.sbatch sample_ids
```
where `sample_ids` is a file where listed all samples id with **path**. 

### File requirments 
`$np_file` 
[modbam2bed](https://github.com/epi2me-labs/modbam2bed) output file with the following **13** fields (from https://github.com/epi2me-labs/modbam2bed):
(**highlighted** fields are required, the other can be filled with any value)

* **1**	reference sequence name
* **2**	0-based start position
* **3**	0-based exclusive end position (invariably start + 1)
* 4	Abbreviated name of modified-base examined
* **5**	"Score" 1000 * (Nmod + Ncanon) / (Nmod + Ncanon + Nno call + Nalt mod + Nfilt + Nsub + Ndel). The quantity reflects the extent to which the calculated modification frequency in Column 11 is confounded by the alternative calls. The denominator here is the total read coverage as given in Column 10.
* **6**	Strand (of reference sequence). Forward "+", or reverse "-".
* 7-9	Ignore, included simply for compatibility.
* **10**	Read coverage at reference position including all canonical, modified, undecided (no calls and filtered), substitutions from reference, and deletions. Nmod + Ncanon + Nno call + Nalt mod + Nfilt + Nsub + Ndel
* 11	Percentage of modified bases, as a proportion of canonical and modified (excluding no calls, filtered, substitutions, and deletions). 100 * Nmod / (Nmod + Nalt mod + Ncanon)
* **12**	Ncanon
* **13**	Nmod

`$epic_file` 
output of methylation calling with **2** columns containing:
* probes id
* beta value
