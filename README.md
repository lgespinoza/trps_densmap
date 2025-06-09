# TRP Densmap Fingerprints Tutorial

## Requirements

Before starting, make sure you have the following tools installed:

- **GROMACS** [cite]
- **Matplotlib** [cite]

## Input Files

To run the analysis, you will need:

- Molecular dynamics **trajectory** files (`.xtc` or `.trr`)  
  *Note: `.trr` files store more detailed information (e.g., velocities and forces) but are significantly larger. For most analyses, `.xtc` files are preferred due to their reduced size.*
- A **run input** file (`.tpr`)
- An **index** file (`.ndx`)
- The system **topology** file (`.top`)
- All associated **force field parameter files** (`.itp`)
