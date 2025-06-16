# TRP Densmap Fingerprints Tutorial

##  Requirements

Before starting, make sure you have the following tools installed:

- **GROMACS** – for trajectory processing and density map generation  
- **Matplotlib** – (optional) for custom visualization of `.xpm` outputs

##  Input Files

To run the analysis, you will need the following:

- Molecular dynamics **trajectory** files (`.xtc` or `.trr`)  
  *Note: `.trr` files store more detailed information (e.g., velocities and forces) but are significantly larger. For most analyses, `.xtc` files are preferred due to their reduced size.*
- A **run input** file (`.tpr`)
- An **index** file (`.ndx`)
- The system **topology** file (`.top`)
- All associated **force field parameter files** (`.itp`)

---

##  Densmap Analysis Workflow

This section describes how to process multiple MD replicas to generate 2D density maps (fingerprints) for the **PVP2** group (e.g., PIP₂ molecules) using GROMACS.

### Input Assumptions

- Trajectories: `rep1.xtc`, `rep2.xtc`, `rep3.xtc` (in `xtc/`)
- Run input file: `tpr/rep1.tpr`
- Index file: `index.ndx`, containing a group named `PVP2`

###  Preprocessing and Densmap Steps

Each replica is:

1. **Unwrapped** (`-pbc whole`)
2. **Centered** (`-center`)
3. **Aligned** in the XY plane (`-fit rotxy`)
4. **Processed** with `gmx densmap` to compute the 2D density projection

Afterward, the replicas are **merged**, and a final combined density map is computed.

---

## ▶️ Bash Script

Save the following as `run_densmap.sh` and make it executable:

```bash
#!/bin/bash

# === VARIABLES ===
reps=(rep1 rep2 rep3)
xtc_dir="xtc"
tpr="tpr/rep1.tpr"
index="index.ndx"
group_pvp2="PVP2"

# === PROCESS EACH REPLICA ===
for rep in "${reps[@]}"; do
    echo ">>> Processing $rep..."

    # Step 1: Remove PBC jumps
    echo 0 | gmx trjconv -f $xtc_dir/${rep}.xtc -s $tpr -pbc whole -o ${rep}_step1.xtc

    # Step 2: Center the system
    echo 1 0 | gmx trjconv -f ${rep}_step1.xtc -s $tpr -pbc mol -center -o ${rep}_step2.xtc

    # Step 3: Fit XY plane
    echo 16 0 | gmx trjconv -f ${rep}_step2.xtc -s $tpr -fit rotxy -o ${rep}_step3.xtc -n $index

    # Step 4: Compute individual densmap
    echo "$group_pvp2" | gmx densmap -f ${rep}_step3.xtc -s $tpr -n $index -bin 0.2 -o ${rep}_pip2_map.xpm
done

# === MERGE REPLICAS ===
echo ">>> Concatenating aligned trajectories..."
gmx trjcat -f rep1_step3.xtc rep2_step3.xtc rep3_step3.xtc -o traj_merged.xtc -cat

# === FINAL DENSITY MAP ===
echo "$group_pvp2" | gmx densmap -f traj_merged.xtc -s $tpr -n $index -bin 0.2 -o merged_pip2_map.xpm


