# TRP Densmap Fingerprints Tutorial

## 📦 Requirements

Before starting, make sure you have the following tools installed:

- **GROMACS** – for trajectory processing and density map generation  
- **Matplotlib** – (optional) for custom visualization of `.xpm` outputs

## 📁 Input Files

To run the analysis, you will need the following:

- Molecular dynamics **trajectory** files (`.xtc` or `.trr`)  
  *Note: `.trr` files store more detailed information (e.g., velocities and forces) but are significantly larger. For most analyses, `.xtc` files are preferred due to their reduced size.*
- A **run input** file (`.tpr`)
- An **index** file (`.ndx`)
- The system **topology** file (`.top`)
- All associated **force field parameter files** (`.itp`)

---

## 🧪 Densmap Analysis Workflow

This section describes how to process multiple MD replicas to generate 2D density maps (fingerprints) for the **PVP2** group (e.g., PIP₂ molecules) using GROMACS.

### 🔖 Input Assumptions

- Trajectories: `rep1.xtc`, `rep2.xtc`, `rep3.xtc` (in `xtc/`)
- Run input file: `tpr/rep1.tpr`
- Index file: `index.ndx`, containing a group named `PVP2`

### ⚙️ Preprocessing and Densmap Steps

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
```

---

##  Output Files

- `repX_pip2_map.xpm` – 2D density map for each replica  
- `merged_pip2_map.xpm` – Combined density map from all replicas  
- Intermediate `.xtc` files (`_step1`, `_step2`, `_step3`) can be deleted after final output

---

##  Python Plotting Script

This script reads GROMACS `.xpm` density map files and converts them into `.svg` images using a consistent colormap (`plasma`) and a fixed value range (`vmin=0.0`, `vmax=0.4`). This allows reliable visual comparison across multiple density maps.

###  Features

- **Colormap**: `plasma`  
- **Scale**: fixed between 0.0 and 0.4  
- **Input**: all `.xpm` files in the current directory  
- **Output**: one `.svg` per `.xpm` file

###  Usage

```bash
python plot_xpm.py
```

### Script

```python
import os
import re
import matplotlib.pyplot as plt

def parse_xpm(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header_line = None
    for line in lines:
        if line.startswith('"') and len(line.split()) == 4:
            header_line = line
            break

    if not header_line:
        raise ValueError("Header line not found")

    header = header_line.strip().strip('"').split()
    width, height, num_colors = int(header[0]), int(header[1]), int(header[2])

    color_map = {}
    for i in range(lines.index(header_line) + 1, lines.index(header_line) + 1 + num_colors):
        parts = lines[i].strip().strip('"').split()
        color_code = parts[0]
        match = re.search(r'/\* "(.*?)" \*/', lines[i])
        if match:
            value = float(match.group(1))
            color_map[color_code] = value
        else:
            raise ValueError(f"Value not found for color code {color_code}")

    vectors = []
    for line in lines[lines.index(header_line) + 1 + num_colors:]:
        if line.strip().startswith('"'):
            vectors.append(line.strip().strip('"'))

    x_axis = []
    y_axis = []
    for line in lines:
        if line.startswith('/* x-axis:'):
            x_axis.extend(map(float, line.split(':')[1].strip().rstrip('*/').split()))
        if line.startswith('/* y-axis:'):
            y_axis.extend(map(float, line.split(':')[1].strip().rstrip('*/').split()))

    return vectors, color_map, x_axis, y_axis

def convert_to_matrix(vectors, color_map):
    matrix = []
    for vector in vectors:
        row = [color_map[char] for char in vector if char in color_map]
        matrix.append(row)
    return matrix

def plot_heatmap(matrix, x_axis, y_axis, output_file):
    plt.figure(figsize=(10, 8))
    heatmap = plt.imshow(
        matrix,
        aspect='auto',
        cmap='plasma',
        origin='lower',
        extent=[x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]],
        vmin=0.0,
        vmax=0.4
    )
    cbar = plt.colorbar(heatmap, label='Density')
    cbar.ax.tick_params(labelsize=12, width=2, colors='black')
    cbar.ax.yaxis.set_tick_params(labelsize=20, width=2, colors='black')
    plt.xlabel('x (nm)', fontsize=20, fontweight='bold')
    plt.ylabel('y (nm)', fontsize=20, fontweight='bold')
    plt.xticks(fontsize=20, fontweight='bold')
    plt.yticks(fontsize=20, fontweight='bold')
    plt.title('Heatmap of PVP2 Number Density', fontsize=16, fontweight='bold')
    plt.savefig(output_file, format='svg')
    plt.close()

# Batch process all .xpm files
for filename in os.listdir('.'):
    if filename.endswith('.xpm') and not filename.endswith('.eps'):
        print(f"Processing {filename}...")
        vectors, color_map, x_axis, y_axis = parse_xpm(filename)
        matrix = convert_to_matrix(vectors, color_map)
        output_svg = filename.replace('.xpm', '.svg')
        plot_heatmap(matrix, x_axis, y_axis, output_svg)
        print(f"SVG saved: {output_svg}")
```
