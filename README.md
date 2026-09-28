# MAGs KEGG Analysis and Visualization

Python scripts for evaluating KEGG module completeness across metagenome-assembled genomes (MAGs) and visualizing KO abundance on KEGG pathway maps.

## Scripts

### `module_completeness.py`

Calculates KEGG module completeness for multiple MAGs from per-sample KO annotation files. It produces:

- a tab-separated module-completeness matrix (`<output>.tsv`); and
- a hierarchically clustered heatmap (`<output>.pdf`).

### `kegg_chart.py`

Maps KO counts onto a KEGG pathway image using the coordinates in the corresponding KEGG map configuration file. Rectangles are colored according to KO counts, and pathway lines associated with detected KOs are highlighted.

## Requirements

- Python 3.6 or later
- NumPy
- pandas
- Matplotlib
- seaborn
- OpenCV for Python

Install the required packages with:

```bash
python -m pip install numpy pandas matplotlib seaborn opencv-python
```

## Usage

### 1. Calculate KEGG module completeness

```bash
python module_completeness.py \
  -d path/to/module_chart_files \
  -m path/to/kegg/module \
  -o results/module_completeness
```

Options:

| Option | Description |
| --- | --- |
| `-d` | Directory containing the per-sample KO annotation files. |
| `-m` | KEGG MODULE flat file containing `ENTRY`, `NAME`, and `DEFINITION` records. Always provide this option unless the default path in the script exists on your system. |
| `-o` | Output prefix. The script writes `<output>.tsv` and `<output>.pdf`. |
| `-s` | Input filename suffix. Default: `_Module_chart.csv`. |
| `-h`, `--help` | Display the built-in help message. |

Each input annotation file must be tab-delimited and contain the KO identifier in its second column. For example:

```text
gene_0001	K00001
gene_0002	K00002
gene_0003	K00003
```

With the default suffix, an input directory may contain files such as:

```text
MAG_001_Module_chart.csv
MAG_002_Module_chart.csv
MAG_003_Module_chart.csv
```

The script interprets KEGG MODULE definitions as follows:

- underscore-separated components are all required;
- comma-separated components are alternatives;
- at least 75% of plus-separated subunits must be detected; and
- nested or referenced modules are evaluated from previously calculated module results.

The TSV file contains module completeness values for all samples. The PDF contains a clustered heatmap of modules whose summed completeness across samples is at least 1.

### 2. Visualize KO counts on a KEGG pathway map

```bash
python kegg_chart.py path/to/map00920.conf pathway_ko_counts.tsv
```

The map configuration file and its PNG image must share the same basename and directory:

```text
maps/map00920.conf
maps/map00920.png
```

The KO-count file must be tab-delimited with three columns:

```text
ko00920	K04091	1
ko00920	K00299	2
ko00920	R07210	3
```

| Column | Description |
| --- | --- |
| 1 | KEGG pathway identifier in `koNNNNN` format. |
| 2 | KO or reaction identifier. |
| 3 | Integer count; values greater than 3 are capped at 3 by the script. |

The script reads the PNG associated with the supplied `.conf` file and writes the annotated PNG to the current working directory using the original map filename. To avoid overwriting an existing file, run the script from a separate output directory when the source PNG is also in the current directory.

Batch example:

```bash
mkdir -p colored_maps
cd colored_maps
for conf in ../maps/*.conf; do
  python ../kegg_chart.py "$conf" ../pathway_ko_counts.tsv
done
```

## Suggested workflow

1. Annotate genes from each MAG with KEGG Orthology identifiers.
2. Store each MAG's KO annotations in a tab-delimited `*_Module_chart.csv` file.
3. Run `module_completeness.py` to generate the completeness matrix and heatmap.
4. Prepare a pathway–KO–count table for pathways of interest.
5. Download the matching KEGG map `.conf` and `.png` files.
6. Run `kegg_chart.py` to generate annotated pathway maps.

## Notes

- Despite the `.csv` suffix used by the default input filenames, both scripts expect tab-delimited data.
- `module_completeness.py` was written for Linux-style paths when extracting sample names from filenames.
- KEGG data and pathway maps are subject to KEGG's licensing and usage terms. Users are responsible for obtaining and using KEGG resources appropriately.

