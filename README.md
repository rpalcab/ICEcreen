# ICEscreen

Pipeline designed to analyze bacterial chromosomes in search of distinctive markers facilitating the identification of ICEs.

## Requirements
**Please, create a specific Conda environment from the available yaml file**
```
cd ICEcreen
conda env create -n icecreen --file icecreen.yml
```

- Conda
- Python3
- Bakta
- macsyfinder
- Mobscan
- blastn
- integron_finder
- hmmsearch
- Databases: PlasmidFinder oriT and Plasmid databases, Bakta database, and integrases_ser_try.hmm (a customized Ser/Tyr database available in this repository). Ensure all paths are accurately specified in ice_characterize.sh

## Usage
The file "samples.txt" (mandatory name) contains absolute or relative paths to the FASTA files to be analyzed (preferably chromosomes), separated by newlines.
```
bash ice_characterize.sh -i samples.txt [-t 30] [-o /path/to/output]
python3 ice_parser.py -p ./
```
