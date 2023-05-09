# Phylogenetic Statistics


### software setup

git clone "VCFToPhylogenetic depository"

cd VCFToPhylogenetic

conda env create -f environment.yml

conda activate divbrowse_dev

for JC Model
python Phylogenetic_Statistics.py example.fasta example.newick JC 1.000 1.000 1.000 1.000 1.000 1.000

for K2P Model
python Phylogenetic_Statistics.py example.fasta example.newick K2P 1.000 K 1.000 1.000 K 1.000 <br>
K value would be transition rate.


this code would be re-modify from Special aspects of Advanced Algorithms exercise 4, "Phylogeny Inference and Application" (2 SWS), Winter 2020/21 of Freie university berlin.
