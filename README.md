# nCov2019_Guide_Design
## 1. Create a set of group-specific guide RNAs ("Horizontal pool")
**Extract a comprehensive set of guide RNAs that target conserved regions of a fasta alignment file**
```
python horizontal/horizontal.py aligned_sequences.fasta -o ./horizontal_pool_output
```

## 2. Create a family-covering guide pool ("Minipool")
**Generate a very small pool of guide RNAs that target all genomes in a fasta file**
```
python minipool/prepare_guide_candidates.py sequences.fasta -o ./minipool_output
python minipool/set_cover.py -o ./minipool_output
```