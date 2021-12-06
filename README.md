# tasean
## Target Sequence Analysis
Tasean is a handy tool for checking a certain sequence in a collection of viral and prokaryotic genomes. 
## Setup (for Linux environment only)
### Dependencies
- Python 3.6
- bwa
- muscle (included)
- samtools
- fastTree 2.1
- biopython
#### Install via conda:
`conda create -n tasean`

`conda activate tasean`

`conda install -y -c conda-forge -c bioconda -c anaconda -c etetoolkit -c defaults --file requirements.txt`

## Usage
`python tasean.py [--gene] -s /path/target/sequence.fasta -g /path/to/genome/genomes.fasta -o /output/folder`

 	-- gene : indicate that target sequence is a gene
 	-s, --seq : target sequence (ref)
 	-g,--genomes: all genomes in fasta format (each samples is a sequence)
 	-o,--output: output folder 
example:
 `python tasean.py --gene -s spike.fna -g Delta1000.fasta -o data`
