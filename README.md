# kmer2probe2detection
software to generate k-mer databases, and molecular (in situ) RNA/DNA detection assays 

INSTALLATION
to install it on a Win10 PC, using conda:
create a fresh virtual envrironment: conda create -n my_shiny_new_env 
activate the just created environment: conda activate my_shiny_new_env
cd to the directory where requirements.txt is located.
run: pip install -r requirements.txt in your shell.
if pip cannot be found, try: conda install -c anaconda pip, then re-run previous pip command
note: requires Python 3.6. or later

RUNNING THE SOFTWARE
set parameters in params.py
type: python main_v4.py makekmerdb (to generate a kmerdb from fasta file(s)
or type: python main_v4.py makeprobes cdnapadlocks (or rnapadlocks or rnasmfish) to create probes


