import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from rnanorm import CPM
from rnanorm import TPM



counts_data = pd.read_csv(snakemake.input.counts, sep='\t', index_col=0).T

matrix = CPM().set_output(transform="pandas").fit_transform(counts_data)

matrix.to_csv(snakemake.output.cpm, sep="\t")

tpm = TPM(gtf=snakemake.input.gtf).set_output(transform="pandas")
matrix = tpm.fit_transform(counts_data)

matrix.to_csv(snakemake.output.tpm, sep="\t")
