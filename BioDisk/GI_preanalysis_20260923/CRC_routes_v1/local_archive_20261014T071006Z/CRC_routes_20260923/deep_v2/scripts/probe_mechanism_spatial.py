import pathlib,pandas as pd
B=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');z=pd.read_csv(B/'results/route2/genomewide_descriptive_effects.tsv.gz',sep='\t');genes=['CGNL1','CREB3L1','PCSK5','SLC9B2','TMEM131L','ZC2HC1A','PLOD2','BTNL8'];print(z[z.gene.isin(genes)][['gene','contrast','n_patients','median_delta','positive','negative']].to_string(index=False))
