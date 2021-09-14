import pandas as pd
import pyBigWig
import numpy as np
nbct_file='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30870/matrix/GSE30870_series_matrix.txt.gz'
md_file='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33233/matrix/GSE33233_series_matrix.txt.gz'
nb_file='http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-Newborn/tracks_hg19/Human_CD4T-Newborn.hmr.bb'
ct_file='http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-100yr/tracks_hg19/Human_CD4T-100yr.hmr.bb'
data='/oak/stanford/groups/smontgom/moqri/data/meth/'
ct_bw='http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-100yr/tracks_hg19/Human_CD4T-100yr.meth.bw'
nb_bw='http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-Newborn/tracks_hg19/Human_CD4T-Newborn.meth.bw'
h_bw='http://smithdata.usc.edu/methbase/data/Xie-Human-2013/Human_H9/tracks_hg19/Human_H9.meth.bw'
i_bw='http://smithdata.usc.edu/methbase/data/Lister-iPSC-2011/Human_ADSiPSC/tracks_hg19/Human_ADSiPSC.meth.bw'
i_bw='http://smithdata.usc.edu/methbase/data/Lister-iPSC-2011/Human_FFiPSC69/tracks_hg19/Human_FFiPSC69.meth.bw'
ezh=pd.read_table('ezh',skiprows=1)
ezh.shape
ggs=[]
ezs=[]
for ch in range(1,22):
    print(ch,end='\n')
    ez=ezh[ezh.chrom=='chr'+str(ch)]
    dfs=[]
    for f in [i_bw,h_bw,nb_bw,ct_bw]:
        df=pyBigWig.open(f)
        df=df.intervals("chr"+str(ch))
        print(len(df))
        df=pd.DataFrame(df)
        df.index=df[0]
        df=df[2]
        dfs.append(df)
    df=pd.concat(dfs,1)
    df.columns=['i','h','n','c']
    ez['r']=ez.apply(lambda x:range(x['chromStart'],x['chromEnd']),1)
    ez=ez['r'].tolist()
    ez=set([item for sublist in ez for item in sublist])
    df['ind']=df.index
    df['ez']=df.ind.apply(lambda x: x in(ez))
    gg=df[~df['ez']].mean()
    ggs.append(gg)
    ez=df[df.h<.5][df['ez']]
    print(ez.shape)
    ezs.append(ez.mean())   
ggs=pd.DataFrame(ggs)
#ggs[['h','i','n','c']].plot(kind='bar')
ezs=pd.DataFrame(ezs)
#ezs[['h','i','n','c']].plot(kind='bar')
ggs.to_csv('ggs')
ezs.to_csv('ezs')