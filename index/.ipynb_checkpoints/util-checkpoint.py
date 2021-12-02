import pandas as pd
import pyBigWig
from bisect import bisect_left
import matplotlib.pyplot as plt
import seaborn as sns

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before
    
def bw2df(bw_url,tss):
    bw=pyBigWig.open(bw_url)
    bws=[]
    count=0
    for i,r in tss.iterrows():
        count+=1
        if (count % 200) == 0:        
            print(r['chrom'],end='.')
        bwi=bw.intervals("chr"+str(r['chrom']),r['chromStart'],r['chromEnd'])
        if bwi is not None:
            for b in bwi:
                bws.append([r['chrom'],b[0],b[2]])
    df=pd.DataFrame(bws)
    df.drop_duplicates(inplace=True)
    df.index=df[0].astype(str)+'_'+df[1].astype(str)
    df.drop([0,1],1,inplace=True)
    df.columns=['beta']
    print(df.shape)
    return df

def get_df(bw_url,tss):

    return df


def run_etl():
    ezh=get_ezh('al')
    dfs=[]
    cells=['h9','ips','nb','ct','sy','so']
    for cell in cells:
        df=etl(cell,ezh)
        df.to_csv(data+cell)
        dfs.append(df)
    for i in range(6):
        dg=dfs[0]
        for df in dfs[1:]:
            dg=dg.merge(df.drop(['g','tss','d','dq'],axis=1),left_index=True,right_index=True,how='left')
        dg.drop(['g','tss','d'],axis=1).mean().plot(kind='bar')     
def merge_dfs():
    dfs=[]
    cells=['h9','ips','nb','ct','sy','so']
    for cell in cells:
        df=pd.read_csv(data+cell,index_col=0)
        #df.to_csv(data+cell)
        dfs.append(df)
    inc=dfs[0].index
    for i in range (1,6):
        ini=dfs[i].index
        inc=inc.intersection(ini)   
    for i in range (6):
        df=dfs[i]
        df=df.loc[inc]
        dfs[i]=df    
    for i in range (6):
        dg=dfs[0]
        dg[cells[i]]=dfs[i].iloc[:,0]   
    return dg
def get_cpg():
    man_='ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
    cpg=pd.read_csv(man_,skiprows=7)
    return cpg

def hg382hg19():
    df=pd.read_csv('ezs.csv',index_col=0)
    df['ch']=df.index.str.split('_').str[0]
    df['b']=df.index.str.split('_').str[1]
    df['e']=df.index.str.split('_').str[1]
    df[['ch','b','e','d']].to_csv('ezh2_38.bed',index=None,header=None,sep='\t')
    #!CrossMap.py bed hg38ToHg19.over.chain.gz ezh2_38.bed ezh2_19    
def rep():
    rep='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM923nnn/GSM923447/suppl/GSM923447_hg19_wgEncodeUwRepliSeqImr90ValleysRep1.bed.gz'
    rep=pd.read_table(rep,header=None)
    rep[[0,1,2]].to_csv(data+'rep19.bed',sep='\t',index=False)
    #!CrossMap.py bed hg19ToHg38.over.chain.gz "{data}rep19.bed" "{data}rep38.bed"
    rep=pd.read_table(data+'rep38.bed',header=None)
    rep=rep[~rep[0].isin(['chrX','chr22_KI270879v1_alt','chr4_GL000008v2_random','chr14_GL000009v2_random'])]
    rep['ch']=rep[0].str[3:].astype(int)
    rep.shape
    rep['g']=rep.ch*10**9+rep[1]
    l=[list(range(g,g+1000)) for g in rep.g]
    l = [item for sublist in l for item in sublist] 
    

    
def get_fig(df,cell,label,pal):
    sns.set(rc={'figure.figsize':(3.5,2.5)})
    pl=df[cell+['dq']].groupby('dq').mean()[cell]
    pl.columns=label
    pl.index=list(range(-2500,2501,500))
    plt.figure()
    ax=sns.lineplot(data=pl,dashes=False,hue_order=label[::-1],
                    palette = sns.color_palette(pal))
    plt.figure()
    dl=df[df.h9<.2]
    pl=dl.groupby('dq').mean()[cell]
    pl.columns=label
    pl.index=list(range(-2500,2501,500))
    ax=sns.lineplot(data=pl,dashes=False,hue_order=label[::-1],
                    palette = sns.color_palette(pal),legend=False)
    plt.figure()
    dle=dl[dl.ezh.notna()].drop(['g','tss','d'],axis=1)
    print(dle.shape)
    pl=dle.groupby('dq').mean()[cell]
    pl.columns=label
    pl.index=list(range(-2500,2501,500))
    ax=sns.lineplot(data=pl,dashes=False,hue_order=label[::-1],
                    palette = sns.color_palette(pal),legend=False)
    plt.figure()
    dh=df[df.h9>.6]
    pl=dh.groupby('dq').mean()[cell]
    pl.columns=label
    pl.index=list(range(-2500,2501,500))
    ax=sns.lineplot(data=pl,dashes=False,hue_order=label[::-1],
                    palette = sns.color_palette(pal),legend=False)
    plt.figure()
    dho=dh[dh.flank.str[0].isin(['A','T'])&dh.flank.str[3].isin(['A','T'])].drop(['g','tss','d'],axis=1)
    pl=dho.groupby('dq').mean()[cell]
    pl.columns=label
    pl.index=list(range(-2500,2501,500))
    ax=sns.lineplot(data=pl,dashes=False,hue_order=label[::-1],
                    palette = sns.color_palette(pal),legend=False)    