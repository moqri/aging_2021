import pandas as pd
import pyBigWig
from bisect import bisect_left

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
            dg=dg.merge(df.drop(['g','tss','d','dq'],1),left_index=True,right_index=True,how='left')
        dg.drop(['g','tss','d'],1).mean().plot(kind='bar')     
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
def get_flank():  
    df['ch']='chr'+df.index.str.split('_').str[0]
    df['b']=df.index.str.split('_').str[1].astype(int)
    df['b1']=df.b-1
    df['b2']=df.b+3
    df[['ch','b1','b2']].to_csv(data+'wg/flank_pos.bed',sep='\t',header=None,index=None)
    #!bedtools getfasta -fi ref/hg38.fa.masked -bed wg/flank_pos.bed -fo wg/flank.fasta
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