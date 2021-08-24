import os, shutil, glob
import pysam
import sys
from Bio import SeqIO
import os
import time
import vcf
import math
from cigar import Cigar
import hashlib
from shutil import copyfile
from datetime import datetime
import argparse
protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"
           }
def translateToProtMap(dna):
    map_protein={}
    count=0
    print(len(dna))
    if len(dna)%3 == 0:
        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]

            count=count+1
            map_protein[i+1]={'codon':str(codon),'p':protein[dna[i:i+3]],'pos':count}
    return map_protein
def translateToProt(dna):
    protein_sequence = ""
    if len(dna)%3 == 0:
        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]
            protein_sequence+= protein[codon]

    return protein_sequence
def GlobalAlignment(v,w,score,sigma):
    if v=='' or w=='':
        return '',0
    s=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    bk=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    s[0][0]=0
    for i in range(1,len(v)+1):
        s[i][0]=s[i-1][0]-sigma
    for j in range(1,len(w)+1):
        s[0][j]=s[0][j-1]-sigma
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            # s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]+score[v[i-1]][w[j-1]])
            if not v[i-1]==w[j-1]:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]-1)
            else:
                s[i][j]=s[i-1][j-1]+score
            if s[i][j]==s[i-1][j]-sigma:
                bk[i][j]=1
            elif s[i][j]==s[i][j-1]-sigma:
                bk[i][j]=2
            else:
                bk[i][j]=3
    AlignV=''
    AlignW=''
    i=len(v)
    j=len(w)
    while True:
        if i>0 and j>0:
            if bk[i][j] ==3:

                AlignV=v[i-1]+AlignV
                AlignW=w[j-1]+AlignW
                i=i-1
                j=j-1
            elif bk[i][j] ==1:

                AlignW='-'+AlignW
                AlignV=v[i-1]+AlignV
                i =i-1
            else:

                AlignW=w[j-1]+AlignW
                AlignV='-'+AlignV
                j =j-1
        elif i>0:
            AlignW='-'+AlignW
            AlignV=v[i-1]+AlignV
            i =i-1
        elif j>0:
            AlignW=w[j-1]+AlignW
            AlignV='-'+AlignV
            j =j-1
        else:
            break
    mm=[]

    for i in range(len(AlignV)):
        if not AlignV[i]==AlignW[i]:

            mm.append({'pos':i+1,'ref':AlignV[i],'alt':AlignW[i]})
    return mm,AlignV,AlignW
def calVariants(ref_fa,query_fa,output):
    '''Index reference genome, call aligment with BWA and convert to vcf'''
    temp_folder=output+"/temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    #alignment
    cmd = 'bwa mem {ref} {query} > {output}/aln-se.sam'.format(
        ref=ref_fa,
        query=query_fa,
        output=output
    )
    os.system(cmd)

    cmd = 'samtools view -S -b {output}/aln-se.sam > {temp}/aln-se.bam'.format(
        output=output,
        temp=temp_folder
    )
    os.system(cmd)
    #sort BAM for SNP calling
    #cmd = 'samtools sort {output}/aln-filtered.bam > {output}/aln-filtered-sorted.bam'.format(
    cmd = 'samtools sort {temp}/aln-se.bam > {temp}/aln-sorted.bam'.format(

        temp=temp_folder
    )
    os.system(cmd)

    #generate raw bcf
    #cmd = 'samtools mpileup -g -f {ref} {output}/aln-filtered-sorted.bam | bcftools call -mv -Ov > {output}/variants.vcf'.format(
    cmd = 'samtools mpileup -g  -f {ref} {temp}/aln-sorted.bam | bcftools call -mv -Ov > {output}/variants.vcf'.format(
        temp=temp_folder,
        ref=ref_fa,
        output=output
    )
    os.system(cmd)


    try:
        shutil.rmtree(temp_folder)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))

    return  output+'/variants.vcf', output+'/aln-se.sam'
def translateDNAMuToProMu(map_prot, pos, ref, aln):
    str_prot_mul=''
    p_start=pos
    order=0
    note=''
    if pos%3==0:
        p_start=p_start-2
        order=2
    elif pos%3==2:
        p_start=p_start-1
        order=1
    #print(map_prot)
    if len(ref)==1 and len(aln)==1:
        ref_p=map_prot[p_start]['p']
        codon=map_prot[p_start]['codon']
        list1 = list(codon)
        list1[order]=aln
        codon=''.join(list1)
        aln_p=protein[codon]

        if not ref_p==aln_p:
            str_prot_mul=ref_p+str(map_prot[p_start]['pos'])+aln_p+','
    else :

        num_diff=abs(len(ref)-len(aln))
        if not  num_diff%3==0:
            str_prot_mul=str(map_prot[p_start]['pos'])+'frameshift'
        else:
            #in-frame
            #number of codon effected
            diff=pos-p_start
            num_b=diff+len(ref)
            num_c=1
            if num_b%3==0:
                num_c=num_b/3
            else:
                num_c=math.floor(num_b/3)+1
            #print('numc:'+str(num_c))
            str_codons=''
            for c in range(int(num_c)):
                str_codons+=map_prot[p_start+c*3]['codon']
            #split and build new str codon
            new_str_codons=str_codons[:diff]+aln  +    str_codons[diff+len(ref):]
            old_p=translateToProt(str_codons)
            new_p=translateToProt(new_str_codons)

            mm,aop,anp=GlobalAlignment(old_p,new_p,1,1)

            p=0
            map_pos_to_pos_align={}
            for i in range(len(aop)):
                if not aop[i]=='-':
                    map_pos_to_pos_align[i]=map_prot[p_start+p*3]
                    p=p+1
            #print(map_pos_to_pos_align)
            str_prot_mul=''
            if len(mm)>0:
                for m in mm:
                    if m['alt']=='-':
                        str_prot_mul=str_prot_mul+m['ref']+str(map_pos_to_pos_align[m['pos']-1]['pos'])+m['alt']+','
                    elif m['ref']=='-':
                        cursor=m['pos']-1
                        lastp=m['pos']-1
                        while cursor>=0:
                            if not aop[cursor]== '-':
                                lastp=cursor
                                break
                            cursor=cursor-1
                        left_p=None
                        if lastp in map_pos_to_pos_align.keys():
                            left_p=map_pos_to_pos_align[lastp]
                        else:
                            left_p=map_prot[p_start-3]
                        cursor=m['pos']-1
                        firstp=m['pos']-1

                        while cursor<len(aop):
                            if not aop[cursor]== '-':
                                firstp=cursor
                                break
                            cursor=cursor+1
                        right_p=None
                        if firstp in map_pos_to_pos_align.keys():
                            right_p=map_pos_to_pos_align[firstp]
                        else:
                            #print(left_p)
                            right_p=map_prot[left_p['pos']*3-2+3]
                        str_prot_mul=str_prot_mul+left_p['p']+str(left_p['pos'])+'_'+right_p['p']+str(right_p['pos'])+'ins'+m['alt']+','
                    else:
                        str_prot_mul=str_prot_mul+m['ref']+str(map_pos_to_pos_align[m['pos']-1]['pos'])+m['alt']+','




    return str_prot_mul[:-1]
def getDNAMuString( pos, ref, aln):
    if len(ref)==len(aln):
        return str(ref+str(pos)+aln)
    elif len(ref)>len(aln):
        num_delete=len(ref)-len(aln)
        if num_delete==1:
            return str(pos+len(aln))+'del'+ref[len(aln):]
        else:
            return str(pos+len(aln))+'_'+str(pos+len(aln)+num_delete)+'del'
    else:
        num_add=len(aln)-len(ref)
        base_add=aln[len(ref):]
        return str(pos+len(ref)-1)+'_'+str(pos+len(ref)-1+num_add)+'ins'+base_add
def readVCF(vcffile):
    arrmul=[]
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    for record in vcf_reader:
        arrmul.append({'pos':record.POS,'ref':record.REF,'alts':record.ALT})
    return arrmul
def readSAM(samfile):
    hit=[]
    with open(samfile) as f:
        lines=f.readlines()
        for i in range(len(lines)):
            if not lines[i].startswith('@'):
                token=lines[i].split('\t')
                if token[1]=='4':
                    continue
                hit.append({'chr':token[0],'ss':token[3],'cigar':token[5],'seq':token[9]})


    return hit
def getGeneLocation(hit,gene):
    #parse CIGAR string
    loc={}
    loc['note']=''
    list_hit=[]
    total_len=0
    for h in hit:
        hit={}
        cigar=Cigar(h['cigar'])
        items=list(cigar.items())
        if items[0][1]=='S'and items[-1][1]=='S':
            hit['seq']=h['seq'][items[0][0]:-items[-1][0]]
        elif items[0][1]=='H'and items[-1][1]=='H':
            hit['seq']=h['seq']
        else:
            hit['seq']=h['seq']
        send=int(h['ss'])
        for item in items:
            if item[1]=='M':
                send=send+int(item[0])
            if item[1]=='D':
                send=send+int(item[0])
            if item[1]=='I':
                send=send-int(item[0])
        hit['pos']= int(items[0][0]  )
        hit['ss']=int(h['ss'])
        hit['send']=send
        list_hit.append(hit)
    list_hit_sorted=sorted(list_hit, key=lambda k: k['ss'])
    if not len(list_hit_sorted)>0:
        loc['consensus']=''
        loc['hit']=[]
        loc['note']='Not found'
        return loc
    scafold=[list_hit_sorted[0]]
    cover_len=0
    for i in range(len(list_hit_sorted)):
        if list_hit_sorted[i]['ss']>scafold[-1]['ss']+len(scafold[-1]['seq']):
            scafold.append(list_hit_sorted[i])

    if len(scafold)<len(list_hit_sorted):
        loc['note']="Multiple sequences found"
    loc['hit']=list_hit_sorted
    loc['pos']=list_hit_sorted[0]['pos']
    for i in range (len(scafold)):
        cover_len=cover_len+len(scafold[i]['seq'])
    loc['consensus']=makeConsensus(scafold,gene)


    loc['coverage']=cover_len/len(gene)
    #make consensus sequence from scafold

    #print(loc)
    return loc
def makeConsensus(scafold,gene):
    list_remain=[]
    first_seq=''
    if scafold[0]['ss']>1:
        first_seq=gene[:scafold[0]['ss']]

    for i in range (len(scafold)-1):
        list_remain.append(gene[scafold[i]['send']-1:scafold[i+1]['ss']-1])

    list_remain.append(gene[scafold[-1]['send']:])
    consensus=''+first_seq
    for i in range (len(scafold)):
        consensus=consensus+scafold[i]['seq']+list_remain[i]


    return consensus
def makePhylo(sequence_type,output,gene_ref):
    if not os.path.exists(output+'/phylo'):
        os.makedirs(output+'/phylo')
    #print(sequence_type)
    f = open(output+'/sequencetypes.fasta', "w")
    for k in sequence_type.keys():

        f.write(">"+sequence_type[k]['type']+'\n')
        f.write(str(sequence_type[k]['seq']+'\n'))
    with open(gene_ref) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            f.write(">"+record.id+'\n')
            f.write(str(record.seq)+'\n')
    f.close()
    cmd = './muscle  -in {}  -out {}'.format(
       output+'/sequencetypes.fasta',output+'/sequencetypes.afa'
        )
    print(cmd)
    os.system(cmd)
    cmd='FastTree {} > {} '.format(output+'/sequencetypes.afa',output+'/sequencetypes.tree')
    os.system(cmd)
    return output+'/sequencetypes.tree'
def doStatistic(report):
    set_prot_mu={}
    set_dna_mu={}
    set_type={}
    for  r in report:
        #split to get mutations
        prots=r['prot'].split(',')
        for p in prots:
            pstr=p.strip()
            if p=='':
                continue
            if not pstr in set_prot_mu:
                set_prot_mu[pstr]={}
                set_prot_mu[pstr]['count']=1
            else:
                set_prot_mu[pstr]['count']=set_prot_mu[pstr]['count']+1
        dnas=r['dna'].split(',')
        for d in dnas:
            dstr=d.strip()
            if dstr =='':
                continue
            if not dstr in set_dna_mu:
                set_dna_mu[dstr]={}
                set_dna_mu[dstr]['count']=1
            else:
                set_dna_mu[dstr]['count']=set_dna_mu[dstr]['count']+1
        # sequence Type
        if not r['seq'] in set_type:
            set_type[r['seq']]={}
            set_type[r['seq']]['count']=1
            set_type[r['seq']]['prot']=r['prot']
            set_type[r['seq']]['dna']=r['dna']
        else:
            set_type[r['seq']]['count']=set_type[r['seq']]['count']+1
    for k in set_prot_mu:
        set_prot_mu[k]['f']=round(set_prot_mu[k]['count']/len(report),3)
    for k in set_dna_mu:
        set_dna_mu[k]['f']=round(set_dna_mu[k]['count']/len(report),3)
    for k in set_type:
        set_type[k]['f']=round(set_type[k]['count']/len(report),3)
    statistics={'prot':set_prot_mu,'dna':set_dna_mu,'type':set_type}
    #print(statistics)
    return statistics

def assignData(line,hit_data,dna_data,prot_data,type_data,phylo_data,seq_name,seq_in,sample_in,dna_len,prot_len,mode):
    if 'f09328d495fa622700aabe5707edf00b' in line:
        line=line.replace('f09328d495fa622700aabe5707edf00b',hit_data)
    if '5427e0bf847d7c379a1d2e09a0b45f8a' in line:
        line=line.replace('5427e0bf847d7c379a1d2e09a0b45f8a',dna_data)
    if '02a42284998e2ce6a3776ad3696dceed' in line:
        line=line.replace('02a42284998e2ce6a3776ad3696dceed',prot_data)
    if '599dcce2998a6b40b1e38e8c6006cb0a' in line:
        line=line.replace('599dcce2998a6b40b1e38e8c6006cb0a',type_data)
    if '55f4634c7ac2899801eac4ced5b8120e' in line:
        line=line.replace('55f4634c7ac2899801eac4ced5b8120e',phylo_data)
    if 'd5d3db1765287eef77d7927cc956f50a' in line:
        line=line.replace('d5d3db1765287eef77d7927cc956f50a',seq_in)
    if 'a43c1b0aa53a0c908810c06ab1ff3967' in line:
        line=line.replace('a43c1b0aa53a0c908810c06ab1ff3967',sample_in)

    now = datetime.now()

    current_time = now.strftime("%H:%M:%S %d, %b, %Y")

    if '07cc694b9b3fc636710fa08b6922c42b' in line:
        line=line.replace('07cc694b9b3fc636710fa08b6922c42b',current_time)
    if 'b068931cc450442b63f5b3d276ea4297' in line:
        line=line.replace('b068931cc450442b63f5b3d276ea4297',seq_name)
    if '12340bf847d7c379a1d2e09a0b2345' in line:
        line=line.replace('12340bf847d7c379a1d2e09a0b2345',str(dna_len))
    if 'abcdbf847d7c379a1d2e09a0b2jkl' in line:
        line=line.replace('abcdbf847d7c379a1d2e09a0b2jkl',str(prot_len))
    if 't6678gghjuyt334567abe5707ediop' in line:
        line=line.replace('t6678gghjuyt334567abe5707ediop',mode)
    return line
def exportReport(report,output,phylo,gene_ref,samles,gene_name,dna_len,prot_len,mode):
    #hit table
    file_out=os.path.join(output,'hit-table.tsv')
    f=open(file_out,'w')
    f.write('Sequence\tPosition\tHits\tCoverage\tProtein mutations\tDNA mutations\tSequence type\tNote\n')

    #statistic_haplotype={}
    statistic_sample={}
    str_hit_data=''

    for r in report:
        print('coverage:'+str(r['coverage']))
        coverage=str(round(r['coverage']*100, 5))+'%'
        f.write(r['sid']+'\t'+str(r['pos'])+'\t'+str(r['hit'])+'\t'+coverage+'\t'+r['prot']+'\t'+r['dna']+'\t'+r['seq']+'\t'+r['note']+'\n')
        str_hit_data=str_hit_data+'["'+r['sid']+'","'+str(r['pos'])+'","'+str(r['hit'])+'","'+coverage+'","'+r['prot'].replace(',',', ')+'","'+r['dna'].replace(',',', ')+'","'+r['seq']+'","'+r['note']+'"],'
    f.close()
    #prepare for html export
    statistics=doStatistic(report)

    str_hit_data='['+str_hit_data[:-1]+']'
    str_dna_data=''
    for p in statistics['dna'].keys():
        str_dna_data=str_dna_data+'["'+p+'","'+str(statistics['dna'][p]['count'])+'","'+str(statistics['dna'][p]['f'])+'"],'
    str_dna_data='['+str_dna_data[:-1]+']'
    str_prot_data=''
    for p in statistics['prot'].keys():
        str_prot_data=str_prot_data+'["'+p+'","'+str(statistics['prot'][p]['count'])+'","'+str(statistics['prot'][p]['f'])+'"],'
    str_prot_data='['+str_prot_data[:-1]+']'
    str_type_data=''
    for t in statistics['type'].keys():
        str_type_data=str_type_data+'["'+t+'","'+str(statistics['type'][t]['count'])+'","'+str(statistics['type'][t]['f'])+'","'+str(statistics['type'][t]['dna'].replace(',',', '))+'","'+str(statistics['type'][t]['prot'].replace(',',', '))+'"],'
    str_type_data='['+str_type_data[:-1]+']'

    phylotree=''
    with open(phylo) as f:
        phylotree = f.read()
    #print(phylotree)
    #copyfile('report_template/report.html', output+'/report.html')
    f_in = open('report_template/report.html', "r")
    f_out = open( output+'/report_'+str(datetime.now().timestamp())+'.html', "w")
    for line_in in f_in:
        f_out.write(assignData(line_in,str_hit_data,str_dna_data,str_prot_data,str_type_data,phylotree,gene_name,gene_ref,samles,dna_len,prot_len,mode))
    f_in.close()
    f_out.close()

def pipeline(args):
    gene_ref=args.seq
    sample_seqs=args.genomes
    output=args.output
    mode='seq'
    if args.gene:
        mode='gene'
    str_dna=''
    gene_name=''
    len_prot=0
    with open(gene_ref) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            gene_name=record.id
            str_dna=record.seq
    m_prot=None
    if mode=='gene':
        if not len(str_dna)%3==0:
            print('Gene length is not approriate!')
            return
        m_prot=translateToProtMap(str_dna)
        len_prot=len(m_prot.keys())
    else:
        m_prot=None
       #indexing reference
    isIndex=False
    if not isIndex:
        cmd='bwa index '+gene_ref
        os.system(cmd)
        cmd='samtools faidx '+gene_ref
        os.system(cmd)
    gene_profiles={}
    sequence_type={}
    report=[]
    stype=1
    if not os.path.exists(output):
        os.makedirs(output)
    with open(sample_seqs) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_f=record.id.replace('/','_')
            if not os.path.exists(os.path.join(output,record_f)):
                os.makedirs(os.path.join(output,record_f))
            with open(os.path.join(output,record_f)+'/'+record_f+".fasta", "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")

            vcf,sam=calVariants(gene_ref,os.path.join(output,record_f)+'/'+record_f+".fasta",os.path.join(output,record_f))
            mutations=readVCF(vcf)
            str_mut=''
            str_dna_mut=''
            for mut in mutations:
                for alt in mut['alts']:
                    if mode=='gene':
                        str_mut=str_mut+translateDNAMuToProMu(m_prot,mut['pos'],mut['ref'],str(alt))+','
                    str_dna_mut=str_dna_mut+getDNAMuString(mut['pos'],mut['ref'],str(alt))+','

            sam=readSAM(sam)
            loc=getGeneLocation(sam,str_dna)

            str_md5 = hashlib.md5(str(loc['consensus']).encode()).hexdigest()
            if len(loc['consensus'])>0:
                if not str_md5 in sequence_type.keys():
                    sequence_type[str_md5]={}
                    sequence_type[str_md5]['type']='Type'+str(stype)
                    sequence_type[str_md5]['seq']=loc['consensus']
                    stype=stype+1

                with open(os.path.join(output,record_f)+"/consensus.fasta", "w") as output_handle:

                    record.seq=loc['consensus']
                    SeqIO.write(record, output_handle, "fasta")

                report.append({'sid':record.id,'prot':str_mut[:-1],'dna':str_dna_mut[:-1],'pos':loc['pos'],'hit':len(loc['hit']),'coverage':loc['coverage'],'seq':sequence_type[str_md5]['type'],'note':loc['note']})
    phylofile=makePhylo(sequence_type,output,gene_ref)
    #print(report)
    exportReport(report,output,phylofile,gene_ref,sample_seqs,gene_name,len(str_dna),len_prot,mode)
    #print(len(sequence_type))

count=0
'''with open("Delta.fasta") as handle:
    with open("Delta1000.fasta", "w") as output_handle:
        for record in SeqIO.parse(handle, "fasta"):
            if count<1000:
                SeqIO.write(record, output_handle, "fasta")
            count=count+1'''

#pipeline('/media/ktht/Store/Quang/bio/check-gene/spike.fna','/media/ktht/Store/Quang/bio/check-gene/Delta1000.fasta','/media/ktht/Store/Quang/bio/check-gene/data4',mode='gene')


def main(arguments=sys.argv[1:]):

    parser = argparse.ArgumentParser(
        prog='Target sequence analysis',
        description='Tool for checking sequence in a genome corpus')

    parser.set_defaults(func=pipeline)
    parser.add_argument('-s','--seq', help='Target sequence file (Fasta format)',type=str)
    parser.add_argument('-g','--genomes', help='Collection of genomes',type=str)
    parser.add_argument('-o','--output', help='Output folder', type=str)
    parser.add_argument('--gene', help='Indicate that sequence is a gene ', action="store_true")

    args = parser.parse_args()
    start = time.time()
    pipeline(args)
    end = time.time()
    print(end - start)

if __name__ == "__main__":
    main()
