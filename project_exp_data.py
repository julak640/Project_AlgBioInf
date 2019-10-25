#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 08:59:30 2019

@author: bioinf
"""
import pandas as pd
from functools import reduce
import mygene
from sklearn.preprocessing import quantile_transform 
import argparse

#os.chdir('/home/bioinf/Documents/ATAC-seq/170619/')

def parse_exp_data(data, samples, key_file, key_samples):
    key = pd.read_csv(key_file, sep='\t').loc[:,['ID','SPOT_ID','GENE_SYMBOL', 'GENE_NAME', 'UNIGENE_ID', 'ENSEMBL_ID']] #'data/RNA_expression/GSE78244/GPL17077-17467.txt'
    key = key.rename(columns={'ID':'ID REF'})
    sample_key = pd.read_csv(key_samples, sep='\t', header=None, index_col=0).T.to_dict('records')[0] #'data/RNA_expression/sample_key.txt'
    
    data_gene_symbol = pd.merge(data, key, how='left', left_on = ['ID REF'], right_on = ['ID REF'])
    data_gene_symbol = data_gene_symbol[['GENE_SYMBOL']+samples].dropna()
    data_gene_symbol = data_gene_symbol.rename(columns = sample_key)
    
    data_gene_symbol.drop(list(data_gene_symbol.filter(regex = 'GSM')), axis = 1, inplace = True) #for tysbari samples
    data_gene_symbol.drop(list(data_gene_symbol.filter(regex = 'MS')), axis = 1, inplace = True) #for tysbari samples

    data_gene_symbol = data_gene_symbol.groupby('GENE_SYMBOL').mean().reset_index() 
    return data_gene_symbol

def parse_cell_report_data(cr_file, key_file, key_samples):
    CR = pd.read_csv(cr_file, sep='\t')#'data/RNA_expression/GSE78244/Matrix med normaliserade data fran Mika.csv'
    CR_samples = ['MS3_unstim', 'MS3_stim',
           'MS4_unstim', 'MS4_stim', 'MS5_unstim', 'MS5_stim', 'MS6_unstim',
           'MS6_stim', 'MS7_unstim', 'MS7_stim', 'MS1_unstim', 'MS1_stim',
           'MS2_unstim', 'MS2_stim', 'MS14_unstim', 'MS14_stim', 'MS11_unstim', 'MS11_stim',
           'MS12_unstim', 'MS12_stim', 'MS13_unstim', 'MS13_stim', 'MS9_unstim',
           'MS9_stim', 'MS10_unstim', 'MS10_stim', 'MS8_unstim', 'MS8_stim']
    CR = parse_exp_data(CR, CR_samples, key_file, key_samples)
    return CR

def parse_tysabri_data(tys1_file, tys2_file, key_file, key_samples):
    tys1 = pd.read_csv(tys1_file, sep='\t') #'data/RNA_expression/GSE44964_Tysabri_HL/GSE44964-GPL14550_series_matrix_2.txt'
    tys1 = tys1.rename(columns={'ID_REF':'ID REF'})
    tys1_samples = ['GSM1094623', 'GSM1094624', 'GSM1094625', 'GSM1094626',
           'GSM1094627', 'GSM1094628', 'GSM1094629', 'GSM1094630', 'GSM1094631',
           'GSM1094632', 'GSM1094633', 'GSM1094634', 'GSM1094635', 'GSM1094636',
           'GSM1094637', 'GSM1094638', 'GSM1094639', 'GSM1094640', 'GSM1094641',
           'GSM1094642', 'GSM1094643', 'GSM1094644', 'GSM1094645', 'GSM1094646']
    tys1 = parse_exp_data(tys1, tys1_samples, key_file, key_samples)
    
    tys2 = pd.read_csv(tys2_file, sep='\t') #'data/RNA_expression/GSE44964_Tysabri_HL/GSE44964-GPL17077_series_matrix_2.txt'
    tys2 = tys2.rename(columns={'ID_REF':'ID REF'})
    tys2_samples = ['GSM1094647', 'GSM1094648', 'GSM1094649', 'GSM1094650', 'GSM1094651', 'GSM1094652', 'GSM1094653', 'GSM1094654']
    tys2 = parse_exp_data(tys2, tys2_samples, key_file, key_samples)
    return tys1, tys2

def gene_symbol_to_entrez(exp, atac_genes):
    mg = mygene.MyGeneInfo()
    entrez = mg.querymany(list(exp.loc[:,'GENE_SYMBOL']), scopes='symbol', fields='entrezgene', as_dataframe=True)
    entrez = entrez.loc[:,'entrezgene'].dropna()
    entrez1 = entrez.reset_index()
    entrez1.index = entrez1.entrezgene.astype(int)
    
    atac_genes = [gene for gene in atac_genes if gene in entrez1.index]
    entrez1 = entrez1.loc[atac_genes,:].dropna() #only genes present in atac-seq
    entrez1 = pd.Series(entrez1.entrezgene.values,index=entrez1.iloc[:,0].values).to_dict()
    
    exp = exp.replace({'GENE_SYMBOL':entrez1})
    exp = exp[~exp['GENE_SYMBOL'].str.contains("[a-zA-Z]").fillna(False)]
    exp.index = exp['GENE_SYMBOL']
    exp = exp.drop('GENE_SYMBOL',1)
    return exp

def normalize_exp(exp):
    exp_qt = quantile_transform(exp, axis=0, n_quantiles=1000, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, random_state=None, copy=False)
    exp = pd.DataFrame(exp_qt, columns=exp.columns, index=exp.index)
    # exp = exp.apply(zscore)
    
    return exp

def main():
    
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--cr_samples', type=str, help='file with samples from the cell report paper.')
    parser.add_argument('--tys1_samples', type=str, help='file with part 1 of the samples from the tysabri paper.')
    parser.add_argument('--tys2_samples', help='file with part 2 of the samples from the tysabri paper.')
    parser.add_argument('--key', help='key with gene annotation')
    parser.add_argument('--key_samples', help='key with sample name translation')
    parser.add_argument('--atac_genes', help='Annotated peaks with CHIPseeker')
    parser.add_argument('--output_file', help='file where results are saved')
    args = parser.parse_args()
    
    cr_samples="/home/bioinf/Documents/ATAC-seq/atac_pipeline/data/RNA_expression/GSE78244/Matrix_med_normaliserade_data_fran_Mika.csv"
    tys1_samples="/home/bioinf/Documents/ATAC-seq/atac_pipeline/data/RNA_expression/GSE44964_Tysabri_HL/GSE44964-GPL14550_series_matrix_2.txt"
    tys2_samples="/home/bioinf/Documents/ATAC-seq/atac_pipeline/data/RNA_expression/GSE44964_Tysabri_HL/GSE44964-GPL17077_series_matrix_2.txt"
    key="/home/bioinf/Documents/ATAC-seq/atac_pipeline/data/RNA_expression/GSE78244/GPL17077-17467.txt"
    key_samples="/home/bioinf/Documents/ATAC-seq/atac_pipeline/data/RNA_expression/sample_key.txt"
    #atac_genes="/home/bioinf/Documents/ATAC-seq/atac_pipeline/results/DAP_DESeq2_AvU_annotated_CHIPseeker.csv"
    #output_file="/home/bioinf/Documents/ATAC-seq/atac_pipeline/results/exp_atac_genes_normalized.csv"
       
    CR = parse_cell_report_data(cr_samples, key, key_samples) #CR = parse_cell_report_data(cr_samples, key, key_samples)

    tys1, tys2 = parse_tysabri_data(tys1_samples, tys2_samples, key, key_samples) #tys1, tys2 = parse_tysabri_data(tys1_samples, tys2_samples, key, key_samples)

    exp = reduce(lambda left,right: pd.merge(left,right,on='GENE_SYMBOL'),[CR, tys1, tys2])
    exp.index = exp.GENE_SYMBOL
    exp = exp.drop('GENE_SYMBOL',1)

    #atac_genes = pd.read_csv(args.atac_genes, sep='\t').geneId.unique().astype(int) #atac_genes = pd.read_csv(atac_genes, sep='\t').geneId.unique().astype(int) 
    
    #exp = gene_symbol_to_entrez(exp, atac_genes)
    
    exp_norm = normalize_exp(exp.T)
    exp_norm.to_csv('/home/bioinf/Documents/Courses/MedBioInfo_courses/Algorithms_bioinf/Algorithms_Bioinformatics/project_data/exp_normalized.csv', sep='\t', index=True)#'results/exp_atac_genes_normalized.csv'

if __name__ == '__main__':
    main()

