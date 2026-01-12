#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def bed_SNP_stat(filename):
    df = pd.read_csv(f"{filename}.bed", sep='\t', names=["Chr", "SNP_start", "SNP_end", "SNP_info", "Chr2", "Alu_start", "Alu_end"])
    df[['Type', 'numbers']]=df['SNP_info'].str.split('-',expand=True).astype(str)
    df[['A1', 'A2', 'B1','B2']]=df['numbers'].str.split('_',expand=True).astype(int)
    AtoG_genome1 = df.loc[df['Type'] == 'A_G', 'A2'].sum()
    AtoG_genome2 = df.loc[df['Type'] == 'T_C', 'A1'].sum()
    AtoG_genome = AtoG_genome1 + AtoG_genome2

    AtoC_genome1 = df.loc[df['Type'] == 'A_C', 'A2'].sum()
    AtoC_genome2 = df.loc[df['Type'] == 'T_G', 'A1'].sum()
    AtoC_genome = AtoC_genome1 + AtoC_genome2

    AtoT_genome1 = df.loc[df['Type'] == 'A_T', 'A2'].sum()
    AtoT_genome2 = df.loc[df['Type'] == 'T_A', 'A1'].sum()
    AtoT_genome = AtoT_genome1 + AtoT_genome2

    TtoG_genome1 = df.loc[df['Type'] == 'T_G', 'A2'].sum()
    TtoG_genome2 = df.loc[df['Type'] == 'A_C', 'A1'].sum()
    TtoG_genome = TtoG_genome1 + TtoG_genome2

    TtoC_genome1 = df.loc[df['Type'] == 'T_C', 'A2'].sum()
    TtoC_genome2 = df.loc[df['Type'] == 'A_G', 'A1'].sum()
    TtoC_genome = TtoC_genome1 + TtoC_genome2

    TtoA_genome1 = df.loc[df['Type'] == 'T_A', 'A2'].sum()
    TtoA_genome2 = df.loc[df['Type'] == 'A_T', 'A1'].sum()
    TtoA_genome = TtoA_genome1 + TtoA_genome2

    CtoG_genome1 = df.loc[df['Type'] == 'C_G', 'A2'].sum()
    CtoG_genome2 = df.loc[df['Type'] == 'G_C', 'A1'].sum()
    CtoG_genome = CtoG_genome1 + CtoG_genome2

    CtoT_genome1 = df.loc[df['Type'] == 'C_T', 'A2'].sum()
    CtoT_genome2 = df.loc[df['Type'] == 'G_A', 'A1'].sum()
    CtoT_genome = CtoT_genome1 + CtoT_genome2

    CtoA_genome1 = df.loc[df['Type'] == 'C_A', 'A2'].sum()
    CtoA_genome2 = df.loc[df['Type'] == 'G_T', 'A1'].sum()
    CtoA_genome = CtoA_genome1 + CtoA_genome2

    GtoA_genome1 = df.loc[df['Type'] == 'G_A', 'A2'].sum()
    GtoA_genome2 = df.loc[df['Type'] == 'C_T', 'A1'].sum()
    GtoA_genome = GtoA_genome1 + GtoA_genome2

    GtoC_genome1 = df.loc[df['Type'] == 'G_C', 'A2'].sum()
    GtoC_genome2 = df.loc[df['Type'] == 'C_G', 'A1'].sum()
    GtoC_genome = GtoC_genome1 + GtoC_genome2

    GtoT_genome1 = df.loc[df['Type'] == 'G_T', 'A2'].sum()
    GtoT_genome2 = df.loc[df['Type'] == 'C_A', 'A1'].sum()
    GtoT_genome = GtoT_genome1 + GtoT_genome2
    
    AtoG_SNP1 = df.loc[df['Type'] == 'A_G', 'B2'].sum()
    AtoG_SNP2 = df.loc[df['Type'] == 'T_C', 'B1'].sum()
    AtoG_SNP = AtoG_SNP1 + AtoG_SNP2

    AtoC_SNP1 = df.loc[df['Type'] == 'A_C', 'B2'].sum()
    AtoC_SNP2 = df.loc[df['Type'] == 'T_G', 'B1'].sum()
    AtoC_SNP = AtoC_SNP1 + AtoC_SNP2

    AtoT_SNP1 = df.loc[df['Type'] == 'A_T', 'B2'].sum()
    AtoT_SNP2 = df.loc[df['Type'] == 'T_A', 'B1'].sum()
    AtoT_SNP = AtoT_SNP1 + AtoT_SNP2

    TtoG_SNP1 = df.loc[df['Type'] == 'T_G', 'B2'].sum()
    TtoG_SNP2 = df.loc[df['Type'] == 'A_C', 'B1'].sum()
    TtoG_SNP = TtoG_SNP1 + TtoG_SNP2

    TtoC_SNP1 = df.loc[df['Type'] == 'T_C', 'B2'].sum()
    TtoC_SNP2 = df.loc[df['Type'] == 'A_G', 'B1'].sum()
    TtoC_SNP = TtoC_SNP1 + TtoC_SNP2

    TtoA_SNP1 = df.loc[df['Type'] == 'T_A', 'B2'].sum()
    TtoA_SNP2 = df.loc[df['Type'] == 'A_T', 'B1'].sum()
    TtoA_SNP = TtoA_SNP1 + TtoA_SNP2

    CtoG_SNP1 = df.loc[df['Type'] == 'C_G', 'B2'].sum()
    CtoG_SNP2 = df.loc[df['Type'] == 'G_C', 'B1'].sum()
    CtoG_SNP = CtoG_SNP1 + CtoG_SNP2

    CtoT_SNP1 = df.loc[df['Type'] == 'C_T', 'B2'].sum()
    CtoT_SNP2 = df.loc[df['Type'] == 'G_A', 'B1'].sum()
    CtoT_SNP = CtoT_SNP1 + CtoT_SNP2

    CtoA_SNP1 = df.loc[df['Type'] == 'C_A', 'B2'].sum()
    CtoA_SNP2 = df.loc[df['Type'] == 'G_T', 'B1'].sum()
    CtoA_SNP = CtoA_SNP1 + CtoA_SNP2

    GtoA_SNP1 = df.loc[df['Type'] == 'G_A', 'B2'].sum()
    GtoA_SNP2 = df.loc[df['Type'] == 'C_T', 'B1'].sum()
    GtoA_SNP = GtoA_SNP1 + GtoA_SNP2

    GtoC_SNP1 = df.loc[df['Type'] == 'G_C', 'B2'].sum()
    GtoC_SNP2 = df.loc[df['Type'] == 'C_G', 'B1'].sum()
    GtoC_SNP = GtoC_SNP1 + GtoC_SNP2

    GtoT_SNP1 = df.loc[df['Type'] == 'G_T', 'B2'].sum()
    GtoT_SNP2 = df.loc[df['Type'] == 'C_A', 'B1'].sum()
    GtoT_SNP = GtoT_SNP1 + GtoT_SNP2
       
    column_names = ['Genome', 'SNP']
    row_names = ['AtoG', 'AtoT', 'AtoC', 'TtoG', 'TtoC', 'TtoA', 'CtoG', 'CtoA', 'CtoT', 'GtoA', 'GtoC', 'GtoT']
    dfresult = pd.DataFrame(index=row_names, columns=column_names)
    
    dfresult.iloc[0, 0] = AtoG_genome
    dfresult.iloc[1, 0] = AtoT_genome
    dfresult.iloc[2, 0] = AtoC_genome
    dfresult.iloc[3, 0] = TtoG_genome
    dfresult.iloc[4, 0] = TtoC_genome
    dfresult.iloc[5, 0] = TtoA_genome
    dfresult.iloc[6, 0] = CtoG_genome
    dfresult.iloc[7, 0] = CtoA_genome
    dfresult.iloc[8, 0] = CtoT_genome
    dfresult.iloc[9, 0] = GtoA_genome
    dfresult.iloc[10, 0] = GtoC_genome
    dfresult.iloc[11, 0] = GtoT_genome
    
    dfresult.iloc[0, 1] = AtoG_SNP
    dfresult.iloc[1, 1] = AtoT_SNP
    dfresult.iloc[2, 1] = AtoC_SNP
    dfresult.iloc[3, 1] = TtoG_SNP
    dfresult.iloc[4, 1] = TtoC_SNP
    dfresult.iloc[5, 1] = TtoA_SNP
    dfresult.iloc[6, 1] = CtoG_SNP
    dfresult.iloc[7, 1] = CtoA_SNP
    dfresult.iloc[8, 1] = CtoT_SNP
    dfresult.iloc[9, 1] = GtoA_SNP
    dfresult.iloc[10, 1] = GtoC_SNP
    dfresult.iloc[11, 1] = GtoT_SNP
    
    dfresult.to_csv(f'{filename}_SNP_stats.csv')    

def main():
    parser = argparse.ArgumentParser(description='Bed file name after bedtools intersect with selected regions.')
    parser.add_argument('filename', type=str, help='bed file name')
    
    args = parser.parse_args()
    
    bed_SNP_stat(args.filename)
    print('Stats calculated.')

if __name__ == '__main__':
    main()
