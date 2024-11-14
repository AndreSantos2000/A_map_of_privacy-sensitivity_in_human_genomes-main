import argparse 
import math
import pandas as pd

filePath = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_nov2014/HumanChr.22.txt"
RefFilePath = "/home/androx/Documents/trabalho/datasets/GenomeAssembly/feb3_2022/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna'"

currentline=0
f = open(filePath, 'r')


str_dataframe = pd.read_csv(filePath, sep="\t")
#print(str_dataframe)
#print(str_dataframe.columns)

#Função que devolve a janela da dimensao desejada à volta do padrão
def getSTRSequence(dataframe, line, size):
    """
    Indices may be specific to single STR file
    """
    
    STRsize = int(dataframe["ArrayLength"][line])
    print(STRsize)
    
    if size < STRsize or size > STRsize + 200:
        print("need a sequence of size at least:", STRsize+2, "and no bigger than: ", STRsize + 200)
    elif size == STRsize:
        sequence = dataframe[line]['ArraySequence']
        lc_sequence =sequence.lower()
        print("returning the STR sequence with no flanks")
        return lc_sequence, len(lc_sequence)
    else:
        flanks_size = size - STRsize
        left = dataframe['FlankingLeft100'][line][-round(flanks_size/2):]
        short_tr = dataframe['ArraySequence'][line]
        rigth = dataframe['FlankingRight100'][line][:math.floor(flanks_size/2)]
        sequence = left + short_tr + rigth
        
        lc_sequence =sequence.lower()
        return lc_sequence, len(lc_sequence)
    

print(getSTRSequence(str_dataframe, 10, 180))
    
def main():

    parser = argparse.ArgumentParser(description='Script to receive gene string sequence from STR txt input file')

    parser.add_argument('-p','--refpath', required=True, type=str, help='Path to the STR file')
    parser.add_argument('-r','--row', required=True, type=str, help='number of the row to obtain the sequence')
    parser.add_argument('-s', '--Sequence_size', required=True, type=str, help='the size of returned sequence')
    parser.add_argument('-c', '--case', required=False, type=str, help='upper or lowercase')

    args = parser.parse_args()

    filePath = args.refpath

    dataset = pd.read_csv(filePath, sep="\t")

    seq = getSTRSequence(dataset, int(args.row), int(args.Sequence_size))

    print(seq)

if __name__ == '__main__':
    main()

#print(str_dataframe["ArrayLength"][10])
#print(str_dataframe["FlankingLeft50"][10])
#print(str_dataframe["FlankingLeft50"][10][:6])
#print(str_dataframe["FlankingLeft50"][10][-6:])