import pandas as pd

def join_dataframes(STR_df1, STR_df2):
    #print("shape of", STR_df1["FastaHeader"], ":", STR_df1.shape)
    #print("shape of", STR_df2["FastaHeader"], ":", STR_df2.shape)
    STR_merged = pd.merge(STR_df1, STR_df2, how="outer")
    #print("shape of STR_merged: ", STR_merged.shape)
    return STR_merged

def CreateFilesArray():
    STRFiles = []
    for i in range (1, 23):
        STRPath = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanChr" + str(i) + ".txt"
        STRFiles.append(STRPath)
    STRFiles.append("/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanChr" + "X" + ".txt")
    STRFiles.append("/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanChr" + "Y" + ".txt")
    return STRFiles

def merge_all(files_list):
    STR_total = pd.read_csv(files_list[0], sep="\t")
    for i in files_list[1:]:
        STR_df = pd.read_csv(i, sep="\t")
        STR_total = join_dataframes(STR_total, STR_df)
        #print(STR_total)
    return STR_total


def writeNewtxt(file_path, df):
    """Writes a pandas dataframe into txt, columns separated by tabs ("\t")
    """
    f = open(file_path, 'w')
    str_cols = ""
    for col in df.columns[:-1]:
        str_cols += col + "\t"
    str_cols += df.columns.values[-1]+ "\n"
    f.write(str(str_cols))
    for row in df.values:
        str_row = ""
        for col in row[:-1]:
            str_row += str(col) + "\t"
        str_row += str(row[-1])+ "\n"
        f.write(str_row)
    f.close


#STRPath_read1 = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanChr.1.txt"
#STRPath_read2 = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanChr.2.txt"
STRwrite = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanGenome.txt"

STRFiles = CreateFilesArray()
STR_all = merge_all(STRFiles)
writeNewtxt(STRwrite, STR_all)
print(STR_all.shape)            #total nº of repeats should be 1.199.628, current is 99?.???, so its missing about 200.000 rows!
print(STR_all.head)
print(STR_all["FastaHeader"][30000])
#merged_df = join_dataframes(STRPath_read1, STRPath_read2)
#print(merged_df.shape)
#print(merged_df.head)
#print(merged_df["FastaHeader"][70000])

