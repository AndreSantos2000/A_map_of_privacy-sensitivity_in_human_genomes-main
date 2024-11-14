import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def add_headersToDf(dataframe):
    return dataframe.rename(columns={0: 'Chromossome', 1: 'First_index', 2: 'Last_index', 3: 'Citoband', 4: 'Unkown'})
    

def subtract_lastindex(dataframe):
    new_row = []
    for row in dataframe['Last_index']:
        new_row.append(row - 1)
    
    dataframe.drop('Last_index', axis = 1, inplace=True)

    dataframe.insert(2, "Last_index", new_row, True)


def add_sizeColumn(dataframe):
    sizes = []
    for row in range(len(dataframe.values)):
        size = dataframe["Last_index"][row] - dataframe["First_index"][row]
        sizes.append(size)
    dataframe.insert(len(dataframe.columns), "Size", sizes, True)


def getChromossome(df, n):
    """Returns pandas dataframe with only the citobands of one chromossome n
    """
    chomossome = str(n)
    chr_array = []
    #print("dataframe columns = ", df.columns)
    #chr_array.append(['Chromossome','First_index', 'Last_index', "Citoband", 'Unkown'])
    for row in df.values:
        cito_array = []
        if row[0] == chomossome:
            for col in row:
                cito_array.append(col)
            chr_array.append(cito_array)
    return pd.DataFrame(chr_array, columns=df.columns)#['Chromossome','First_index', 'Last_index', 'Citoband', 'Unkown', 'Size'] )


######################################## calc densities ############################################
def calc_strDensity(citobands_df, STR_df):
    densities = []
    #print(citobands_df["Chromossome"].unique())
    #print(citobands_df.columns)
    for chr in citobands_df["Chromossome"].unique():                #ordem esta definida aqui
        print("searching citobands in chromossome: ", chr)
        citobands_chr = getChromossome(citobands_df, chr)
        #print(citobands_chr)
        citobands_size = citobands_chr["Last_index"][-1:]+1
        citobands_strs = np.zeros((citobands_size), dtype=int)
        #print(len(citobands_strs))
        for row in range(len(STR_df.values)):
            if STR_df["FastaHeader"][row] == chr:
                for i in range (STR_df["FirstIndex"][row], STR_df["LastIndex"][row]):
                    citobands_strs[i]+=1
                    #print("inside cicle: ", np.count_nonzero(citobands_strs))
        #print("total STR length in chromossome ", chr, ": ", np.count_nonzero(citobands_strs))
        for citoband in range(len(citobands_chr.values)):
            citoband_array = citobands_strs[int(citobands_chr["First_index"][citoband]):int(citobands_chr["Last_index"][citoband])]
            #print(citoband_array)
            str_length = np.count_nonzero(citoband_array)
            citoband_length = citobands_chr["Size"][citoband]
            density = str_length/citoband_length
            #print("length of citoband ", citobands_df["Citoband"][citoband] ,": ", citoband_length)
            #print("length of strs in citoband ", citobands_df["Citoband"][citoband] ,": ", str_length)
            #print("resulting str density in citoband ", citobands_df["Citoband"][citoband] ,": ", density)
            densities.append(round(density, 3))
        
        #print(np.where(citoband_strs > 0))
    citobands_df.insert(len(citobands_df.columns), "STR_Density", densities, True)



def writetxt(txt_path, df):
    with open(txt_path, 'w') as f:
        df_string = df.to_string(header=False, index=False)
        f.write(df_string)
        f.close

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

citoPath_read = "/home/androx/Documents/trabalho/citobands/cytobandFiltered.txt"
citoPath_write = "/home/androx/Documents/trabalho/citobands/cytobandFiltered_processed_complete.txt"
#STRPath_read = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_nov2014/HomosapiensHG38_FullGenome_processed.txt"
STRPath_read = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanGenome_processed.txt"

citobands_df = add_headersToDf(pd.read_csv(citoPath_read, sep="\t", header=None))
str_df = pd.read_csv(STRPath_read, sep="\t")
#print("columns: ", cito_dataframe.columns)

subtract_lastindex(citobands_df)
add_sizeColumn(citobands_df)
#print(citobands_df)

calc_strDensity(citobands_df, str_df)
print(citobands_df)
#print(str_df)

#writetxt(citoPath_write, citobands_df)
writeNewtxt(citoPath_write, citobands_df)