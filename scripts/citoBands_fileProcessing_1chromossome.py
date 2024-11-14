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
    chomossome = 'chr' + str(n)
    chr_array = []
    #chr_array.append(['Chromossome','First_index', 'Last_index', "Citoband", 'Unkown'])
    for row in df.values:
        cito_array = []
        if row[0] == chomossome:
            for col in row:
                cito_array.append(col)
            chr_array.append(cito_array)
    return pd.DataFrame(chr_array, columns=['Chromossome','First_index', 'Last_index', 'Citoband', 'Unkown', 'Size'] )


######################################## calc densities ############################################
def Somatorio_STR_lengths(citoband, STR_df):
    """
    Requires: STR_df is pandas dataframe with the according citoband column
    """
    totalSTRlength = 0
    for row in range(len(STR_df.values)):
        if STR_df["Citoband"][row] == citoband[3]:# or citoband[3] in STR_df["Citoband"][row]:
            totalSTRlength += int(STR_df["ArrayLength"][row])
        elif citoband[3] in STR_df["Citoband"][row]:
            if citoband[3] + " &" in STR_df["Citoband"][row]:
                print("primeira citoband do par")
                totalSTRlength += int(citoband[2]) - int(STR_df["FirstIndex"][row])
            elif "& " + citoband[3] in STR_df["Citoband"][row]:
                print("segunda citoband do par")
                totalSTRlength += int(STR_df["LastIndex"][row]) - int(citoband[1])

    return totalSTRlength

#Esta função só considera STRs em linhas seguidas
def STR_overlap(citoband, STR_df):
    to_subtract = 0
    for row in range(len(STR_df)-1):
        if str(citoband) in STR_df["Citoband"][row] and str(citoband) in STR_df["Citoband"][row+1]:
            last_index = STR_df["LastIndex"][row]
            next_first_index = STR_df["FirstIndex"][row + 1]
            if next_first_index<last_index:
                print("overlap found")
                to_subtract += last_index - next_first_index
    return to_subtract


def STR_dens_toCitoband(citobands_df, STR_df):
    """
    """
    densities = []
    for citoband in citobands_df.values:
        density = (Somatorio_STR_lengths(citoband, STR_df) - STR_overlap(citoband, STR_df)) /int(citoband[-1])
        densities.append(density)
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

def makeaplot(x, y):
    plt.figure(figsize=(10,10))
    plt.plot(x, y)
    plt.xlabel("Citoband")
    plt.ylabel("STR_density")
    plt.grid()
    plt.show()

def calc_strDensity(citobands_df, STR_df):
    densities = []
    citobands_size = citobands_df["Last_index"][-1:]+1
    citobands_strs = np.zeros((citobands_size), dtype=int)
    #print(len(citobands_strs))
    for row in range(len(STR_df.values)):
        for i in range (STR_df["FirstIndex"][row], STR_df["LastIndex"][row]):
            citobands_strs[i]+=1
    print("total STR length in chromossome ", chr, ": ", np.count_nonzero(citobands_strs))
    for citoband in range(len(citobands_df.values)):
        citoband_array = citobands_strs[int(citobands_df["First_index"][citoband]):int(citobands_df["Last_index"][citoband])]
        #print(citoband_array)
        str_length = np.count_nonzero(citoband_array)
        citoband_length = citobands_df["Size"][citoband]
        density = str_length/citoband_length
        densities.append(density)
        
        #print(np.where(citoband_strs > 0))
    citobands_df.insert(len(citobands_df.columns), "STR_Density", densities, True)




citoPath_read = "/home/androx/Documents/trabalho/citobands/cytobandFiltered.txt"
citoPath_write = "/home/androx/Documents/trabalho/citobands/cytobandFiltered_Processed.txt"
STRPath_read = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_nov2014/HumanChr.18_processed.txt"

citobands_df = add_headersToDf(pd.read_csv(citoPath_read, sep="\t", header=None))
str_df = pd.read_csv(STRPath_read, sep="\t")
#print("columns: ", cito_dataframe.columns)

subtract_lastindex(citobands_df)
add_sizeColumn(citobands_df)
#print(citobands_df)

citobands_df_18 = getChromossome(citobands_df, 18)

calc_strDensity(citobands_df_18, str_df)
print(citobands_df_18)

makeaplot(citobands_df_18["Citoband"], citobands_df_18["STR_Density"])