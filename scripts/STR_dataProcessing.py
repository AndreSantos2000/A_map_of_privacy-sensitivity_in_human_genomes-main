import pandas as pd

def addCitobandTo_STR(str_df, cb_df):
    """Adds the citoband the STR belongs to the STR dataframe as a new column
    """
    citobands = []
    for row in range(len(str_df.values)):
        citoband_found = False
        f_index = str_df["FirstIndex"][row]
        l_index = str_df["LastIndex"][row]
        chr = str_df["FastaHeader"][row]
        for citoband in cb_df.values:
            #print(citoband[0])
            if citoband[0] == chr:
                if int(citoband[1]) < int(f_index) and int(citoband[2])-1 > int(l_index) and not citoband_found:
                    citobands.append(citoband[3])
                    citoband_found = True
        if not citoband_found:
            for citoband_index in range (len(cb_df.values) - 1):
                citoband = cb_df.values[citoband_index]
                next_citoband = cb_df.values[citoband_index + 1]
                if citoband[0] == chr and next_citoband[0] == chr:
                    if int(citoband[1]) < int(f_index) and int(next_citoband[2])-1 > int(l_index) and not citoband_found:
                        citobands.append(citoband[3] + " & " + next_citoband[3])
                        citoband_found = True
            if not citoband_found:
                    citobands.append("citoband not found")
    
    str_df.insert(len(str_df.columns), "Citoband", citobands, True)


def calc_size1(dataframe):
    """Calculate size of STR by subtractong last index with first index
    More accurate so far, still the "ArrayLength" column (calc_size1 + 1) is better
    """
    sizes = []
    for row in dataframe.values:
        size = int(row[2]) - int(row[1])
        sizes.append(size)
    dataframe.insert(len(str_df.columns), "size", sizes, True)

def calc_size2(dataframe):
    """Calculate size of STR by multiplying pattern_size and copy_number
    Its less optimal
    """
    sizes = []
    for row in dataframe.values:
        size = float(row[3]) * float(row[4])
        sizes.append(size)
    dataframe.insert(len(str_df.columns), "size", sizes, True)


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

CitobandPath = "/home/androx/Documents/trabalho/citobands/cytobandFiltered.txt"
#STRPath_read = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_nov2014/HomosapiensHG38_FullGenome.txt"
#STRPath_write = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_nov2014/HomosapiensHG38_FullGenome_processed.txt"
STRPath_read = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanGenome.txt"
STRPath_write = "/home/androx/Documents/trabalho/datasets/STR/HomosapiensHG38_jan2016/HumanGenome_processed.txt"

citobands_df = pd.read_csv(CitobandPath, sep="\t", header=None)
str_df = str_dataframe = pd.read_csv(STRPath_read, sep="\t")


addCitobandTo_STR(str_df, citobands_df)

writeNewtxt(STRPath_write, str_df)

#print(str_df.columns)
#print(str_df["FastaHeader"][1])
