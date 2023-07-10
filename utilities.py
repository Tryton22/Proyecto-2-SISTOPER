def colored(seq):
    # Diccionario para colorear los nucleotidos
    bcolors = {
        'A': '\033[92m',  # Verde
        'C': '\033[94m',  # Azul
        'G': '\033[93m',  # Amarillo
        'T': '\033[91m',  # Rojo
        'U': '\033[91m',  # Rojo (para ARN)
        'reset': '\033[0;0m'  # Resetear el color
    }

    tmpStr = ""

    # Colorear cada nucleotido segun corresponda
    for nuc in seq:
        if nuc in bcolors:
            tmpStr += bcolors[nuc] + nuc
        else:
            tmpStr += bcolors['reset'] + nuc

    return tmpStr + '\033[0;0m'


def read_FASTA(filePath):
    # Leer el archivo FASTA y almacenar las lineas en una lista
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    # Recorrer las lineas del archivo FASTA
    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict

