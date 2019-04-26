#!usr/bin/env python
from collections import Counter

"""Functions for dealing with DNA sequences."""



CODONS = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

"""Class DNA """
class dna:
    def __init__(self):
        print("Nucleotide counter:")
        self.outputData=[]
        self.proteinSeq=''

    def __iter__(self):
        self.n=0
        return self

    def __next__(self):
        if self.n<len(self.outputData):
            self.n += 1
            return (self.outputData[self.n-1])

        else:
            raise StopIteration

    def countNucl(self):
        """The function counts nucleotides in sequence and print it in format number_of_nucleotide:amount"""
        self.outputData=[]
        with open(readFile,"r") as gene:
            next (gene)
            for line in gene:
                    nucleotide_counts = Counter(char for line in gene for char in line.strip())
            #print (nucleotide_counts) litery to klucze a ilosc danych nukleotydow to wartosci
            for (nucleotide, count) in nucleotide_counts.most_common():
                    print ("Number of %s 's : %d" % (nucleotide, count))
                    self.outputData.append((nucleotide,count))
                    #s.write("Number of %s 's : %d" % (nucleotide, count) + '\n')

    def translate(self):
        """The function translates each three of nucleotides to protein"""
        with open(readFile, "r") as gene:
            next(gene)
            for line in gene:
                for i in range(0, len(line)-3):
                    if i % 3 == 0:
                        self.proteinSeq += CODONS[line[i] + line[i + 1] + line[i + 2]]
            print(self.proteinSeq)

    def saveOutput(self):
        """The function saves a record of function countNucl to file"""
        with open(saveFile, "w") as outputFile:
        #iterator:
            for (nucleotide,count) in self:
                outputFile.write("Number of %s 's : %d \n" % (nucleotide, count))
            outputFile.write(self.proteinSeq)


if __name__ == '__main__':
    readFile = input("Input a filename to read:")
    saveFile = input("Input a filename to write:")
    code = dna ()
    code.countNucl()
    #for (nucleotide, count) in code.outputData:
        #print (nucleotide, count)
    code.translate()
    code.saveOutput()