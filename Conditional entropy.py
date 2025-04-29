#This script takes in a .txt file with information about morphological paradigms
#and outputs a .txt file with information about conditional entropy.
#The input file is a tab-separated table.
#It has names of morphological categories going along row 1,
#and names of declension classes going along column 1.
#(The top-left cell is ignored.)
#All other cells contain the realisation of a given morphological category
#in a given declension class.
#The output file is also a tab-separated table.
#It contains conditional entropies for each morphological category
#given every other morphological category.
#It gives the average conditional entropy for each row and column,
#as well as the average conditional entropy for the whole paradigm.
#The script was written in May 2024 by Samuel Andersson, based entirely on:
#Ackerman, Farrell & Robert Malouf (2013) Morphological Organization:
#The low conditional entropy conjecture. In Language 89(3): 429-464.
#https://doi.org/10.1353/lan.2013.0054

import math
from collections import defaultdict

#Some of the code is easier to write if we're allowed to take the logarithm
#of 0. The result always gets multiplied by 0, so it doesn't actually matter
#what the function returns, as long as it doesn't throw an error, which
#Python's native math.log2() does
def log2(x):

    if x == 0:

        return 0

    else:

        return math.log2(x)

#Read in data
paradigm = []

with open("ot_paradigm_input_nofreq.txt", encoding = "utf-8") as f:

    paradigm = f.read().split("\n")

#The morphological categories, e.g. "nominative singular"
cells = paradigm[0].split("\t")[1:]

#The declension classes
classes = [x.split("\t")[0] for x in paradigm][1:]

#realisations[cell] = a list of how cell is realised in each declension class
realisations = defaultdict(lambda: [])

for i in range(len(cells)):

    realisations[cells[i]] = [x.split("\t")[i + 1] for x in paradigm][1:]

#pcr[cell][realisation] is the number of declension classes in which
#cell is realised as realisation.
#This implements equation (7) from Ackerman & Malouf (2013: 439)
#This definition depends on the assumption that the declension classes are
#equiprobable. In other words, we're using Ackerman & Malouf's (2013: 438)
#equation (5) rather than (6) on p. 439.
pcr = defaultdict(lambda: defaultdict(lambda: 0))

for cell in cells:

    for realisation in realisations[cell]:

        pcr[cell][realisation] = realisations[cell].count(realisation) / len(classes)

#H[cell] is the entropy of cell
#This produces the results in (8) from Ackerman & Malouf (2013: 439)
#Note that these calculations are optional in the sense that Hc isn't used
#by any later parts of this script.
Hc = defaultdict(lambda:  0)

for cell in cells:

    Hc[cell] = -sum([pcr[cell][realisation] * log2(pcr[cell][realisation]) for realisation in list(set(realisations[cell]))])

#pcrcr[cell1][realisation1][cell2][realisation2] is the number
#of declension classes in which
#cell1 is realised as realisation1, AND
#cell2 is realised as realisation2
#This implements equation (10) from Ackerman & Malouf (2013: 440)
#This definition depends on the assumption that the declension classes are
#equiprobable. In other words, we're using equation (5) from
#Ackerman & Malouf (2013: 438) rather than equation (6) from
#Ackerman & Malouf (2013: 439)
pcrcr = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))))

for cell1 in cells:

    for realisation1 in list(set(realisations[cell1])):

        for cell2 in cells:

            if cell2 == cell1:

                continue

            for realisation2 in list(set(realisations[cell2])):

                for i in range(len(classes)):

                    if realisations[cell1][i] == realisation1 and realisations[cell2][i] == realisation2:

                        pcrcr[cell1][realisation1][cell2][realisation2] += 1

                pcrcr[cell1][realisation1][cell2][realisation2] /= len(classes)

#cpcrcr[cell1][realisation1][cell2][realisation2] is the conditional probability
#that
#cell1 is realised as realisation1, GIVEN THAT
#cell2 is realised as realisation2,
#p(cell1 == realisation1|cell2 == realisation2)
#This implements equation (11) from Ackerman & Malouf (2013: 440)
cpcrcr = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))))

for cell1 in cells:

    for realisation1 in list(set(realisations[cell1])):

        for cell2 in cells:

            if cell2 == cell1:

                continue

            for realisation2 in list(set(realisations[cell2])):

                cpcrcr[cell1][realisation1][cell2][realisation2] = pcrcr[cell1][realisation1][cell2][realisation2] / pcr[cell2][realisation2]

#Hcc[cell1][cell2] is the conditional entropy of cell1 given knowledge of cell2,
#H(cell1|cell2)
#This implements equation (12) from Ackerman & Malouf (2013: 441)
Hcc = defaultdict(lambda: defaultdict(lambda: 0))

for cell1 in cells:

    for cell2 in cells:

        if cell2 == cell1:

            continue

        for realisation2 in list(set(realisations[cell2])):

            Hcc[cell1][cell2] += pcr[cell2][realisation2] * -sum([cpcrcr[cell1][realisation1][cell2][realisation2] * log2(cpcrcr[cell1][realisation1][cell2][realisation2]) for realisation1 in list(set(realisations[cell1]))])

#All code below this point assumes that all paradigm cells are equiprobable.
#In other words, we're using equation (15) from Ackerman & Malouf (2013: 441),
#instead of equation (16) on the same page.

#Ecol[cell1] is the average uncertainty in guessing cell1 based on another
#randomly chosen paradigm cell
#This implements part 1 of equation (17) from Ackerman & Malouf (2013: 442)
Ecol = defaultdict(lambda: 0)

for cell1 in cells:

    Ecol[cell1] = sum([Hcc[cell1][cell2] for cell2 in cells if not cell2 == cell1]) / (len(cells) - 1)

#Erow[cell1] is the average uncertainty in guessing a randomly chosen
#paradigm cell based on cell1
#This implements part 2 of equation (17) from Ackerman & Malouf (2013: 442)
Erow = defaultdict(lambda: 0)
    
for cell2 in cells:

    Erow[cell2] = sum(Hcc[cell1][cell2] for cell1 in cells if not cell1 == cell2) / (len(cells) - 1)

#Hp is the average conditional entropy of the whole paradigm
#It can be defined in terms of Ecol or Erow, both yield the same results
#This implements the Ecol-based part of
#equation (18) from Ackerman & Malouf (2013: 442)
Hp = sum([Ecol[cell] for cell in cells]) / len(cells)

#The code below creates a tab-separated table which reproduces
#Table 2 from Ackerman & Malouf (2013: 441)
output = "H(col|row)\t"

for cell in cells:

    output += f"{cell}\t"

output += "E[row]\n"

for row in cells:

    output += f"{row}\t"

    for col in cells:

        if col == row:

            output += "--\t"

        else:

            output += f"{round(Hcc[col][row], 3)}\t"

    output += f"{round(Erow[row], 3)}\n"

output += "E[col]\t"

for cell in cells:

    output += f"{round(Ecol[cell], 3)}\t"

output += str(round(Hp, 3))

#Save the output table to a .txt file
with open("ot_paradigm_output_nofreq.txt", mode = "w", encoding = "utf-8") as f:

    f.write(output)
