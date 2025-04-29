#This script takes in a .txt file with information about morphological paradigms
#and outputs a .txt file with information about conditional entropy.
#The input file is a tab-separated table.
#It has names of morphological categories going along row 1,
#and names of declension classes going along column 1.
#The final row contains token frequencies (raw counts) of each morphological
#category, and
#the final column contains type frequencies (raw counts) of each declension.
#(The corners of the table are ignored.)
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

with open("ot_paradigm_input_frequency.txt", encoding = "utf-8") as f:

    paradigm = f.read().split("\n")

#The morphological categories, e.g. "nominative singular"
cells = paradigm[0].split("\t")[1:-1]

#The declension classes
classes = [x.split("\t")[0] for x in paradigm][1:-1]

#realisations[cell] = a list of how cell is realised in each declension class
realisations = defaultdict(lambda: [])

for i in range(len(cells)):

    realisations[cells[i]] = [x.split("\t")[i + 1] for x in paradigm][1:-1]

#The token frequencies for each paradigm cell
tokenC = [int(x) for x in paradigm[-1].split("\t")[1:-1]]

#The type frequencies for each declension class
typeD = [int(y) for y in [x.split("\t")[-1] for x in paradigm][1:-1]]

#pcr[cell][realisation] is the probabilitiy that a given lexeme
#has paradigm cell cell realised as realisation.
#This implements equation (7) from Ackerman & Malouf (2013: 439).
pcr = defaultdict(lambda: defaultdict(lambda: 0))

for cell in cells:

    for realisation in list(set(realisations[cell])):

        for i in range(len(classes)):

                if  realisations[cell][i] == realisation:

                    pcr[cell][realisation] += typeD[i] / sum(typeD)

#pcrcr[cell1][realisation1][cell2][realisation2] is the number
#of declension classes in which
#cell1 is realised as realisation1, AND
#cell2 is realised as realisation2
#This implements equation (10) from Ackerman & Malouf (2013: 440)
pcrcr = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))))

for cell1 in cells:

    for realisation1 in list(set(realisations[cell1])):

        for cell2 in cells:

            if cell2 == cell1:

                continue

            for realisation2 in list(set(realisations[cell2])):

                for i in range(len(classes)):

                    if realisations[cell1][i] == realisation1 and realisations[cell2][i] == realisation2:

                        pcrcr[cell1][realisation1][cell2][realisation2] += typeD[i] / sum(typeD)

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

#Ecol[cell1] is the average uncertainty in guessing cell1 based on another
#randomly chosen paradigm cell
#This implements part 1 of equation (17) from Ackerman & Malouf (2013: 442)
#Although they don't make this explicit, the numbers in Ackerman & Malouf (2013)
#can only be arrived at if it is crucially assumed that c1 and c2 in (17) are
#different cells. In other words, we do not consider the entropy of a cell
#given itself, which is definitionally 0. This affects the probability
#calculations below. Instead of using the probability of a paradigm cell
#directly, the calculations below use the conditional probability of each cell,
#given that we're ignoring the cell whose Ecol value we're calculating.
Ecol = defaultdict(lambda: 0)

for cell1 in cells:

    for i in range(len(cells)):

        if cells[i] == cell1:

            continue

        Ecol[cell1] += Hcc[cell1][cells[i]] * (tokenC[i] / sum([tokenC[j] for j in range(len(tokenC)) if not cells[j] == cell1]))

#Erow[cell1] is the average uncertainty in guessing a randomly chosen
#paradigm cell based on cell1
#This implements part 2 of equation (17) from Ackerman & Malouf (2013: 442)
#Although they don't make this explicit, the numbers in Ackerman & Malouf (2013)
#can only be arrived at if it is crucially assumed that c1 and c2 in (17) are
#different cells. In other words, we do not consider the entropy of a cell
#given itself, which is definitionally 0. This affects the probability
#calculations below. Instead of using the probability of a paradigm cell
#directly, the calculations below use the conditional probability of each cell,
#given that we're ignoring the cell whose Erow value we're calculating.
Erow = defaultdict(lambda: 0)
    
for cell2 in cells:

    for i in range(len(cells)):

        if cells[i] == cell2:

            continue

        Erow[cell2] += Hcc[cells[i]][cell2] * (tokenC[i] / sum([tokenC[j] for j in range(len(tokenC)) if not cells[j] == cell2]))

#Hp is the average conditional entropy of the whole paradigm
#It can be defined in terms of Ecol or Erow, both yield the same results
#This implements the Ecol-based part of
#equation (18) from Ackerman & Malouf (2013: 442)
Hp = 0

for i in range(len(cells)):

    Hp += Ecol[cells[i]] * (tokenC[i] / sum(tokenC))

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
with open("ot_paradigm_output_frequency.txt", mode = "w", encoding = "utf-8") as f:

    f.write(output)
