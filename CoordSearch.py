import gzip
import os
import math
import xlwt

rangeList = list()
LanthAtom = b"C1 " or b"Ln "

molNum = len(os.listdir()) - 9
for i in range(1, molNum):

    obsList = list()
    bondList = list()

    if i < 10:
        global x_file
        x_file = gzip.open(os.path.join('M000' + str(i), "output.gz"), 'r')
    else:
        x_file = gzip.open(os.path.join('M00' + str(i), "output.gz"), 'r')


    lines = x_file.readlines()
    length = len(lines)
    NeededAtoms = list()

    for line in range(length):
        if lines[line].startswith(b"                     Cartesian Coordinates", 0):  # ДОБЫВАЕМ НУЖНЫЕ КУСКИ OUTPUT'А
            i = line + 4
            while len(lines[i]) > 30:
                obsList.append(lines[i])
                i = i + 1
        if lines[line].startswith(b" Bond Orders              Mulliken", 0):
            j = line + 1
            while len(lines[j]) > 30:
                bondList.append(lines[j])
                j = j + 1

    datList = list()
                                            # НАХОДИМ ДЕКАРТОВЫЕ КООРДИНАТЫ АТОМА ЛАНТАНИДА
    for line in bondList:
        if line.count(LanthAtom) > 0:
            datList.append(str(line.decode().split(" ")[11]) + " ")

    lanthCoordList = list()
    FinallyCoords = list()
    lantAtomCoord = list()

    for line in obsList:
        if line.count(LanthAtom) > 0:
            lantAtomCoord = line.decode().split(" ")
            for i in lantAtomCoord:
                if i != '':
                    FinallyCoords.append(i)

    for line in obsList:
        for i in range(len(datList)):
            if line.count(datList[i].encode()) > 0:
                lanthCoordList = line.decode().split(" ")
                for i in lanthCoordList:
                    if i != '':
                        FinallyCoords.append(i)


    for elem in FinallyCoords:
        if len(elem) < 5:
            FinallyCoords.remove(elem)

    for elem in FinallyCoords:              #Дублируем цикл с удалением, чтобы не оставалось ненужных элементов
        if len(elem) < 5:
            FinallyCoords.remove(elem)



    for i in range(3, len(FinallyCoords), 3):
        rangeList.append(math.sqrt((float(FinallyCoords[0]) - float(FinallyCoords[i]))**2 + (float(FinallyCoords[1]) -
        float(FinallyCoords[i + 1]))**2 + (float(FinallyCoords[2]) - float(FinallyCoords[i + 2])*2)))

    rangeList.append("\n")

print(rangeList)
wb = xlwt.Workbook()
ws = wb.add_sheet("Bonds")


for i in range(len(rangeList)):
    ws.write(i, 0, rangeList[i])

wb.save("spartanOutput.xls")
