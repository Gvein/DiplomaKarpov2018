import os
import gzip
import xlwt

FinalList = list()
PropList = list()
SplitList = list()

filename = "EbindTest_EDF-2_6-31G%2A.spardir"
molNum = len(os.listdir(filename)) - 7
print(len(os.listdir()))
print(os.listdir())

for i in range(1, molNum):
# ОТКРЫВАЕМ И РАСПАКОВЫВАЕМ ФАЙЛ СО СВОЙСТВАМИ МОЛЕКУЛЫ
    if i < 10:
        global x_file
        x_file = gzip.open(os.path.join(filename,'M000' + str(i), "proparc.gz"), 'r')
    else:
        x_file = gzip.open(os.path.join(filename, 'M00' + str(i), "proparc.gz"), 'r')

    lines = x_file.readlines()
    length = len(lines)

# ВЫРЕЗАЕМ НУЖНЫЕ СТРОКИ ИЗ OUTPUT SPARTAN14
    for line in range(length):

        if lines[line].startswith(b"PROP VALUE E_HOMO", 0):
            PropList.append(lines[line])
        if lines[line].startswith(b"PROP VALUE E_LUMO", 0):
            PropList.append(lines[line])
        # if lines[line].startswith(b"PROP VALUE TOT_ELECTRONEG", 0):
        #     PropList.append(lines[line])
        # if lines[line].startswith(b"PROP VALUE TOT_HARDNESS", 0):
        #     PropList.append(lines[line])
        # if lines[line].startswith(b"PROP VALUE EST_POLARIZ", 0):
        #     PropList.append(lines[line])
        if lines[line].startswith(b"PROP VALUE DIPOLE_MAG =", 0):
            PropList.append(lines[line])
        if lines[line].startswith(b"PROP VALUE CORR_ENERGIES", 0):
            PropList.append(lines[line + 1])

    for j in PropList:
        j.decode()

for q in range(0, len(PropList)):
    SplitList.append(PropList[q].split())
print(SplitList)

# МУТИМ МАТРИЦУ

StNum = molNum
ColNum = 4

Matrix = [[0 for x in range(StNum)] for y in range(ColNum)]


Matrix[0][0] = "DIPOLE_MOMENTUM, Debye"
Matrix[1][0] = "TOTAL ENERGY, Eh"
Matrix[2][0] = "E_HOMO,eV"
Matrix[3][0] = "E_LUMO,eV"

# ЗАПОЛНЯЕМ МАТРИЦУ

q1 = 1
for i in range(0, len(SplitList), 4):
    Matrix[0][q1] = SplitList[i][4]

    q1 += 1

q2 = 1
for i in range(1, len(SplitList), 4):
    Matrix[1][q2] = SplitList[i][0]

    q2 += 1

q3 = 1
for i in range(2, len(SplitList), 4):
    Matrix[2][q3] = SplitList[i][4]
    print(SplitList[i][4])
    q3 += 1

q4 = 1
for i in range(3, len(SplitList), 4):
    Matrix[3][q4] = SplitList[i][4]
    q4 += 1

for i in range(ColNum):
    for j in range(1, StNum):
        for k in str(Matrix[i][j]):
            print(Matrix[i][j])
            k.replace("b" or "'", " ")

wb = xlwt.Workbook()
ws = wb.add_sheet("Properties")
ws.write(0, 0, "DIPOLE_MOMENTUM, Debye")
ws.write(0, 1, "TOTAL ENERGY, Eh")
ws.write(0, 2, "E_HOMO,eV")
ws.write(0, 3, "E_LUMO,eV")

for i in range(1, StNum):
    ws.write(i, 0, Matrix[0][i].decode())
    ws.write(i, 1, Matrix[1][i].decode())
    ws.write(i, 2, Matrix[2][i].decode())
    ws.write(i, 3, Matrix[3][i].decode())


wb.save("spartanOutput.xls")
