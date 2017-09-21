__author__ = 'VikiQiu'
import csv
import numpy as np

# filename = '../data/L-I-20C.csv'

def getData(filename):
    I = []
    P = []
    U = []
    f = open(filename)
    f_csv = csv.reader(f)
    headers = next(f_csv)
    for row in f_csv:
        I.append(float(row[0]))
        P.append(float(row[1]))
        U.append(float(row[2]))
    f.close()
    return np.array(I), np.array(P), np.array(U)

def getData2(filename):
    f = []
    s21 = []
    csv_file = open(filename)
    f_csv = csv.reader(csv_file)
    headers = next(f_csv)
    for row in f_csv:
        f.append(float(row[0]))
        s21.append(float(row[1]))
    csv_file.close()
    f = np.array(f)*1e9
    s21 = np.array(s21)
    h0 = 10**(s21[0]/20)
    Hf = (10**(s21/20)/h0)**2
    return f, s21, Hf

# I, P, U = getData(filename)
# print(I)
