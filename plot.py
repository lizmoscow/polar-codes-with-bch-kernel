import math
import pandas as pd
import matplotlib.pyplot as plt
#не забываем отредактировать working directory
data1 = pd.read_csv("out/15_5_7.csv")
X = data1['ratio'] #первая колонка с абсциссами (отношение сигнал / шум)
Y = data1['errors'] #вторая колонка с ординатами (отношение ошиок к общему числу декодирований)
with open('out/15_5_7.txt', 'w') as file:
    X = [x for x in X if Y[int(x * 10 - 1)] != 0]
    Y = [y for y in Y if y != 0]
    for j in range(len(Y)):
        Y[j] = math.log(Y[j], 10)
        file.write("%f %f \n" % (X[j], Y[j]))


