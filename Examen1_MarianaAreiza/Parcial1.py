#!/bin/bash
# Script para resolver los ejercicios del parcial 1 de Mariana Areiza Sánchez

from Bio import SeqIO

WT = SeqIO.read("Ecoli_phoA.fasta", "fasta")
Alpha = SeqIO.read("Ecoli_phoA_mutalpha.fasta", "fasta")
Beta = SeqIO.read("Ecoli_phoA_mutbeta.fasta", "fasta")

print(WT)
print(Alpha)
print(Beta)

# Ahora procedo a calcular la longitud de los nucleótidos
WT_nt = len(WT.seq)
Alpha_nt = len(Alpha.seq)
Beta_nt = len(Beta.seq)

print("La longitud de los nucleótidos WT es:", WT_nt)
print("La longitud de los nucleótidos Alpha es:", Alpha_nt)
print("La longitud de los nucleótidos Beta es:", Beta_nt)

# Ahora traducimos esos nucleótidos a proteína
WT_prot = WT.seq.translate()
Alpha_prot = Alpha.seq.translate()
Beta_prot = Beta.seq.translate()

# Luego, calculamos la longitud de las proteínas
WT_aa = len(WT_prot)
Alpha_aa = len(Alpha_prot)
Beta_aa = len(Beta_prot)

print("La longitud de la proteína de WT es:", WT_aa)
print("La longitud de la proteína Alpha es:", Alpha_aa)
print("La longitud de la proteína Beta es:", Beta_aa)

# Finalmente, guardamos los resultados obtenidos en un dataframe
dataframe = [
["Tipo", "Longitud_nt", "Longitud_aa"],
["WT", WT_nt, WT_aa],
["Alpha", Alpha_nt, Alpha_aa],
["Beta", Beta_nt, Beta_aa]
]

print("\nDataFrame:")
dataframe = [
    ["Tipo", "Longitud_nt", "Longitud_aa"],
    ["WT", WT_nt, WT_aa],
    ["Alpha", Alpha_nt, Alpha_aa],
    ["Beta", Beta_nt, Beta_aa]
]

print("\nDataFrame:")
for fila in dataframe:
    print(fila)

# Y ahora realizamos el gráfico
import matplotlib.pyplot as plt

labels = ["WT", "Alpha", "Beta"]
nt = [WT_nt, Alpha_nt, Beta_nt]
aa = [WT_aa, Alpha_aa, Beta_aa]

x = range(len(labels))

plt.bar(x, nt)
plt.xticks(x, labels)
plt.ylabel("Longitud nucleótidos")
plt.title("Longitud de secuencias")
plt.savefig("nucleotidos.png")

plt.clf()

plt.bar(x, aa)
plt.xticks(x, labels)
plt.ylabel("Longitud proteínas")
plt.title("Longitud proteinas")
plt.savefig("proteinas.png")

# Conclusión

# Primero calculé la longitud de las secuencias de nucleótidos para cada cepa y luego traduje las secuencias para obtener la longitud de las 
# proteínas. Observé que la cepa WT tiene 1416 nucleótidos y produce una proteína de 472 aminoácidos, lo que indica que corresponde a la 
# secuencia completa y normal del gen y, por lo tanto, debería tener la actividad normal de la enzima fosfatasa alcalina. 

# La cepa Alpha también tiene 1416 nucleótidos y 472 aminoácidos, es decir, la longitud es igual a la de la cepa silvestre. Esto sugiere que la 
# mutación no cambió el tamaño de la proteína, por lo que probablemente su actividad sea parecida a la de la WT o solo un poco diferente. Pero por
# otro lado, como se puede observar en el DataFrame y en las gráficas, la cepa Beta tiene 1359 nucleótidos y una proteína de 453 aminoácidos, lo 
# que indica que la secuencia es más corta y que la proteína quedó truncada. Debido a esto, es probable que la enzima no funcione correctamente y 
# que su actividad de fosfatasa alcalina sea menor.

# Por esta razón, esperaría que la línea con mayor actividad en la gráfica corresponda a la WT, la actividad intermedia o similar a Alpha, y la 
# actividad más baja a Beta, ya que su proteína es más corta y probablemente menos funcional.
