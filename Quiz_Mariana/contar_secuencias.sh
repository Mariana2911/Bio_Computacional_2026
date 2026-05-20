# Quiz 1 - Fundamentos Biología Computacional
# Script para contar el número total de secuencias

#!/bin/bash

# Luego buscamos el patrón ">" ya que como lo vimos en clase, todas las líneas en un archivo FASTA comienzan con este
grep -ic ">" CytBDNA.txt

