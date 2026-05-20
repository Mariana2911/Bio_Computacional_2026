from Bio import SeqIO

# Leer secuencia desde FASTA
record = SeqIO.read("mystery_seq.fasta", "fasta")

seq = record.seq

# Transcripción a ARN
print("ARN:", mi_secuencia.transcribe())

# Traducción a proteína
print("Proteína:", mi_secuencia.translate())

# Traducir (frame 1 por defecto)
#protein = seq.translate(to_stop=False)

#print("Secuencia proteica:")
#print(protein)

#print("\nLongitud en aminoácidos:")
#print(len(protein))
