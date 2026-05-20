# Clase filogenia

# install.packages("ape")
# install.packages("BiocManager")
# BiocManager::install("msa")
# install.packages("phangorn")
# install.packages("seqinr")

#-----------------------------------------------------------
# 1. Instalar y cargar librerías
#-----------------------------------------------------------
# if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
# if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")
# if (!requireNamespace("msa", quietly = TRUE)) install.packages("msa")a


library(ape)
library(phangorn)
library(msa)   # Para alinear las secuencias (multiple sequence aln)

#-----------------------------------------------------------
# 2. Leer las secuencias originales (sin alinear) desde GitHub
#-----------------------------------------------------------

# ADN
cytb_dna_raw <- readDNAStringSet(
  "https://raw.githubusercontent.com/lauraalazar/BiologiaComputacional/main/CytBDNA.txt",
  format = "fasta"
)

# Proteína
cytb_prot_raw <- readAAStringSet(
  "https://raw.githubusercontent.com/lauraalazar/BiologiaComputacional/main/CytBProt.txt",
  format = "fasta"
)

#-----------------------------------------------------------
# 3. Alinear las secuencias con MSA (ClustalW)
#-----------------------------------------------------------
alignment_dna <- msa(cytb_dna_raw, method = "ClustalW")
alignment_prot <- msa(cytb_prot_raw, method = "ClustalW")

# Convertir los alineamientos a formato DNAbin y AAbin
alignment_dna_bin <- as.DNAbin(alignment_dna)
alignment_prot_bin <- as.AAbin(alignment_prot)

#-----------------------------------------------------------
# 4. Convertir a formato phyDat (para phangorn)
#-----------------------------------------------------------
cytb_dna_phy <- phyDat(alignment_dna_bin, type = "DNA")
cytb_prot_phy <- phyDat(alignment_prot_bin, type = "AA")

#-----------------------------------------------------------
# 5. Calcular matrices de distancia evolutiva
#-----------------------------------------------------------

# Modelo Jukes-Cantor (ADN)
dist_dna <- dist.ml(cytb_dna_phy, model = "JC69")

# Modelo JTT (Proteína)
dist_prot <- dist.ml(cytb_prot_phy, model = "JTT")

#-----------------------------------------------------------
# 6. Construir árboles filogenéticos
#-----------------------------------------------------------

# Neighbor-Joining (NJ)
tree_nj_dna <- NJ(dist_dna)
tree_nj_prot <- NJ(dist_prot)

# UPGMA
tree_upgma_dna <- upgma(dist_dna)
tree_upgma_prot <- upgma(dist_prot)

#-----------------------------------------------------------
# 8. Visualización
#-----------------------------------------------------------
par(mfrow = c(2, 2))
plot(tree_nj_dna, main = "CytB ADN - Neighbor Joining (JC69)", cex = 0.7)
plot(tree_upgma_dna, main = "CytB ADN - Upgma", cex = 0.7)
plot(tree_nj_prot, main = "CytB Proteína - Neighbor Joining (JTT)", cex = 0.7)
plot(tree_upgma_prot, main = "CytB Proteína - Upgma", cex = 0.7)