   

# awk 'BEGIN{FS=""} !/^>/{print NF}' /Users/pauline/Desktop/Stage_MNHN/PF07966.15.B |uniq

#

import struct
import sys
# on enregistre un entier, un réel et 4 caractères
i = 10
nom = "seq1"
seq = "ATGGBEGDYHD"


with open("uneseq.ascii", "w") as f:
     for i in range(5):
         f.write(f">{nom}\n")
         f.write(f"{seq}\n")

# écriture
with open("uneseq.bin", "ab") as fb:
     fb.write(struct.pack("i", 20))
     for i in range(20):
         fb.write(struct.pack("i", len(nom)))
         fb.write(struct.pack("i", len(seq)))
         #fb.write(struct.pack("d", x))
         # il faut convertir les caractères en bytes
         nom2 = nom.encode("ascii")
         seq2 = seq.encode("ascii")
         monFormat=str(len(nom))+"s"
         fb.write(struct.pack(monFormat, nom2))
         monFormat=str(len(seq))+"s"
         fb.write(struct.pack(monFormat, seq2))


# lecture
with open("uneseq.bin", "rb") as fb:
     nbSeq = struct.unpack("i", fb.read(4))[0]
     lesSeq=[]
     for i in range(nbSeq):
         lgNom = struct.unpack("i", fb.read(4))[0]
         print(lgNom)
         lgSeq = struct.unpack("i", fb.read(4))[0]
         monFormat=str(lgNom)+"s"
         nomLu = struct.unpack(monFormat, fb.read(lgNom))[0]
         monFormat=str(lgSeq)+"s"
         seqLu = struct.unpack(monFormat, fb.read(lgSeq))[0]
         lesSeq.append((nomLu, seqLu))
# affichage pour vérifier que les données ont été bien lues
print(lgNom, nomLu)
print(lgSeq, seqLu)
for n,s in lesSeq:
     print(n,s)
