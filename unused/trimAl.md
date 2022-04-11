# Part1: General Information
### Please cite:
Salvador Capella-Gutierrez, Jose M. Silla-Martinez and Toni Gabaldon. trimAl: a tool for automated alignment trimming (2009).

### The general command line used
trimal -in \<inputfile\> -out \<outputfile\> -\<trimming type\>
            
### Documentation link
https://vicfero.github.io/trimal/whatcanido.html            
        
           

# Part2: User Guide: 
## To trim the alignments in the folder Pfam_fasta 
### - with nogaps option 
1. cd X/trimAl/source   # replace X by your path to trimAl 

2. mkdir Y/Pfam_fasta_trimAl_nogaps # replace Y by your path to the BLOSUM folder downloaded from this repository
            
3. for file in Y/Pfam_fasta/*; do ./trimal -in "$file" -out Y/Pfam_fasta_trimAl_nogaps/"$(basename "$file.trimAl.nogaps")" -nogaps; done

                 
### - with gappyout option            
1. cd X/trimAl/source   # replace X by your path to trimAl 
            
2. mkdir Y/Pfam_fasta_trimAl_gappyout # replace Y by your path to the BLOSUM folder downloaded from this repository
            
3. for file in Y/Pfam_fasta/*; do ./trimal -in "$file" -out Y/Pfam_fasta_trimAl_gappyout/"$(basename "$file.trimAl.gappyout")" -gappyout; done

### Other TrimAl options that could be used:
-noallgaps 
            
-strict  
            
-strictplus
            
-automated1         



#### test 2:
1. cd X/trimAl/source   # replace X by your path to trimAl 
            
2. mkdir Y/Pfam_fasta_trimAl_gappyout # replace Y by your path to the BLOSUM folder downloaded from this repository
            
3. for file in Y/Pfam_fasta/*; do ./trimal -in "$file" -out Y/Pfam_fasta_trimAl_gappyout/"$(basename "$file.trimAl.gappyout")" -gappyout; done

trimal -in example1 -out output2 -htmlout output2.html -gt 0.8 -st 0.001



cd /Users/pauline/Desktop/trimAl/source
mkdir /Users/pauline/Desktop/data/PfamFastaTrimAl_50
for file in /Users/pauline/Desktop/data/Pfam_fasta/*; do ./trimal -in "$file" -out /Users/pauline/Desktop/data/PfamFastaTrimAl_50/"$(basename "$file.trimAl_50")" -htmlout output2.html -gt 0.5 -st 0.001; done

# Error: the symbol 'B' accesing the matrix is not defined in this object