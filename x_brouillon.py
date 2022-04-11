import os

def countFile(folder):
    path, dirs, files = next(os.walk(folder))
    file_count = len(files)
    return file_count



if __name__ == '__main__': 
    print("Pfam_fasta_99:", countFile("/Users/pauline/Desktop/data/Pfam_fasta_99"))  
    print("Pfam_fasta_99_trimmed:", countFile("/Users/pauline/Desktop/data/Pfam_fasta_99_trimmed")) 