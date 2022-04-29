def readFastaMul(nomFi):
    seq=""
    nom=""
    lesSeq=[]
    with open(nomFi,"r") as f:  
        for l in f:
            if l[0] == '>':
                if seq != "":
                    tmp=(nom,seq)
                    lesSeq.append(tmp)
                nom=l[1:-1]
                seq=""
            else:
                seq = seq + l.strip()
    if seq != "":
        tmp=(nom,seq)
        lesSeq.append(tmp)
    return lesSeq





def readFasta(path_fasta_file):
    """
    :return: tuples of (title, seq)
    """
    title = None
    data = None
    with open(path_fasta_file, "r") as f:
        for line in f:
            if line[0] == ">":
                if title:
                    yield (title.strip(), data)   # USE GENERATOR INSTEAD OF for and in ????????????
                title = line[1:]
                data = ''
            else:
                data += line.strip()
        if not title:
            yield None
    yield (title.strip(), data) 