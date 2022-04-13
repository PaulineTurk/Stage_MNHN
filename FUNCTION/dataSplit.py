import os, shutil
from pathlib import Path
import numpy as np
from timer import Timer
from fileNumber import countFile
from sklearn.model_selection import train_test_split




def addFileFromFolder(file_name, folder_name_source, folder_path_target, extension_file_name_target):
    path_file_source = folder_name_source + "/" + file_name
    accession_num = os.path.basename(file_name).split(".")[0] + '.' + os.path.basename(file_name).split(".")[1]
    path_file_target = folder_path_target + "/" + accession_num + extension_file_name_target
    shutil.copy2(path_file_source, path_file_target)


def trainTestSplit(name_folder_fasta_total, name_folder_data_train_test, name_folder_data_train, name_folder_data_test, percentage_train):
    t = Timer()
    t.start()
    path_folder_data_train_test = name_folder_data_train_test + '/'
    if os.path.isdir(path_folder_data_train_test):
        shutil.rmtree(path_folder_data_train_test) 
    os.mkdir(path_folder_data_train_test)

    path_folder_data_train = path_folder_data_train_test + '/' + name_folder_data_train + '/'
    if os.path.isdir(path_folder_data_train):
        shutil.rmtree(path_folder_data_train) 
    os.mkdir(path_folder_data_train)

    path_folder_data_test = path_folder_data_train_test + '/' + name_folder_data_test + '/'
    if os.path.isdir(path_folder_data_test):
        shutil.rmtree(path_folder_data_test) 
    os.mkdir(path_folder_data_test)
    

    path_folder_fasta_total = Path(name_folder_fasta_total+ '/')
    files_in_path_folder_fasta_total = path_folder_fasta_total.iterdir()
    data_name = []

    for file_path in files_in_path_folder_fasta_total:
        file_name = str(file_path).split("/")[-1]
        data_name.append(file_name)

    fraction_train = percentage_train/100
    x_train ,x_test = train_test_split(data_name, train_size = fraction_train)  
    # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
    # train_sizefloat or int, default=None
    # If float, should be between 0.0 and 1.0 and represent the proportion of the dataset to include in the train split. 
    # If int, represents the absolute number of train samples. 
    # If None, the value is automatically set to the complement of the test size.

    # Rq.
    # shufflebool, default=True
    # Whether or not to shuffle the data before splitting. If shuffle=False then stratify must be None.

    for file_name in x_train:
        addFileFromFolder(file_name, name_folder_fasta_total, path_folder_data_train, ".train")

    countFile(path_folder_data_train)
    
    for file_name in x_test:
        addFileFromFolder(file_name, name_folder_fasta_total, path_folder_data_test, ".test")
    countFile(path_folder_data_test)

    t.stop("Split data_total in data_train and data_test")






if __name__ == '__main__': 

    # data split train/test
    # not trimmed version

    path_folder_data = "/Users/pauline/Desktop/data/"    # dossier à créer ou l'on veut mettre toutes les data entrées/sorties
    name_folder_data_train = "PfamTrain"
    name_folder_data_test = "PfamTest"
    name_folder_fasta_total = "Pfam_fasta_99"
    path_folder_fasta_total = path_folder_data + name_folder_fasta_total

    list_percentage = [0.05, 0.5, 5, 50]   # respectively 9, 98, 981, 9816 seeds for the data_train 
    for percentage_train in list_percentage:
        name_folder_data_train_test = "PfamSplit" + "_" + str(percentage_train)
        path_folder_data_train_test = path_folder_data + name_folder_data_train_test
        trainTestSplit(path_folder_fasta_total, path_folder_data_train_test, name_folder_data_train, name_folder_data_test, percentage_train) 
    
    
    
    
    # trimmed version

    #path_folder_data = "/Users/pauline/Desktop/data/"    # dossier à créer ou l'on veut mettre toutes les data entrées/sorties
    #name_folder_data_train = "PfamTrain"
    #name_folder_data_test = "PfamTest"
    #name_folder_fasta_total = "Pfam_fasta_99_trimmed"
    #path_folder_fasta_total = path_folder_data + name_folder_fasta_total

    #list_percentage = [0.05, 0.5, 5, 50]   # respectively 9, 91, 910, 9101 seeds for the data_train
    #for percentage_train in list_percentage:
    #    name_folder_data_train_test = "PfamSplit_trimmed" + "_" + str(percentage_train)
    #    path_folder_data_train_test = path_folder_data + name_folder_data_train_test
    #    trainTestSplit(path_folder_fasta_total, path_folder_data_train_test, name_folder_data_train, name_folder_data_test, percentage_train) 
    pass