import os

def countFile(folder):
    path, dirs, files = next(os.walk(folder))
    file_count = len(files)
    print(file_count)
    return file_count