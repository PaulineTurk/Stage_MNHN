import os
import shutil


def copyFileRenameComplete(path_folder, path_original, list_delay_number):

    if os.path.isdir(path_folder):
        shutil.rmtree(path_folder) 
    os.mkdir(path_folder)

    for delay_num in list_delay_number:
        for kp_SeqChoice in ["k", "p"]:
            path_target = f"{path_folder}/cube{delay_num}_{kp_SeqChoice}.sh"
            shutil.copy(path_original, path_target)


        
            file = open(path_target, "r")
            list_line = file.readlines()
            list_line[9] = f"delay_num = {delay_num}\n"
            list_line[10] = f"kp_SeqChoice = {kp_SeqChoice}\n"

            file = open(path_target, "w")
            file.writelines(list_line)


            file.close()



if __name__ == '__main__': 
    path_folder = '/Users/pauline/Desktop/Stage_MNHN/CubesGenerator/ScriptCubes'
    path_original = '/Users/pauline/Desktop/Stage_MNHN/CubesGenerator/cube_original.sh'
    list_delay_number = [k for k in range(-10, 11) if k!=0]

    copyFileRenameComplete(path_folder, path_original, list_delay_number)
