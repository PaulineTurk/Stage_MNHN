def characterList_prob(total_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "-"], 
                   excluded_character = ["B", "Z", "X", "-"]):
    for elem in excluded_character:
        total_list.remove(elem)
    return total_list



def characterList():
    return ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


