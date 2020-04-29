import os

'''
For the given path, get the List of all files in the directory tree 
'''
def getListOfFiles(dirName):
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(dirName):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    return listOfFiles

def main():
    dirName="./Beta40_rs4.0_lambda0.2301"

    for elem in getListOfFiles(dirName):
        print(elem)    
                
if __name__ == '__main__':
    main()
