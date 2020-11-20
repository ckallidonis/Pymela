import json

def parse(ini_file):
    with open(ini_file) as f:
        data = json.load(f)
    return data
#-------------------------------

def dumpDictObject(dictObj,printMessage=''):
    if printMessage != '':
        print(printMessage)
    print(json.dumps(dictObj, indent = 2, sort_keys=True))
#-------------------------------
