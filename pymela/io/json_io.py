import json

def parse(ini_file):
    with open(ini_file) as f:
        data = json.load(f)
    return data
#-------------------------------

def dump(iniDict):
    print(iniDict)
#-------------------------------
