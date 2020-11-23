import json
import io_conventions as ioConv

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

def makeInputChecks(infoTag,infoDict,runType):
    if infoTag not in ioConv.inputInfoTags[runType]:
        raise ValueError('Expected object %s in JSON input file' % (infoTag))

    for input in ioConv.expectedInput[infoTag]:
        if input not in infoDict:
            raise ValueError('Expected entry "%s" in object "%s" of JSON input file' % (input,infoTag))
#-------------------------------
