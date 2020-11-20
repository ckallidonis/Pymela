import sys, os
sys.path.insert(len(sys.path), '../pymela/')

import pymela.io.json_io as JSONio

iniDict = JSONio.parse('/Users/kallidoc/work_local/pdf_distillation/analysis/scratch/test_json.py')

JSONio.dump(iniDict)

if iniDict["languages"][1] == "English":
    print('OK')
else:
    print('NOT OK')
