# File of homebrew debugging tools.
import numpy as np

verbose_debugging_mode = False
def vprint(msg, outFile = None, printToTerminal = True): 
    if verbose_debugging_mode == True: 
        if outFile == None: 
            outFile = 'verbose.log' 
        with open(outFile, "a") as f: 
            f.write(msg + '\n') 
        if printToTerminal: 
            print(msg) 
 
def vdiag(func, preface = '', printToTerminal = True, outFile = None, **kwargs):
    if verbose_debugging_mode == True: 
        vprint(preface, outFile, printToTerminal) 
        for key, value in kwargs.iteritems(): 
            if isinstance(value,np.ndarray): 
                vprint('[' + func + ']: max(' + key + ') = ' + str(max(value)),outFile,printToTerminal) 
                vprint('[' + func + ']: min(' + key + ') = ' + str(min(value)),outFile,printToTerminal) 
                vprint('[' + func + ']: size(' + key + ') = ' + str(np.size(value)),outFile,printToTerminal) 
                vprint('[' + func + ']: memory usage of ' + key + ' = ' + str(value.nbytes),outFile,printToTerminal) 
            else: 
                print('[' + func + ']: ' + key + ' = ' + str(value)) 




