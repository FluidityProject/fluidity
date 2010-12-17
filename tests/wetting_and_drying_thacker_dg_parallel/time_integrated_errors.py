#!/usr/bin/python

import time_integrated_errors_mod
from fluidity_tools import stat_parser
import getopt
import sys
from math import ceil

def usage():
        print 'Usage:'
        print 'time_integrate_l2error.py --file=stat_filename\n'


################# Main ###########################
def main(argv=None):

        filename=''
        try:                                
                opts, args = getopt.getopt(sys.argv[1:], "", ['file=',''])
        except getopt.GetoptError:  
                usage()                     
                sys.exit(2)                     
        for opt, arg in opts:                
                if opt == '--file':      
                        filename=arg
                elif opt == '-h' or opt == '--help':
                    usage()                     
                    sys.exit(2) 
        if filename=='':
                print 'No stat filename specified. You have to give the detectors filename.'
                usage()   
                sys.exit(2) 

        
        ####################### Print time plot  ###########################

      
        # one wave lasts 1/((8*9.81*50/(430620*430620))^0.5)*2*pi=43193 seconds. 
        
        get_time_integrated_fs_error(filename)

       
        
if __name__ == "__main__":
    main()


