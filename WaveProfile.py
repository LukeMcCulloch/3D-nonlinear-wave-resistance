## TLM routine to read and graph the wave profile
## 11-2013
import csv
import os.path
import numpy as np
import matplotlib.pyplot as plt

# Open a file
with open("out", "rw+") as fo:
    print "Name of the file: ", fo.name


    for i,l in enumerate(fo):
        pass
    totlines=i+1

fo.close()


with open("out", "rw+") as fo:
    print "Name of the file: ", fo.name

    flag=0

    count = 0

    profile = []
    for line in fo:

        count =count+1
        if  line.startswith(' now'):
            print "Read Line: %s" % (line)
            savel = line
            flag=1
            size = totlines-count-1
            profile = np.zeros((2,size),float)
            loc=0
        elif flag==1 and count<totlines:
            
            thing=line.split()
            #print thing
            profile[0,loc]=float(thing[0])
            profile[1,loc]=float(thing[1])
            loc+=1

        
        
            
            
        
    #line = fo.readlines(2)
    #print "Read Line: %s" % (line)

    # Close opend file
fo.close()

plt.plot(profile[0,:],profile[1,:])
plt.show()
