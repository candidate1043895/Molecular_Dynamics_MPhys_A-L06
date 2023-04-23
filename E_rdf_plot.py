import numpy as np
import matplotlib.pyplot as plt
import re

#file location
where_your_rdf_file_lives='/Users/candidate123456/MPhys_MD/lammps-rel/examples/USER/rel/magneto-statics/18_runs_e/ewald'
rdf_file='{}/whatever_your_rdf_file_name_is.txt'
file=rdf_file.format(where_your_rdf_file_lives)


#Depends on how many times you dumped your data
Nbegin=20
Nend=100
# your setting for curoff radius, this is used as the max plot range( it doesn't have to be r_cutoff, you can write whatever you are happy with it)
r_cut=1.2


#how your rdf file is generated in lammps:
#  compute myrdf all rdf 1000 1 1
#  fix 345 all ave/time 1 100 10000 c_myrdf[*] file E_rdf.txt mode vector

# function rdf_Nend collects data in E_rdf.txt from the 0th dump(the starting point) to the {Nend}_th dump.
#~ from the 0th timestep to the {Nend}_th timestep.
def rdf_Nend(testFile, Nend):
    #the first part reads E_rdf.txt file as a list
    data=[]
    with open(testFile, 'r') as f:
        lines=f.readlines()[3:]
        
        for row in lines:
            new_row=re.split(r'\s', row)
            if len(new_row)==5:
                new_row2=list(filter(None, new_row))
                new_row3=[float(n) for n in new_row2]
                data.append(new_row3)
    
    #the second part extracts data in E_rdf.txt from the 0th dump to the {Nend}_th dump.
    RUN =  [[0] * 4] * 1000
    
    for n in range(1, Nend+1, 1):
        run=np.array(data[1000*(n-1): 1000*n])
        RUN=RUN+run
        
    
    
    RUN=RUN/(Nend-0)
    
    return RUN

# function rdf_Nend collects data in E_rdf.txt from the {Nbegin}_th dump to the {Nend}_th dump.
#~ from the {Nend}_th timestep to the {Nend}_th timestep.
def rdf_Nbegin_Nend(testFile, Nbegin, Nend):
    #the first part reads E_rdf.txt file as a list
    data=[]
    with open(testFile, 'r') as f:
        lines=f.readlines()[3:]
        
        for row in lines:
            new_row=re.split(r'\s', row)
            if len(new_row)==5:
                new_row2=list(filter(None, new_row))
                new_row3=[float(n) for n in new_row2]
                data.append(new_row3)
    
    #the second part extracts data in E_rdf.txt from the {Nbegin}_th dump to the {Nend}_th dump.
    RUN =  [[0] * 4] * 1000
    
    for n in range(Nbegin+1, Nend+1, 1):
        run=np.array(data[1000*(n-1): 1000*n])
        RUN=RUN+run
        
    
    
    RUN=RUN/(Nend-Nbegin)
    
    #radial position is RUN[:, 1], rdf value is RUN[:, 2]
    return RUN


rdfs=rdf_Nbegin_Nend(file, Nbegin, Nend)


#now you can plot rdf
plt.plot(rdfs[:, 1],rdfs[:, 2], '.')
plt.xlim(0, r_cut)
plt.xlabel('Distance(r)/Angstrom')
plt.ylabel('g(r)')
plt.title('what_ever_you_want_to_name_your_plot')

#in case you want to save your beautiful plot
plt.savefig('what_ever_you_want_to_name_your_plot.jpg', format = 'jpg',bbox_inches='tight', dpi=400)

