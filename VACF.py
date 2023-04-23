import numpy as np
import matplotlib.pyplot as plt
import re
from operator import itemgetter

#file location, coord.txt contains
#           'some nonsense strings we want to get rid of'
#           ATOMS_id type x y z vx vy vz c_px c_py c_pz
#           data
file_loc='/home/zhang/lammps-rel/examples/USER/rel/magneto-statics/VACF'
save_loc='/home/zhang/lammps-rel/examples/USER/rel/magneto-statics/VACF'


#calculating Velocity auto-correlation function
def VACF(File):
    #first part:extracting velocities from coord.txt (get rid of all annoying strings, x y z an other stuff
    data=[]
    with open(File, 'r') as f:
        lines=f.readlines()[int(0):]
        
    nop=int(float(lines[3])) #number_of_particles
    dump_interval=float(lines[9+nop+1]) #dump_interval
    
    for row in lines:
            if 'ITEM' not in row:
                new_row=re.split(r'\s', row)
            
                if len(new_row)>=9:
                    new_row2=list(filter(None, new_row))
                    new_row3=[float(n) for n in new_row2]
                    data.append(new_row3)
    datalist=data
    time_length=int(len(datalist)/nop)
    datalist1=np.reshape(datalist, (time_length,nop,11))
    
    datalist_good=[]
    for sublist in datalist1:
        sublist_sorted=sorted(sublist,key=itemgetter(0))
        datalist_good.append(sublist_sorted)
    
    #calculate vacf
    VACF_matrix=[]
    
    #for d in range (0, len(datalist_good), 1): #d is laptime
    for d in range (0, len(datalist_good), 1): #d is laptime
        delta=d
        time_tot=len(datalist_good)-delta
        v=datalist_good
            
        vj=0
        vi=0
        for j in range(0, nop, 1): #j is particle index
                
            for i in range(0, time_tot, 1 ):  #i is timestep index
                    
                vij=v[i][j]
                vijd=v[i+delta][j]
                VCi=vij[5]*vijd[5]+vij[6]*vijd[6]+vij[7]*vijd[7]
                vi=vi+VCi
                    #print(VCi)
                i=i+1
                        
                          
            VCj=vi/time_tot
                
            vj=vj+VCj
            j=j+1
            
        NVCF_del=vj/(3*nop)
            
            
        VACF_matrix.append([delta*dump_interval, NVCF_del])
            
        delta=delta+1
        
    VACF_matrix=np.array(VACF_matrix)
        
    return VACF_matrix

def VACF_c(File):
    #first part:extracting velocities from coord.txt (get rid of all annoying strings, x y z an other stuff
    data=[]
    with open(File, 'r') as f:
        lines=f.readlines()[int(0):]
        
    nop=int(float(lines[3])) #number_of_particles
    dump_interval=float(lines[9+nop+1]) #dump_interval
    
    for row in lines:
            if 'ITEM' not in row:
                new_row=re.split(r'\s', row)
            
                if len(new_row)>=9:
                    new_row2=list(filter(None, new_row))
                    new_row3=[float(n) for n in new_row2]
                    data.append(new_row3)
    datalist=data
    time_length=int(len(datalist)/nop)
    datalist1=np.reshape(datalist, (time_length,nop,8))
    
    datalist_good=[]
    for sublist in datalist1:
        sublist_sorted=sorted(sublist,key=itemgetter(0))
        datalist_good.append(sublist_sorted)
    
    #calculate vacf
    VACF_matrix=[]
    for d in range (0, len(datalist_good), 1): #d is laptime
        delta=d
        time_tot=len(datalist_good)-delta
        v=datalist_good
            
        vj=0
        vi=0
        for j in range(0, nop, 1): #j is particle index
                
            for i in range(0, time_tot, 1 ):  #i is timestep index
                    
                vij=v[i][j]
                vijd=v[i+delta][j]
                VCi=vij[5]*vijd[5]+vij[6]*vijd[6]+vij[7]*vijd[7]
                vi=vi+VCi
                    #print(VCi)
                i=i+1
                        
                          
            VCj=vi/time_tot
                
            vj=vj+VCj
            j=j+1
            
        NVCF_del=vj/(3*nop)
            
            
        VACF_matrix.append([delta*dump_interval, NVCF_del])
            
        delta=delta+1
        
    VACF_matrix=np.array(VACF_matrix)
        
    return VACF_matrix

#file_name='e_6new_3_mag'
#coord_file='coord_{}.txt'.format(file_name)
#file='{}/{}'.format(file_loc, coord_file)
#VACF_data=VACF(file)
#
##normalisation
#VACF_data_norm=VACF_data[:,1]/VACF_data[0][1]

#plt.figure()
#x=np.linspace(0, len(VACF_data)-1, len(VACF_data))
#y=np.zeros(len(VACF_data))
#plt.plot(x, y, '--')
#plt.plot(x, VACF_data)
#plt.title('what_ever_you_want_to_name_your_plot')
#plt.xlabel('time/(number of laps)')
#plt.ylabel('normalised VACF')
#
#
#plt.savefig('VACF_{}.jpg'.format(file_name), format = 'jpg',bbox_inches='tight', dpi=400)
#
##np.savetxt('{}/VACF_{}.txt'.format(save_loc, file_name), VACF_data)
#
#file_name='e_6new_6'
#coord_file='coord_{}.txt'.format(file_name)
#file='{}/{}'.format(file_loc, coord_file)
#VACF_data=VACF(file)
#np.savetxt('{}/VACF_{}.txt'.format(save_loc, file_name), VACF_data)
#
#file_name='e_6new_6_mag'
#coord_file='coord_{}.txt'.format(file_name)
#file='{}/{}'.format(file_loc, coord_file)
#VACF_data=VACF(file)
#np.savetxt('{}/VACF_{}.txt'.format(save_loc, file_name), VACF_data)
#
#
file_name='e_6new_6_1e7'
coord_file='coord_{}.txt'.format(file_name)
file='{}/{}'.format(file_loc, coord_file)
VACF_data=VACF(file)
np.savetxt('{}/VACF_{}.txt'.format(save_loc, file_name), VACF_data)
print('{} finished haha'.format(file_name))

file_name='e_6new_6_mag_1e7'
coord_file='coord_{}.txt'.format(file_name)
file='{}/{}'.format(file_loc, coord_file)
VACF_data=VACF(file)
np.savetxt('{}/VACF_{}.txt'.format(save_loc, file_name), VACF_data)
print('{} finished haha'.format(file_name))


file_name='e_6new_3_mag_1e7'
coord_file='coord_{}.txt'.format(file_name)
file='{}/{}'.format(file_loc, coord_file)
VACF_data=VACF(file)
np.savetxt('{}/VACF_{}.txt'.format(save_loc, file_name), VACF_data)
print('{} finished haha'.format(file_name))
