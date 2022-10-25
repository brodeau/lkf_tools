import os,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#----- INPUT -----
model='creg025' # creg12 ou creg025
#col1='b'
#col2='r'
#col3='c'
label1='$P_{ref}$'
label2='$P_{ref}$*1.25'
label3='$P_{ref}$*1.5'
PATH_OUT='/home/jfl001/data/Lemieux2022/LKF_diag/Histo_width'
#-----------------

fileout=PATH_OUT + '/mean_width_Feb2005_'+model+'.png'

y1=np.zeros(4) # Pstar
y2=np.zeros(4) # 1.25Pstar
y3=np.zeros(4) # 1.5Pstar
x=np.zeros(4)

x[0]=1.25
x[1]=1.5
x[2]=1.75
x[3]=2.0

y1[0]=4.34
y1[1]=4.59
y1[2]=4.89
y1[3]=4.98

y2[0]=4.21
y2[1]=4.55
y2[2]=4.87
y2[3]=4.94

y3[0]=4.21
y3[1]=4.45
y3[2]=4.76
y3[3]=4.91

#foutKEa=os.path.join(PATH_OUT + 'tsKEam_' +model +'_' +EXP1s +'_' +EXP2s +'_' +EXP3s +'_' +EXP4s +'.png')
plt.figure(1)
plt.plot(x, y1, marker='o')
plt.plot(x, y2, marker='o')
plt.plot(x, y3, marker='o')
plt.ylim(4, 5)
plt.xlim(1.1, 2.1)
plt.legend([label1, label2, label3], loc ="lower right")
plt.ylabel('Mean LKF width [nb of grid cells]')
plt.xlabel('e', fontsize=14)
#plt.xticks(x, mylabels)
#plt.locator_params(axis='x', nbins=7)
#plt.xticks(rotation=45)
plt.savefig(fileout)
