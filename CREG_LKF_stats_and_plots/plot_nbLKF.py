import os,sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#----- INPUT -----
#label1='$P_{ref}$'
#label2='$P_{ref}$*1.25'
#label3='$P_{ref}$*1.5'
#PATH_OUT='/home/jfl001/data/Lemieux2022/LKF_diag/Histo_width'
main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
exp1='eg1p0_ef1p5'
exp2='eg1p5_ef1p5'
exp3='eg2p25_ef1p5'
suffix='20050101_20050531.npy'
#-----------------

file1=os.path.join(main_dir+'/'+exp1+'/nbLKFs/number_lkf_'+suffix)
file2=os.path.join(main_dir+'/'+exp2+'/nbLKFs/number_lkf_'+suffix)
file3=os.path.join(main_dir+'/'+exp3+'/nbLKFs/number_lkf_'+suffix)

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
df3 = pd.read_csv(file3)

print('mean nb1:', df1['nb_of_LKFS'].mean())
print('mean nb2:', df2['nb_of_LKFS'].mean())
print('mean nb3:', df3['nb_of_LKFS'].mean())

plt.figure(1)
plt.plot(df1['nb_of_LKFS'])
plt.plot(df2['nb_of_LKFS'])
plt.plot(df3['nb_of_LKFS'])
#plt.plot(x, y2, marker='o')
#plt.plot(x, y3, marker='o')
#plt.ylim(4, 5)
#plt.xlim(1.1, 2.1)
#plt.legend([label1, label2, label3], loc ="lower right")
#plt.ylabel('Mean LKF width [nb of grid cells]')
#plt.xlabel('e', fontsize=14)
#plt.xticks(x, mylabels)
#plt.locator_params(axis='x', nbins=7)
#plt.xticks(rotation=45)
#plt.savefig(fileout)
plt.show()
