import os,sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#----- INPUT -----
#PATH_OUT='/home/jfl001/data/Lemieux2022/LKF_diag/Histo_width'
main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
EXP1='run_eg1p16_ef1p75'
label1='e_f=1.75, e_g=1.16'
EXP2='run_eg1p75_ef1p75'
label2='e_f=1.75, e_g=1.75'
EXP3='run_eg2p63_ef1p75'
label3='e_f=1.75, e_g=2.63'

suffix='20050101_20050531.npy'
#-----------------

file1=os.path.join(main_dir+'/'+EXP1+'/nbLKFs/number_lkf_'+suffix)
file2=os.path.join(main_dir+'/'+EXP2+'/nbLKFs/number_lkf_'+suffix)
file3=os.path.join(main_dir+'/'+EXP3+'/nbLKFs/number_lkf_'+suffix)

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
df3 = pd.read_csv(file3)

print('mean nb1:', df1['nb_of_LKFS'].mean())
print('mean nb2:', df2['nb_of_LKFS'].mean())
print('mean nb3:', df3['nb_of_LKFS'].mean())

ax = df1.plot(x = 'date', y = 'nb_of_LKFS', color = "dodgerblue")
ax.set_prop_cycle(None)
df2.plot(ax=ax, x = 'date', y = 'nb_of_LKFS', color = "orange")
df3.plot(ax=ax, x = 'date', y = 'nb_of_LKFS', color = "darkviolet")
ax.legend([label1, label2, label3])
ax.set_xlabel("", fontsize='2')
ax.set_ylabel("Nb of LKFs")

plt.savefig('FIGS/Nb_LKFs_ef1p75_2005.png')
#plt.show()
