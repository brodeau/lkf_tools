import os,sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#----- INPUT -----
#PATH_OUT='/home/jfl001/data/Lemieux2022/LKF_diag/Histo_width'
main_dir='/home/jfl001/data/LKF_rips_analysis'
exp1='exp1'
exp2='control'
#exp3='eg2p25_ef1p5'
suffix='20200101_20201230.npy'
#-----------------

file1=os.path.join(main_dir+'/'+exp1+'/nbLKFs/number_lkf_'+suffix)
file2=os.path.join(main_dir+'/'+exp2+'/nbLKFs/number_lkf_'+suffix)
#file3=os.path.join(main_dir+'/'+exp3+'/nbLKFs/number_lkf_'+suffix)

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
print(df1)
#df3 = pd.read_csv(file3)

print('mean nb1:', df1['nb_of_LKFS'].mean())
print('mean nb2:', df2['nb_of_LKFS'].mean())
#print('mean nb3:', df3['nb_of_LKFS'].mean())

ax = df1.plot(x = 'date', y = 'nb_of_LKFS', color = "dodgerblue")
ax.set_prop_cycle(None)
df2.plot(ax=ax, x = 'date', y = 'nb_of_LKFS', color = "orange")
#df3.plot(ax=ax, x = 'date', y = 'nb_of_LKFS', color = "darkviolet")
ax.legend([exp1, exp2])
ax.set_xlabel("", fontsize='2')
ax.set_ylabel("Nb of LKFs")

plt.savefig('Nb_LKFs_2020.png')
#plt.show()
