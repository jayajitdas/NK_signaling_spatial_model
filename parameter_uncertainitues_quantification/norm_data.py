#!/usr/bin/env python3

import numpy as np
import math

def calc_para_cluster(f,prob,ran):
    #bact_data = np.zeros((Nx,Ny,Nz))
    
    i = -1
    x = []
    x1 = []
    x2 = []
    x3 = []
    x4 = []
    x5 = []
    x6 = []
    x7 = []
    x8 = []
    x9 = []
    x10 = []
    x11 = []
    x12 = []
    x13 = []
    x14 = []
    x15 = []
    x16 = []
    x17 = []
    x18 = []
    x19 = []
    x20 = []
    x_max1 = 0.0
    x_max2 = 0.0
    x_max3 = 0.0
    x_max4 = 0.0
    x_max5 = 0.0
    x_max6 = 0.0
    x_max7 = 0.0
    x_max8 = 0.0
    x_max9 = 0.0
    x_max10 = 0.0
    x_max11 = 0.0
    x_max12 = 0.0
    x_max13 = 0.0
    x_max14 = 0.0
    x_max15 = 0.0
    x_max16 = 0.0
    x_max17 = 0.0
    x_max18 = 0.0
    x_max19 = 0.0
    for lines in f:
        i = i+1
        #print((f[i].strip('\n').split(' ')[1]))
        x1.append(float(f[i].strip('\n').split('\t')[0]))
        #print i #,file_line
        x2.append(float(f[i].strip('\n').split('\t')[1]))
        x3.append(float(f[i].strip('\n').split('\t')[2]))
        x4.append(float(f[i].strip('\n').split('\t')[3]))
        x5.append(float(f[i].strip('\n').split('\t')[4]))
        x6.append(float(f[i].strip('\n').split('\t')[5]))
        x7.append(float(f[i].strip('\n').split('\t')[6]))
        x8.append(float(f[i].strip('\n').split('\t')[7]))
        x9.append(float(f[i].strip('\n').split('\t')[8]))
        x10.append(float(f[i].strip('\n').split('\t')[9]))
        x11.append(float(f[i].strip('\n').split('\t')[10]))
        x12.append(float(f[i].strip('\n').split('\t')[11]))
        x13.append(float(f[i].strip('\n').split('\t')[12]))
        x14.append(float(f[i].strip('\n').split('\t')[13]))
        x15.append(float(f[i].strip('\n').split('\t')[14]))
        x16.append(float(f[i].strip('\n').split('\t')[15]))
        x17.append(float(f[i].strip('\n').split('\t')[16]))
        x18.append(float(f[i].strip('\n').split('\t')[17]))
        #print(i,x1[i],x2[i],x3[i],x4[i],x5[i],x6[i],x7[i],x8[i],x9[i],x10[i],x11[i],x12[i],x13[i],x14[i],x15[i],x16[i],x17[i],x18[i],x19[i]) #,file_line

    # making the data points dimensionless
    x_min1= min(x1) 
    x_min2= min(x2) 
    x_min3= min(x3) 
    x_min4= min(x4) 
    x_min5= min(x5) 
    x_min6= min(x6) 
    x_min7= min(x7) 
    x_min8= min(x8)
    x_min9= min(x9) 
    x_min10= min(x10) 
    x_min11= min(x11) 
    x_min12= min(x12)
    x_min13= min(x13) 
    x_min14= min(x14)
    x_min15= min(x15) 
    x_min16= min(x16)
    x_min17= min(x17) 
    
   
    x_max1= max(x1)  # rate kon
    x_max2= max(x2)# beta - local clustering pvav dependent factor upper limit =0
    x_max3= max(x3) # AR-AL +SFk kon
    x_max4= max(x4) # AR-AL +SFk koff
    x_max5= max(x5)  # ARp kcat
    x_max6= max(x6) # ARp+Vav kon
    x_max7= max(x7)  # ARp+Vav koff
    x_max8= max(x8)  # ARp-Vav +SFK kon
    x_max9= max(x9) # ARp-Vav +SFK koff
    x_max10= max(x10) # ARp-Vavp kcat
    x_max11= max(x11) # ARp-Vav dp kcat
    x_max12= max(x12) # AR dp kcat
    x_max13= max(x13)  # number of activating ligands
    x_max14= max(x14) # periphery centripatal force component upper limit =0
    x_max15= max(x15)  # AR-AL influx
    x_max16= max(x16) # k - pvav dependent centripetal force component

    x1_1 = []
    x1_2 = []
    x1_3 = []
    x1_4 = []
    x1_5 = []
    x1_6 = []
    x1_7 = []
    x1_8 = []
    x1_9 = []
    x1_10 = []
    x1_11 = []
    x1_12 = []
    x1_13 = []
    x1_14 = []
    x1_15 = []
    x1_16 = []
    x1_17 = []
    all_dis = []
    
    fo = open("0.29cut_off_param_pso_norm.txt","w")
    for i in range(len(x1)):
        x1_1.append((x1[i] - x_min1)/ ((x_max1 - x_min1)*1.0) )
        x1_2.append((x2[i] - x_min2)/ ((x_max2 - x_min2)*1.0) )
        x1_3.append((x3[i] - x_min3)/ ((x_max3 - x_min3)*1.0) )
        x1_4.append((x4[i] - x_min4)/ ((x_max4 - x_min4)*1.0) )
        x1_5.append((x5[i] - x_min5)/ ((x_max5 - x_min5)*1.0) )
        x1_6.append((x6[i] - x_min6)/ ((x_max6 - x_min6)*1.0) )
        x1_7.append((x7[i] - x_min7)/ ((x_max7 - x_min7)*1.0) )
        x1_8.append((x8[i] - x_min8)/ ((x_max8 - x_min8)*1.0) )
        x1_9.append((x9[i] - x_min9)/ ((x_max9 - x_min9)*1.0)) 
        x1_10.append((x10[i] - x_min10)/ ((x_max10 - x_min10)*1.0) )
        x1_11.append((x11[i] - x_min11)/ ((x_max11 - x_min11)*1.0) )
        x1_12.append((x12[i] - x_min12)/ ((x_max12 - x_min12)*1.0) )
        x1_13.append((x13[i] - x_min13)/ ((x_max13 - x_min13)*1.0) )
        x1_14.append((x14[i] - x_min14)/ ((x_max14 - x_min14)*1.0) )
        x1_15.append((x15[i] - x_min15)/ ((x_max15 - x_min15)*1.0) )
        x1_16.append((x16[i] - x_min16)/ ((x_max16 - x_min16)*1.0) )
        x1_17.append((x17[i] - x_min17)/ ((x_max17 - x_min17)*1.0) )
        fo.write("%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\n"%(x1_1[i],x1_2[i],x1_3[i],x1_4[i],x1_5[i],x1_6[i],x1_7[i],x1_8[i],x1_9[i],x1_10[i],x1_11[i],x1_12[i],x1_13[i],x1_14[i],x1_15[i],x1_16[i],x1_17[i],x18[i]))

    fo.close()

    return None

prob1 = 0.02 #2%

f6 = open("0.29cut_off_param_pso.txt").readlines()
ran = 0.29
calc_para_cluster(f6, prob1, 0.29)
