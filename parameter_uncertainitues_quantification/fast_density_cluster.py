#!/usr/bin/env python3

import numpy as np
import math

def calc_para_cluster(f,prob,ran):
    
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
        x1.append(float(f[i].strip('\n').split('\t')[0]))
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
    
    # Method-1 : divide the each column of variables with the max value set in the PSO range  of that column; 
    x_max1= 1.0 # rate kon
    x_max2= 1.0 # beta - local clustering pvav dependent factor upper limit =0
    x_max3= 1.0  # AR-AL +SFk kon
    x_max4= 1.0 # AR-AL +SFk koff
    x_max5= 1.0  # ARp kcat
    x_max6= 1.0 # ARp+Vav kon
    x_max7= 1.0 # ARp+Vav koff
    x_max8= 1.0  # ARp-Vav +SFK kon
    x_max9= 1.0 # ARp-Vav +SFK koff
    x_max10= 1.0  # ARp-Vavp kcat
    x_max11= 1.0  # ARp-Vav dp kcat
    x_max12= 1.0  # AR dp kcat
    x_max13= 1.0  # number of activating ligands
    x_max14= 1.0 # periphery centripatal force component upper limit =0
    x_max15= 1.0 # AR-AL influx
    x_max16= 1.0  # k - pvav dependent centripetal force component
    x_max17= 1.0  # periphery rate constant
    

    all_dis = []
    for i in range(len(x1)):
        for j in range(i+1,len(x1)):
                dis1 = ((x1[i]/x_max1) - (x1[j]/x_max1))**2 
                dis2 = ((x2[i]/x_max2) - (x2[j]/x_max2))**2 
                dis3 = ((x3[i]/x_max3) - (x3[j]/x_max3))**2 
                dis4 = ((x4[i]/x_max4) - (x4[j]/x_max4))**2 
                dis5 = ((x5[i]/x_max5) - (x5[j]/x_max5))**2 
                dis6 = ((x6[i]/x_max6) - (x6[j]/x_max6))**2 
                dis7 = ((x7[i]/x_max7) - (x7[j]/x_max7))**2 
                dis8 = ((x8[i]/x_max8) - (x8[j]/x_max8))**2 
                dis9 = ((x9[i]/x_max9) - (x9[j]/x_max9))**2 
                dis10 = ((x10[i]/x_max10) - (x10[j]/x_max10))**2 
                dis11 = ((x11[i]/x_max11) - (x11[j]/x_max11))**2 
                dis12 = ((x12[i]/x_max12) - (x12[j]/x_max12))**2 
                dis13 = ((x13[i]/x_max13) - (x13[j]/x_max13))**2 
                dis14 = ((x14[i]/x_max14) - (x14[j]/x_max14))**2 
                dis15 = ((x15[i]/x_max15) - (x15[j]/x_max15))**2 
                dis16 = ((x16[i]/x_max16) - (x16[j]/x_max16))**2 
                dis17 = ((x17[i]/x_max17) - (x17[j]/x_max17))**2 
                dis = round(np.sqrt(dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10+dis11+dis12+dis13+dis14+dis15+dis16+dis17),4) 
                all_dis.append(dis)
    
    # calulate d_c - distance cut-off
    sort_all_dis = sorted(all_dis)
    max_all_dis = max(all_dis)
    len_dis = (len(x1) * (len(x1)-1)*1.0)/2.0
    percentage_nbh = prob*len_dis #prob = % of the nbh, N C 2 = 1540, where N = total number of points in the data set; 1540 - N = 1484
    #dc = max_all_dis*percentage_nbh
    dc = percentage_nbh
    index_dis_cut_off = round(dc)-1
    print("\n","len_x1",len(x1),"len_sort_all_dis", len(sort_all_dis),"len_dis",len_dis,"percentage_nbh",percentage_nbh, "max_all_dis",max_all_dis)
    dis_cut_off = round(sort_all_dis[index_dis_cut_off],4)
    print("\n","dis_cut_off",dis_cut_off,"index_dis_cut_off",index_dis_cut_off,"dc",dc,"\n")
    
    fo=open(str(prob)+"_percen_"+str(dis_cut_off)+"_dc_"+str(ran)+"cluster_param_cost_m1_norm.txt","w")
    fo1=open(str(prob)+"_percen_"+str(dis_cut_off)+"_dc_"+str(ran)+"interval_param__m1_norm.txt","w")

    # calculate rho for each data point
    
    roh = []
    for i in range(len(x1)):
        roh.append(0)

    rho_dict = {}
    for i in range(len(x1)):
        for j in range(i+1, len(x1)):
            dis1 = ((x1[i]/x_max1) - (x1[j]/x_max1))**2 
            dis2 = ((x2[i]/x_max2) - (x2[j]/x_max2))**2 
            dis3 = ((x3[i]/x_max3) - (x3[j]/x_max3))**2 
            dis4 = ((x4[i]/x_max4) - (x4[j]/x_max4))**2 
            dis5 = ((x5[i]/x_max5) - (x5[j]/x_max5))**2 
            dis6 = ((x6[i]/x_max6) - (x6[j]/x_max6))**2 
            dis7 = ((x7[i]/x_max7) - (x7[j]/x_max7))**2 
            dis8 = ((x8[i]/x_max8) - (x8[j]/x_max8))**2 
            dis9 = ((x9[i]/x_max9) - (x9[j]/x_max9))**2 
            dis10 = ((x10[i]/x_max10) - (x10[j]/x_max10))**2 
            dis11 = ((x11[i]/x_max11) - (x11[j]/x_max11))**2 
            dis12 = ((x12[i]/x_max12) - (x12[j]/x_max12))**2 
            dis13 = ((x13[i]/x_max13) - (x13[j]/x_max13))**2 
            dis14 = ((x14[i]/x_max14) - (x14[j]/x_max14))**2 
            dis15 = ((x15[i]/x_max15) - (x15[j]/x_max15))**2 
            dis16 = ((x16[i]/x_max16) - (x16[j]/x_max16))**2 
            dis17 = ((x17[i]/x_max17) - (x17[j]/x_max17))**2 
            dis =round( np.sqrt(dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10+dis11+dis12+dis13+dis14+dis15+dis16+dis17),4) 
           
            if ( dis <dis_cut_off):
                roh[i] = roh[i] + 1.0
                roh[j] = roh[j] + 1.0

        rho_dict.update({i:roh[i]})
        rho_dict.update({j:roh[j]})
    roh_dict = {}
    for i in range(5):
        roh_dict.update({i:roh[i]})
        print(rho_dict[i],roh_dict[i], roh[i])
        
    # calculate delta
    delt = np.zeros(len(x1))
    delta_dict = {}
    sigma_dict = {}
    sigma_i = np.zeros(len(x1))
    max_roh = max(roh)
    print("max_roh",max_roh)
    for i in range(len(x1)):
        dis_del = []
        dis_del_dict = {}
        cal_sigma = []
        if (roh[i] != max_roh):
            for j in range(len(x1)):
                if ((i != j) and (roh[j] > roh[i])):
                    dis1 = ((x1[i]/x_max1) - (x1[j]/x_max1))**2 
                    dis2 = ((x2[i]/x_max2) - (x2[j]/x_max2))**2 
                    dis3 = ((x3[i]/x_max3) - (x3[j]/x_max3))**2 
                    dis4 = ((x4[i]/x_max4) - (x4[j]/x_max4))**2 
                    dis5 = ((x5[i]/x_max5) - (x5[j]/x_max5))**2 
                    dis6 = ((x6[i]/x_max6) - (x6[j]/x_max6))**2 
                    dis7 = ((x7[i]/x_max7) - (x7[j]/x_max7))**2 
                    dis8 = ((x8[i]/x_max8) - (x8[j]/x_max8))**2 
                    dis9 = ((x9[i]/x_max9) - (x9[j]/x_max9))**2 
                    dis10 = ((x10[i]/x_max10) - (x10[j]/x_max10))**2 
                    dis11 = ((x11[i]/x_max11) - (x11[j]/x_max11))**2 
                    dis12 = ((x12[i]/x_max12) - (x12[j]/x_max12))**2 
                    dis13 = ((x13[i]/x_max13) - (x13[j]/x_max13))**2 
                    dis14 = ((x14[i]/x_max14) - (x14[j]/x_max14))**2 
                    dis15 = ((x15[i]/x_max15) - (x15[j]/x_max15))**2 
                    dis16 = ((x16[i]/x_max16) - (x16[j]/x_max16))**2 
                    dis17 = ((x17[i]/x_max17) - (x17[j]/x_max17))**2 
                    dis = 0.0
                    dis = round(np.sqrt(dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10+dis11+dis12+dis13+dis14+dis15+dis16+dis17),4)
                    dis_del.append(round(np.sqrt(dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10+dis11+dis12+dis13+dis14+dis15+dis16+dis17),4))
                    dis_del_dict.update({j:dis})
                    cal_sigma.append(j)
                   
            delt[i] = min(dis_del)
            delta_dict.update({i:delt[i]})
           
            
        else:
            for j in range(len(x1)):
                    dis1 = ((x1[i]/x_max1) - (x1[j]/x_max1))**2 
                    dis2 = ((x2[i]/x_max2) - (x2[j]/x_max2))**2 
                    dis3 = ((x3[i]/x_max3) - (x3[j]/x_max3))**2 
                    dis4 = ((x4[i]/x_max4) - (x4[j]/x_max4))**2 
                    dis5 = ((x5[i]/x_max5) - (x5[j]/x_max5))**2 
                    dis6 = ((x6[i]/x_max6) - (x6[j]/x_max6))**2 
                    dis7 = ((x7[i]/x_max7) - (x7[j]/x_max7))**2 
                    dis8 = ((x8[i]/x_max8) - (x8[j]/x_max8))**2 
                    dis9 = ((x9[i]/x_max9) - (x9[j]/x_max9))**2 
                    dis10 = ((x10[i]/x_max10) - (x10[j]/x_max10))**2 
                    dis11 = ((x11[i]/x_max11) - (x11[j]/x_max11))**2 
                    dis12 = ((x12[i]/x_max12) - (x12[j]/x_max12))**2 
                    dis13 = ((x13[i]/x_max13) - (x13[j]/x_max13))**2 
                    dis14 = ((x14[i]/x_max14) - (x14[j]/x_max14))**2 
                    dis15 = ((x15[i]/x_max15) - (x15[j]/x_max15))**2 
                    dis16 = ((x16[i]/x_max16) - (x16[j]/x_max16))**2 
                    dis17 = ((x17[i]/x_max17) - (x17[j]/x_max17))**2 
                    dis = 0.0
                    dis = round(np.sqrt(dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10+dis11+dis12+dis13+dis14+dis15+dis16+dis17),4)
                    dis_del.append(round(np.sqrt(dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10+dis11+dis12+dis13+dis14+dis15+dis16+dis17),4))
            delt[i] = max(dis_del)
            delta_dict.update({i:delt[i]})
            

    gamma_i = []
    gamma_dict = {}
    for i in range(len(x1)):
        gamma_i.append(roh[i]*delt[i])
        gamma_dict.update({i:rho_dict[i]*delta_dict[i]})
    for i in range(len(gamma_i)):
        fo.write("%d %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %f %f %f %f\n"%(i,x1[i],x2[i],x3[i],x4[i],x5[i],x6[i],x7[i],x8[i],x9[i],x10[i],x11[i],x12[i],x13[i],x14[i],x15[i],x16[i],x17[i],roh[i],delt[i],gamma_i[i],x18[i]))

    fo.close()

    return None

prob1 = 0.02 #2%
f6 = open("29cut_off_param_pso_norm.txt").readlines()
ran = 0.29
calc_para_cluster(f6, prob1, 0.29)
