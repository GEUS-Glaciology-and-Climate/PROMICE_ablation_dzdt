#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:34:22 2021

code to 1.) ID ablation time series outliers and 2.) (optional cell) output monthly dzdt

@author: jeb


"""

import datetime
import time
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from numpy.polynomial.polynomial import polyfit

os.chdir('/Users/jason/Dropbox/AWS/PROMICE/PROMICE_ablation_dzdt/')

do_plot=1

fs=22 # font size
th=1
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.5
plt.rcParams['grid.color'] = "grey"
plt.rcParams["font.size"] = fs

PROMICE_stations = [
    #('EGP',(75.6247,-35.9748), 2660),  #OK
                    # ('KAN_B',(67.1252,-50.1832), 350), 
                    ('KAN_L',(67.0955,-35.9748), 670), #Height sensor boom unusable, Height stake do not capture winter accumulation
                    ('KAN_M',(67.0670,-48.8355), 1270), # minor adjustment left
                    ('KAN_U',(67.0003,-47.0253), 1840), # OK
                    ('KPC_L',(79.9108,-24.0828), 370),
                    # ('KPC_U',(79.8347,-25.1662), 870), # pressure transducer not working
                    ('MIT',(65.6922,-37.8280), 440), # ok minor adjustment left
                    ('NUK_K',(64.1623,-51.3587), 710), 
                    ('NUK_L',(64.4822,-49.5358), 530),
                    ('NUK_U',(64.5108,-49.2692), 1120),
                    ('QAS_L',(61.0308,-46.8493), 280),
                    ('QAS_M',(61.0998,-46.8330), 630), 
                    ('QAS_U',(61.1753,-46.8195), 900), 
                    ('SCO_L',(72.2230,-26.8182), 460),
                    ('SCO_U',(72.3933,-27.2333), 970),
                    ('TAS_A',(65.7790,-38.8995), 890),
                    ('TAS_L',(65.6402,-38.8987), 250),
                    ('THU_L',(76.3998,-68.2665), 570),
                    ('THU_U',(76.4197,-68.1463), 760),
                    ('UPE_L',(72.8932,-54.2955), 220), 
                    ('UPE_U',(72.8878,-53.5783), 940)
                    ]

for st_index,ws in enumerate(PROMICE_stations):
    if st_index>=3:      
        site=ws[0]
        print('--------------------------------------  '+site)
        outliers_file='./outliers/'+site+'_outliers.csv'
        out=open(outliers_file,'w+')
        out.write('date,dzdt\n')
        
        fn="/Users/jason/0_dat/v03_L3/"+site+"_hour_v03_L3.txt"
        df=pd.read_csv(fn,sep='\t')
        # print(df.columns)

        nams=["DepthPressureTransducer_Cor_adj(m)"]
        
        df["date"]=pd.to_datetime(df.time)
        df['year'] = df['date'].dt.year
        df['month'] = df['date'].dt.month
        
        df['jt']=df.DayOfYear+df['HourOfDay(UTC)']/24.
        
        iyear=df['year'][0]
        fyear=2020
        n_years=fyear-iyear+1
        
        if n_years>0:
            
            plt.close()
            fig, ax = plt.subplots(figsize=(14,7))
            
            x = df.date
            for nam in nams:
                df[nam][df[nam]<-998]=np.nan
                y = df[nam]
                y2 = df[nam]*np.nan
                if do_plot:
                    plt.plot_date(x, y,'-',lw=th)

                for k in range (len(df)-1):
                    dzdt=y[k+1]-y[k]
                    if abs(dzdt)>0.15:
                        out.write(str(df["date"][k])+",{:.2f}".format(dzdt)+'\n')

                        print(df["date"][k],"{:.2f}".format(dzdt))
                        y2[df["date"]==df["date"][k]]=y[k]
                        y2[df["date"]==df["date"][k+1]]=y[k+1]
                        # plt.text(x[k],y[k],"{:.2f}".format(dzdt),fontsize=fs,rotation=45)

                if do_plot:
                    plt.plot_date(x, y,'-b',lw=th,label=nam)
                    # plt.plot_date(x, y2,'or',lw=th)              
                    plt.scatter(x, y2, s=120, facecolors='none', edgecolors='r',label='abs(dzdt)>0.15')


            if do_plot:
                ax.xaxis.set_tick_params(rotation=30, labelsize=fs)
                ax.set_ylabel('m')
                ax.set_title(site+' hourly L3 ')
                plt.legend()
                # plt.show()
                ly='p'
        
                figpath='./figures/'
                if ly == 'x':plt.show()
                
                if ly == 'p':
                    figname=figpath+site+'_hourly.png'
                    plt.savefig(figname, bbox_inches='tight', dpi=150)
            out.close()

            #%% 
            # compute monthly SMB
            dzdt=np.zeros(((2,n_years,12)))
            counts=np.zeros(((2,n_years,12)))
            years=np.zeros(n_years*12)
            dec_year=np.zeros(n_years*12)
            months=np.zeros(n_years*12)
            t1mean=np.zeros(n_years*12)
            t1count=np.zeros(n_years*12)
            
            cc=0
            for yy in range(n_years):
                for mm in range(12):
                    if mm>=0:
                        v=((df.year==yy+iyear)&(df.month==mm+1)&(np.isfinite(df[nams[0]])))
                        c=np.sum(v)
                        # print(yy+iyear,mm+1,c)
                        
                        # if c>710:
                        #     dzdt[0,yy,mm]=np.nanmean(df[nams[0]][v])
                        #     counts[0,yy,mm]=np.sum(v)
                        if c>710:
                            y=df[nams[0]][v].values
                            dz=y[c-1]-y[0]
                            x=df['jt'][v].values
                            dt=(x[c-1]-x[0])
                            b, m = polyfit(x,y, 1)
                            dz=m*dt
                            # plt.scatter(x,y)
                            if ((dz>-3)&(dz<=0.02)):
                                dzdt[0,yy,mm]=dz
                                counts[0,yy,mm]=np.sum(v)
                                print(yy+iyear,mm+1,c,dz,dt)
                        else:
                            dzdt[0,yy,mm]=np.nan
                            counts[0,yy,mm]=np.nan            
                        years[cc]=yy+iyear
                        months[cc]=mm+1
                        dec_year[cc]=yy+iyear+(mm/12)
                        t1mean[cc]=dzdt[0,yy,mm]
                        t1count[cc]=counts[0,yy,mm]
                        cc+=1
            
            # data=[[years],[dzdt],[stds]]
            df2 = pd.DataFrame(columns = ['year', 'month','dz','count']) 
            df2.index.name = 'index'
            df2["year"]=pd.Series(years)
            df2["month"]=pd.Series(months)
            df2["dz"]=pd.Series(t1mean)
            df2["count"]=pd.Series(t1count)
            ofile="./output/"+site
            df2.to_csv(ofile+'.csv')
            
            if do_plot:
                plt.close()
                fig, ax = plt.subplots(figsize=(14,7))
                plt.plot(dec_year,t1mean,'s-')
                plt.title(site)
                plt.ylabel('m per month ice equiv')
                if ly == 'x':plt.show()
                
                if ly == 'p':
                    figname=figpath+site+'_monthly.png'
                    plt.savefig(figname, bbox_inches='tight', dpi=150)
