import math
from typing import List
import numpy as np
from astropy.io import fits
import scipy.special as scisp
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import sunpy.coordinates
from datetime import datetime
import sunpy.util.net
import requests
from bs4 import BeautifulSoup
from urllib.request import urlopen
import os

from numpy import ndarray
## User defined
date = list()
# this is the time at which the magnetic field is needed.
date.append('2007-10-07T20:00:00') 
# date in UTC format
# These two days ensure the existence of two moments, one earlier than the one to be interpolated and one later than the one to be interpolated, in the reconstruction stencil.
date.append('2007-10-06T11:00:00') 
date.append('2007-10-08T11:00:00') 
map_type = 'ADAPT' #'GONG' #  'HMI' # 'WSO', 'GONG', 'ADAPT', 'HMI, 'KPVT', 'MDI', 'SOLIS', 'MWO'
# NB: for the moment, only GONG, ADAPT and HMI are operational
output_dir = './MapData'

if(map_type == 'ADAPT'):
  out_dir = './MapData'
elif (map_type == 'GONG'):
  out_dir = './test_maps/GONGzqs/'
lmax = 15 #20 #9 #8 #15 #25 #7 #  numbers of modes
shift_adapt0 = 320.0
adapt_map = 0 # for ADAPT maps, between 0 and 11
amp: float = 1.0 #0.75 # amplitude factor of the map
r_st = 1.0 # radius at which the magnetic field is computed
smooth_surfmap='no' # Whether to smooth the magnetograms or not
smooth_Num=4
flux_correct = 'no' # keep flux to zero
B_scale = 'no' # Limit |Br_max| under threshold
Brmax_const='no' # Whether to keep the max magnetic field strength of magnetograms const or not
candence=6 #1 #2 #3
hour_ave=int(24/candence) #12  # keep the max magnetic field strength of magnetograms const as daily average value
sinlat2lat = 'no' #'yes' #
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
corotation = 'yes' #'yes' # only for GONGzqs
FR_insert = "no" # ”yes“ for insertion of flux rope field
#########！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
write_map = 'yes' # 'yes' or 'no' for COCONUT code
show_map = 'yes' # 'yes' or 'no'
plot_relative = 'no'
visu_type = 'sinlat' # 'lat' or 'sinlat'
write_map_correct = 'yes'  # for SIP code
TransCorr = 'no'
Interp_order=2
alpha=3*10**(-6)
second_filter_factor=0.0
FRfactor=2.0

Brmin = -20.0
Brmax = 20.0
days_Feb=28
Frame_num = 31*hour_ave+hour_ave
ii_begin= 0 # default set to 0; for a restart, set to the last (ii-hour_ave)+hour_ave+1.


Bref_Scale=6.6/2.2
#Bref_Scale = 42.0 / 2.2
Bref_Scale0 = Bref_Scale
dailyBr_max = np.zeros(Frame_num)
## Name of the ouput file
year = list()
month = list()
day = list()
hour = list()
minute = list()
second = list()
for i in range(3):
  year.append(int(date[i].split('-')[0]))
  month.append(int(date[i].split('-')[1]))
  day.append(int(date[i].split('-')[2].split('T')[0]))
  hour.append(int(date[i].split('T')[1].split(':')[0]))
  minute.append(int(date[i].split(':')[1]))
  second.append(int(date[i].split(':')[2]))
days_tore: ndarray = np.zeros(31-day[0], int)


for ii in range(0,ii_begin):
  FR_factor=1.0/(FRfactor**(min(ii*candence,10)))
  print("FR_factor=",FR_factor)
  if(ii>0):
    hour[0]+=candence
  if(hour[0]>=24):
    day_add = hour[0] // 24
    hour[0]=hour[0]%24
    for jj in range(3):
      day[jj]+=day_add
      if(month[jj] in (1, 3, 5, 7, 8, 10, 12)):
        if(day[jj]>31):
          day[jj] -= 31
          month[jj]+=1
      elif(month[jj] == 2):
        if (day[jj] > days_Feb):
          day[jj] -= days_Feb
          month[jj]+=1
      else:
        if (day[jj] > 30):
          day[jj] -= 30
          month[jj] += 1
      if(month[jj]>12):
        month[jj]=1
        year[jj]+=1
        if(year[jj]%4==0):
          if(year[jj]%100==0):
            if(year[jj]%400==0):
              days_Feb=29
          else:
            days_Feb = 29

for ii in range(ii_begin,Frame_num):
  FR_factor=1.0/(FRfactor**(min(ii*candence,10)))
  print("FR_factor=",FR_factor)
  if(ii>0):
    hour[0]+=candence
  if(hour[0]>=24):
    day_add = hour[0] // 24
    hour[0]=hour[0]%24
    for jj in range(3):
      day[jj]+=day_add
      if(month[jj] in (1, 3, 5, 7, 8, 10, 12)):
        if(day[jj]>31):
          day[jj] -= 31
          month[jj]+=1
      elif(month[jj] == 2):
        if (day[jj] > days_Feb):
          day[jj] -= days_Feb
          month[jj]+=1
      else:
        if (day[jj] > 30):
          day[jj] -= 30
          month[jj] += 1
      if(month[jj]>12):
        month[jj]=1
        year[jj]+=1
        if(year[jj]%4==0):
          if(year[jj]%100==0):
            if(year[jj]%400==0):
              days_Feb=29
          else:
            days_Feb = 29


  if(ii>=hour_ave):
    Bref_Scale=0.0
    for jj in range(0,hour_ave):
      Bref_Scale+=dailyBr_max[ii-1-jj]
    Bref_Scale/=hour_ave
    Bref_Scale = 6.6 / 2.2

  date_datetime = datetime(year[0], month[0], day[0], hour[0], minute[0], second[0])
  print('The data is', date_datetime)
  cr_number = int(sunpy.coordinates.sun.carrington_rotation_number(date_datetime))
  print('crNo=', cr_number)
  #exit(0)
  data_out = str(year[0]).rjust(4,'0') + str(month[0]).rjust(2,'0') + str(day[0]).rjust(2,'0') +\
           str(hour[0]).rjust(2,'0') + str(minute[0]).rjust(2,'0') + str(second[0]).rjust(2,'0')
  data_outSix=str(month[0]) + str(day[0]).rjust(2,'0') + str(hour[0]).rjust(2,'0')
  # WSO
  if (map_type == 'WSO'):
    output_name = output_dir + 'map_wso_lmax' + str(lmax) + '_cr' + str(cr_number) + '.dat'
  # GONG
  # NB: only mrzqs for now
  elif (map_type == 'GONG'):
    output_name = out_dir + 'COCONUT/201906/' + 'map_gong_lmax' + str(lmax) + '_' + data_out + '.dat'
    output_Figre_name = out_dir + 'Figure/201906/' + 'map_gong_lmax' + str(lmax) + '_' + data_out + '.png'
  # ADAPT
  elif (map_type == 'ADAPT'):
    output_name = out_dir  + 'map_adapt_lmax' + str(lmax) + '_' + data_out + '.dat'
    output_Figre_name = out_dir  + 'map_adapt_lmax' + str(lmax) + '_' + data_out
  ##2023-10-12 data-->data_out
  # HMI
  # NB: only HMI carrington for now
  elif (map_type == 'HMI'):
    output_name = output_dir + 'map_hmi_lmax' + str(lmax) + '_cr' + str(cr_number) + '.dat'
  #elif (map_type == 'MDI'):
  #  output_name = output_dir + 'map_mdi_lmax' + str(lmax) + '_cr' + str(cr_number) + '.dat'
  # KPVT
  #output_name = '/Users/bperri/Documents/Data/Leuven/map_studies/cr1914/map_kpvt_cr1914_lmax30.dat'
  # MDI
  #output_name = '/Users/bperri/Documents/Data/Leuven/map_studies/cr2071/map_mdi_010808_lmax30.dat'
  # SOLIS
  #output_name = '/Users/bperri/Documents/Data/Leuven/map_studies/cr2193/map_solis_2193.dat'

  ## Name of the input magnetogram
  # WSO
  if (map_type == 'WSO'):
    map_name = 'wso_cr2194.txt'
  #map_name = 'wso_cr2192.txt'
  #map_name = 'wso_cr2219.txt'
  #map_name = 'wso_cr1902.txt'
  #map_name = 'wso_cr2072.txt'
  # GONG
  elif (map_type == 'GONG'):
    year_str=list()
    month_str=list()
    day_str=list()
    for i in range(3):
      year_str.append(str(year[i])[2:])
      if (month[i] < 10):
        month_str.append('0' + str(month[i]))
      else:
        month_str.append(str(month[i]))
      if (day[i] < 10):
        day_str.append('0' + str(day[i]))
      else:
        day_str.append(str(day[i]))
    file_id = 'mrzqs'
    file_id_str = file_id[2:]
    remote_dir = 'https://gong.nso.edu/data/magmap/QR/' + file_id_str + '/' + str(year[0]) + \
               month_str[0] + '/' + file_id + year_str[0] + month_str[0] + day_str[0] + '/'
    remote_dir_m = 'https://gong.nso.edu/data/magmap/QR/' + file_id_str + '/' + str(year[1]) + \
               month_str[1] + '/' + file_id + year_str[1] + month_str[1] + day_str[1] + '/'
    remote_dir_p = 'https://gong.nso.edu/data/magmap/QR/' + file_id_str + '/' + str(year[2]) + \
               month_str[2] + '/' + file_id + year_str[2] + month_str[2] + day_str[2] + '/'
    # Find closest maps
    page_text = requests.get(remote_dir).text
    # print('page_text=',page_text)
    soup = BeautifulSoup(page_text, "html.parser")
    file_names=[node.get("href") for node in soup.find_all("a") if file_id in node.get("href")]
    n0=len(file_names)
    #print('n0=',n0)
    page_text1 = requests.get(remote_dir_m).text
    soup = BeautifulSoup(page_text1, "html.parser")
    file_names_m=[node.get("href") for node in soup.find_all("a") if file_id in node.get("href")]
    n1=len(file_names_m)
    page_text2 = requests.get(remote_dir_p).text
    soup = BeautifulSoup(page_text2, "html.parser")
    file_names_p = [node.get("href") for node in soup.find_all("a") if file_id in node.get("href")]
    n2=len(file_names_p)
    file_namefile =  file_names_m + file_names + file_names_p
    print('file_namefile=',file_namefile)
    #time_deltas = list()
    file_date = list()
    m=0
    for file_name in file_namefile:
      file_date.append(datetime.strptime(file_name.split("c")[0], file_id + "%y%m%dt%H%M"))
      m=m+1
    print('m=',m)
    for i in range(1, m-2):
      temp1=(date_datetime-file_date[i]).total_seconds()
      temp2=(file_date[i+1]-date_datetime).total_seconds()
      if (temp1*temp2>=0):
        coef1 = temp2 / (temp1 + temp2)
        coef2 = temp1 / (temp1 + temp2)
        temp=temp1+temp2
        tempm = (file_date[i] - file_date[i-1]).total_seconds()
        tempp = (file_date[i+2] - file_date[i + 1]).total_seconds()
        timem=(file_date[i - 1]-file_date[i - 1]).total_seconds()/3600.0
        time=(file_date[i]-file_date[i - 1]).total_seconds()/3600.0
        time1 = (file_date[i+1]-file_date[i - 1]).total_seconds()/3600.0
        timep = (file_date[i+2]-file_date[i - 1]).total_seconds()/3600.0
        time_interp=(date_datetime-file_date[i-1]).total_seconds()/3600.0
        print('temp1=',temp1,'temp2=',temp2,'tempm=',tempm,'tempp=',tempp)
        break
      # print('file_date=', file_date)
    if(temp<=0):
      print('Error! Check the data!')
      exit(0)
      #time_deltas.append((file_date - date_datetime).total_seconds())
    #map_name = file_names[abs(np.array(time_deltas)).argmin()]
    print('i=',i,'timem=',timem,'time=',time,'time1=',time1,'timep=',timep,'time_interp=',time_interp)
    print('n1=', n1, 'n0=', n0, 'n2=', n2)
    map_name = file_namefile[i]
    map_name1 = file_namefile[i+1]
    map_namem = file_namefile[i-1]
    map_namep = file_namefile[i + 2]
    if(i<n1):
      remote_filem = remote_dir_m + map_namem
      remote_file = remote_dir_m + map_name
      if((i+1)<n1):
        remote_file1 = remote_dir_m + map_name1
        if((i+2)<n1):
          remote_filep = remote_dir_m + map_namep
        elif((i+2)>=n1 and (i+2)<(n1+n0)):
          remote_filep = remote_dir + map_namep
        else:
          remote_filep = remote_dir_p + map_namep
      elif((i+1)>=n1 and (i+1)<(n1+n0)):
        remote_file1 = remote_dir + map_name1
        if ((i + 2) >= n1 and (i + 2) < (n1 + n0)):
          remote_filep = remote_dir + map_namep
        else:
          remote_filep = remote_dir_p + map_namep
      else:
        remote_file1 = remote_dir_p + map_name1
        remote_filep = remote_dir_p + map_namep
    elif(i>=n1 and i<(n1+n0)):
      remote_file = remote_dir + map_name
      if((i-1)<n1):
        remote_filem = remote_dir_m + map_namem
      else:
        remote_filem = remote_dir + map_namem
      if ((i + 1) >= n1 and (i + 1) < (n1 + n0)):
        remote_file1 = remote_dir + map_name1
        if ((i + 2) >= n1 and (i + 2) < (n1 + n0)):
          remote_filep = remote_dir + map_namep
        else:
          remote_filep = remote_dir_p + map_namep
      else:
        remote_file1 = remote_dir_p + map_name1
        remote_filep = remote_dir_p + map_namep
    else:
      remote_file = remote_dir_p + map_name
      remote_file1 = remote_dir_p + map_name1
      remote_filep = remote_dir_p + map_namep
      if ((i - 1) < n1):
        remote_filem = remote_dir_m + map_namem
      elif((i - 1) >= n1 and (i - 1) < (n1 + n0)):
        remote_filem = remote_dir + map_namem
      else:
        remote_filem = remote_dir_p + map_namem

  # ADAPT
  elif (map_type == 'ADAPT'):
    remote_dir = 'https://gong.nso.edu/adapt/maps/gong/' + str(year[0]) + '/'
    page_text = requests.get(remote_dir).text
    soup = BeautifulSoup(page_text, "html.parser")
    file_id = 'adapt40311'
    file_namefile = [node.get("href") for node in soup.find_all("a") if file_id in node.get("href")]
    #time_deltas = list()
    file_date = list()
    m=0
    for file_name in file_namefile:
      file_date.append(datetime.strptime(file_name.split("_")[2], "%Y%m%d%H%M"))
      m=m+1
    print('m=', m)
    for i in range(m):
      temp1 = (date_datetime - file_date[i]).total_seconds()
      temp2 = (file_date[i + 1] - date_datetime).total_seconds()
      if (temp1 * temp2 >= 0):
        coef1 = temp2 / (temp1 + temp2)
        coef2 = temp1 / (temp1 + temp2)
        temp = temp1 + temp2
        tempm = (file_date[i] - file_date[i - 1]).total_seconds()
        tempp = (file_date[i + 2] - file_date[i + 1]).total_seconds()
        print('temp1=', temp1, 'temp2=', temp2)
        break
      # print('file_date=', file_date)
    if (temp1 + temp2 <= 0):
      print('Error! Check the data!')
      exit(0)
    print('i=', i)
    map_name = file_namefile[i]
    map_name1 = file_namefile[i + 1]
    map_namem = file_namefile[i-1]
    map_namep = file_namefile[i + 2]
    remote_file = remote_dir + map_name
    remote_file1 = remote_dir + map_name1
    remote_filem = remote_dir + map_namem
    remote_filep = remote_dir + map_namep

      #time_deltas.append((file_date - date_datetime).total_seconds())
    #map_name = file_names[abs(np.array(time_deltas)).argmin()]
    #remote_file = remote_dir + map_name
  # HMI
  elif (map_type == 'HMI'):
    map_name = 'hmi.Synoptic_Mr_small.' + str(cr_number) + '.fits'
    remote_file = 'http://jsoc.stanford.edu/data/hmi/synoptic/' + map_name
  # Download file
  file0_exist = False
  file1_exist = False
  filem_exist = False
  filep_exist = False
  path='./test_maps/'
  # 20231204 #
  for files in os.walk(path):
    for filename in files:
      if map_namem in filename:
        print('filem already exist')
        filem_exist = True
        print('map_namem=', map_namem)
      if map_name in filename:
        print('file0 already exist')
        file0_exist = True
        print('map_name=',map_name)
      if map_name1 in filename:
        print('file1 already exist')
        file1_exist = True
        print('map_name1=',map_name1)
      if map_namep in filename:
        print('filep already exist')
        filep_exist = True
        print('map_namep=', map_namep)
  if (filem_exist == False):
    local_filem = sunpy.util.net.download_file(remote_filem, directory = output_dir, overwrite = True)
  if (file0_exist == False):
    local_file = sunpy.util.net.download_file(remote_file, directory = output_dir, overwrite = True)
  if (file1_exist == False):
    local_file1 = sunpy.util.net.download_file(remote_file1, directory = output_dir, overwrite = True)
  if (filep_exist == False):
    local_filep = sunpy.util.net.download_file(remote_filep, directory = output_dir, overwrite = True)

  # MDI
  #map_name = 'synop_Mr_0.2047.fits'
  #map_name = 'synop_Mr_0.2071.fits'
  # SOLIS
  #map_name = 'kbv7g170802t2030c2193_000_int-mas_dim-180.fits'
  #map_name = 'kbv7g170802t2030c2193_000_int-mas_dim-900.fits'
  # KPVT
  #map_name = 'kbv7g060907t1457c2047_000_int-mas_dim-180.fits'
  #map_name = 'm1914f.fits'

  # Opening fits
  #Br_max= 3.729889179298422 Scale= 1.0
  print('Reading file')
  input_file = output_dir + map_name
  input_file1 = output_dir + map_name1
  input_filem = output_dir + map_namem
  input_filep = output_dir + map_namep
  nb_modes_tot = int((lmax+1)*(lmax+2)/2 - 1)
  # ADAPT
  if (map_type == 'ADAPT'):
    input_data = fits.getdata(input_file, ext=0)
    input_data1 = fits.getdata(input_file1, ext=0)
    input_datam = fits.getdata(input_filem, ext=0)
    input_datap = fits.getdata(input_filep, ext=0)
    shape = np.shape(input_data)
    nb_maps = shape[0]
    nb_th = shape[1]
    nb_phi = shape[2]
    Br_data = input_data[adapt_map,::-1,:]
    Br_data1 = input_data1[adapt_map,::-1,:]
    Br_datam = input_datam[adapt_map,::-1, :]
    Br_datap = input_datap[adapt_map,::-1, :]
    #Br_data = coef1 * Br_data + coef2 * Br_data1
    # Linerar temporal interpolation
    if(Interp_order==1):
      Br_Linearinterp = np.zeros((nb_th, nb_phi))
      Br_Linearinterp=coef1 * Br_data + coef2 * Br_data1
      Br_data=Br_Linearinterp
    # Cubic Hermit temporal interpolation
    elif(Interp_order==2):
      Br_Linearinterp = np.zeros((nb_th, nb_phi))
      Br_Linearinterp = coef1 * Br_data + coef2 * Br_data1
      time_norm = coef2
      h00 = 2.0 * time_norm**3.0 - 3.0 * time_norm**2.0 + 1.0
      h10 = time_norm**3.0 - 2.0 * time_norm**2.0 + time_norm
      h01 = -2.0 * time_norm**3.0 + 3.0 * time_norm**2.0
      h11 = time_norm**3.0 - time_norm**2.0
      derivative1 = np.zeros((nb_th, nb_phi))
      derivative1 = 0.5 * ((Br_data - Br_datam) / tempm + (Br_data1 - Br_data) / temp)
      derivative2 = np.zeros((nb_th, nb_phi))
      derivative2 = 0.5 * ((Br_datap - Br_data1) / tempp + (Br_data1 - Br_data) / temp)
      Br_data = Br_data * h00 + derivative1 * temp * h10 + Br_data1 * h01 + derivative2 * temp * h11

    print('Br_data[0][0]=', Br_data[0][0], 'Br_mode[179][0]=', Br_data[179][0],'Br_data[12][0]=', Br_data[12][0], 'Br_mode[167][0]=', Br_data[167][0])
    i=20
    j=30
    print('Br_Linearinterp[i][j]-Br_data[i][j])=', Br_Linearinterp[i][j]-Br_data[i][j], 'i, j=',i,j)
    print('(Br_Linearinterp[i][j]-Br_data[i][j])/(0.5*(Br_Linearinterp[i][j]+Br_data[i][j]))=',(Br_Linearinterp[i][j]-Br_data[i][j])/(0.5*(Br_Linearinterp[i][j]+Br_data[i][j])), 'i, j=',i,j)
    Br_data = np.nan_to_num(Br_data)
    if (corotation == 'no'):
      # 将读出的磁图数据向左偏移回 dLong.
      # Reoffset the downloaded magnetograph data dLong degrees to the left.
      shift_adapt = shift_adapt0-ii*candence*360.0/653.0
      while shift_adapt < 0.0:
        shift_adapt+=360.0
      while shift_adapt >= 360.0:
        shift_adapt-=360.0
      dLong_adapt=round(shift_adapt/360.0*nb_phi)
      R_dLong=nb_phi-dLong_adapt
      Br_data = np.hstack((Br_data[:, -R_dLong:], Br_data[:, :-R_dLong]))
      Br_Linearinterp = np.hstack((Br_Linearinterp[:, -R_dLong:], Br_Linearinterp[:, :-R_dLong]))
    Br_map = Br_data
    Br_Linearinterp = np.nan_to_num(Br_Linearinterp)
    Br_mapLinear = Br_Linearinterp
    d1 = np.pi / nb_th
    d2 = 2.0*np.pi / nb_phi
    theta = np.linspace(0.5*d1,np.pi-0.5*d1,nb_th)  # 0:pi
    phi = np.linspace(d2,2.0*np.pi,nb_phi)
    Theta = np.tile(theta, (nb_phi,1)).T
    Phi = np.tile(phi, (nb_th,1))
    print('End of reading file')
  # WSO
  elif (map_type == 'wso'):
    fwso = open(input_file,'r')
    line = fwso.readline().split()
    if ('sine' in line):
      lat_type = 'sinlat'
    else:
      lat_type = 'lat'
    nb_th = int(line[1])
    nb_phi = int(360/5+1)
    nb_lines = 4*nb_phi
    nb_thplus = 4
    nb_th2 = nb_th + 2*nb_thplus
    Br_read = np.zeros((nb_th, nb_phi))
    fwso.readline()
    idx_th = nb_thplus
    idx_ph = nb_phi
    for ll in range(nb_lines):
      line = fwso.readline()
      if (line.split()[0][0] == 'C'):
        #idx_th = nb_thplus
        idx_th = 0
        idx_ph = idx_ph - 1
        for k in range(len(line.split())-1):
          Br_read[idx_th,idx_ph] = float(line.split()[k+1])
          idx_th = idx_th + 1
      else:
        for k in range(len(line.split())):
          Br_read[idx_th,idx_ph] = float(line.split()[k])
          idx_th = idx_th + 1
    fwso.close()
    if (lat_type == 'lat'):
      print('Extending Br')
      Br_ext = np.zeros((nb_th2, nb_phi))
      Br_ext[nb_thplus:nb_th2-nb_thplus,:] = Br_read
      idx_th = 0
      for k in range(nb_thplus):
        Br_ext[idx_th,:] = Br_read[0,:]
        idx_th = idx_th + 1
      idx_th = 0
      for k in range(nb_thplus):
        Br_ext[nb_th2-1 - idx_th,:] = Br_read[-1,:]
        idx_th = idx_th + 1
      Br_map = Br_ext*0.01 # from micro-tesla to gauss
      theta = (np.linspace(-90.,90.,nb_th2)+90.)*np.pi/180.
      phi = np.linspace(0.,360.,nb_phi)*np.pi/180.
      #theta = np.linspace(90.,-90.,nb_th2)
      #phi = np.linspace(0.,360.,nb_phi)
      Theta = np.tile(theta, (nb_phi,1)).T
      Phi = np.tile(phi, (nb_th2,1))
      nb_th = nb_th2
    else:
      Br_data = Br_read[::-1,:]*0.01 # from micro-tesla to gauss
      sinlat = np.linspace(-14.5/15.,14.5/15.,nb_th)
      theta_map = np.arcsin(sinlat) + np.pi/2.
      theta = np.linspace(0.,np.pi,nb_th)
      phi = np.linspace(0.,360.,nb_phi)*np.pi/180.
      Theta = np.tile(theta, (nb_phi,1)).T
      Theta_map = np.tile(theta_map, (nb_phi,1)).T
      Phi = np.tile(phi, (nb_th,1))
      #Br_data = Br_read[::-1,:]*0.01/np.cos(Theta_map) # from micro-tesla to gauss + from LOS to Br
      fbr = interpolate.RectBivariateSpline(theta_map,phi,Br_data)
      Br_map = fbr(theta,phi)
      Br_map = Br_map[::-1,:]
      #Br_map = Br_map/np.cos(Theta) # from LOS to Br
    print('End of reading file')
  # GONG, HMI, MDI, KPVT, SOLIS, MWO
  else:
    input_data = fits.getdata(input_file, ext=0)
    input_data1 = fits.getdata(input_file1, ext=0)
    input_datam = fits.getdata(input_filem, ext=0)
    input_datap = fits.getdata(input_filep, ext=0)
    shape = np.shape(input_data)
    nb_th = shape[0]
    nb_phi = shape[1]
    print('shape[0]=,',shape[0],'shape[1]=',shape[1])
    Br_data = input_data[::-1,:]
    Br_data1 = input_data1[::-1, :]
    Br_datam = input_datam[::-1, :]
    Br_datap = input_datap[::-1, :]
    dLong = int(map_name.split("_")[-1].split(".")[0])-1 # GONG磁图在经度上的偏移
    dLong1 = int(map_name1.split("_")[-1].split(".")[0])-1
    dLongm = int(map_namem.split("_")[-1].split(".")[0])-1
    dLongp = int(map_namep.split("_")[-1].split(".")[0])-1
    print('dLong=', dLong, 'dLong1=', dLong1,'dLongm=', dLongm, 'dLongp=', dLongp)
    # 将读出的磁图数据向右偏移 dLong.
    # Offset the downloaded magnetograph data dLong degrees to the right.
    Br_data = np.hstack((Br_data[:, -dLong:], Br_data[:, :-dLong]))
    Br_data1 = np.hstack((Br_data1[:, -dLong1:], Br_data1[:, :-dLong1]))
    Br_datam = np.hstack((Br_datam[:, -dLongm:], Br_datam[:, :-dLongm]))
    Br_datap = np.hstack((Br_datap[:, -dLongp:], Br_datap[:, :-dLongp]))

    # 20231204 #
    # corrected version
    d1=2.0/nb_th
    d2=2.0*np.pi/nb_phi
    sinlat = [-1.0 + (0.5+i)*d1 for i in range(nb_th)]
    #sinlat = np.linspace(-1.+0.5*d1, 1.-0.5*d1, nb_th)
    # print('sinlat[0]=',sinlat[0],'sinlat[nb_th-1]=',sinlat[nb_th-1])
    # original version
    #sinlat = np.linspace(-1.,1.,nb_th)
    theta = np.arcsin(sinlat) + np.pi/2.  # 0:pi
    print('theta[0]=',theta[0],'theta[90]=',theta[90],'theta[179]=',theta[179])
    phi = np.linspace(d2,2.0*np.pi,nb_phi)
    #print('phi[1]=', phi[1], 'phi[nb_phi-1]=', phi[nb_phi - 1])
    Theta = np.tile(theta, (nb_phi,1)).T
    Phi = np.tile(phi, (nb_th,1))

    #Br_data = coef1 * Br_data + coef2 * Br_data1
    # Linerar temporal interpolation
    Br_Linearinterp = np.zeros((nb_th, nb_phi))
    Br_Linearinterp=coef1 * Br_data + coef2 * Br_data1
    if (TransCorr == 'yes' or Interp_order == 1):
      Br_data = Br_Linearinterp
    # Cubic Hermit temporal interpolation
    else:
      time_norm = coef2
      h00 = 2.0 * time_norm**3.0 - 3.0 * time_norm**2.0 + 1.0
      h10 = time_norm**3.0 - 2.0 * time_norm**2.0 + time_norm
      h01 = -2.0 * time_norm**3.0 + 3.0 * time_norm**2.0
      h11 = time_norm**3.0 - time_norm**2.0
      derivative1 = np.zeros((nb_th, nb_phi))
      derivative1 = 0.5 * ((Br_data - Br_datam) / tempm + (Br_data1 - Br_data) / temp)
      derivative2 = np.zeros((nb_th, nb_phi))
      derivative2 = 0.5 * ((Br_datap - Br_data1) / tempp + (Br_data1 - Br_data) / temp)
      Br_data = Br_data * h00 + derivative1 * temp * h10 + Br_data1 * h01 + derivative2 * temp * h11
    i=20
    j=30
    print('Br_Linearinterp[i][j]-Br_data[i][j])=', Br_Linearinterp[i][j]-Br_data[i][j], 'i, j=',i,j)
    print('(Br_Linearinterp[i][j]-Br_data[i][j])/(0.5*(Br_Linearinterp[i][j]+Br_data[i][j]))=',(Br_Linearinterp[i][j]-Br_data[i][j])/(0.5*(Br_Linearinterp[i][j]+Br_data[i][j])), 'i, j=',i,j)
    Br_data = np.nan_to_num(Br_data)
    if (corotation == 'no'):
      # 将读出的磁图数据向左偏移回 dLong.
      # Reoffset the downloaded magnetograph data dLong degrees to the left.
      R_dLong=nb_phi-dLong
      Br_data = np.hstack((Br_data[:, -R_dLong:], Br_data[:, :-R_dLong]))
      Br_Linearinterp = np.hstack((Br_Linearinterp[:, -R_dLong:], Br_Linearinterp[:, :-R_dLong]))

    Br_map = Br_data
    Br_Linearinterp = np.nan_to_num(Br_Linearinterp)
    Br_mapLinear = Br_Linearinterp
    # 20231204 #
    # don't require invertion
    print('Br_data[0][0]=', Br_data[0][0], 'Br_data[179][0]=', Br_data[179][0])
    print('Br_data[0][359]=', Br_data[0][359], 'Br_data[179][359]=', Br_data[179][359])
    with fits.open(input_file) as hdu:
      img_data_eit = hdu[0].data
      img_data_eit = np.nan_to_num(img_data_eit)
      if (corotation == 'yes'):
        img_data_eit = np.hstack((img_data_eit[:, -dLong:], img_data_eit[:, :-dLong]))
      # require invertion
      print('img_data_eit[0][0]=', img_data_eit[0][0], 'img_data_eit[179][0]=', img_data_eit[179][0])
    print('End of reading file')
  Br = Br_map
  BrLinear = Br_mapLinear
  print('nb_th=',nb_th,'nb_phi=',nb_phi)
  #raise SystemExit()

  if(flux_correct=='yes'):
    print('Beginning of correction-flux')
    # calculate magnetic flux
    N = nb_phi  # 将每个纬线圈均分成360份
    R = 6371393  # 地球半径，单位m
    S = 4 * np.pi * R ** 2  # 球表面积
    flux = 0
    # 计算总的磁通量
    for i, line in enumerate(Br):
      if(i==0):
        line_flux = (abs(4. * np.pi * R ** 2. / N * (\
          math.cos(0.0) - math.cos(0.5*(theta[i]+theta[i+1])))) * line).sum()
      elif(i==nb_th-1):
        line_flux = (abs(4. * np.pi * R ** 2. / N * ( \
              math.cos(0.5*(theta[i]+theta[i-1])) - math.cos(np.pi))) * line).sum()
      else:
        line_flux = (abs(4. * np.pi * R ** 2. / N * ( \
              math.cos(0.5 * (theta[i] + theta[i-1])) - math.cos(0.5 * (theta[i] + theta[i+1])))) * line).sum()
      flux += line_flux
    print('flux / S=', flux / S)
    # print(flux / S)
    # print(img_data_eit2.flatten().sum())
    #print('Br[0][0]=',Br[0][0],'Br[179][359]=',Br[179][359],'Br[90][180]=',Br[90][180])
    Br -= flux / S
    #print('Br[0][0]=', Br[0][0], 'Br[179][359]=', Br[179][359], 'Br[90][180]=', Br[90][180])
    # print(img_data_eit2.flatten().sum())
    # exit()
    print('End of correction-flux')

  # Decomposition of Br on spherical harmonics
  print('Beginning of projection')
  theta_mid = np.zeros(nb_th+1) #list()
  theta_mid[0]=0.0
  for i in range(1, nb_th):
    theta_mid[i] = 0.5 * (theta[i-1]+theta[i])
  theta_mid[nb_th] = np.pi
  #shape = np.shape(theta_mid)
  #print('shape=,', shape, 'shape[0]=', shape[0])![](test_maps/ADAPT/Figure/correct/map_adapt_lmax15_20110618020000.dat.png)
  #dtheta = np.tile(np.concatenate([np.diff(theta), [theta[-1] - theta[-2]]]), (nb_phi, 1)).T
  dtheta = np.tile(np.concatenate([np.diff(theta_mid)]),(nb_phi,1)).T
  dphi = np.tile(np.concatenate([np.diff(phi),[phi[1]-phi[0]]]),(nb_th,1))
  ##2023_10_12-->coefbr = np.zeros(nb_modes_tot, dtype=np.complex)
  coefbr = np.zeros(nb_modes_tot, dtype=complex)
  coefbrLinear = np.zeros(nb_modes_tot, dtype=complex)
  ylm = np.zeros((nb_modes_tot, nb_th, nb_phi), dtype=complex)
  mod = 0
  l=0
  for l in range(1, lmax+1):
    # 20231204 #
    #print('{}'.format(l))
    for m in range(0, l+1):
      ylm[mod] = scisp.sph_harm(m, l, Phi, Theta)*1/(1+second_filter_factor*l**2*(l+1)**2)
      ylm_c = np.conj(ylm[mod])
      integrand_a = Br*ylm_c
      integrand_a = integrand_a*np.sin(Theta)*dphi*dtheta
      coefbr[mod] = np.sum(integrand_a)
      integrand_aLinear = BrLinear * ylm_c
      integrand_aLinear = integrand_aLinear * np.sin(Theta) * dphi * dtheta
      coefbrLinear[mod] = np.sum(integrand_aLinear)
      mod = mod+1
  print('l_max=','{}'.format(l))
  print('End of projection')
  if (map_type == 'ADAPT'):
    sinlat2lat='no'
  if(sinlat2lat=='yes'):
    d1 = np.pi / nb_th
    d2 = 2.0 * np.pi / nb_phi
    theta = np.linspace(0.5 * d1, np.pi - 0.5 * d1, nb_th)  # 0:pi
    phi = np.linspace(d2, 2.0 * np.pi, nb_phi)
    Theta = np.tile(theta, (nb_phi, 1)).T
    Phi = np.tile(phi, (nb_th, 1))

  # Reconstruction of the field
  print('Reconstructing Br')
  Br_mode = np.zeros((nb_th,nb_phi))
  Br_modeLinear = np.zeros((nb_th,nb_phi))
  mod = 0
  # Reconstruct field up to lmax
  for l in range(1,lmax+1):
    # 20231204 #
    #print('{}'.format(l))
    for m in range(0,l+1):
      ylm = scisp.sph_harm(m, l, Phi, Theta)
      Br_mode = Br_mode + np.real(coefbr[mod]*ylm)
      Br_modeLinear = Br_modeLinear + np.real(coefbrLinear[mod] * ylm)
      mod = mod+1

  Br_mode = Br_mode /2.2 # normalization of CF
  Br_mode = Br_mode*amp # amplitude factor
  Br_modeLinear = Br_modeLinear /2.2 # normalization of CF
  Br_modeLinear = Br_modeLinear*amp # amplitude factor
  print('End of reconstructing Br')


  print('Beginning of scaleing Br')
  Br_max=0.00000001

  for j in range(nb_th):
    for k in range(nb_phi):
      Br_max = max(np.abs(Br_mode[j, k]), Br_max)
  dailyBr_max[ii] = Br_max
  Scale=1.0
  if(ii>=hour_ave):
    if(B_scale=='yes'):
      if(Brmax_const=='yes'):
        Scale = Bref_Scale / Br_max
      else:
        Scale = min(Bref_Scale0, Br_max) / Br_max
  print('ii-hour_ave=',ii-hour_ave,'Br_max=',Br_max*2.2,'Bref_Scale=',Bref_Scale*2.2,'Scale=',Scale)
  print('End of scaleing Br and Br_max=',Br_max*2.2)
  # <- by HP Wang

  # Show maps
  if (show_map == 'yes'):
    if (plot_relative == 'yes'):
      fig, (ax1, ax2, ax4, ax3) = plt.subplots(4, 1, figsize=(8, 7), sharex=True)
    else:
      fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4), sharex=True)
    if (map_type == 'ADAPT'):
      visu_type = 'lat'
    # Latitudes and longitudes
    if (visu_type == 'sinlat'):
      if (map_type == 'wso'):
        if (lat_type == 'lat'):
          print('Careful, input file in latitudes, switching to lat plot!')
        else:
          longi = 180.*phi/np.pi
          Sinlat, Sinlong = np.meshgrid(sinlat, longi, indexing='ij')
      else:
        longi = 180.*phi/np.pi
        Sinlat, Sinlong = np.meshgrid(sinlat, longi, indexing='ij')
    lat = 90. - 180.*theta/np.pi
    longi = 180.*phi/np.pi
    Lat, Long = np.meshgrid(lat, longi, indexing='ij')
    # Plot original map
    if (visu_type == 'lat'):
      if (map_type == 'wso'):
        im1 = ax1.pcolormesh(Long,Lat,Br,cmap='seismic',vmin=-np.max(Br)/5.,vmax=np.max(Br)/5.)
      else:
        im1 = ax1.pcolormesh(Long,Lat,Br,cmap='seismic',vmin=-np.max(Br)/10.,vmax=np.max(Br)/10.)
      ax1.set_ylabel('Latitude')
    else:
      if (map_type == 'wso'):
        if (lat_type == 'lat'):
          im1 = ax1.pcolormesh(Long,Lat,Br,cmap='seismic',vmin=-np.max(Br)/5.,vmax=np.max(Br)/5.)
        else:
          im1 = ax1.pcolormesh(Sinlong,Sinlat,Br_data,cmap='seismic',vmin=-np.max(Br_data)/5.,vmax=np.max(Br_data)/5.)
      else:
        im1 = ax1.pcolormesh(Sinlong,Sinlat,Br_data[::-1],cmap='seismic',vmin=-np.max(Br_data)/10.,vmax=np.max(Br_data)/10.)
      longi_pos = np.arange(0.,360.,60.)
      ax1.set_xticks(longi_pos)
      #ax1.grid(b=True, which='major', color='k', linestyle='-')
      ax1.set_ylabel('Sine Latitude')
    ax1.set_title('Original map')
    plt.colorbar(im1,ax=ax1)
    # Plot lmax reconstruction
    #im2 = ax2.pcolormesh(Long,Lat,Br_mode,cmap='seismic',vmin=-np.max(Br_mode),vmax=np.max(Br_mode))
    im2 = ax2.pcolormesh(Long,Lat,Br_mode,cmap='seismic',vmin=Brmin,vmax=Brmax)
    ax2.set_title('BC file')
    ax2.set_ylabel('Latitude')
    #ax2.set_xlabel('Longitude')
    plt.colorbar(im2,ax=ax2)
    if(plot_relative == 'yes'):
      # Plot Absolute difference
      im4 = ax4.pcolormesh(Long,Lat,Br_mode-Br_modeLinear,cmap='seismic',vmin=-np.max(Br_mode-Br_modeLinear),vmax=np.max(Br_mode-Br_modeLinear))
      #im2 = ax2.pcolormesh(Long,Lat,Br_mode,cmap='seismic',vmin=-5,vmax=5)
      ax4.set_title('Absolute difference')
      ax4.set_ylabel('Latitude')
      #ax2.set_xlabel('Longitude')
      plt.colorbar(im4,ax=ax4)
      # Plot Relative difference
      im3 = ax3.pcolormesh(Long,Lat,(Br_mode-Br_modeLinear)/(0.5*(abs(Br_mode)+abs(Br_modeLinear))),cmap='seismic',\
                       vmin=-np.max((Br_mode-Br_modeLinear)/(0.5*(abs(Br_mode)+abs(Br_modeLinear)))),\
                       vmax=np.max((Br_mode-Br_modeLinear)/(0.5*(abs(Br_mode)+abs(Br_modeLinear)))))
      #im2 = ax2.pcolormesh(Long,Lat,Br_mode,cmap='seismic',vmin=-5,vmax=5)
      ax3.set_title('Relative difference')
      ax3.set_ylabel('Latitude')
      ax3.set_xlabel('Longitude')
      plt.colorbar(im3,ax=ax3)
    plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.2)
    plt.savefig(output_Figre_name +'.png', bbox_inches='tight', pad_inches=0.05)
    #plt.savefig(output_Figre_name +'.png')
    plt.close('all')

    #plt.show()
    #plt.draw()
    #plt.pause(4)

    # Write boundary conditions file
  if (write_map == 'yes'):
    print('Writing BC file')
    F = open(output_name, 'w')
    F.write('1 \n')
    F.write('!PHOTOSPHERE {} \n'.format((nb_th - 2) * nb_phi + 2))
    #[0.5d1,pi-0.5d1],0:179
    for j in range(nb_th):
      # [d2,2*pi],0:359
      for k in range(nb_phi):
        xcoord = r_st * np.sin(theta[j]) * np.cos(phi[k])
        ycoord = r_st * np.sin(theta[j]) * np.sin(phi[k])
        zcoord = r_st * np.cos(theta[j])
        if ((j == 0) & (k != 0)):
          break
        if ((j == nb_th - 1) & (k != 0)):
          break
        F.write('{:.16e} {:.16e} {:.16e} {:.16e} \n'.format(xcoord, ycoord, zcoord, Br_mode[j, k] * Scale))
    F.close()
    print('End of writing BC file')

  