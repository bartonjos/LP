#Calculate campaig integrated parameters for Metal Rings
#Campaign shot range: 167064 - 167622


start = time.clock()

#d3dRDB access example
"""
d3dRDB_summary = OMFITrdb(query="select * from summaries where shot="+str(shot), db='d3drdb', server='d3drdb', by_column=True)
d3dRDB_shot = OMFITrdb(query="select * from shots where shot="+str(shot), db='d3drdb', server='d3drdb', by_column=True)

root['SCRIPTS']['20160624_WR']['167334']['OUTPUTS']['d3dRDB_summary'] = d3dRDB_summary
root['SCRIPTS']['20160624_WR']['167334']['OUTPUTS']['d3dRDB_shot'] = d3dRDB_summary
"""



#Step through all shots and get RVSOUT statistics

time_ShelfRing = 0 #Time spent on shelf ring
R_shelf_ring = (1.404, 1.454) #[m] shelf ring location

time_Ip = 0 #Time with Ip greater than 0
Ip_thres = 0.5e6 #[A] Threshold for Ip to be counted





"""
#Collect statistics on OSP_R
shot_succ = 0 #Counter of good shots

#for shot in range(167450,167450+2):
for shot in range(167064,167622+1):
    try:
        d_OSP_R = OMFITmdsValue(server='DIII-D', shot = shot, TDI='RVSOUT')
        
        #Find times on shelf ring
        temp = (d_OSP_R.data()>R_ShelfRing[0]) *  (d_OSP_R.data()<R_ShelfRing[1])
    except MDSplus.MdsException: continue

    if np.sum(temp)<5: continue #Check that the bool array is not too short


    #Update total time
    time_ShelfRing+= np.sum( np.diff(d_OSP_R.dim_of(0))[temp[:-1]]  )
        
    #Update counter
    shot_succ+=1

print('Number of shots with RVSOUT1 data: {:d}'.format(int(shot_succ)))
print('Time on the Shelf Ring: {:.1f} s'.format(time_ShelfRing/1000.))


"""


#Collect statistics on OSP_R

#Counter of successful shots
shot_succ=0

for shot in range(167064,167622+1):
    try:
        d_Ip = OMFITmdsValue(server='DIII-D', shot = shot, TDI='IP')
        
        #Find times on shelf ring
        temp = np.abs(d_Ip.data())>Ip_thres
    except MDSplus.MdsException: continue

    if np.sum(temp)<5: continue #Check that the bool array is not too short

    #Update total time
    time_Ip+= np.sum( np.diff(d_Ip.dim_of(0))[temp[:-1]]  )
        
    #Update counter
    shot_succ+=1



print("\n\n")
print('Number of shots with Ip> threshold: {:d}'.format(int(shot_succ)))
print('Time with Ip>Ip_thres: {:.1f} s'.format(time_Ip/1000.))


print('Time spent: {:2.2f} s'.format(time.clock() - start))
