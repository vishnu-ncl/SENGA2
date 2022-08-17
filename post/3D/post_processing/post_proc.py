import numpy as np

flag_main = 0

def read_files(s1,s2,s,num_spc,lenx,leny,lenz):
    import h5py

    for i in range(s1,s2+1,s):
        
        print("reading file:",i)
        
        if i<10:
            filename = "../3B000{}.h5".format(i)
        elif i<100:
            filename = "../3B00{}.h5".format(i)
        else:
            filename = "../3B0{}.h5".format(i)
        
        f = h5py.File(filename, "r")

        global varname
        global fld_num
        global vardata
        global xyzdata
        global xdata
        global ydata
        global zdata
        varname = []
        fld_num = 2*num_spc+6
        vardata = np.zeros((fld_num,lenz,leny,lenx),dtype=float) ##Need to check if it is z, y, x or y, z, x
        xyzdata = np.zeros((lenx,leny,lenz),dtype=float)
        xdata = np.zeros(lenx,dtype=float)
        ydata = np.zeros(leny,dtype=float)
        zdata = np.zeros(lenz,dtype=float)

        for key in f.keys():
            ##print(key) #Names of the root level object names in HDF5 file - can be groups or datasets.
            ##print(type(f[key])) # get the object type: usually group or dataset
        
            #Get the HDF5 group; key needs to be a group name from above
            group = f[key]
            
            knum=0
        
            #Checkout what keys are inside that group.
            for key in group.keys():
                ##print(key)
                varname.append(key)
            
                # This assumes group[some_key_inside_the_group] is a dataset,
                # and returns a np.array:
                data = group[key][()]
                ##print(data)
                #Do whatever you want with data

                ddim = data.ndim
                if ddim == 3:
                    dlen1 = len(data)
                    dlen2 = len(data[0])
                    dlen3 = len(data[0][0])
                    for c1 in range(dlen1):
                        for c2 in range(dlen2):
                            for c3 in range(dlen3):
                                vardata[knum][c1][c2][c3]=data[c1][c2][c3]
                else:
                    dlen = len(data)
                    for c1 in range(dlen):
                        if knum==0:
                            xdata[c1]=data[c1]
                        if knum==1:
                            ydata[c1]=data[c1]
                        if knum==2:
                            zdata[c1]=data[c1]

                knum = knum+1
        
        #After you are done
        f.close()

def flame_speed():
    if flag_main==0:
        global num_spc
        global lenx
        global leny
        global lenz
        fn = int(input("file number to be read:"))
        num_spc = int(input("Input number of species: "))
        lenx = int(input("Number of grid points in x direction :"))
        leny = int(input("Number of grid points in y direction :"))
        lenz = int(input("Number of grid points in z direction :"))
        read_files(fn,fn,1,num_spc,lenx,leny,lenz)
    rctnt = input("Species symbol for reactant :")
    rrt_rt = "RRTE_"+rctnt
    for i in range(fld_num):
        if varname[i] == rrt_rt:
            flag1 = i
        if varname[i] == rctnt:
            flag2 = i
        if varname[i] == "DRUN":
            flag3 = i
    dx = xdata[1]-xdata[0]
    s_rrt = 0.5*(vardata[flag1][0][0][0]*dx+vardata[flag1][0][0][lenx-1]*dx)
    for c1 in range(1,lenx-1):
        s_rrt = s_rrt+(vardata[flag1][0][0][c1]*dx)

    fs = -(1.0/(vardata[flag3][0][0][0]*vardata[flag2][0][0][0]))*s_rrt

    return fs

def flame_thick():
    if flag_main==0:
        global num_spc
        global lenx
        global leny
        global lenz
        fn = int(input("file number to be read:"))
        num_spc = int(input("Input number of species: "))
        lenx = int(input("Number of grid points in x direction :"))
        leny = int(input("Number of grid points in y direction :"))
        lenz = int(input("Number of grid points in z direction :"))
        read_files(fn,fn,1,num_spc,lenx,leny,lenz)
    for i in range(fld_num):
        if varname[i] == "TRUN":
            flag1 = i
    dx = xdata[1]-xdata[0]
    eps = 0.0
    dT = 0.0
    for i in range(1,lenx-1):
        dT = (vardata[flag1][0][0][i+1]-vardata[flag1][0][0][i-1])/(2.0*dx)
        if dT>eps:
            eps=dT
    ft = (vardata[flag1][0][0][lenx-1]-vardata[flag1][0][0][0])/eps
    
    return ft

def xtrapol():
    import matplotlib.pyplot as plt
    from scipy import interpolate
    if flag_main==0:
        global num_spc
        global lenx
        global leny
        global lenz
        fn = int(input("file number to be read :"))
        num_spc = int(input("Input number of species :"))
        lenx = int(input("Number of grid points in x direction :"))
        leny = int(input("Number of grid points in y direction :"))
        lenz = int(input("Number of grid points in z direction :"))
        read_files(fn,fn,1,num_spc,lenx,leny,lenz)
    len3x=int(input("Number of grid points in x direction in 3D domain :"))
    xstart=float(input("Start location :"))
    xstop=float(input("Stop location :"))
    x3data=np.linspace(xstart,xstop,len3x,endpoint=False,dtype=float)
    print(x3data)
    flagS = np.zeros(num_spc,dtype=int)
    for i in range(fld_num):
        if varname[i]=="DRUN":
            flagD = i
        elif varname[i]=="TRUN":
            flagT = i
        elif varname[i]=="URUN":
            flagU = i
    for i in range(num_spc):
        flagS[i]=flagD+i+1
    x31data = []
    i=0
    while x3data[i]<xdata[lenx-1]:
        x31data.append(x3data[i])
        i+=1
    f = interpolate.interp1d(xdata,vardata[flagD][0][0])
    drun3 = f(x31data)
    f = interpolate.interp1d(xdata,vardata[flagT][0][0])
    trun3 = f(x31data)
    f = interpolate.interp1d(xdata,vardata[flagU][0][0])
    urun3 = f(x31data)
    len31x=len(x31data)
    sp32 = np.zeros((num_spc,len31x))
    for i in range(num_spc):
        f = interpolate.interp1d(xdata,vardata[flagS[i]][0][0])
        sp32[i] = f(x31data)

    sp31=np.zeros((num_spc,(len3x-len31x)))
    sp3=np.concatenate((sp32,sp31),axis=1)
    for i in range(len31x,len3x):
        temp = drun3[len31x-1]
        drun3=np.append(drun3,temp)#drun3.append(temp)
        temp = trun3[len31x-1]
        trun3=np.append(trun3,temp)#trun3.append(temp)
        temp = urun3[len31x-1]
        urun3=np.append(urun3,temp)#urun3.append(temp)
        for j in range(num_spc):
            temp=sp3[j][len31x-1]
            sp3[j][i]=temp
    #plt.plot(xdata, vardata[flagS[1]][0][0], 'o', x3data, sp3[1], '-')
    #plt.show()

    t=open('lamflame_py.dat','w')
    t.write('%4i \n'%(len3x))
    t.write("DRUN \n")
    for i in range(len3x):
        t.write('%21.16e \n'%(drun3[i]))
    t.write("TRUN \n")
    for i in range(len3x):
        t.write('%21.16e \n'%(trun3[i]))
    t.write("URUN \n")
    for i in range(len3x):
        t.write('%21.16e \n'%(urun3[i]))
    for j in range(num_spc):
        t.write("{} \n".format(varname[flagS[j]]))
        for i in range(len3x):
            t.write('%21.16e \n'%(sp3[j][i]))
    t.close()
    print(varname)

def main():

    start = int(input("Start file number :"))
    stop = int(input("Stop file number :"))
    step = int(input("Step file size :"))
    global num_spc
    global lenx
    global leny
    global lenz
    num_spc = int(input("Input number of species: "))
    lenx = int(input("Number of grid points in x direction :"))
    leny = int(input("Number of grid points in y direction :"))
    lenz = int(input("Number of grid points in z direction :"))

    read_files(start,stop,step,num_spc,lenx,leny,lenz)

    menu_flag=1
    while menu_flag == 1:
        print("Menu:")
        print("1) Flame speed calculation")
        print("2) Flame thickness calculation")
        print("3) Extrapolation")
        print("99) Exit Menu")
        opt = int(input("Option :"))

        if opt==1:
            sL = flame_speed()
            print("flame speed is :",sL)
        elif opt==2:
            df = flame_thick()
            print("flame thickness:", df)
        elif opt==3:
            xtrapol()
            print("Executed the extrapolation")
        elif opt==99:
            print("Exiting!")
            menu_flag=0
        else:
            print("Invalid!")

if __name__ == "__main__":
    flag_main=1
    main()
