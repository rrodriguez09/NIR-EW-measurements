#Complete Code - Equivalent Width Calculations
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab



def romy_temp(filename):
    #filename = "/Users/Romy/Desktop/20160122/proc/287-290_tc.fits"
    #filename = "fake_spec82.fits"
    #filename = "/Users/Romy/Desktop/20160121/proc/269-272_tc.fits"

    hdu =pyfits.open(filename)
    tabledata= hdu[0].data



    datafile= "/Users/Romy/desktop/line_positions"
    names= ["line","line_start","line_end","blue_start","blue_end", "red_start", "red_end"]
    sim=np.genfromtxt(datafile,dtype=None,delimiter=' ',autostrip=True,names=names)

    line=sim[names[0]]
    line_start = sim[names[1]]
    line_end=sim[names[2]]
    blue_start = sim[names[3]]
    blue_end = sim[names[4]]
    red_start = sim[names[5]]
    red_end = sim[names[6]]


    #for-loop implementation
    Equiv_widths = []
    Sigmas = []
    for i in range(9):

        my_feature = (line_start[i] + line_end[i])/2
        if i <= 6:
            wave = tabledata[1,0,:]
            flux = tabledata[1,1,:]
            err = tabledata[1,2,:]

        elif i > 6:
            wave = tabledata[0,0,:]
            flux = tabledata[0,1,:]
            err = tabledata[0,2,:]

        goodind = (wave > blue_start[i]) & (wave < red_end[i])


        my_lam = wave[goodind]
        my_data= flux[goodind]
        my_sigmas = err[goodind]


        line_ind1 = (my_lam > blue_start[i]) & (my_lam < blue_end[i]) 
        line_ind2 = (my_lam > red_start[i]) & (my_lam < red_end[i])
        line_ind = line_ind1 + line_ind2

    #Define our models
        def my_model(lam,center,height,width):

            return height * (width/2.)**2 / ((lam - center)**2 + (width/2.)**2)
        def my_Gaussian_model(lam,center,height,width):
            return height * np.exp(-0.5 * (lam - center)**2/width**2)
        def lin_model(lam,m,b):
            return m*lam + b
        def metric_chi2(theta,lam,data,sigmas):
            center, height, width = theta
            return np.sum((my_model(lam,center,height,width)-data)**2/(sigmas)**2)

        p0 = [0,0]
        popt,pcov = curve_fit(lin_model,my_lam[line_ind],my_data[line_ind],sigma=my_sigmas[line_ind],p0=p0)
        my_data = my_data/(lin_model(my_lam,popt[0],popt[1]))

        #p0 = [my_feature,-0.2,0.0005]
        #popt,pcov = curve_fit(my_Gaussian_model,my_lam,my_data -1,sigma=my_sigmas,p0=p0)

        #equivalent width with Trapz

        goodind2 = (my_lam > line_start[i]) & (my_lam < line_end[i])
        #matrix = np.zeros((100,len(goodind2)))
        #orig_data = my_data[goodind2]

        #print len(goodind2)
        #print len(my_data)
        #for i in range(100):
            #noise_scale = np.random.uniform(low = -1, high = 1, size = 118)*my_sigmas[goodind2]

            #fake_data = orig_data + noise
        ew = np.trapz((1 - my_data[goodind2]), x=my_lam[goodind2])
        EW = ew*10000    
        Equiv_widths.append(EW)

        #Error on the Equivalent Widths of each line:

        #sigma_center = np.sqrt(pcov[0,0])
        #sigma_height = np.sqrt(pcov[1,1])
        #sigma_width = np.sqrt(pcov[2,2])

        #sigma_EW = np.sqrt((np.sqrt(2*np.pi)*popt[2]*sigma_height*10000)**2 + (np.sqrt(2*np.pi)*popt[1]*sigma_width*10000)**2)
        #Sigmas.append(sigma_EW)



    #print Equiv_widths
    #print Sigmas

    #-----------------------------------------------------------------------
    # Stellar Parameters:
    #Line indices 
    Mg1 = Equiv_widths[0]
    Mg2 = Equiv_widths[1]
    Mg3 = Equiv_widths[2]
    Al = Equiv_widths[3]
    Ala = Equiv_widths[4]
    Alb = Equiv_widths[5]
    Mg4 = Equiv_widths[6]
    Na = Equiv_widths[7]
    Ca = Equiv_widths[8]

    Teff = 271.4*(Ala) + 392.7 * (Mg2/Alb) + 2427
    R_Rsolar = -0.0489*(Mg3) + 0.275*(Ala) + 0.201*(Mg3/Ala) - 0.216
    L_Lsolar = 0.832*(Mg4) - 0.176*((Mg4)**2) + 0.266*(Mg2) - 3.491

    #print "The equivalent widths are: %.4f for Mg(1.50), %.4f for Al-a, %.4f for Al-b, %.4f for Mg (1.57), %.4f for Mg (1.71)" % (EW_1, EW_2, EW_3, EW_4, EW_5)
    #print "Temperature = %.2f K" % (Teff) 
    #print "R/R_sun =  is %.2f solar radii" % (R_Rsolar)
    #print "Log(L/L_sun) = %.2f : " % (L_Lsolar)
    #print "Mg(1.50) = %.4f, Mg(1.57) = %.4f, Mg(1.71) = %.4f, Al-a = %.4f, Al-b = %.4f " % (Mg2, Mg3, Mg4, Ala, Alb)


    #Error/Uncertainty on the Temperature, Radius and Luminosity 

    #Error = np.sqrt(((271.4)*Sigmas[4])**2 + (((392.7)*(1/Alb)*Sigmas[1]))**2 + ( (392.7)*(Mg2/(Alb**2))  *Sigmas[5])**2)
    #print Error

    #Error_gaussian_T = np.sqrt(((271.4)*Sigmas[4])**2 + (((392.7)*(1/Alb)*Sigmas[1]))**2 + ( (392.7)*(Mg2/(Alb**2))  *Sigmas[5])**2)
    #print "The error on the temperature is: %.2f " % (Error_gaussian_T)

    #Error_gaussian_R = np.sqrt(((-0.0489 + 0.201/Ala)*Sigmas[2])**2 + ( (0.275 - 0.201*Mg3/(Ala**2))*Sigmas[4])**2)
    #print "The error on R/Rsun is: %.2f" % (Error_gaussian_R)

    #Ala_dist = np.random.normal(Ala,Sigmas[4],200)
    #plt.hist(Ala_dist)
    #plt.show()
    #print Ala
    #Alb_dist = np.random.normal(Alb, Sigmas[5],200)
    #plt.hist(Alb_dist)
    #plt.show()
    #print Alb
    #Mg2_dist = np.random.normal(Mg2, Sigmas[1],200)
    #plt.hist(Mg2_dist)
    #plt.show()
    #print Mg2

    #Teffs = []
    #for i in range(200):
        #Teff = 271.4*(Ala_dist[i]) + 392.7*(Mg2_dist[i])/(Alb_dist[i]) + 2427
        #Teffs.append(Teff)

    #print Teffs

    #plt.hist(Teffs)
    #plt.show()
    #xmin = 2000
    #xmax = 7000

    #plt.xlabel("Temperature (K)")
    #plt.ylabel("Frequency of temperature")
    #plt.xlim((min(Teffs), max(Teffs)))
    #plt.title("EPIC 211336288")
    #plt.hist(Teffs, normed=True)
    #x = np.linspace(0, 7000, 10)
    #plt.plot(x,mlab.normpdf(x, np.mean(Teffs), np.var(Teffs)))
    #plt.show()

    #mean = np.mean(Teffs)
    #variance = np.var(Teffs)
    #sigma = np.sqrt(variance)
    #x = np.linspace(min(Teffs), max(Teffs),100)
    #dx = result[1][1] - result[1][0]
    #scale = len(Teffs)*dx
    #plt.plot(x, mlab.normpdf(x,mean,sigma)*scale)

    #plt.show()

    #-----------------------------------------------------------------------
    # H2O-K2 Water Index 

    wave_ind1 = np.where((tabledata[0,0,:] > 2.070) & (tabledata[0,0,:] < 2.090))
    wave_ind2 = np.where((tabledata[0,0,:] > 2.235) & (tabledata[0,0,:] < 2.255))
    wave_ind3 = np.where((tabledata[0,0,:] > 2.360) & (tabledata[0,0,:] < 2.380))

    median1 = np.median(tabledata[0,1,:][wave_ind1])
    median2 = np.median(tabledata[0,1,:][wave_ind2])
    median3 = np.median(tabledata[0,1,:][wave_ind3])

    #H2O-K2 = (median1/median2)/(median2/median3)
    n1 = median1/median2
    n2 = median2/median3
    H2OK2 = n1/n2


    #-----------------------------------------------------------------------
    # Mann's Lines for Metallicity Calculation


    datafile= "/Users/Romy/Desktop/metal_lines"
    names= ["line","line_start","line_end","blue_start","blue_end", "red_start", "red_end"]
    sim=np.genfromtxt(datafile,dtype=None,delimiter=' ',autostrip=True,names=names)

    line=sim[names[0]]
    line_start = sim[names[1]]
    line_end=sim[names[2]]
    blue_start = sim[names[3]]
    blue_end = sim[names[4]]
    red_start = sim[names[5]]
    red_end = sim[names[6]]

    #for-loop implementation
    Equiv_widths_2 = []
    for i in range(3):

        wave = tabledata[0,0,:]
        flux = tabledata[0,1,:]
        err = tabledata[0,2,:]

        goodind = (wave > blue_start[i]) & (wave < red_end[i])


        my_lam = wave[goodind]
        my_data= flux[goodind]
        my_sigmas = err[goodind]


        line_ind1 = (my_lam > blue_start[i]) & (my_lam < blue_end[i]) 
        line_ind2 = (my_lam > red_start[i]) & (my_lam < red_end[i])
        line_ind = line_ind1 + line_ind2

    #Define our models
        def my_model(lam,center,height,width):

            return height * (width/2.)**2 / ((lam - center)**2 + (width/2.)**2)
        def my_Gaussian_model(lam,center,height,width):
            return height * np.exp(-0.5 * (lam - center)**2/width**2)
        def lin_model(lam,m,b):
            return m*lam + b
        def metric_chi2(theta,lam,data,sigmas):
            center, height, width = theta
            return np.sum((my_model(lam,center,height,width)-data)**2/(sigmas)**2)

        p0 = [0,0]
        popt,pcov = curve_fit(lin_model,my_lam[line_ind],my_data[line_ind],sigma=my_sigmas[line_ind],p0=p0)
        my_data = my_data/(lin_model(my_lam,popt[0],popt[1]))

        #p0 = [my_feature,-0.2,0.0005]
        #popt,pcov = curve_fit(my_Gaussian_model,my_lam,my_data -1,sigma=my_sigmas,p0=p0)

        #equivalent width with Trapz

        goodind2 = (my_lam > line_start[i]) & (my_lam < line_end[i])
        ew = np.trapz((1 - my_data[goodind2]), x=my_lam[goodind2])
        EW = ew*10000    
        Equiv_widths_2.append(EW)


    #print Equiv_widths_2

    F19 = Equiv_widths_2[0]
    F20 = Equiv_widths_2[1]
    F22 = Equiv_widths_2[2]

    metallicity = 0.19*(F19) + 0.069*(F22) + 0.083*(F20) + 0.218*(H2OK2) - 1.55
    #print "The metallicity is: %.4f"  % (metallicity)
    #print "H2O-K2 index = %.4f" % (H2OK2)
    #print "F19 = %.4f, F20= %.4f, F22 = %.4f" % (F19, F20, F22)
    #print Sigmas[1],Sigmas[4],Sigmas[5]


    #print noise_scale
    #print Equiv_widths[1], Equiv_widths[2], Equiv_widths[6], Equiv_widths[4], Equiv_widths[5]
    #print Equiv_widths_2
    #print H2OK2  
    #print metallicity
    return Teff, R_Rsolar, L_Lsolar, metallicity

R_Rsolars = np.zeros(100)
Teffs = np.zeros(100)
L_Lsolars = np.zeros(100)
metallicities = np.zeros(100)

names = ["spectra", "realfilename"]
a=np.genfromtxt("file1.txt",dtype=None,autostrip=True,names=names) 
list1 = a[names[0]]
names1 = a[names[1]]
directory = list2[0]
name = names1[0]

for j in range(len(list2)):
#for j in range(24,25):
    directory = list1[j]
    name = names1[j]
    #print j
    for i in range(100):
        #print i
        Teff, R_Rsolar, L_Lsolar, metallicity = romy_temp("/Users/Romy/desktop/fakespec1/" + directory+ "/fakespec" + str(i) +".fits")
        Teffs[i] = Teff
        R_Rsolars[i] = R_Rsolar
        L_Lsolars[i] = L_Lsolar
        metallicities[i] = metallicity
        Mass = (-0.6063 + np.sqrt((0.6063)**2 - 4*0.32*(0.0906 - R_Rsolar)))/(2*0.32)

    #print Teffs
    #plt.show()
    plt.subplot(221)
    plt.hist(Teffs)
    plt.suptitle(name)
    plt.xlabel("Temperature")
    plt.ylabel("Count")
    #plt.show()
    plt.subplot(222)
    plt.hist(R_Rsolars)
    plt.xlabel("Stellar Radius")
    plt.ylabel("Count")
    #plt.show()
    plt.subplot(223)
    plt.hist(L_Lsolars)
    plt.xlabel("Luminosity")
    plt.ylabel("Count")
    #plt.show()
    plt.subplot(224)
    plt.hist(metallicities)
    plt.xlabel("Metallicity")
    plt.ylabel("Count")
    #plt.show()
    plt.savefig("/Users/Romy/desktop/fakespec1/" + directory+ "/Stellar parameters.jpeg",dpi = 300)
    plt.close()
    print directory
    print np.mean(Teffs)
    print np.mean(R_Rsolars)
    print np.std(Teffs)
    print np.mean(L_Lsolars)
    print np.mean(metallicities)

    text_file = open("/Users/Romy/desktop/fakespec1/" + directory+ "/Stellar parameters.txt", "w")
    text_file.write("Temperature ="+ str(np.mean(Teffs))+ "\n")
    text_file.write("Temperature error ="+str(np.std(Teffs))+ "\n")
    text_file.write("Rsolar ="+str(np.mean(R_Rsolars))+ "\n")
    text_file.write("Radius error =" +str(np.std(R_Rsolars))+ "\n")
    text_file.write("Luminosity ="+str(np.mean(L_Lsolars))+ "\n")
    text_file.write("Luminosity error=" +str(np.std(L_Lsolars))+ "\n")
    text_file.write("Metallicity =" +str(np.mean(metallicities))+ "\n")
    text_file.write("Metallicity error =" + str(np.std(metallicities))+ "\n")
    text_file.write("Mass ="+ str(Mass))

    text_file.close()
