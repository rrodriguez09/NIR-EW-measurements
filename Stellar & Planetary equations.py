import numpy as np
import math
from fractions import Fraction

#Planet Radius


R_earth = 6371000
transit_depth = 0.18044
R_sun = 6.957*10**(8)
R_Rsun = 0.494977203684
R_star = R_sun*R_Rsun
Rp = transit_depth*R_star
Rp_earthunits = Rp/R_earth


#Planet radius error

transit_depth = 0.1957
t_depth_error = 0.007398


R_Rsolar = 0.639818482237

Rstar_error = 0.023039884392

Rp_error = np.sqrt((transit_depth*Rstar_error)**2 + (0.5*t_depth_error*(R_Rsolar/transit_depth))**2)

#Planet temperature

T = 32.5
period = T*86400

expo = Fraction('1/3')
semimajor = ((G*Mass_star*(period**2))/39.4784)**expo
semi_AU = (semimajor) / (1.496*10**(11))

l_lsun =-0.415459581621
luminosity = L_solar*math.exp(l_lsun)

planet_Teff = (np.sqrt(np.sqrt(luminosity/(16*np.pi*SB_factor*semimajor**2))))

#Stellar Mass

Sun_mass = 1.989*10**(30)
M_star = 0.436*Sun_mass
#Boyajian et al 2012 Mass-radius relation
R_Rsolar = 0.412953450707
Mass = (-0.6063 + np.sqrt((0.6063)**2 - 4*0.32*(0.0906 - R_Rsolar)))/(2*0.32)

#Stellar Luminosity
L_solar = 3.846*10**(26)
l_lsun = -0.41
luminosity = L_solar*math.exp(l_lsun)
Lum = luminosity/L_solar
#luminosity = 0.3
#distance = 9*10**12

SB_factor = 5.670373*10**(-8)
T = 3.4
period = T*86400
Teff = 3621
G = 6.67408*10**(-11)
expo = Fraction('1/3')
distance = ((G*M_star*(period**2))/39.47)**expo
Teff_planet = np.sqrt(np.sqrt(luminosity/(50.265*SB_factor*distance**2)))


#Stellar Luminosity
L_solar = 3.846*10**(26)
l_lsun = -0.415459581621
luminosity = L_solar*math.exp(l_lsun)


M_solar = 1.9*10**(30) print Teff_planet
