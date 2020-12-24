import numpy as np
import scipy.integrate as spi
import scipy.special as sps
import matplotlib.pyplot as plt


# Function for J1(x)
def J1(x, n):
    '''
    Approximation for values of first order first kind Bessel function.

    Takes arguments:
    x -- array of x values for the range of the Bessel function.
    n -- number of points of integration used in Simpson's rule.

    Uses Simpson's rule to approximate the integral between 0 and pi in 
    the first order first kind Bessel function.
    scipy.integrate.simps is used for Simpson's rule calulation.
    '''

    max_t = np.pi
    min_t = 0
    t = np.linspace(start = min_t, stop = max_t, num = n, endpoint = True)
    
    y = np.cos(t - x * np.sin(t))

    J1_simps.append((1/np.pi) * spi.simps(y, x=None, dx=max_t/n, axis=- 1, 
					  even='first'))
    
    return

# Number of intgration points to be used in function J1.
n = 1000
# Range of x values for Bessel function to be approximated over.
min_x = 0
max_x = 20

# Array of x values in the range above.
x = np.linspace(start = min_x, stop = max_x, num = n, endpoint = True)

J1_simps = []

# Using function J1 to approximate Bessel function values at x values.
for i in range(n):
    J1(x[i], n)
    
# Scipy Bessel functions for x values.
J1_scipy = sps.jv(1,x)

# Difference in the Bessel function values from the two methods.
diff = J1_scipy - J1_simps

# Plotting the values of the Bessel functions with n = 0, 1, 2, for x 
# values in the array "x".
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(7,7))

p1 = plt.subplot(2, 1, 1)
plt.plot(x, J1_simps, color = "blue", linestyle = "solid", 
	 label = "Simpson's Rule Approximation (n=%s)" %(n))
plt.plot(x, J1_scipy, color = "black", linestyle = "solid", 
   	 label = "Scipy Bessel Function")

# Adding axis labels and a title.
plt.ylabel("Spherical Bessel Function Value")
plt.xlabel("x")
plt.title("First Order of the First Kind Spherical Bessel Function", 
          fontsize=15)

# Adding axis lines and a legend.
plt.axhline(0,color='black')
plt.axvline(0,color='black') 
plt.legend()

# Plotting the difference between the two methods value's.
p2 = plt.subplot(2, 1, 2)
plt.plot(x, diff, color = "red", linestyle = "solid")

# Adding axis labels and title.
plt.ylabel("Difference")
plt.xlabel("x")
plt.title("Difference Between the Two Calulations", fontsize=15)

# Adding axis lines
plt.axhline(0,color='black')
plt.axvline(0,color='black') 

fig.tight_layout()
plt.show()

# Intensity of diffration pattern
# Wavelength and wavenumber of light to be used.
wl = 500 * 10**-9    # metres
k = (2 * np.pi) / wl
# Radius of the pattern to be computed.
radius = 1 * 10**-6    #metres

# Creating arrays of x and y points to produce a square region of
# lengths = 2 * radius.
x = np.linspace(start = -radius, stop = radius, num = 100, 
                endpoint = True)
y = np.linspace(start = -radius, stop = radius, num = 100, 
                endpoint = True)

# Intensities within a column and intensities of the region.
column_inten = []
inten = []

# Calulating the intensities of the points in each column. Intensities 
# of region given in "inten".
for i in range(100):
    for j in range(100):
        
        r = (x[i]**2 + y[j]**2)**(1/2)
        intensity = ((2 * sps.jv(1,(k * r))) / (k * r))**2
        column_inten.append(intensity)
        
    inten.append(np.array(column_inten))
    column_inten = []

# Plotting the intensities of the region as a density plot. Using 
# "interpolation="bilinear"" to smooth out the data.
plt.imshow(inten, interpolation ="bilinear", vmax=0.05, cmap = "gray", 
           extent=(-radius, radius, -radius, radius))
plt.colorbar()

# Adding axis labels and title.
plt.ylabel("y (m)")
plt.xlabel("x (m)")
plt.title('''Diffraction Pattern Produced by a Point
Source When Viewed Through a Telescope
''', fontsize=12)

# Adding formating to the labels and ticks of the axis.
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0), 
                     useOffset=None, useLocale=None, useMathText=None)
plt.xticks(np.linspace(-radius, radius, 5, endpoint = True))
plt.yticks(np.linspace(-radius, radius, 5, endpoint = True))

plt.show()