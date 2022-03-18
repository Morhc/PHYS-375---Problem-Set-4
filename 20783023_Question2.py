import os
import numpy as np
from matplotlib import pyplot as plt, rc

rc('text', usetex=True)

def compute_rho_and_mass(r,spiral=1):
    """Compute the value of rho and mass using basic ODE solving.
    INPUTS:
        r - The radius to compute rho and mass at.
        spiral - The initial conditions for rho. I accidentally found out by messing with it
                 that if you set it > 1 then you get the spiral bit in Bonnor 1956!
    OUTPUTS:
        rho - The computed densities.
        mass - The computed masses.
    """

    rho = np.zeros(r.size, dtype=float)
    drho = np.zeros(r.size, dtype=float)
    mass = np.zeros(r.size, dtype=float)

    #initial conditions
    rho[0] = spiral
    drho[0] = 0

    dr = r[1] - r[0] #get the r-spacing
    for k in range(1, rho.size):
        #this is from equation (10) in my solutions
        ddrho = -(4*np.pi*rho[k-1]**2 + ((2./r[k]) - (1.0/rho[k-1])*drho[k-1])*drho[k-1])
        #this is just basic ODE solving technique
        drho[k] = drho[k-1] + ddrho * dr
        rho[k] = rho[k-1] + drho[k] * dr

        #this is from equation (8) in my solutions
        mass[k] = mass[k-1] + 4*np.pi*r[k]**2*rho[k]*dr



    return rho, mass

def part_b(savepath=""):
    """Solve the reduced dM/dr and dp/dr equations found.
    INPUTS:
        savepath - Path to save the plot to. If nothing is given, it will just display.
    OUTPUTS:
        Outputs two plots, one for the rho as a function of radius and the other mass as a function of radius.
    """


    r = np.linspace(0, 5, 1000)

    rho, mass = compute_rho_and_mass(r)

    #exclude the r=0 solution
    plt.plot(r[1:], mass[1:])

    plt.title(r'$\tilde{M}$($\tilde{r}$) vs $\tilde{r}$')
    plt.ylabel(r'$\tilde{M}$($\tilde{r}$)')
    plt.xlabel(r'$\tilde{r}$')
    plt.xlim(r[0], r[-1])
    plt.ylim(0, max(mass)*1.1)

    if savepath == "":
        plt.show()
    else:
        plt.savefig(savepath.replace('.png', '_mass.png'))

    plt.close('all')

    #exclude the r=0 solution
    plt.plot(r[1:], rho[1:])

    plt.title(r'$\tilde{\rho}$($\tilde{r}$) vs $\tilde{r}$')
    plt.ylabel(r'$\tilde{\rho}$($\tilde{r}$)')
    plt.xlabel(r'$\tilde{r}$')
    plt.xlim(r[0], r[-1])

    if savepath == "":
        plt.show()
    else:
        plt.savefig(savepath.replace('.png', '_rho.png'))

    plt.close('all')

def part_c(savepath=""):
    """Regenerate at least the important part of Bonnor 1956 Figure 1
    INPUTS:
        savepath - Path to save the plot to. If nothing is given, it will just display.
    OUTPUTS:
        Outputs the plot from Bonnor 1956, which is surface pressure vs radius
    """

    #constants
    mu = 2.4
    mp = 1.67e-27
    G = 6.67e-11
    k = 1.38e-23
    AU = 1.496e11

    #our params
    T = 10
    Msun = 1.99e30

    r_tilde = np.linspace(0,5,100000)
    rho_tilde, mass_tilde = compute_rho_and_mass(r_tilde, spiral=30000)
    #cutting the first index to avoid an error
    rho_tilde, mass_tilde = rho_tilde[1:], mass_tilde[1:]

    #equation (11) in my solutions
    Ps = np.power(k*T/(mu*mp), 4) * 1./(G**3) * np.power(mass_tilde/Msun, 2) * rho_tilde

    #converts from unscaled to meters to AU, based on equation (5) in my solutions
    r = r_tilde[1:]*G*mu*mp/(k*T)*(Msun/mass_tilde) / AU

    plt.plot(r, Ps)
    #plt.xlim(0, 1e5) #if no spiral
    plt.xlim(0, 4e4) #if spiral
    plt.ylim(0, max(Ps)*1.1)
    plt.title(r'Surface Pressure $P_s$ vs Radius $r$')
    plt.ylabel(r'$P_s$ (N/m$^2$)')
    plt.xlabel(r'$r$ (AU)')
    if savepath == "":
        plt.show()
    else:
        plt.savefig(savepath)

    plt.close('all')


def part_d(savepath=""):
    """Plots the Jeans mass as a function of radius
    INPUTS:
        savepath - Path to save the plot to. If nothing is given, it will just display.
    OUTPUTS:
        Outputs a plot of Jeans mass as a function of radius
    """
    mu = 2.4
    mp = 1.67e-27
    G = 6.67e-11
    k = 1.38e-23
    AU = 1.496e11

    #our params
    T = 10
    Msun = 1.99e30

    r_tilde = np.linspace(0,5,100000)

    _, mass_tilde = compute_rho_and_mass(r_tilde, spiral=1)
    #cutting the first index to avoid an error
    mass_tilde = mass_tilde[1:]

    #converts from unscaled to meters to AU, based on equation (5) in my solutions
    r = r_tilde[1:]*G*mu*mp/(k*T)*(Msun/mass_tilde) / AU

    #from equations (12) and (14) in my solutions
    mean_rho = 3/(4*np.pi) * Msun/np.power(r, 3)
    Mj = 0.2*Msun*np.power(T/10, 3/2)*np.power(mean_rho/3e-15, -1/2)

    plt.plot(r, Mj)
    plt.xlim(0, 1e5) #without spiral
    plt.ylim(0, 1e15) #without spiral
    plt.title(r'Jeans Mass $M_J$ vs Radius $r$')
    plt.ylabel(r'$M_J$ (kg)')
    plt.xlabel(r'$r$ (AU)')
    if savepath == "":
        plt.show()
    else:
        plt.savefig(savepath)

    plt.close('all')


def main():

    here = os.path.dirname(os.path.realpath(__file__))

    b_path = os.path.join(here, 'PS4-Q2b.png')
    part_b(b_path)

    c_path = os.path.join(here, 'PS4-Q2c.png')
    part_c(c_path)

    d_path = os.path.join(here, 'PS4-Q2d.png')
    part_d(d_path)

if __name__ == "__main__":
    main()
