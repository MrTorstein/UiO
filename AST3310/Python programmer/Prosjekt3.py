from sys import exit
import fvis3 as FVis    # Visualiser
from numpy import zeros, roll, max, abs, rot90, exp, linspace, mgrid
from numpy.ma import masked_greater_equal, masked_less, masked_where

class TwoDConvection():
    """
    Class to model a two dimensional box as a cutout of the Sun.
    The class doesn't take any parameters, but you need to run the Initialise function before calling the hydrosolver.
    """
    def __init__(self):
        """ Define variables """
        
        totx            = 12 * 10 ** 6                      # Width of the box [m]
        toty            = 4 * 10 ** 6                       # Length of the box [m]
        self.m_u        = 1.660539040 * 10 ** (-27)         # Atomic mass [kg]
        self.T_sp       = 5770                              # Temperature of photosphere [K]
        self.P_sp       = 1.8 * 10 ** 8                     # Pressure of photosphere [Pa]
        self.G          = 6.674 * 10 ** (-11)               # Gravitational constant [m^3 / kg^1 s^2]
        self.M          = 1.989 * 10 ** 30                  # Mass of the Sun [kg]
        self.R          = 6.96 * 10 ** 8                    # Radius of the Sun [m]
        self.g          = - self.G * self.M / self.R ** 2   # Gravitational acceleration [m / s^2]
        self.k          = 1.38064852 * 10 ** (-23)          # Boltzmann's constant[J / K]
        self.mu         = 0.61                              # Average atomic weight
        self.nx         = 300                               # Number of boxes in horizontal direction
        self.ny         = 100                               # Number of boxes in vertical direction
        self.Dx         = totx / self.nx                    # Width of one of the boxes [m]
        self.Dy         = toty / self.ny                    # Height of one of the boxes [m]
        self.T          = zeros([self.nx, self.ny])         # Temperature matrix [K]
        self.P          = zeros([self.nx, self.ny])         # Pressure matrix [Pa]
        self.rho        = zeros([self.nx, self.ny])         # Density matrix [kg / m^3]
        self.u          = zeros([self.nx, self.ny])         # Horizontal velosity matrix [m / s]
        self.w          = zeros([self.nx, self.ny])         # Vertical velosity matrix [m / s]
        self.e          = zeros([self.nx, self.ny])         # Inner energy matrix [J / m^3]
        self.rho_dt     = zeros([self.nx, self.ny])         # Density time derivative matrix
        self.rhou_dt    = zeros([self.nx, self.ny])         # Horizontal momentum time derivative matrix
        self.rhow_dt    = zeros([self.nx, self.ny])         # Vertical momentum time derivative matrix
        self.e_dt       = zeros([self.nx, self.ny])         # Inner energy time derivative matrix

    def Initialise(self, init_par = False):
        """ Initialise temperature, pressure, density and internal energy """

        if init_par == False:
            add_gauss   = input("Turn on Gaussian distribution? [Y/N] ")

            gauss = zeros([self.nx, self.ny])
            y       = mgrid[0:self.nx, 0:self.ny][-1]
            if add_gauss == "Y" or add_gauss == "y":
                n_gauss     = int(input("Number of Gaussian distributions? [Int] "))

                A       = 30000
                x       = mgrid[0:self.nx, 0:self.ny][0]
                x_      = (self.nx - 1) / (n_gauss + 1)
                sig_x   = 50
                y_      = (self.ny - 1)
                sig_y   = 50
                for i in range(n_gauss):
                    gauss   += A * exp( - ( (x - (i + 1) * x_) ** 2 / (2 * sig_x ** 2) + (y - y_) ** 2 / (2 * sig_y ** 2) ) )
        
            Nabla           = 2 / 5 + 0.01

            self.T[:]       = self.T_sp - self.mu * self.m_u * self.g * Nabla / self.k * y * self.Dy + gauss
            self.P[:]       = self.P_sp * (self.T / self.T_sp) ** (1 / Nabla)
            self.e[:]       = self.P / ( 2 / 3 )
            self.rho[:]     = self.P * self.mu * self.m_u / ( self.k * self.T )
        
        elif type(init_par) == dict:
            self.T[:]       = rot90(init_par["T"], -1)
            self.P[:]       = rot90(init_par["P"], -1)
            self.e[:]       = rot90(init_par["e"], -1)
            self.rho[:]     = rot90(init_par["rho"], -1)
            self.u[:]       = rot90(init_par["u"], -1)
            self.w[:]       = rot90(init_par["w"], -1)
        
        else:
            print("Initial parameters have to be given as a dictionary!")
            exit()

    def _timestep(self):
        """ Calculate timestep """

        p       = 0.1
        rho_dt  = self.rho_dt
        rhou_dt = self.rhou_dt
        rhow_dt = self.rhow_dt
        e_dt    = self.e_dt
        rho     = self.rho
        u       = self.u + masked_where(self.u == 0, self.u).mask
        w       = self.w + masked_where(self.w == 0, self.w).mask
        e       = self.e

        rel_rho = max(abs( rho_dt / rho ) )                                                                 # relative change in density
        rel_u   = max(abs( (rhou_dt - u * rho_dt) / (rho * u) * masked_where(self.u != 0, self.u).mask ) )  # relative change in velocity, x direction
        rel_w   = max(abs( (rhow_dt - w * rho_dt) / (rho * w) * masked_where(self.w != 0, self.w).mask ) )  # relative change in velocity, y direction
        rel_x   = max(abs( u / self.Dx ) )                                                                  # relative change in x direction
        rel_y   = max(abs( w / self.Dy ) )                                                                  # relative change in y direction
        rel_e   = max(abs( e_dt / e ) )                                                                     # relative change in energy

        d = max([rel_rho, rel_u, rel_w, rel_x, rel_y, rel_e])
        if d < 0.01:
            d = 0.01
        if d > 10:
            d = 10
        
        self.Dt = p / d

    def _central_x(self, func):
        """ Central difference scheme in x-direction """
        return ( roll(func, 1, axis = 0) - roll(func, -1, axis = 0) ) / ( 2 * self.Dx )

    def _central_y(self, func):
        """ Central difference scheme in y-direction """
        return ( roll(func, 1, axis = 1) - roll(func, -1, axis = 1) ) / ( 2 * self.Dy )

    def _upwind_x(self, func, u):
        """ Upwind difference scheme in x-direction """
        return masked_greater_equal(u, 0).mask * ( func - roll(func, -1, axis = 0) ) / self.Dx + masked_less(u, 0).mask * ( roll(func, 1, axis = 0) - func ) / self.Dx

    def _upwind_y(self, func, w):
        """ Upwind difference scheme in y-direction """
        return masked_greater_equal(w, 0).mask * ( func - roll(func, -1, axis = 1) ) / self.Dy + masked_less(w, 0).mask * ( roll(func, 1, axis = 1) - func ) / self.Dy

    def _boundary_conditions(self):
        """ Boundary conditions for energy, density and velocity """
        self.w[:, 0]    = 0
        self.w[:, -1]   = 0

        self.u[:, 0]    = ( 4 * self.u[:, 1] - self.u[:, 2] ) / 3
        self.u[:, -1]   = ( 4 * self.u[:, -2] - self.u[:, -3] ) / 3

        self.e[:, 0]    = ( - self.e[:, 2] + 4 * self.e[:, 1] ) / ( 3 - 2 * self.mu * self.m_u * self.g / ( self.k * self.T[:, 0] ) * self.Dy )
        self.e[:, -1]   = ( - self.e[:, -3] + 4 * self.e[:, -2] ) / ( 3 + 2 * self.mu * self.m_u * self.g / ( self.k * self.T[:, -1] ) * self.Dy )
        
        self.rho[:, 0]  = self.e[:, 0] * 2 / 3 * self.mu * self.m_u / (self.k * self.T[:, 0])
        self.rho[:, -1] = self.e[:, -1] * 2 / 3 * self.mu * self.m_u / (self.k * self.T[:, -1])
    
    def hydro_solver(self):
        """ Hydrodynamic equations solver """

        rho = self.rho
        P   = self.P
        u   = self.u
        w   = self.w
        e   = self.e

        self.rho_dt  = - rho * (self._central_x(u) + self._central_y(w)) - u * self._upwind_x(rho, u) - w * self._upwind_y(rho, w)
        self.rhou_dt = rho * u * (self._upwind_x(u, u) + self._upwind_y(w, u)) - u * self._upwind_x(rho * u, u) - w * self._upwind_y(rho * u, w) - self._central_x(P)
        self.rhow_dt = rho * w * (self._upwind_y(w, w) + self._upwind_x(u, w)) - w * self._upwind_y(rho * w, w) - u * self._upwind_x(rho * w, u) - self._central_y(P) \
                       + rho * self.g
        self.e_dt = - e * (self._central_x(u) + self._central_y(w)) - u * self._upwind_x(e, u) - w * self._upwind_y(e, w) - P * (self._central_x(u) + self._central_y(w))

        self._timestep()
                
        self.rho[:]    = rho + self.rho_dt * self.Dt
        self.u[:]      = ( rho * u + self.rhou_dt * self.Dt ) / self.rho
        self.w[:]      = ( rho * w + self.rhow_dt * self.Dt ) / self.rho
        self.e[:]      = e + self.e_dt * self.Dt

        self._boundary_conditions()

        self.P[:] = 2 / 3 * self.e
        self.T[:] = self.P * self.m_u * self.mu / (self.rho * self.k)

        return self.Dt


if __name__ == '__main__':
    A = TwoDConvection()
    dictionary = {"Arrows": "Velocity"}
    #A.Initialise()
    vis = FVis.FluidVisualiser(fontsize = 30)
    #vis.save_data(150, A.hydro_solver, rho = rot90(A.rho, 1), u = rot90(A.u, 1), w = rot90(A.w, 1), e = rot90(A.e, 1), P = rot90(A.P, 1), T = rot90(A.T, 1), \
    #              sim_fps = 1.0, sim_params = dictionary)
    #
    #answer = " "
    #while type(answer) == str:
    #    answer = input("What do you want to do with the animation? [Visualise/Save] ")
    #    if answer == "visualise" or answer == "Visualise":
    #        vis.animate_2D("T", save = False)
    #        answer = True
    #    elif answer == "save" or answer == "Save":
    #        vis.animate_2D("T", save = True)
    #        answer = True
    #    elif answer == "exit" or answer == "Exit":
    #        exit()
    #    elif answer == "pass" or answer == "Pass":
    #        answer = True
    #    else:
    #        print("'%s' is not one of the options."%answer)
    #
    arrs, params = vis.get_last_data("FVis_output_2020-05-20_22-54")
    A.Initialise(arrs)
    vis.save_data(20, A.hydro_solver, rho = rot90(A.rho, 1), u = rot90(A.u, 1), w = rot90(A.w, 1), e = rot90(A.e, 1), P = rot90(A.P, 1), T = rot90(A.T, 1), \
                  sim_fps = 1.0, appendMode = True, folder = "FVis_output_2020-05-20_22-54")
    vis.animate_2D("T", save = True, folder = "FVis_output_2020-05-20_22-54")