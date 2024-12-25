#Embedded file name: AST2000SolarSystem.py
"""
This is a module for use with the project work in course AST2000 (previously AST1100)
at the University of Oslo. It is designed to work with 3D visualization MCast.
There have been many contributors. Thanks to all of them!
Currently, the responsible person is Robert Hagala (robert.hagala@astro.uio.no)
"""
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import sys
import os.path
import time
import traceback
import logging
import scipy.integrate as integrate
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy import interpolate
from math import pi, floor, ceil
from numpy import cos, sin, sqrt, tan, arctan, log, arctan2, exp, arccos
from PIL import Image
from lxml import etree
from six.moves import range
from six.moves import zip
from six.moves import input
c = 299792.458
years_to_seconds = 60 * 60 * 24 * 365.24
seconds_to_years = 1.0 / years_to_seconds
AU_to_m = 149597870700L
AU_to_km = 149597870.7
m_to_AU = 1 / AU_to_m
km_to_AU = 1 / AU_to_km
c_AU = c * km_to_AU
AU_pr_year_to_m_pr_s = 4740.57172
sun_mass_to_kg = 1.9891e+30
H2_g_pr_mol = 2.01588
N_a = 6.02214179e+23
mass_1_particle = H2_g_pr_mol / N_a * 0.001
__path__ = os.path.abspath(os.path.dirname(__file__))

def uniform(start, stop):
    return start + random.random() * (stop - start)


def randint(start, stop):
    return int(uniform(start, stop))


class AST2000SolarSystem():

    class Moons():
        number_of_moons = 0
        incli = 0.0
        radius = []
        phi0 = []
        omegao = []
        omegam = []
        distance = []

        def __init__(selff, self, i):
            number_of_moons = randint(0, max(1, int(self.radius[i] / 2500)))
            if number_of_moons > 98:
                number_of_moons = 98
            selff.number_of_moons = number_of_moons
            selff.incli = uniform(-0.1, 0.1)
            selff.radius = np.zeros(number_of_moons)
            selff.omegao = np.zeros(number_of_moons)
            selff.omegam = np.zeros(number_of_moons)
            selff.phi0 = np.zeros(number_of_moons)
            selff.distance = np.zeros(number_of_moons)
            for j in range(number_of_moons):
                selff.radius[j] = uniform(400.0, min(4000.0, 0.4 * self.radius[i]))
                selff.phi0[j] = uniform(0.0, 2.0 * pi)
                selff.distance[j] = uniform(3.0, 100.0) * self.radius[i] * km_to_AU
                omegao = np.sqrt(self.G * self.mass[i] / selff.distance[j] ** 3)
                selff.omegao[j] = omegao
                selff.omegam[j] = uniform(omegao, 2 * omegao)

        def findPosRot(selff, t):
            phi = selff.phi0 + selff.omegam * t
            theta = selff.phi0 + selff.omegao * t
            x = selff.distance * np.cos(theta)
            y = selff.distance * np.sin(theta)
            z = x * selff.incli
            return (np.transpose(np.array([x, y, z])), phi)

    a = []
    e = []
    radius = []
    omega = []
    psi = []
    mass = []
    period = []
    x0 = []
    y0 = []
    vx0 = []
    vx0 = []
    mu = []
    rho0 = []
    label = []
    moons = []
    muList = [34.5,
     27.555556,
     29.0,
     31.0,
     29.25,
     22.0,
     25.2,
     27.142857,
     29.4,
     30.8,
     26.6,
     27.666667,
     23.333333,
     29.636364,
     24.5,
     27.090909,
     24.857143,
     26.6,
     22.0,
     26.0,
     27.142857,
     31.333333,
     26.0,
     28.0,
     28.571429,
     29.555556,
     22.666667,
     27.25,
     31.333333,
     25.2,
     28.571429,
     34.666667,
     29.454545,
     21.333333,
     30.666667,
     30.666667,
     29.0,
     23.2,
     29.454545,
     29.4,
     27.090909,
     31.428571,
     27.2,
     29.454545,
     30.5,
     28.5,
     28.0,
     27.090909,
     24.333333,
     30.666667,
     27.090909,
     28.285714,
     28.5,
     28.5,
     29.0,
     24.857143,
     26.0,
     27.090909,
     29.4,
     39.2,
     29.555556,
     25.5,
     31.0,
     31.5,
     28.2,
     26.6,
     27.0,
     21.666667,
     28.5,
     29.333333,
     28.5,
     28.181818,
     27.142857,
     33.5,
     29.142857,
     24.25,
     30.333333,
     28.666667,
     31.2,
     27.555556,
     27.090909,
     33.0,
     24.666667,
     27.0,
     31.5,
     27.555556,
     25.0,
     31.2,
     30.571429,
     24.666667,
     26.666667,
     28.5,
     28.181818,
     29.333333,
     30.444444,
     25.142857,
     28.5,
     29.454545,
     27.090909,
     27.75]
    seed = 0
    G = 4 * pi * pi
    eps = 0.25
    number_of_planets = 0
    temperature = 0
    star_radius = 0
    star_mass = 0
    sFile = 's.bin'
    RSV = [-0.4, 1.2]
    RSA = [102.5, 4.1]
    dpdt1 = 0.0
    numBox = 0.0
    partPerSec = 0.0
    total_fuel_calculated = 0.0
    fuel_left = 0.0
    inittime = 0.0
    fuelPrint = False

    def __init__(self, seed, resolution = 64, has_moons = True, fuel_out = False, debug = False, data_path = None):
        self.seed = seed
        self.debug = debug
        self.resolution = resolution
        self.prevdir = [0, 1, 0]
        self.fuelPrint = fuel_out
        self.mass_sat = 1100
        self.engine_checked = False
        self.launch_accomplished = False
        if data_path:
            self.path = data_path.strip()
        else:
            try:
                pfl = open('path.txt')
                self.path = pfl.read().strip()
                pfl.close()
            except:
                self.path = ''

        if len(self.path) > 0:
            slash = '/'
            if '\\' in list(self.path):
                slash = '\\'
            if self.path[-1] != slash:
                self.path = self.path + slash
        if self.debug:
            print('Debug prints enabled.')
            print("path set to '%s'" % self.path)
        try:
            random.seed(seed, version=1)
        except:
            random.seed(seed)

        if seed < 100:
            forceMu = self.muList[seed]
        else:
            ii = int(str(seed)[-2:])
            forceMu = self.muList[ii]
        if os.path.exists(self.sFile):
            statusFile = open(self.sFile, 'rb')
            s = np.load(statusFile)
            statusFile.close()
            if s[0][0] != self.seed:
                print('The random seed has changed.')
                statusFile = open(self.sFile, 'wb')
                np.save(statusFile, [[self.seed], [0], []])
                statusFile.close()
        else:
            statusFile = open(self.sFile, 'wb')
            np.save(statusFile, [[self.seed], [0], []])
            statusFile.close()
        self.number_of_planets = randint(7, 9)
        self.temperature = uniform(3000, 12000)
        self.star_mass = np.amax([self.temperature ** 2 / 33640000.0 * random.gauss(1, 0.1), 0.1])
        self.star_radius = abs(self.star_mass ** 0.724 * 700000.0 * random.gauss(1, 0.1))
        self.a = np.zeros(self.number_of_planets)
        self.e = np.zeros(self.number_of_planets)
        self.mass = np.zeros(self.number_of_planets)
        self.radius = np.zeros(self.number_of_planets)
        self.omega = np.zeros(self.number_of_planets)
        self.psi = np.zeros(self.number_of_planets)
        self.period = np.zeros(self.number_of_planets)
        self.x0 = np.zeros(self.number_of_planets)
        self.y0 = np.zeros(self.number_of_planets)
        self.vx0 = np.zeros(self.number_of_planets)
        self.vy0 = np.zeros(self.number_of_planets)
        self.mu = np.zeros(self.number_of_planets)
        self.rho0 = np.zeros(self.number_of_planets)
        pi = np.pi
        for i in range(self.number_of_planets):
            if i == 0:
                self.label.append('Home')
                self.a[i] = self.star_radius * km_to_AU / 2 * (self.temperature / uniform(320, 350)) ** 2
                self.e[i] = abs(random.gauss(0, 0.02))
                self.omega[i] = 0
                self.psi[i] = 0
                densi = uniform(4500.0, 6500)
                self.period[i] = uniform(0.8, 1.3)
                self.mass[i] = uniform(1.0, 11) * 1e-06
                self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                self.mu[i] = 28.3
                self.rho0[i] = uniform(1, 1.5)
                R0 = self.a[i] * (1 - self.e[i] * self.e[i]) / (1 - self.e[i] * cos(self.omega[i] - self.psi[i]))
                v = sqrt(self.star_mass * self.G * (2 / R0 - 1 / self.a[i]))
                h = sqrt(self.star_mass * self.G * self.a[i] * (1 - self.e[i] * self.e[i]))
                vPerp = h / R0
                vParr = 0
                self.vx0[i] = -vPerp * sin(self.omega[i]) + vParr * cos(self.omega[i])
                self.vy0[i] = vPerp * cos(self.omega[i]) + vParr * sin(self.omega[i])
                self.x0[i] = R0 * cos(self.omega[i])
                self.y0[i] = R0 * sin(self.omega[i])
            elif i == 1:
                self.label.append('Habitable')
                self.a[i] = self.star_radius * km_to_AU / 2 * (self.temperature / uniform(270, 300)) ** 2
                self.e[i] = abs(random.gauss(0, 0.02))
                self.omega[i] = uniform(0, 2 * pi)
                self.psi[i] = uniform(0, 2 * pi)
                self.period[i] = uniform(0.6, 2.1)
                densi = uniform(3000.0, 6500)
                self.mass[i] = 10 ** uniform(-5, -7)
                self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                self.mu[i] = uniform(20, 35)
                self.rho0[i] = uniform(0.5, 10)
                R0 = self.a[i] * (1 - self.e[i] * self.e[i]) / (1.0 - self.e[i] * cos(self.omega[i] - self.psi[i]))
                v = sqrt(self.star_mass * self.G * (2 / R0 - 1 / self.a[i]))
                h = sqrt(self.star_mass * self.G * self.a[i] * (1 - self.e[i] * self.e[i]))
                vPerp = h / R0
                if self.omega[i] > self.psi[i]:
                    if self.omega[i] > self.psi[i] + pi:
                        vParr = sqrt(v ** 2 - vPerp ** 2)
                    else:
                        vParr = -sqrt(v ** 2 - vPerp ** 2)
                elif self.omega[i] < self.psi[i] - pi:
                    vParr = -sqrt(v ** 2 - vPerp ** 2)
                else:
                    vParr = sqrt(v ** 2 - vPerp ** 2)
                self.vx0[i] = -vPerp * sin(self.omega[i]) + vParr * cos(self.omega[i])
                self.vy0[i] = vPerp * cos(self.omega[i]) + vParr * sin(self.omega[i])
                self.x0[i] = R0 * cos(self.omega[i])
                self.y0[i] = R0 * sin(self.omega[i])
            else:
                self.a[i] = uniform(0.2 * self.a[0], 8 * self.a[1])
                ok = True
                count = 1.0
                while ok:
                    ok = False
                    count = count + 0.01
                    for k in range(i):
                        while abs(self.a[i] - self.a[k]) < self.eps * (self.a[i] + self.a[k]) * 0.5:
                            self.a[i] = uniform(0.2 * self.a[0], 8 * count * self.a[1])
                            ok = True

                self.e[i] = abs(random.gauss(0, 0.05))
                self.omega[i] = uniform(0, 2 * pi)
                self.psi[i] = uniform(0, 2 * pi)
                self.period[i] = uniform(0.6, 40)
                R0 = self.a[i] * (1 - self.e[i] * self.e[i]) / (1.0 - self.e[i] * cos(self.omega[i] - self.psi[i]))
                v = sqrt(self.star_mass * self.G * (2 / R0 - 1 / self.a[i]))
                h = sqrt(self.star_mass * self.G * self.a[i] * (1 - self.e[i] * self.e[i]))
                vPerp = h / R0
                self.mu[i] = uniform(10, 30)
                if self.omega[i] > self.psi[i]:
                    if self.omega[i] > self.psi[i] + pi:
                        vParr = sqrt(v ** 2 - vPerp ** 2)
                    else:
                        vParr = -sqrt(v ** 2 - vPerp ** 2)
                elif self.omega[i] < self.psi[i] - pi:
                    vParr = -sqrt(v ** 2 - vPerp ** 2)
                else:
                    vParr = sqrt(v ** 2 - vPerp ** 2)
                self.vx0[i] = -vPerp * sin(self.omega[i]) + vParr * cos(self.omega[i])
                self.vy0[i] = vPerp * cos(self.omega[i]) + vParr * sin(self.omega[i])
                self.x0[i] = R0 * cos(self.omega[i])
                self.y0[i] = R0 * sin(self.omega[i])
                r_vect = (self.x0[i] * AU_to_m, self.y0[i] * AU_to_m)
                T_planet = self.temperature * np.sqrt(self.star_radius * 1000 / (2 * np.linalg.norm(r_vect)))
                earth_radius = 6371 * 1000
                if T_planet < 170:
                    R = (self.temperature / T_planet) ** 2 * (self.star_radius / 2)
                    R_max = (self.temperature / 170) ** 2 * (self.star_radius / 2) * 4
                    if R < R_max:
                        prob = uniform(0, 1)
                        if prob < 0.3:
                            self.label.append('StoDw')
                            self.mass[i] = 10 ** uniform(-7.5, -8.5)
                            densi = uniform(2000.0, 6500)
                            self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                            self.rho0[i] = uniform(1, 1.5)
                        elif prob < 0.31:
                            self.label.append('GasDw')
                            self.mass[i] = 10 ** uniform(-5, -6.5)
                            densi = uniform(400.0, 2000)
                            self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                            self.rho0[i] = uniform(18, 25)
                        else:
                            self.label.append('GasGi')
                            self.mass[i] = 10 ** uniform(np.log10(self.star_mass * 0.01), -5)
                            densi = uniform(500.0, 1000) + 1000 * max(np.log10(3500 * self.mass[i]), 0)
                            self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                            self.rho0[i] = uniform(18, 25)
                    else:
                        prob = uniform(0, 1)
                        if prob < 0.4:
                            self.label.append('IceGi')
                            self.mass[i] = 10 ** uniform(-4, -6)
                            densi = uniform(900.0, 2500)
                            self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                            self.rho0[i] = uniform(1.3, 3)
                        else:
                            self.label.append('IceDw')
                            self.mass[i] = 10 ** uniform(-6.5, -8.5)
                            densi = uniform(1500.0, 2500)
                            self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                            self.rho0[i] = uniform(1, 1.5)
                else:
                    prob = uniform(0, 1)
                    if prob < 0.96:
                        self.label.append('HotRoc')
                        self.mass[i] = 10 ** uniform(-5.5, -7.5)
                        densi = uniform(4000.0, 6500)
                        self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                        self.rho0[i] = uniform(1, 1.5)
                    else:
                        self.label.append('HotGas')
                        self.mass[i] = 10 ** uniform(np.log10(self.star_mass * 0.01), -4.5)
                        densi = uniform(600.0, 2500)
                        self.radius[i] = (3 * self.mass[i] * sun_mass_to_kg / (4 * pi * densi)) ** 0.3333333333333333 * 0.001
                        self.rho0[i] = uniform(18, 25)
                if self.radius[i] > 15000:
                    self.period[i] = uniform(0.2, 0.9)
        #print self.label
        if has_moons:
            for i in range(self.number_of_planets):
                self.moons.append(self.Moons(self, i))

        self.mu[:] = forceMu
        self.RSV = [uniform(-2, 2), uniform(-2, 2)]
        a1 = uniform(0, 360)
        a2 = (a1 + uniform(40, 140)) % 360
        self.RSA = [a2, a1]

    def print_info(self):
        print('Information about solar system with seed', self.seed, '\n')
        print('Number of planets:', self.number_of_planets)
        print('Star surface temperature:', self.temperature, 'K')
        print('Star radius:', self.star_radius, 'km')
        print('Star mass:', self.star_mass, 'solar masses')
        print('\nIndividual planet information. Masses in units of M_sol, radius in km,\nmean molecular weight mu in multiples of proton mass,\natmospheric density rho_0 in kg/m^3, day length in earth days.\nPlanet |      Mass      |     Radius     |     rho_0      |   Day length   |')
        for pla in range(self.number_of_planets):
            print('%6d |%15.9g |%15.7f |%15.9f |%15.8f |' % (pla,
             self.mass[pla],
             self.radius[pla],
             self.rho0[pla],
             self.period[pla]))

        print('\nIndividual planet initial positions (in AU) and velocoties (in AU/year).\nPlnt |       x_0       |       y_0       |       vx0       |       vy0       |')
        for pla in range(self.number_of_planets):
            print('%4d |%16.10f |%16.10f |%16.10f |%16.10f |' % (pla,
             self.x0[pla],
             self.y0[pla],
             self.vx0[pla],
             self.vy0[pla]))

    def findPositions(self, time):
        positions = np.zeros(shape=(2, self.number_of_planets))
        for i in range(self.number_of_planets):
            th = self.theta(time, i)
            R = self.a[i] * (1 - self.e[i] ** 2) / (1 - self.e[i] * cos(th - self.psi[i]))
            positions[:, i] = [R * cos(th), R * sin(th)]

        return positions

    def findPositionsArray(self, times):
        positions = np.zeros((2, self.number_of_planets, len(times)))
        for i in range(self.number_of_planets):
            sys.stdout.write('.')
            sys.stdout.flush()
            th = self.thetaArray(times, i)
            R = self.a[i] * (1 - self.e[i] ** 2) / (1 - self.e[i] * np.cos(th - self.psi[i]))
            positions[0, i] = R * np.cos(th)
            positions[1, i] = R * np.sin(th)

        sys.stdout.write(' Done.\n')
        sys.stdout.flush()
        return positions

    def theta(self, time, i):
        dx = 1e-05
        tolerance = 1e-07
        maxIterations = 50
        h = sqrt(self.star_mass * self.G * self.a[i] * (1 - self.e[i] * self.e[i]))
        sqt = sqrt(self.star_mass * self.G / self.a[i] ** 3)
        iterations = 1
        x0 = time * sqt + self.omega[i]
        x00 = x0
        ff = self.f(x0 - self.psi[i], self.e[i])
        ffo = self.f(self.omega[i] - self.psi[i], self.e[i])
        try:
            while abs(ff - time * h / (self.a[i] * self.a[i] * (1 - self.e[i] * self.e[i])) - ffo) > tolerance and iterations < maxIterations:
                iterations = iterations + 1
                dfdx = (self.f(x0 + dx - self.psi[i], self.e[i]) - ff) / dx
                x0 = x0 - (ff - h * time / (self.a[i] * self.a[i] * (1 - self.e[i] * self.e[i])) - ffo) / dfdx
                ff = self.f(x0 - self.psi[i], self.e[i])
                if iterations == maxIterations:
                    print(dfdx, ff, ffo)
                    raise StopIteration()
                    exit()

        except StopIteration as error:
            traceback.print_exc(limit=1, file=sys.stdout)
            print("\nError: Newton's method did not converge after %i iterations." % maxIterations)
            print('time          = ', time)
            print('planet number = ', i)

        return x0

    def thetaArray(self, times, i):
        dx = 1e-05
        tolerance = 1e-07
        maxIterations = 50
        h = sqrt(self.star_mass * self.G * self.a[i] * (1 - self.e[i] * self.e[i]))
        sqt = sqrt(self.star_mass * self.G / self.a[i] ** 3)
        iterations = 1
        x0 = times * sqt + self.omega[i]
        x00 = x0
        ff = self.f(x0 - self.psi[i], self.e[i])
        ffo = self.f(self.omega[i] - self.psi[i], self.e[i])
        lenmask = 1
        try:
            while lenmask > 0 and iterations < maxIterations:
                mask = abs(ff - times * h / (self.a[i] * self.a[i] * (1 - self.e[i] * self.e[i])) - ffo) > tolerance
                lenmask = np.sum(mask)
                iterations = iterations + 1
                dfdx = (self.f(x0[mask] + dx - self.psi[i], self.e[i]) - ff[mask]) / dx
                x0[mask] = x0[mask] - (ff[mask] - h * times[mask] / (self.a[i] * self.a[i] * (1 - self.e[i] * self.e[i])) - ffo) / dfdx
                ff[mask] = self.f(x0[mask] - self.psi[i], self.e[i])
                if iterations == maxIterations:
                    raise StopIteration()
                    exit()

        except StopIteration as error:
            traceback.print_exc(limit=1, file=sys.stdout)
            print("\nError: Newton's method did not converge after %i iterations." % maxIterations)
            print('time          = ', time)
            print('planet number = ', i)

        return x0

    def getAcceleration(self, r, m):
        R = np.linalg.norm(r)
        return -self.G * m / R ** 3 * r

    def f(self, x, e):
        sq = np.sqrt(1.0 - e * e)
        y = (e + 1) * np.tan(x / 2.0) / sq
        n = np.floor((x + pi) / (2 * pi))
        return -e * np.sin(x) / (e * np.cos(x) - 1) + 2 * (np.arctan(y) + n * pi) / sq

    def check_planet_positions(self, planet_pos, T, N, writeFile = True):
		print ("herpdiderp")
		pos = planet_pos
		TN = len(planet_pos[0, 0, :])
		print(TN)
		home_period = sqrt(self.a[0] ** 3 / self.star_mass)
		if abs(TN - T * N) > 1:
		    print('The size of T*N does not match the size of (the last index of) the position array.')
		    return 3
		if T < 20 * home_period:
		    print('Need to include at least 20 orbits of the home planet!')
		    return 2
		if TN < 1000:
		    print('You should definitely include more than 1000 time steps total.')
		    return 1
		timeslong = np.linspace(0, T, TN)
		i = 0
		while timeslong[i] < 20 * home_period:
		    i += 1

		times = np.copy(timeslong[:i])
		TN = len(times)
		pos = pos[:, :, :TN]
		if TN < 300:
			print('You should have more time steps. \nYou might be simulating for too long time as well.')
			return 1
		skips = TN // 75
		pos = pos[:, :, ::skips]
		times = times[::skips]
		TN = len(pos[0, 0, :])
		print('TN', TN)
		print(len(times))
		max_deviation = 0.0
		planet_deviation = 0
		index_deviation = 0
		deviation_tol = 1e10
		print('Checking the deviations of the planets. \nThis might take a couple of minutes!...')
		poses = np.zeros((2, self.number_of_planets, TN))
		toolbar_width = 75
		toolbar_fin = 0
		sys.stdout.write('[%s ]' % (' ' * toolbar_width))
		sys.stdout.flush()
		sys.stdout.write('\x08' * (toolbar_width + 2))
		for i in range(TN):
			t = times[i]
			pos_computed = self.findPositions(t)
			poses[:, :, i] = pos_computed
			if toolbar_fin < 76:
			    sys.stdout.write('-')
			    sys.stdout.flush()
			    toolbar_fin += 1
			for j in xrange(self.number_of_planets):
			    try:
			        calcrel = abs(np.linalg.norm(pos[:, j, i] - pos_computed[:, j]) / np.linalg.norm(pos_computed[:, j]))
			    except:
			        print(np.shape(pos), np.shape(pos_computed), j, i)
		
			    if calcrel > max_deviation:
					print ("DUUUUUUUU", i, j)
					max_deviation = calcrel
					planet_deviation = j
					index_deviation = i

			if t > 20 * home_period:
			    break

		sys.stdout.write('\n')
		print('The biggest relative deviation was detected at planet %i,\nwhich drifted %.4f %% from its actual position, "index: %i"' % (planet_deviation, 100 * max_deviation, index_deviation))
		if max_deviation * 100 > deviation_tol:
			print('Your planets are not where they should be after 20 orbits of planet 0.')
			print('Check your program for flaws or experiment with your time step for more precise trajectories.')
			return 5
		print('Your planet trajectories after %i years are satisfyingly calculated.' % np.rint(t))
		print()
		if writeFile:
			print('Calculating exact binary array file planet_positions.npy. Congratulations!')
			plapo = self.findPositionsArray(timeslong)
			print("[planet_positions, times] = np.load('planet_positions.npy') #how to load")
			print('planet_positions will have indexes [xy, planet, time]')
			outFile = open('planet_positions.npy', 'wb')
			np.save(outFile, [plapo, timeslong])
			outFile.close()
		return 0

    def send_satellite(self, fileName = 'satCommands.txt', planet_pos = False, dt = 0.001):
        if self.engine_checked == False:
            print('You need a working engine before you can send out the satellite')
            print('Activate the function engine_settings before proceding ')
            sys.exit()
        if self.launch_accomplished == False:
            print('You need to finish the launching sequence before continuing.')
            print('Make sure your final position after launch is correct.')
            sys.exit()

        def satelliteRHS(xv, t):
            position = self.findPositions(t)
            acceleration = np.zeros(2)
            for j in range(self.number_of_planets):
                rij = xv[0:2] - position[:, j]
                acceleration += self.getAcceleration(rij, self.mass[j])
                if np.linalg.norm(rij) < 0.99 * self.radius[j] * km_to_AU:
                    print('CRASH!! You collided with planet nr %i at time %g. Too bad :(' % (j, t))
                    sys.exit()

            if np.linalg.norm(xv[0:2]) < 0.99 * self.star_radius * km_to_AU:
                print('CRASH!! You collided with with the sun at time %g. Too bad :(' % t)
                sys.exit()
            acceleration += self.getAcceleration(xv[0:2], self.star_mass)
            dxv = np.zeros(4)
            dxv[0] = xv[2]
            dxv[1] = xv[3]
            dxv[2] = acceleration[0]
            dxv[3] = acceleration[1]
            return dxv

        statusFile = open(self.sFile, 'rb')
        status = np.load(statusFile)
        orientBool = status[1][0]
        statusFile.close()
        eventType, eventTime, eventCommand = self.read_satellite_input_file(fileName, dt)
        xv = np.zeros(shape=(1, 4))
        nIntegrate = 1000

        def intforward(t1, t2, xv):
            tarr = np.linspace(t1, t2, nIntegrate)
            return (tarr, integrate.odeint(satelliteRHS, xv, tarr, mxstep=5000))

        sameTime = np.zeros(len(eventTime))
        ii = 0
        while ii < len(eventTime) - 1:
            if abs(eventTime[ii] - eventTime[ii + 1]) < 1e-10:
                sameTime[ii] = 1
                ii += 2
            else:
                ii += 1

        tFullArray = np.zeros((len(eventTime) - 1 - int(sum(sameTime))) * nIntegrate)
        pos = np.zeros(shape=(2, self.number_of_planets + 1, (len(eventTime) - 1 - int(sum(sameTime))) * nIntegrate))
        lenthBar = 40
        ii = -1
        orientNumber = 0
        videoTimes = []
        videoPos = []
        videoCommand = []
        for i in range(len(eventTime)):
            if eventType[i] == 0:
                self.dump_to_xml(eventTime[i], xv[-1, 0:2], [eventCommand[i][0],
                 eventCommand[i][1],
                 0,
                 0,
                 1], i)
                print('Picture taken! selfie_%04d.xml' % i)
            elif eventType[i] == 90:
                videoTimes.append(eventTime[i])
                videoPos.append(xv[-1, 0:2])
                videoCommand.append(eventCommand[i])
            elif eventType[i] == 1:
                fueluse = self.mass_needed(xv[-1, 2:4], np.array(eventCommand[i]))
                if self.fuel_left < -0.15 * self.total_fuel_calculated:
                    print('You ran out of fuel. Too bad!')
                    sys.exit()
                xv[-1, 2:4] += np.array(eventCommand[i])
            elif eventType[i] == -1:
                xv[-1] = np.array(eventCommand[i])
            elif eventType[i] == 5:
                if orientBool:
                    print('AU_to_matic orientation.', end=' ')
                    orientFile = open('orient' + str(orientNumber) + '.npy', 'wb')
                    orientPos = self.findPositions(eventTime[i])
                    np.save(orientFile, [[eventTime[i]], orientPos, xv[-1, :]])
                    orientFile.close()
                    orientNumber = orientNumber + 1
                    print('Time:', eventTime[i])
                    print('Pos: (', xv[(-1, 0)], ',', xv[(-1, 1)], ')')
                    print('Vel: (', xv[(-1, 2)], ',', xv[(-1, 3)], ')')
                else:
                    print('Manual orientation.')
                    inputSat = xv[-1, :]
                    if self.orientForRealz(eventTime[i], inputSat):
                        orientFile = open('orient' + str(orientNumber) + '.npy', 'wb')
                        orientPos = self.findPositions(eventTime[i])
                        np.save(orientFile, [[eventTime[i]], orientPos, xv[-1, :]])
                        orientFile.close()
                        orientNumber = orientNumber + 1
                        print('Time:', eventTime[i])
                        print('Pos: (', xv[(-1, 0)], ',', xv[(-1, 1)], ')')
                        print('Vel: (', xv[(-1, 2)], ',', xv[(-1, 3)], ')')
                    else:
                        print('Orientation NOT successful. Exiting.')
                        sys.exit()
            if i < len(eventTime) - 1 - int(sum(sameTime)):
                ii += 1 + int(sameTime[ii + 1])
                t, xv = intforward(eventTime[ii], eventTime[ii + 1], xv[-1, :])
                pos[:, -1, i * nIntegrate:(i + 1) * nIntegrate] = xv[:, 0:2].transpose()
                tFullArray[i * nIntegrate:(i + 1) * nIntegrate] = t
            if len(eventTime) > 1:
                progress = float(i) / (len(eventTime) - 1)
            else:
                progress = 1.0
            finished = len(eventTime) - 1

        if len(videoTimes) > 0:
            if len(videoTimes) > 2:
                print('You can only make 1 video at a time, try again with 2 video calls, 1 to start the video and 1 to end it')
                sys.exit()
            sat_pos = interpolate.interp1d(tFullArray, pos[:, -1, :])
            if len(videoCommand[0]) + len(videoCommand[1]) == 4:
                theta0, phi0, theta1, phi1 = (videoCommand[0][0],
                 videoCommand[0][1],
                 videoCommand[1][0],
                 videoCommand[1][1])
                t0, t1 = videoTimes[0], videoTimes[1]
                nFrames = 10000
                video_time_array = np.linspace(t0, t1, nFrames)
                videoPos = [ sat_pos(i) for i in video_time_array ]
                phi = np.linspace(phi0, phi1, nFrames)
                theta = np.linspace(theta0, theta1, nFrames)
                videoAngles = [theta, phi]
            elif len(videoCommand[0] + videoCommand[1]) == 2:
                t0, t1 = videoTimes[0], videoTimes[1]
                nFrames = 10000
                video_time_array = np.linspace(t0, t1, nFrames)
                videoPos = [ sat_pos(i) for i in video_time_array ]
                camDirec = [ sat_pos(i) - self.findPositions(i)[:, videoCommand[0][0]] for i in video_time_array ]
                phi = np.zeros(nFrames)
                theta = np.zeros(nFrames)
                for i in range(nFrames):
                    phi[i] = np.arctan2(camDirec[i][1], camDirec[i][0])

                phi = interpolate.interp1d(video_time_array, phi)
                videoAngles = [theta + np.pi / 2, [ phi(i) + np.pi for i in video_time_array ]]
            else:
                print('Something went wrong in your videocommand, please make sure that you have the right amount of commands and try again')
                sys.exit()
            videoAngles = np.transpose(np.array(videoAngles))
            self.dump_to_video(video_time_array, videoAngles, videoPos, [0.0, 0.0, 1.0], fname='satVideo.xml')
            print('Video recorded to satVideo.xml')
        print(' * Progress: [%s] %4.1f%s' % (lenthBar * '=', 100.0, '%'))
        if planet_pos:
            if planet_pos:
                for i in range(len(tFullArray)):
                    pos[:, 0:-1, i] = self.findPositions(tFullArray[i])
                    if i % 15 == 0:
                        progress = float(i) / len(tFullArray)
                        print(' * Computing planet positions:     [%s%s] %4.1f%s \r' % (int(ceil(lenthBar * progress)) * '=',
                         int(floor(lenthBar * (1 - progress))) * ' ',
                         progress * 100,
                         '%'), end=' ')

                print(' * Computing planet positions:     [%s] %4.1f%s' % (lenthBar * '=', 100.0, '%'))
        statusFile = open(self.sFile, 'rb')
        status = np.load(statusFile)
        orientBool = status[1][0]
        statusFile.close()
        if orientBool:
            statusFile = open(self.sFile, 'rb')
            status = np.load(statusFile)
            statusFile.close()
            status[-1] = [eventTime[-1], xv[-1, :]]
            statusFile = open(self.sFile, 'wb')
            np.save(statusFile, status)
            statusFile.close()
            print('Recorded current satellite position and velocity (at time ', eventTime[-1], '):')
            print('Position = (', xv[(-1, 0)], ', ', xv[(-1, 1)], ')')
            print('Velocity = (', xv[(-1, 2)], ', ', xv[(-1, 3)], ')')
            if self.fuelPrint:
                print()
                print('Total fuel usage: %.2f kg.' % (self.total_fuel_calculated - self.fuel_left))
        if planet_pos:
            return (pos, tFullArray)
        else:
            return 0

    def land_on_planet(self, p, fileName = 'landCommands.txt', dt = 10):
        global landingEngineUsed
        global onGround
        global mafri
        global initialTime
        global landingEngine
        global velvel
        global cameraIsOK
        xv = np.zeros(6)
        if len(self.mu) < p + 1:
            print('You are trying to land on a planet that does not exist, there is no planet nr %i!' % p)
            sys.exit()
        cd = 1.0
        A = 16.0
        mu = self.mu[p]
        atmosphereDensity = self.rho0[p]
        landerMass = 1100.0
        parachuteArea = 42.0
        velvel = 0
        mafri = -1
        initialTime = 0
        landingEngine = [0, 0]
        onGround = False
        landingEngineUsed = False
        cameraIsOK = -1
        M_sun = 1.989e+30
        G_SI = 6.674 * 1e-11
        T0 = self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[p]))
        kB = 1.381 * 1e-23
        protonMass = 1.673 * 1e-27
        gc = self.mass[p] * M_sun * G_SI
        g = gc / (1000 * self.radius[p]) ** 2
        c2 = mu * protonMass * g / (kB * T0)
        cc = mu * protonMass / kB
        h0 = 1.0 / c2
        period = self.period[p] * 24 * 3600
        escVel = sqrt(2 * gc / (1000 * self.radius[p]))
        densteps = 300000
        pivot = 1000 * h0 / densteps
        hArray = np.logspace(np.log10(pivot), np.log10(1000 * h0), densteps) - pivot
        hArray[0] = 0
        density = np.zeros(densteps)
        density[0] = self.rho0[p]
        T = T0
        P = self.rho0[p] * kB * T0 / (mu * protonMass)
        adiab = P ** (-0.4) * T ** 1.4
        for i in range(0, densteps - 1):
            dh = hArray[i + 1] - hArray[i]
            r = 1000 * self.radius[p] + hArray[i]
            dP = -density[i] * gc / (r * r)
            P += dP * dh
            T = (adiab * P ** 0.4) ** 0.7142857142857143
            T = max(T, 0.5 * T0)
            density[i + 1] = cc * P / T

        interpolatedDensity = interpolate.interp1d(hArray, density, bounds_error=False, fill_value=0.0)

        def FD(r):
            return 0.5 * A * cd * interpolatedDensity(r)

        def landingRHS(xv, t):
            global landingEngineUsed
            global onGround
            global velvel
            global cameraIsOK
            global mafri
            vatm = 2 * pi / period * np.array([-xv[1], xv[0], 0])
            vel = xv[3:6] - vatm
            velvel = abs(np.dot(xv[3:6], xv[0:3] / np.linalg.norm(xv[0:3])))
            r = np.linalg.norm(xv[0:3]) - self.radius[p] * 1000
            acc = 0
            if r < landingEngine[0]:
                acc += landingEngine[1] * xv[0:3] / np.linalg.norm(xv[0:3]) / landerMass
                if not landingEngineUsed and r > 0:
                    print('Landing engine engaged at time %g s.' % t)
                    landingEngineUsed = True
            if r < 0 and not onGround:
                onGround = True
                if velvel > 3.0:
                    if cameraIsOK < 0:
                        print('CRASH!! You collided with planet nr %i at time %g s with velocity %g m/s. Too bad :(' % (p, t, np.linalg.norm(vel)))
                        cameraIsOK = t
                        return np.zeros(6)
                elif velvel < escVel:
                    print('Successful landing! Landed on planet nr %i at time %g s with velocity %g m/s. Congratulations!' % (p, t, velvel))
                    rot = self.computeRotation(t * seconds_to_years + initialTime)
                    rot -= self.computeRotation(initialTime)
                    phi = arctan2(xv[1], xv[0]) - rot[p]
                    while phi < 0:
                        phi += 2 * pi

                    theta = arccos(xv[2] / sqrt(xv[0] * xv[0] + xv[1] * xv[1] + xv[2] * xv[2]))
                    rad2deg = 180 / pi
                    print('Landing site: phi = ', phi * rad2deg, 'deg,  theta = ', theta * rad2deg, 'deg.')
                else:
                    print('Tried to land on "planet" nr %i at time %g s... \nBut bounced off because %g m/s is larger than escape velocity!!' % (p, t, np.linalg.norm(vel)))
                    dxv = np.zeros(6)
                    dxv[0] = -xv[3]
                    dxv[1] = -xv[4]
                    dxv[2] = -xv[5]
                    dxv[3] = 0
                    dxv[4] = 0
                    dxv[5] = 0
                    return dxv
            acc -= FD(r) * vel * np.linalg.norm(vel) / landerMass
            fric = abs(FD(r) * np.linalg.norm(vel) ** 2)
            if fric > mafri:
                mafri = fric
            R = np.linalg.norm(xv[0:3])
            acc += -G_SI * self.mass[p] * M_sun / R ** 3 * xv[0:3]
            if onGround:
                xv[3:6] = vatm
                acc = np.zeros(3)
            dxv = np.zeros(6)
            dxv[0] = xv[3]
            dxv[1] = xv[4]
            dxv[2] = xv[5]
            dxv[3] = acc[0]
            dxv[4] = acc[1]
            dxv[5] = acc[2]
            return dxv

        eventType, eventTime, eventCommand, landingEngine = self.read_satellite_input_file(fileName, dt)
        nIntegrate = 1000

        def intforward(t1, t2, xv, maxStep = 0):
            tarr = np.linspace(t1, t2, nIntegrate)
            return (tarr, integrate.odeint(landingRHS, xv, tarr, mxstep=1000000, hmax=maxStep))

        tFullArray = np.zeros((len(eventTime) - 1) * nIntegrate)
        pos = np.zeros(shape=(3, 1, (len(eventTime) - 1) * nIntegrate))
        maxStep = 0
        lenthBar = 40
        videoTimes = []
        videoPos = []
        videoCommand = []
        landerIsLaunched = False
        for i in range(len(eventTime)):
            rr = np.linalg.norm(xv[0:3]) - self.radius[p] * 1000
            if velvel > 0:
                maxStep = 0.001 * rr / velvel
            if onGround:
                xv[0:3] *= self.radius[p] * 1000 / np.linalg.norm(xv[0:3])
                maxStep = 1.0
            if mafri > 1000.0 and not landerIsLaunched:
                print('Orbital decay due to air resistance at time %g!' % eventTime[i - 1])
                sys.exit(9)
            if mafri / 0.3 > 1000000.0 and cameraIsOK < 0:
                if A > 0.31:
                    print('Parachute ropes broke due to extreme pull!')
                print('Lander module burned up due to air friction at time %g!' % eventTime[i - 1])
                cameraIsOK = eventTime[i - 1]
            mafri = 0
            if eventType[i] == 0 and cameraIsOK < 0:
                planet_pos = self.findPositions(initialTime + eventTime[i] * seconds_to_years)
                pp = np.array([planet_pos[0, p], planet_pos[1, p], 0.0])
                if len(eventCommand[i]) == 5:
                    upp = [eventCommand[i][2], eventCommand[i][3], eventCommand[i][4]]
                else:
                    upp = km_to_AU * xv[0:3] / 1000.0
                    upp /= np.linalg.norm(upp)
                self.dump_to_xml(initialTime + eventTime[i] * seconds_to_years, km_to_AU * xv[0:3] / 1000.0 + pp, [eventCommand[i][0],
                 eventCommand[i][1],
                 upp[0],
                 upp[1],
                 upp[2]], i + 100)
                print('Picture taken! selfie_%04d.xml' % (i + 100))
            elif eventType[i] == 90:
                planet_pos = self.findPositions(initialTime + eventTime[i] * seconds_to_years)
                pp = np.array([planet_pos[0, p], planet_pos[1, p], 0.0])
                videoTimes.append(initialTime + eventTime[i] * seconds_to_years)
                videoCommand.append(eventCommand[i])
            elif eventType[i] == 1:
                if landerIsLaunched:
                    print('Cannot boost, lander is launched.')
                elif len(eventCommand[i]) == 3:
                    xv[3:6] += np.array(eventCommand[i])
                else:
                    print('Not a proper boost command on line %i' % i)
                    sys.exit()
            elif eventType[i] == 2:
                if landerIsLaunched and cameraIsOK < 0:
                    if len(eventCommand[i]) > 0:
                        A = float(eventCommand[i][0])
                    else:
                        A = parachuteArea
                    A = max(A, 0.3)
                    print('Parachute with area', A, 'm^2 deployed.')
                else:
                    print('Cannot deploy parachute.')
            elif eventType[i] == 4:
                xv[3:6] += np.array(eventCommand[i])
                A = 0.3
                landerMass = 90
                maxStep = 0.1
                landerIsLaunched = True
                print('Lander is launched at time', eventTime[i])
            elif eventType[i] == -1:
                initialTime = eventTime[i]
                dt2 = 1e-08
                initPos = self.findPositions(eventTime[i])
                initVel = (initPos - self.findPositions(eventTime[i] - dt2)) / dt2
                udlim = self.a[p] * np.sqrt(self.mass[p] / self.star_mass / 10.0)
                udlim = udlim / km_to_AU * 1000
                xv[0:2] = (eventCommand[i][0:2] - initPos[:, p]) / km_to_AU * 1000
                xv[3:5] = (eventCommand[i][2:4] - initVel[:, p]) / km_to_AU * 1000 * seconds_to_years
                eventTime[i] = 0
                if np.linalg.norm(xv[0:2]) > udlim:
                    print('The satellite is not close enough to planet ', p, ' to engage the landing module land_on_planet():')
                    print('Distance             (SI units) = %g' % np.linalg.norm(xv[0:2]))
                    print('Upper distance limit (SI units) = %g' % udlim)
                    print('Exiting.')
                    sys.exit(1)
                else:
                    print('Starting landing module, land_on_planet():')
                    print('Initial relative position = (', eventCommand[i][0] - initPos[0, p], ', ', eventCommand[i][1] - initPos[1, p], ', ', 0, ')')
                    print('Initial relative velocity = (', eventCommand[i][2] - initVel[0, p], ', ', eventCommand[i][3] - initVel[1, p], ', ', 0, ')')
                    print('----------------------------------------------------------')
                    print('Initial relative position (SI units) = (', xv[0], ', ', xv[1], ', ', 0, ')')
                    print('Initial relative velocity (SI units) = (', xv[3], ', ', xv[4], ', ', 0, ')')
                print('Satellite is now receiving commands. This might take some minutes...')
                if eventTime[-1] > 200000.0:
                    print('Especially if you want to command it for several days!')
                if eventTime[-1] > 10000000.0:
                    print('WARNING: are you sure you started at T=0 when calling land_on_planet?')
            elif eventType[i] == 5:
                print('AU_to_matic orientation for landing.', end=' ')
                print('Time:', eventTime[i])
                print('Pos: ', xv[0:3])
                print('Vel: ', xv[3:6])
            if i < len(eventTime) - 1:
                if cameraIsOK < 0:
                    t, xv1 = intforward(eventTime[i], eventTime[i + 1], xv[:], maxStep)
                    xv = xv1[-1, :]
                    pos[:, -1, i * nIntegrate:(i + 1) * nIntegrate] = xv1[:, 0:3].transpose()
                    tFullArray[i * nIntegrate:(i + 1) * nIntegrate] = t
                else:
                    pos[:, -1, i * nIntegrate:(i + 1) * nIntegrate] = pos[:, -1, i * nIntegrate - 1][:, None]
                    tFullArray[i * nIntegrate:(i + 1) * nIntegrate] = np.linspace(eventTime[i], eventTime[i + 1], nIntegrate)

        if len(videoTimes) > 0:
            if len(videoTimes) == 1:
                print('You need 2 video calls, first to start the video and second to end it.')
                sys.exit()
            if len(videoTimes) > 2:
                print('You can only make 1 video at a time, try again with 2 video calls, 1 to start the video and 1 to end it')
                sys.exit()
            t0, t1 = videoTimes[0], videoTimes[1]
            nFrames = 10000
            video_time_array = np.linspace(t0, t1, nFrames)
            sat_pos = interpolate.interp1d(initialTime + tFullArray * seconds_to_years, pos[:, -1, :])
            min_dt = 1000000.0
            if len(videoCommand[0]) + len(videoCommand[1]) == 4:
                theta0, phi0, theta1, phi1 = (videoCommand[0][0],
                 videoCommand[0][1],
                 videoCommand[1][0],
                 videoCommand[1][1])
                videoPos = [ list(sat_pos(i) * km_to_AU / 1000.0 + np.array(list(self.findPositions(i)[:, p]) + [0])) for i in video_time_array ]
                phi = np.linspace(phi0, phi1, nFrames)
                phi = interpolate.interp1d(video_time_array, phi)
                theta = np.linspace(theta0, theta1, nFrames) + np.pi / 2
                theta = interpolate.interp1d(video_time_array, theta)
                videoAngles = [[ theta(i) for i in video_time_array ], [ phi(i) + np.pi for i in video_time_array ]]
                videoCommand = [p]
            elif len(videoCommand[0] + videoCommand[1]) == 2:
                videoPos = [ list(sat_pos(i) * km_to_AU / 1000.0 + np.array(list(self.findPositions(i)[:, p]) + [0])) for i in video_time_array ]
                camDirec = [ -sat_pos(i) for i in video_time_array ]
                phi = np.zeros(nFrames)
                theta = np.zeros(nFrames)
                for i in range(nFrames):
                    phi[i] = np.arctan2(camDirec[i][1], camDirec[i][0])
                    theta[i] = np.arccos(camDirec[i][2] / np.linalg.norm(camDirec[i]))

                videoAngles = [theta, phi]
                videoCommand = [0.0, 0.0, 1.0]
            else:
                print('The video commands were not entered correctly, try again')
            videoAngles = np.transpose(np.array(videoAngles))
            if cameraIsOK > 0:
                cameraIsOK = initialTime + cameraIsOK * seconds_to_years
            self.dump_to_video(video_time_array, videoAngles, videoPos, videoCommand, p, isOK=cameraIsOK, fname='landerVideo.xml')
            print('Video recorded to landerVideo.xml')
        if not onGround:
            print('No landing detected at time = ', eventTime[-1], 's. Exiting.')

    def temperature_to_RGB(self):
        temp = self.temperature / 100
        if temp <= 66:
            red = 255
        else:
            red = temp - 60
            red = 329.698727446 * red ** (-0.1332047592)
        if temp <= 66:
            green = temp
            green = 99.4708025861 * log(green) - 161.1195681661
        else:
            green = temp - 60
            green = 288.1221695283 * green ** (-0.0755148492)
        if temp >= 66:
            blue = 255
        elif temp <= 19:
            blue = 0
        else:
            blue = temp - 10
            blue = 138.5177312231 * log(blue) - 305.0447927307
        rgb = [red, green, blue]
        for c in rgb:
            if c < 0:
                c = 0
            elif c > 255:
                c = 255

        return rgb

    def dump_to_xml(self, time, pos, command, ID):
        planet_pos = self.findPositions(time)
        rotation = self.computeRotation(time)
        try:
            widtheight = int(command[5])
        except:
            widtheight = 900

        planets = etree.Element('Planets')
        objects = etree.Element('Objects')
        star = etree.SubElement(objects, 'SerializedMCAstObject')
        etree.SubElement(star, 'category').text = str('star')
        etree.SubElement(star, 'pos_x').text = str(0)
        etree.SubElement(star, 'pos_z').text = str(0)
        etree.SubElement(star, 'pos_y').text = str(0)
        etree.SubElement(star, 'rot_y').text = str(0)
        etree.SubElement(star, 'radius').text = str(100000)
        etree.SubElement(star, 'temperature').text = str(4000)
        etree.SubElement(star, 'seed').text = str(int(self.seed * 1000 + 990))
        etree.SubElement(star, 'atmosphereDensity').text = str(10)
        etree.SubElement(star, 'atmosphereHeight').text = str(1.025)
        etree.SubElement(star, 'outerRadiusScale').text = str(1.0025)
        etree.SubElement(star, 'name').text = str('The star')
        star_objects = etree.SubElement(star, 'Objects')
        for i in range(self.number_of_planets):
            planet = etree.SubElement(star_objects, 'SerializedMCAstObject')
            etree.SubElement(planet, 'radius').text = str(self.radius[i])
            for j, s in zip(range(2), ['x', 'z']):
                etree.SubElement(planet, 'pos_%s' % s).text = str(planet_pos[j, i])

            etree.SubElement(planet, 'pos_y').text = str(0)
            etree.SubElement(planet, 'rotation').text = str(-rotation[i])
            etree.SubElement(planet, 'temperature').text = str(self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[i])))
            etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + i))
            etree.SubElement(planet, 'parentPlanet').text = str(int(self.seed * 1000 + 990))
            etree.SubElement(planet, 'atmosphereDensity').text = str(np.log(self.rho0[i]) / np.log(25))
            etree.SubElement(planet, 'atmosphereHeight').text = str(1.025)
            etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
            etree.SubElement(planet, 'category').text = str('planet')
            etree.SubElement(planet, 'name').text = 'Planet %i' % i
            if len(self.moons) == self.number_of_planets:
                moonPos, moonRot = self.moons[i].findPosRot(time)
                moonPos = np.transpose(moonPos)
                for ii in range(self.moons[i].number_of_moons):
                    planet = etree.SubElement(planets, 'SerializedPlanet')
                    for j, s in zip(range(3), ['x', 'z', 'y']):
                        etree.SubElement(planet, 'pos_%s' % s).text = str(moonPos[j, ii])

                    etree.SubElement(planet, 'rotation').text = str(-moonRot[ii])
                    etree.SubElement(planet, 'radius').text = str(self.moons[i].radius[ii])
                    etree.SubElement(planet, 'temperature').text = str(self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[i])))
                    etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + i + 10 * (ii + 1)))
                    etree.SubElement(planet, 'parentPlanet').text = str(int(self.seed * 1000 + i))
                    etree.SubElement(planet, 'atmosphereDensity').text = str(0)
                    etree.SubElement(planet, 'atmosphereHeight').text = str(0)
                    etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
                    etree.SubElement(planet, 'category').text = str('moon')
                    etree.SubElement(planet, 'name').text = 'Moon %i of Planet %i' % (ii + 1, i)

        cameras = etree.Element('Cameras')
        for i in range(1):
            camera = etree.SubElement(cameras, 'SerializedCamera')
            for j, s in zip(range(2), ['x', 'z']):
                etree.SubElement(camera, 'cam_%s' % s).text = str(pos[j])

            if len(pos) == 3:
                etree.SubElement(camera, 'cam_y').text = str(pos[2])
            else:
                etree.SubElement(camera, 'cam_y').text = str(0)
            thetai = float(command[0])
            phii = float(command[1])
            dirvec = [np.sin(thetai) * np.cos(phii), np.sin(thetai) * np.sin(phii), np.cos(thetai)]
            upvec = [command[2], command[3], command[4]]
            downvec = [-command[2], -command[3], -command[4]]
            if dirvec == upvec or dirvec == downvec:
                upvec == self.prevdir
            else:
                self.prevdir = upvec
            for j, s in zip(range(3), ['x', 'z', 'y']):
                etree.SubElement(camera, 'dir_%s' % s).text = str(dirvec[j])

            for j, s in zip(range(3), ['x', 'z', 'y']):
                etree.SubElement(camera, 'up_%s' % s).text = str(upvec[j])

            etree.SubElement(camera, 'fov').text = str(70)
            etree.SubElement(camera, 'time').text = str(time)
            etree.SubElement(camera, 'frame').text = str(0)

        tree1 = etree.ElementTree(planets)
        tree2 = etree.ElementTree(cameras)
        outFile = open(self.path + 'selfie_%04d.xml' % ID, 'w')
        outFile.write('<?xml version="1.0" encoding="utf-8"?>\n')
        outFile.write('<SerializedWorld xmlns:xsi="http://www.w3.org/2001/')
        outFile.write('XMLSchema-instance"\n xmlns:xsd="http://www.w3.org/2001/XMLSchema">\n')
        rgb = self.temperature_to_RGB()
        outFile.write('<sun_col_r>%.3f</sun_col_r>\n' % (rgb[0] / 255.0))
        outFile.write('<sun_col_g>%.3f</sun_col_g>\n' % (rgb[1] / 255.0))
        outFile.write('<sun_col_b>%.3f</sun_col_b>\n' % (rgb[2] / 255.0))
        outFile.write('<sun_intensity>%.3f</sun_intensity>\n' % 0.1)
        outFile.write('<screenshot_width>%d</screenshot_width>\n' % widtheight)
        outFile.write('<screenshot_height>%d</screenshot_height>\n' % widtheight)
        outFile.write('<global_radius_scale>0.985</global_radius_scale>\n')
        outFile.write('<resolution>%.3d</resolution>\n' % 48)
        outFile.write('<skybox>%.3d</skybox>\n' % 0)
        outFile.write(str(etree.tostring(objects, pretty_print=True)))
        outFile.write(str(etree.tostring(cameras, pretty_print=True)))
        outFile.write('</SerializedWorld>')
        outFile.close()
        outFile = open(self.path + 'selfie_%04d.xml' % ID, 'r')
        returnString = outFile.read()
        outFile.close()
        return returnString

    def dump_to_video(self, times, angles, pos, command, p = -1, withCam = True, isOK = -1, fname = 'video.xml'):
        times = np.array(times)
        angles = np.array(angles)
        pos = np.array(pos)
        numframes = len(times)
        if isOK < 0:
            isOK = max(times) + 1
        velies = np.zeros(len(times))
        accies = np.zeros(len(times))
        for i in range(numframes - 1):
            dt = times[i + 1] - times[i]
            if times[i + 1] > isOK:
                break
            if dt == 0:
                print('Uh, (at least) one of your time steps is zero!')
                continue
            velies[i] = np.linalg.norm(pos[i] - pos[i + 1]) / dt
            if i > 0:
                accies[i] = abs(velies[i] - velies[i - 1]) / dt

        accies.sort()
        accies = accies[::-1]
        if len(accies) < 10:
            print('You should have more than 10 time steps!')
            avgcci = np.average(accies)
        else:
            avgcci = np.average(accies[0:10])
        if avgcci <= 0:
            print('Your spacecraft is not accelerating during this video...')
            avgcci = 1e-11
        newdt = sqrt(1e-05 / avgcci)
        skipframes = min([int(newdt / dt), numframes // 100])
        skipframes = max([skipframes, 1])
        times = times[::skipframes]
        angles = angles[::skipframes]
        pos = pos[::skipframes]
        numframes = min([len(times), len(angles), len(pos)])
        if numframes < 10:
            print('Warning! Very few frames, are there enough time points in the arrays?')
        print('Writing video, skipping', skipframes - 1, 'frames for each output frame.')
        planet_poses = []
        rotations = []
        for j in range(numframes):
            planet_poses.append(self.findPositions(times[j]))
            rotations.append(self.computeRotation(times[j], False))

        widtheight = 900
        nulltid = times[0]
        planets = etree.Element('Planets')
        planz = list(range(self.number_of_planets))
        if p >= 0:
            planz = [p]
        objects = etree.Element('Objects')
        star = etree.SubElement(objects, 'SerializedMCAstObject')
        etree.SubElement(star, 'category').text = str('star')
        etree.SubElement(star, 'pos_x').text = str(0)
        etree.SubElement(star, 'pos_z').text = str(0)
        etree.SubElement(star, 'pos_y').text = str(0)
        etree.SubElement(star, 'rot_y').text = str(0)
        etree.SubElement(star, 'radius').text = str(100000)
        etree.SubElement(star, 'temperature').text = str(4000)
        etree.SubElement(star, 'seed').text = str(int(self.seed * 1000 + 990))
        etree.SubElement(star, 'atmosphereDensity').text = str(10)
        etree.SubElement(star, 'atmosphereHeight').text = str(1.025)
        etree.SubElement(star, 'outerRadiusScale').text = str(1.0025)
        etree.SubElement(star, 'name').text = str('The star')
        star_objects = etree.SubElement(star, 'Objects')
        for i in planz:
            planet = etree.SubElement(star_objects, 'SerializedMCAstObject')
            for j, s in zip(range(2), ['x', 'z']):
                etree.SubElement(planet, 'pos_%s' % s).text = str(planet_poses[0][j, i])

            etree.SubElement(planet, 'pos_y').text = str(0)
            etree.SubElement(planet, 'rotation').text = str(-rotations[0][i])
            etree.SubElement(planet, 'radius').text = str(self.radius[i])
            etree.SubElement(planet, 'temperature').text = str(self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[i])))
            etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + i))
            etree.SubElement(planet, 'atmosphereDensity').text = str(np.log(self.rho0[i]) / np.log(25))
            etree.SubElement(planet, 'atmosphereHeight').text = str(1.025)
            etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
            etree.SubElement(planet, 'category').text = str('planet')
            etree.SubElement(planet, 'name').text = str('planet ' + str(i))
            frames = etree.SubElement(planet, 'Frames')
            for k in range(0, numframes):
                frame = etree.SubElement(frames, 'Frame')
                etree.SubElement(frame, 'id').text = str(k)
                for j, s in zip(range(2), ['x', 'z']):
                    etree.SubElement(frame, 'pos_%s' % s).text = str(planet_poses[max([k, 0])][j, i])

                etree.SubElement(frame, 'pos_y').text = str(0)
                etree.SubElement(frame, 'rotation').text = str(-rotations[max([k, 0])][i])

            if len(self.moons) == self.number_of_planets:
                moonPoses = np.zeros([numframes, 3, self.moons[i].number_of_moons])
                moonRots = np.zeros([numframes, self.moons[i].number_of_moons])
                for k in range(0, numframes):
                    moonPos, moonRot = self.moons[i].findPosRot(times[k])
                    moonPoses[k] = np.transpose(moonPos)
                    moonRots[k] = moonRot

                for ii in range(self.moons[i].number_of_moons):
                    planet = etree.SubElement(planets, 'SerializedPlanet')
                    for j, s in zip(range(3), ['x', 'z', 'y']):
                        etree.SubElement(planet, 'pos_%s' % s).text = str(moonPoses[0][j, ii])

                    etree.SubElement(planet, 'rotation').text = str(-moonRots[0][ii])
                    frames = etree.SubElement(planet, 'Frames')
                    for k in range(0, numframes):
                        frame = etree.SubElement(frames, 'Frame')
                        etree.SubElement(frame, 'id').text = str(k)
                        for j, s in zip(range(2), ['x', 'z', 'y']):
                            etree.SubElement(frame, 'pos_%s' % s).text = str(moonPoses[max([k, 0])][j, ii])

                        etree.SubElement(frame, 'rotation').text = str(-moonRots[max([k, 0])][ii])

                    etree.SubElement(planet, 'radius').text = str(self.moons[i].radius[ii])
                    etree.SubElement(planet, 'temperature').text = str(self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[i])))
                    etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + i + 10 * (ii + 1)))
                    etree.SubElement(planet, 'parentPlanet').text = str(int(self.seed * 1000 + i))
                    etree.SubElement(planet, 'atmosphereDensity').text = str(0)
                    etree.SubElement(planet, 'atmosphereHeight').text = str(0)
                    etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
                    etree.SubElement(planet, 'category').text = str('moon')
                    etree.SubElement(planet, 'name').text = 'Moon %i of Planet %i' % (ii + 1, i)

        if withCam == False:
            planet = etree.SubElement(planets, 'SerializedPlanet')
            for j, s in zip(range(2), ['x', 'z']):
                etree.SubElement(planet, 'pos_%s' % s).text = str(pos[0][j])

            if len(pos[0]) == 3:
                etree.SubElement(planet, 'pos_y').text = str(pos[0][2])
            else:
                etree.SubElement(planet, 'pos_y').text = str(0)
            etree.SubElement(planet, 'rotation').text = str(0)
            frames = etree.SubElement(planet, 'Frames')
            for i in range(numframes):
                frame = etree.SubElement(frames, 'Frame')
                etree.SubElement(frame, 'id').text = str(i)
                for j, s in zip(range(2), ['x', 'z']):
                    etree.SubElement(frame, 'pos_%s' % s).text = str(pos[i][j])

                if len(pos[i]) == 3:
                    etree.SubElement(frame, 'pos_y').text = str(pos[i][2])
                else:
                    etree.SubElement(frame, 'pos_y').text = str(0)
                etree.SubElement(frame, 'rotation').text = str(0)

            etree.SubElement(planet, 'radius').text = str(4000)
            etree.SubElement(planet, 'temperature').text = str(300)
            etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + 950))
            etree.SubElement(planet, 'parentPlanet').text = str(int(self.seed * 1000 + 990))
            etree.SubElement(planet, 'atmosphereDensity').text = str(1)
            etree.SubElement(planet, 'atmosphereHeight').text = str(1.025)
            etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
            etree.SubElement(planet, 'name').text = 'Spacecraft'
            etree.SubElement(planet, 'category').text = str('spacecraft')
            cameras = etree.Element('Cameras')
        else:
            cameras = etree.Element('Cameras')
            for i in range(0, numframes):
                camera = etree.SubElement(cameras, 'SerializedCamera')
                for j, s in zip(range(2), ['x', 'z']):
                    etree.SubElement(camera, 'cam_%s' % s).text = str(pos[max([i, 0])][j])

                if len(pos[i]) == 3:
                    etree.SubElement(camera, 'cam_y').text = str(pos[max([i, 0])][2])
                else:
                    etree.SubElement(camera, 'cam_y').text = str(0)
                thetai = float(angles[max([i, 0])][0])
                phii = float(angles[max([i, 0])][1])
                dirvec = [np.sin(thetai) * np.cos(phii), np.sin(thetai) * np.sin(phii), np.cos(thetai)]
                if len(command) == 1:
                    pp = command[0]
                    if pp not in planz:
                        print('ERROR, planet', pp, 'not visible!')
                        sys.exit(10)
                    upvec = pos[max([i, 0])] - np.array(list(planet_poses[max([i, 0])][:, pp]) + [0])
                else:
                    upvec = [0.0, 0.0, 1.0]
                for j, s in zip(range(3), ['x', 'z', 'y']):
                    etree.SubElement(camera, 'dir_%s' % s).text = str(dirvec[j])

                for j, s in zip(list(range(3)), ['x', 'z', 'y']):
                    etree.SubElement(camera, 'up_%s' % s).text = str(upvec[j])

                etree.SubElement(camera, 'fov').text = str(70)
                etree.SubElement(camera, 'time').text = str(times[i] - nulltid)
                etree.SubElement(camera, 'frame').text = str(i)
                if times[i] > isOK:
                    etree.SubElement(camera, 'status').text = str(1)

        tree1 = etree.ElementTree(planets)
        tree2 = etree.ElementTree(cameras)
        outFile = open(self.path + fname, 'w')
        outFile.write('<?xml version="1.0" encoding="utf-8"?>\n')
        outFile.write('<SerializedWorld xmlns:xsi="http://www.w3.org/2001/')
        outFile.write('XMLSchema-instance"\n xmlns:xsd="http://www.w3.org/2001/XMLSchema">\n')
        rgb = self.temperature_to_RGB()
        outFile.write('<sun_col_r>%.3f</sun_col_r>\n' % (rgb[0] / 255.0))
        outFile.write('<sun_col_g>%.3f</sun_col_g>\n' % (rgb[1] / 255.0))
        outFile.write('<sun_col_b>%.3f</sun_col_b>\n' % (rgb[2] / 255.0))
        outFile.write('<sun_intensity>%.3f</sun_intensity>\n' % 0.1)
        outFile.write('<screenshot_width>%d</screenshot_width>\n' % widtheight)
        outFile.write('<screenshot_height>%d</screenshot_height>\n' % widtheight)
        outFile.write('<global_radius_scale>0.985</global_radius_scale>\n')
        outFile.write('<resolution>%.3d</resolution>\n' % self.resolution)
        outFile.write('<uuid>5acbd644-37c7-11e6-ac61-9e71128cae77</uuid> \n')
        outFile.write('<skybox>%.3d</skybox>\n' % 0)
        outFile.write(str(etree.tostring(objects, pretty_print=True)))
        outFile.write(str(etree.tostring(cameras, pretty_print=True)))
        outFile.write('</SerializedWorld>')
        outFile.close()

    def orbit_xml(self, planet_pos, times):
        pos = planet_pos
        numframes = len(times)
        if len(times) < 100:
            print('Please give more time steps! 100 would be nice, but more is better.')
            return
        if len(planet_pos[0, 0, :]) != numframes:
            print('Shape of position array does not match length of time array!')
            print('Shape should be [2, nPlanets, len(times)]')
            return
        maxdrr = 0.0
        for j in range(numframes - 1):
            drr = sqrt((pos[0, :, j] - pos[0, :, j + 1]) ** 2 + (pos[1, :, j] - pos[1, :, j + 1]) ** 2) / self.a
            if max(drr) > maxdrr:
                maxdrr = max(drr)

        maxdrr = max(10.0 / numframes, maxdrr)
        skipstep = max([int(0.01 / maxdrr), 1])
        dtt = skipstep * (times[2] - times[1])
        rotationdt = self.computeRotation(dtt, False) - self.computeRotation(0.0, False)
        rotationdt[rotationdt == 0] = 1e-12
        rotationfactor = min([3.0 / max(rotationdt), 1.0])
        numfra = int(numframes // skipstep)
        print('Writing xml with', numfra, 'frames. (Skipping', skipstep - 1, 'for each).')
        print('Rotation factor is', rotationfactor, '(IRL planets rotate', int(1.0 / rotationfactor), 'times faster).')
        rotations = []
        for j in range(numframes):
            rotfram = self.computeRotation(times[j], False)
            rotations.append(rotfram * rotationfactor)

        widtheight = 900
        objects = etree.Element('Objects')
        star = etree.SubElement(objects, 'SerializedMCAstObject')
        etree.SubElement(star, 'category').text = str('star')
        etree.SubElement(star, 'pos_x').text = str(0)
        etree.SubElement(star, 'pos_z').text = str(0)
        etree.SubElement(star, 'pos_y').text = str(0)
        etree.SubElement(star, 'rot_y').text = str(0)
        etree.SubElement(star, 'radius').text = str(100000)
        etree.SubElement(star, 'temperature').text = str(4000)
        etree.SubElement(star, 'seed').text = str(int(self.seed * 1000 + 990))
        etree.SubElement(star, 'atmosphereDensity').text = str(10)
        etree.SubElement(star, 'atmosphereHeight').text = str(1.025)
        etree.SubElement(star, 'outerRadiusScale').text = str(1.0025)
        etree.SubElement(star, 'name').text = str('The star')
        star_objects = etree.SubElement(star, 'Objects')
        planz = list(range(self.number_of_planets))
        for i in planz:
            planet = etree.SubElement(star_objects, 'SerializedMCAstObject')
            for j, s in zip(range(2), ['x', 'z']):
                etree.SubElement(planet, 'pos_%s' % s).text = str(pos[j][i][0])

            etree.SubElement(planet, 'pos_y').text = str(0)
            etree.SubElement(planet, 'rotation').text = str(-rotations[0][i])
            etree.SubElement(planet, 'radius').text = str(self.radius[i])
            etree.SubElement(planet, 'temperature').text = str(self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[i])))
            etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + i))
            etree.SubElement(planet, 'atmosphereDensity').text = str(np.log(self.rho0[i]) / np.log(25))
            etree.SubElement(planet, 'atmosphereHeight').text = str(1.025)
            etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
            etree.SubElement(planet, 'category').text = str('planet')
            etree.SubElement(planet, 'name').text = 'Planet %i' % i
            frames = etree.SubElement(planet, 'Frames')
            for k in range(0, numframes, skipstep):
                frame = etree.SubElement(frames, 'Frame')
                etree.SubElement(frame, 'id').text = str(k)
                for j, s in zip(range(2), ['x', 'z']):
                    etree.SubElement(frame, 'pos_%s' % s).text = str(pos[j][i][k])

                etree.SubElement(frame, 'pos_y').text = str(0)
                etree.SubElement(frame, 'rotation').text = str(-rotations[k][i])

            if len(self.moons) == self.number_of_planets:
                moonPoses = np.zeros([numframes, 3, self.moons[i].number_of_moons])
                moonRots = np.zeros([numframes, self.moons[i].number_of_moons])
                for k in range(0, numframes):
                    moonPos, moonRot = self.moons[i].findPosRot(times[k] * rotationfactor)
                    moonPoses[k] = np.transpose(moonPos)
                    moonRots[k] = moonRot

                for ii in range(self.moons[i].number_of_moons):
                    planet = etree.SubElement(planets, 'SerializedPlanet')
                    for j, s in zip(range(3), ['x', 'z', 'y']):
                        etree.SubElement(planet, 'pos_%s' % s).text = str(moonPoses[0][j, ii])

                    etree.SubElement(planet, 'rotation').text = str(-moonRots[0][ii])
                    frames = etree.SubElement(planet, 'Frames')
                    for k in range(0, numframes, skipstep):
                        frame = etree.SubElement(frames, 'Frame')
                        etree.SubElement(frame, 'id').text = str(k)
                        for j, s in zip(range(2), ['x', 'z', 'y']):
                            etree.SubElement(frame, 'pos_%s' % s).text = str(moonPoses[max([k, 0])][j, ii])

                        etree.SubElement(frame, 'rotation').text = str(-moonRots[max([k, 0])][ii])

                    etree.SubElement(planet, 'radius').text = str(self.moons[i].radius[ii])
                    etree.SubElement(planet, 'temperature').text = str(self.temperature * sqrt(self.star_radius * km_to_AU / (2 * self.a[i])))
                    etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + i + 10 * (ii + 1)))
                    etree.SubElement(planet, 'atmosphereDensity').text = str(0)
                    etree.SubElement(planet, 'atmosphereHeight').text = str(0)
                    etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
                    etree.SubElement(planet, 'category').text = str('moon')
                    etree.SubElement(planet, 'name').text = 'Moon %i of Planet %i' % (ii + 1, i)

        tree1 = etree.ElementTree(planet)
        outFile = open(self.path + 'video_%i.xml' % self.seed, 'w')
        outFile.write('<?xml version="1.0" encoding="utf-8"?>\n')
        outFile.write('<SerializedWorld xmlns:xsi="http://www.w3.org/2001/')
        outFile.write('XMLSchema-instance"\n xmlns:xsd="http://www.w3.org/2001/XMLSchema">\n')
        rgb = self.temperature_to_RGB()
        outFile.write('<sun_col_r>%.3f</sun_col_r>\n' % (rgb[0] / 255.0))
        outFile.write('<sun_col_g>%.3f</sun_col_g>\n' % (rgb[1] / 255.0))
        outFile.write('<sun_col_b>%.3f</sun_col_b>\n' % (rgb[2] / 255.0))
        outFile.write('<sun_intensity>%.3f</sun_intensity>\n' % 0.1)
        outFile.write('<screenshot_width>%d</screenshot_width>\n' % widtheight)
        outFile.write('<screenshot_height>%d</screenshot_height>\n' % widtheight)
        outFile.write('<global_radius_scale>0.985</global_radius_scale>\n')
        outFile.write('<resolution>%.3d</resolution>\n' % self.resolution)
        outFile.write('<uuid>5acbd644-37c7-11e6-ac61-9e71128cae77</uuid> \n')
        outFile.write('<skybox>%.3d</skybox>\n' % 0)
        outFile.write(str(etree.tostring(objects, pretty_print=True)))
        outFile.write('</SerializedWorld>')
        outFile.close()
        returnString = ''
        return returnString

    def read_satellite_input_file(self, inputFile, bigdt = 0.001):
        fileReader = open(inputFile, 'r')
        lines = fileReader.readlines()
        numberOfLines = len(lines) - 1
        fourout = False
        eventType = []
        eventTime = []
        eventCommand = []
        thistime = 1e+70
        mindt = bigdt * 0.0001
        landingEngine = np.zeros(2)
        lastTime = 0
        line = numberOfLines
        while lastTime == 0:
            try:
                lastLine = lines[line].split()
                lastTime = float(lastLine[1])
            except:
                pass

            line -= 1

        try:
            for line in lines:
                line = line.split()
                if len(line) == 0 or line[0][0] == '#':
                    continue
                prevtime = thistime
                if len(line) > 1:
                    thistime = float(line[1])
                if prevtime < thistime - mindt:
                    if (thistime - prevtime) / bigdt > 20000000000.0:
                        print('WARNING: Are you landing? Did you remember init on first line?')
                    prevtime += mindt
                    eventType.append(404)
                    eventTime.append(prevtime)
                    eventCommand.append([])
                    while prevtime < thistime - bigdt:
                        prevtime += bigdt
                        eventType.append(404)
                        eventTime.append(prevtime)
                        eventCommand.append([])

                    eventType.append(404)
                    eventTime.append(thistime - mindt)
                    eventCommand.append([])
                if line[0].lower() == 'picture':
                    eventType.append(0)
                    eventTime.append(float(line[1]))
                    if len(line) > 6:
                        eventCommand.append([float(line[2]),
                         float(line[3]),
                         float(line[4]),
                         float(line[5]),
                         float(line[6])])
                    elif len(line) > 3:
                        eventCommand.append([float(line[2]), float(line[3]), float(line[4])])
                    else:
                        print('Too few arguments for picture!')
                elif line[0].lower() == 'video':
                    eventType.append(90)
                    eventTime.append(float(line[1]))
                    if len(line) > 3:
                        eventCommand.append([float(line[2]), float(line[3])])
                    elif len(line) > 2:
                        eventCommand.append([int(line[2])])
                    else:
                        print('Too few arguments for video!')
                elif line[0].lower() == 'boost':
                    eventType.append(1)
                    eventTime.append(float(line[1]))
                    if len(line) == 4:
                        eventCommand.append([float(line[2]), float(line[3])])
                    elif len(line) == 5:
                        eventCommand.append([float(line[2]), float(line[3]), float(line[4])])
                elif line[0].lower() == 'parachute':
                    eventType.append(2)
                    eventTime.append(float(line[1]))
                    if len(line) > 2:
                        eventCommand.append([float(line[2])])
                    else:
                        eventCommand.append([])
                elif line[0].lower() == 'launch':
                    if len(line) != 4:
                        print('Not the right amount of arguments to launch.')
                        sys.exit()
                    if self.t0_launch != float(line[1]):
                        print('Warning: launch times not equal. Using t = %f years from engine_settings' % self.t0_launch)
                    fueluse_launch, launch_time, final_launch_pos, final_vel = self.mass_needed_launch()
                    self.fuel_left -= fueluse_launch
                    self.launch_pos = final_launch_pos
                    self.launch_vel = final_vel
                    eventType.append(-1)
                    eventTime.append(launchTime)
                    launchComm = [final_launch_pos[0],
                     final_launch_pos[1],
                     float(line[2]),
                     float(line[3])]
                    boost = np.array((float(line[2]), float(line[3]))) - final_vel
                    boost_fuel = self.mass_needed(final_vel, boost)
                    print('Used %f kg fuel for correction boost after escape velocity launch.' % boost_fuel)
                    if self.fuel_left < -0.15 * self.total_fuel_calculated:
                        print('You ran out of fuel. Too bad!')
                        sys.exit()
                    eventCommand.append(launchComm)
                elif line[0].lower() == 'init':
                    fourout = True
                    statusFile = open(self.sFile, 'rb')
                    status = np.load(statusFile)
                    t, xv = status[2]
                    statusFile.close()
                    eventType.append(-1)
                    eventTime.append(t)
                    thistime = 0
                    if bigdt < 0.01:
                        mindt = 0.1
                        bigdt = 10
                        print('For some reason dt was small during landing. We are in SI now.')
                    eventCommand.append([xv[0],
                     xv[1],
                     xv[2],
                     xv[3]])
                elif line[0].lower() == 'launchlander':
                    if len(line) != 5:
                        print('Not the right amount of arguments to launchLander.')
                        sys.exit()
                    eventType.append(4)
                    eventTime.append(float(line[1]))
                    eventCommand.append([float(line[2]), float(line[3]), float(line[4])])
                elif line[0].lower() == 'landing':
                    if len(line) != 4:
                        print('Not the right amount of arguments to landing.')
                        sys.exit()
                    eventType.append(3)
                    eventTime.append(float(line[1]))
                    eventCommand.append([])
                    landingEngine[0] = float(line[2])
                    landingEngine[1] = float(line[3])
                elif line[0].lower() == 'orient':
                    eventType.append(5)
                    eventTime.append(float(line[1]))
                    eventCommand.append([])
                else:
                    raise ValueError(line[0])
                    exit()

        except ValueError as error:
            traceback.print_exc(limit=1, file=sys.stdout)
            print("\nError: Unable to interpret keyword '%s' in file %s" % (error.message, inputFile))
        finally:
            fileReader.close()

        if fourout:
            return (eventType,
             eventTime,
             eventCommand,
             landingEngine)
        else:
            return (eventType, eventTime, eventCommand)

    def engine_settings(self, dpdt, numBox, partPerSec, total_fuel_calculated, launch_T, init_sat_pos, t0_launch):
        """ Method which sets the variables to calculate the needed fuel for the voyage as well
        as the positions and velocity of satellite and planet before launching. MOST units should
        be in SI-units... 
        dpdt is impulse from one box [kg*m/s\xc2\xb2]
        numBox is number of boxes,
        partPerSec is the amout of particles escaping from one box,
        total_fuel_calculated is the initial amount of fuel in the rocket [kg],
        launch_T is the time expected to reach escape velocity [s],
        init_sat_pos is the desired initial position of the rocket. NB, must be close to the surface of the planet! (x[AU], y[AU])
        t0_launch is the time after the initial time (t=0) when the launch starts. [years]"""
        self.engine_checked = True
        self.dpdt1 = dpdt
        self.numBox = numBox
        self.partPerSec = partPerSec
        self.total_fuel_calculated = total_fuel_calculated
        self.fuel_left = self.total_fuel_calculated
        self.init_sat_pos = np.array(init_sat_pos) * AU_to_m
        dt2 = launch_T * seconds_to_years
        initPos = self.findPositions(t0_launch)[:, 0]
        initVel = (self.findPositions(t0_launch + dt2)[:, 0] - initPos) / dt2
        self.init_plan_pos = initPos * AU_to_m
        self.init_plan_vel = initVel * AU_pr_year_to_m_pr_s
        self.launch_T = launch_T
        self.launch_N = int(launch_T * 100)
        self.t0_launch = t0_launch

    def escape_velocity_movie(self):
        if not self.engine_checked:
            print('You must first set engine settings!')
            sys.exit()
        mass_needed_launch(movie=True)

    def mass_needed_launch(self, final_pos = (0, 0), test = False, movie = False):
		"""
		Calculates how much boost you need.
		The launch direction is always radially upwards from the initial satellite position
		Variables:
		final_pos:  The final position calculated by the students in AU-units
		test:       Boolean to check the final calculated position after launch
		returns:    How much fuel you have used up in SI-units,
		            The time at the end of the launch, in years,
		            the final position, and final velocity in AU-units.
		"""
		T = self.launch_T
		N = self.launch_N
		t0 = self.t0_launch * years_to_seconds
		sat_pos = self.init_sat_pos
		final_pos = np.array(final_pos)
		planet_pos = self.init_plan_pos
		planet_vel = self.init_plan_vel
		boost_successful = False
		TakeOff = False
		G = 6.67408e-11
		A = self.dpdt1 * self.numBox
		B = self.partPerSec * mass_1_particle * self.numBox
		time = np.linspace(t0, t0 + T, N)
		dt = float(T) / (N - 1)
		e_b_vect = (sat_pos - planet_pos) / np.linalg.norm(sat_pos - planet_pos)
		v_0_boost = 0.0
		boost_vect = []
		boost_vect.append(v_0_boost)
		boost_check = 0
		initial_fuel_mass = self.fuel_left
		M = self.mass_sat + initial_fuel_mass
		dist_from_center = np.linalg.norm(sat_pos - planet_pos)
		planet_radius = self.radius[0] * 1000
		if T > 3600:
			print('WARNING: do you really need a launch lasting more than one hour?')
		if T > 86400:
			print("ERROR: Doesn't allow >24 hour launches.")
			sys.exit()
		if abs(dist_from_center - planet_radius) < planet_radius * 0.01:
			diff = planet_radius - dist_from_center
			dist_from_center = planet_radius
			if diff > 0:
				print('Rocket moved up by %f m to stand on planet surface.' % diff)
			elif diff < 0:
				print('Rocket moved down by %f m to stand on planet surface.' % -diff)
			else:
				print('Rocket starting exactly at planet surface!!')
		else:
			print('Launch position not at surface of planet! Do you know where the planet is?')
			print('Dist from center: %f m. Planet radius %f m.' % (dist_from_center, planet_radius))
			sys.exit()
		for i in range(N - 1):
			grav_acc = -G * self.mass[0] * sun_mass_to_kg / dist_from_center ** 2
			boost_check = 1.1*A / M * dt + grav_acc * dt
			if boost_check > 0:
				v_0_boost += boost_check
				dist_from_center += v_0_boost * dt
				boost_vect.append(v_0_boost)
				TakeOff = True
				M -= B * dt
				if M < self.mass_sat:
					print('You ran out of fuel. Try a larger amount')
					sys.exit()
				if TakeOff == True:
					kin_energy = 0.5 * M * v_0_boost ** 2
					pot_energy = G * self.mass[0] * sun_mass_to_kg * M / np.linalg.norm(dist_from_center)
					if kin_energy >= pot_energy:
						print(dist_from_center)
						fuel_needed = self.mass_sat + initial_fuel_mass - M
						boost_successful = True
						T = dt * (i + 1)
						break

				if TakeOff == False:
					print('You have a too weak engine. Gravity pulled the satellite back down.')
					sys.exit()
		if boost_successful == True:
			if T / 60.0 < 5:
			    print('The Launch should not take less than 5 minutes. Be realistic!')
			    sys.exit()
			print('Launch completed, reached escape velocity in %f seconds.' % T)
		else:
			print('Rocket never reached escape velocity, try launching over a longer time.')
			sys.exit()
		radial_unit_vector = (sat_pos - planet_pos) / np.linalg.norm(sat_pos - planet_pos)
		sat_pos = dist_from_center * radial_unit_vector
		sat_pos_copy = sat_pos.copy()
		w = 2.0 * np.pi * self.period[0] / 86400.0
		theta = w * T
		sat_pos[0] = sat_pos_copy[0] * np.cos(theta) - sat_pos_copy[1] * np.sin(theta)
		sat_pos[1] = sat_pos_copy[0] * np.sin(theta) + sat_pos_copy[1] * np.cos(theta)
		init_rot = np.array((-sat_pos[1], sat_pos[0]))
		final_rot_unit = init_rot / np.linalg.norm(init_rot)
		final_rot_vel = final_rot_unit * w * dist_from_center
		final_vel = v_0_boost * sat_pos / np.linalg.norm(sat_pos) + final_rot_vel + self.init_plan_vel
		print("HHHHHHHHHHHHHHHHHHHHEEEEEEEEEEEERRRRRRRRRRRR", final_vel)
		final_vel /= AU_pr_year_to_m_pr_s
		sat_pos += planet_vel * T
		sat_pos += planet_pos
		sat_pos *= m_to_AU
		if test == True:
			check_vector = final_pos - sat_pos
			print ("satpos=", sat_pos)
			print ("finalpos", final_pos)
			check_dist = np.linalg.norm(check_vector)
			norm_dist = planet_radius * m_to_AU
			if check_dist > 0.05 * norm_dist:
				print('Your final launch position deviates too much from the actual position. Try again')
				print('TIP: Make sure you have included the rotation and orbital velocity of your home planet.')
				print('Coriolis effect should NOT be included; Rocket is always directly above launch point.')
				print(check_dist)
			else:
				print('Your final launch position was satisfyingly calculated. Well done!')
				self.launch_accomplished = True
				return (fuel_needed,
				time[i] * seconds_to_years,
				sat_pos,
				final_vel)

    def mass_needed(self, v0, boost):
        """
        Calculates the amount of fuel needed for 1 boost.
        Variables:
        v0: initial velocity before boost, given as vector
        boost: the change in velocity, given as vector
        no_boxes: the number of fuel boxes you need
        A: dp/dt for 1 fuel box
        B: dm/dt for 1 fuel box
        no_boxes: number of boxes needed
        T = time in seconds for the boost to possibly take
        N = points in the time array
        mass_1_particle = calculated mass of one particle [kg]
        ###    n = total no. of particles in a single combustion chamber
        return: The amount of fuel used to achieve the boost
        """
        boost *= AU_pr_year_to_m_pr_s
        v0 *= AU_pr_year_to_m_pr_s
        A = self.dpdt1 * self.numBox
        B = self.partPerSec * mass_1_particle * self.numBox
        boost_ms = np.linalg.norm(boost)
        e_b_vect = boost / np.linalg.norm(boost)
        v_0_boost = 0.0
        M = self.mass_sat + self.fuel_left
        N = 100000.0
        dt1 = boost_ms * M / N / A
        dt2 = M / N / B
        dt = min(dt1, dt2)
        dt = max(dt, 1e-08)
        fuel_needed = 0.0
        T = 0.0
        while v_0_boost < boost_ms:
            v_0_boost += A / M * dt
            M -= B * dt
            fuel_needed += B * dt
            T += dt
            if M < 0:
                print('You ran very much out of fuel. Try a larger amount')
                sys.exit()

        self.fuel_left -= fuel_needed
        if self.fuelPrint:
            print('Fuel needed: %.2f kg. (Boost time: %.1f s)' % (fuel_needed, T))
        v0 /= AU_pr_year_to_m_pr_s
        return fuel_needed

    def test_random_number_generator(self):
        seed = 1000
        n = 10
        failure = False
        r1 = np.zeros(n, dtype=int)
        r2 = np.zeros(n)
        r3 = np.zeros(n)
        R1 = [777,
         669,
         99,
         352,
         467,
         534,
         978,
         130,
         671,
         364]
        R2 = [488.83570716,
         203.01221073,
         666.19837557,
         227.66303121,
         458.0640583,
         40.72239755,
         974.28979538,
         487.47607427,
         461.61386364,
         714.14715581]
        R3 = [998.19400883,
         1001.05693777,
         1001.87332251,
         1000.27614852,
         998.20890161,
         1000.35925887,
         1002.69118956,
         999.04239509,
         1000.69290386,
         1001.38523942]
        try:
            random.seed(seed, version=1)
        except:
            random.seed(seed)

        for i in range(n):
            r1[i] = randint(0, seed)

        for i in range(n):
            r2[i] = uniform(0, seed)

        for i in range(n):
            r3[i] = random.gauss(seed, 1.0)

        for i in range(n):
            if R1[i] != r1[i] or abs(R2[i] - r2[i]) > 1e-07 or abs(R3[i] - r3[i]) > 1e-07:
                failure = True

        random.seed(self.seed)
        if failure:
            print('Random number generator test failed!! Contact group teacher immediately!')
            return 1
        else:
            print('Random number generator test passed.')
            return 0

    def computeRotation(self, time, flooring = False):
        rotation = np.zeros(self.number_of_planets)
        for i in range(self.number_of_planets):
            R = self.period[i] / 365.24
            rotation[i] = time / R
            rotation[i] += 0.01618033988749895 * (i + 1)
            if rotation[i] > 1.0 and flooring:
                rotation[i] = rotation[i] - floor(rotation[i])
            rotation[i] = 2 * np.pi * rotation[i]

        return rotation

    @staticmethod
    def ang2pix(theta, phi, nside = 512):
        if nside < 1 or nside > 2050:
            print('nside out of range')
            sys.exit()
        if theta < 0.0 or theta > pi:
            print('theta', theta, 'out of range')
            sys.exit()
        z = cos(theta)
        za = abs(z)
        tt = phi % (2.0 * pi) / (pi / 2.0)
        ipix = 0
        if za <= 0.6666666666666666:
            temp1 = nside * (0.5 + tt)
            temp2 = nside * 0.75 * z
            jp = int(floor(temp1 - temp2))
            jm = int(floor(temp1 + temp2))
            ir = nside + 1 + jp - jm
            kshift = 1 - ir % 2
            nl4 = 4 * nside
            ip = int(floor((jp + jm - nside + kshift + 1) // 2))
            if ip >= nl4:
                ip = ip - nl4
            ipix = 2 * nside * (nside - 1) + nl4 * (ir - 1) + ip
        else:
            tp = tt - int(floor(tt))
            tmp = nside * sqrt(3.0 * (1.0 - za))
            jp = int(floor(tp * tmp))
            jm = int(floor((1.0 - tp) * tmp))
            ir = jp + jm + 1
            ip = int(floor(tt * ir))
            if ip >= 4 * ir:
                ip = ip - 4 * ir
            if z > 0.0:
                ipix = 2 * ir * (ir - 1) + ip
            else:
                ipix = 12 * nside ** 2 - 2 * ir * (ir + 1) + ip
        return ipix

    def generateImage(self, phii, filename = 'find_orient.png'):
        FOV = pi * 70.0 / 180.0
        xymaxmin = 2 * sin(FOV / 2.0) / (1.0 + cos(FOV / 2.0))
        xresolution = 640
        yresolution = 480
        outpic = np.zeros([yresolution, xresolution, 3], dtype=np.uint8)
        try:
            inFile = open('himmelkule.npy', 'rb')
            himmelkulen = np.load(inFile)
            inFile.close()
        except:
            print('Could not open himmelkule.npy. Please download from course webpage.')
            return

        xarray = np.linspace(-xymaxmin, xymaxmin, xresolution)
        yarray = np.linspace(xymaxmin, -xymaxmin, yresolution)
        for i in range(xresolution):
            x = xarray[i]
            for j in range(yresolution):
                y = yarray[j]
                rho = np.sqrt(x ** 2 + y ** 2)
                c = 2 * np.arctan(rho / 2.0)
                theta = pi / 2.0 - np.arcsin(y * np.sin(c) / rho)
                phi0 = phii * pi / 180.0
                phi = phi0 + np.arctan(x * np.sin(c) / (rho * np.cos(c)))
                pixnum = self.ang2pix(theta, phi)
                pix = himmelkulen[pixnum]
                rgb = [pix[2], pix[3], pix[4]]
                outpic[j, i] = rgb

        img2 = Image.fromarray(outpic)
        img2.save(filename)

    def getRefStars(self):
        print('Reference star 1, at phi = %f degrees,\nhas delta lambda of %16.12f nm in the home star rest frame.' % (self.RSA[0], self.make_shift(self.RSV[0] / (km_to_AU / 1000) * seconds_to_years)))
        print('Reference star 2, at phi = %f degrees,\nhas delta lambda of %16.12f nm in the home star rest frame.' % (self.RSA[1], self.make_shift(self.RSV[1] / (km_to_AU / 1000) * seconds_to_years)))

    def orientForRealz(self, time, xv):
        rad2deg = 180.0 / pi
        pos = self.findPositions(time)
        r = np.zeros(self.number_of_planets)
        planetBeta = 0
        planetTheta = 0
        for i in range(self.number_of_planets):
            r[i] = np.linalg.norm(pos[:, i] - xv[0:2])
            if r[i] < 0.01:
                relPos = np.array(pos[:, i]) - np.array(xv[0:2])
                planetBeta = arctan2(relPos[1], relPos[0])
                planetTheta = arctan(self.radius[i] / (r[i] / km_to_AU - self.radius[i]))

        sunTheta = (pi + arctan2(xv[1], xv[0])) * rad2deg
        a = (planetBeta + planetTheta) * rad2deg
        b = sunTheta - 19
        difference = b - a
        while difference < 0:
            difference += 360

        while difference > 360:
            difference -= 360

        if difference > 360 - 2 * planetTheta:
            difference = 360 - 2 * planetTheta
        if difference > 70:
            beta = 1.41 * self.seed * time
            amod = a + 35
            bmod = b - 35
            if amod > bmod:
                bmod += 360
            delta = bmod - amod
            beta = (beta - floor(beta)) * delta + a + 35
            while beta < 0:
                beta += 360

            while beta >= 360:
                beta -= 360

        else:
            beta = 1.41 * self.seed * time
            a = (planetBeta - planetTheta) * rad2deg
            b = sunTheta + 19
            amod = a - 35
            bmod = b + 35
            if bmod > amod:
                amod += 360
            delta = amod - bmod
            beta = (beta - floor(beta)) * delta + b + 35
            while beta < 0:
                beta += 360

            while beta >= 360:
                beta -= 360

        beta = randint(0, 360)
        self.generateImage(beta)
        print("Picture generated ('find_orient.png'). Compute the orientation angle.")
        ang = input('Orientation angle (degrees): ')
        if ang == 'juksepave':
            print('Okay, but only this time!')
        else:
            ang = float(ang)
            d = abs(ang - beta) % 360
            if d > 180:
                r0 = 360 - d
            else:
                r0 = d
            if abs(r0) < 2:
                print('Orientation angle correctly calculated.')
            else:
                print('Orientation angle incorrectly calculated. Exiting.')
                sys.exit(1)
        v = -xv[2] * cos(self.RSA[0] / rad2deg) - xv[3] * sin(self.RSA[0] / rad2deg) + self.RSV[0]
        v = v / (km_to_AU / 1000) * seconds_to_years
        dl = self.make_shift(v)
        print('Reference star 1, measured delta lambda of 656.3 nm line: %16.12f nm.' % dl)
        v = -xv[2] * cos(self.RSA[1] / rad2deg) - xv[3] * sin(self.RSA[1] / rad2deg) + self.RSV[1]
        v = v / (km_to_AU / 1000) * seconds_to_years
        dl = self.make_shift(v)
        print('Reference star 2, measured delta lambda of 656.3 nm line: %16.12f nm.' % dl)
        vx = eval(input('Input x-velocity component: '))
        vy = eval(input('Input y-velocity component: '))
        if abs(vx - xv[2]) < 0.01 and abs(vy - xv[3]) < 0.01:
            print('Velocity correctly calculated.')
        else:
            print('Velocity incorrectly calculated. Exiting.')
            sys.exit(1)
        posFile = open('pos.npy', 'wb')
        posF = np.zeros(self.number_of_planets + 1)
        posF[0:-1] = r
        posF[-1] = np.linalg.norm(xv[0:2])
        np.save(posFile, posF)
        posFile.close()
        print("Distance between the satellite and the planets (including the home star) dumped to file 'pos.npy', calculate the position of the satellite (the distance from the home star is the last element of the array).")
        satPosX = eval(input('Input x-position: '))
        satPosY = eval(input('Input y-position: '))
        if np.linalg.norm(np.array([satPosX, satPosY] - xv[0:2])) < 0.01:
            print('Position correctly calculated.')
            print('**********************************************')
            print('*Achievement unlocked: ORIENTATE THIS, BITCH!*')
            print('*  AU_to_matic orientation software unlocked.  *')
            print('**********************************************')
            statusFile = open(self.sFile, 'rb')
            s = np.load(statusFile)
            statusFile.close()
            s[1] = [1]
            statusFile = open(self.sFile, 'wb')
            np.save(statusFile, s)
            statusFile.close()
            return True
        else:
            print('Position incorrectly calcualted. Exiting.')
            sys.exit()
            return False

    def make_shift(self, v):
        lambda0 = 656.3
        lambda_center = v / 300000000.0 * lambda0 + lambda0
        return lambda_center - lambda0

    def make_spectrum(self, v, plot = False):
        lambda0 = 656.3
        lambda_center = v / 300000000.0 * lambda0 + lambda0
        lambda_min = 654.0
        lambda_max = 658.0
        nlambda = 100000
        sigma = 0.01
        fmin = 0.5
        seed = 34868
        rms = 0.01
        lambdas = np.linspace(lambda_min, lambda_max, nlambda)
        obs = (fmin - 1.0) * exp(-(lambdas - lambda_center) ** 2 / (2.0 * sigma ** 2)) + 1.0 - np.random.rand(nlambda) * rms
        if plot:
            plt.plot(lambdas, obs)
            plt.show()
        return (obs, lambdas)

    def make_junk(self, plot = False):
        lambda0 = 656.3
        lambda_min = 654.0
        lambda_max = 658.0
        nlambda = 100000
        rms = 0.01
        lambdas = np.linspace(lambda_min, lambda_max, nlambda)
        obs = np.random.rand(nlambda) * rms
        if plot:
            plt.plot(lambdas, obs)
            plt.show()
        return (obs, lambdas)

    def sneak_peek(self):
        xv = np.array([300, 0])
        command = np.array([pi / 2.0,
         0,
         0,
         0,
         1])
        self.dump_to_xml(0, xv, command, 0)
        print("XML for viewing the different planets generated ('selfie_0000.xml').")

    def write_to_xml(self, time_array, cam_pos, cam_dir, planet_pos = np.array([0, 0, 0]), other_objects = [], camera_messages = None, planet_messages = None, ruler = None, up_vec = np.array([0, 0, 1]), FOV = 90, filename = 'test_data.xml', chosen_planet = 0, c = c_AU, laser_scale = 1, use_obj_scaling = 1, cheat_light_speed = False, origo_location = np.array([1, 0, 0]), black_hole = False, play_speed = None):
        """
        All positions are given relative to a virtual origo, which is then moved to [1,0,0]AU before sent to MCAst.
        @ time_array  =  Array of the timepoints to simulate.
        @ cam_pos  =  Array of camera positions at given timepoints. Shape (nr_frames, 3)
                      or (3,) for static camera position.
        @ cam_dir  =  Array of camera directions. Either direction vector
                      with shape (nr_frames, 3) or direction angles
                      with shape (nr_frames, 2). Can also be (3,) and (2,) for static direction.
                      First angle is upwards/downwards (0 - pi), second angle is in plane (0 - 2pi).
        @ planet_pos  =  Static position of planet. Shape (nr_frames, 3) or (3,). Defaults to [0,0,0].
        @ other_objects  =  Nested list with info about all non-planet/rocket objects. Shape (nr_objects, 7)
                     [Object name, Object string, Object positions array, Size scale, Color, Message list, Sound list, Visibility list, Orientation]
                     @ Object name      =  Name of object (only shows up in SSView)
                     @ Object string    =  Object type. Valid input: "Satellite", "Sphere01", "Explosion", Laser.
                     @ Object positions array  =  Shape (nr_frames, 3)
                     @ Size scale       =  Shape scaling of object, scalar.
                     @ Color            =  RGB. Values between 0 and 1. Shape (3,)
                     @ Message list     =  List of object-specific messages each frame, Shape = (nr_frames). Send None for no messages.
                     @ Sound list       =  List of object specific sounds. Shape = (nr_frames). Send None for no sounds.
                     @ Visibility list  =  Visibility of object each frame, 0 or 1. Shape = (nr_frames). Send None for always visible.
                     @ Orientation      =  3-vector orientation of object. Optional, Auto-orient by default (Value None).
        @ camera_messages  =  messages displayed on screen. Shape = (nr_frames).
        @ planet_messages  =  messages displayed on planet. Shape = (nr_frames).
        @ ruler    = Information about ruler. List with [rulerStart, rulerEnd, ticks, rulerUnit].
                     If no argument is provided there will be no ruler.
        @ up_vec   = Upward direction vector of camera, Shape = (3,) or (nr_frames, 3)
        @ c        = Light speed in your system. If you want no componentwise
                     scaling of moving objects, simply set this to a reasonably large value. By default c in AU/s.
        @ laser_scale = Laser scale in the direction it is moving in.
        @ use_obj_scaling = 1 by default, if u dont want any objects scaled change to not 1
        @ origo_location = Location of events are sent relative to a local origo (e.g. center of planet).
                           This variable transitions that origo to a position relative to the sun.
        """
        self.cheat_light_speed = cheat_light_speed
        if self.debug:
            print('Entering SR XML writer with seed: %d and chosen planet: %d' % (self.seed, chosen_planet))
        nr_frames = len(time_array)
        time_array = np.copy(time_array)
        cam_pos = np.copy(cam_pos)
        cam_dir = np.copy(cam_dir)
        planet_pos = np.copy(planet_pos)
        other_objects = np.copy(other_objects)
        if time_array[-1] - time_array[0] >= 1:
            clock_time_array = np.copy(time_array)
            clockUnit = 'seconds'
        elif time_array[-1] - time_array[0] >= 0.001:
            clock_time_array = np.copy(time_array) * 1000.0
            clockUnit = 'milli seconds'
        elif time_array[-1] - time_array[0] >= 1e-06:
            clock_time_array = np.copy(time_array) * 1000000.0
            clockUnit = 'micro seconds'
        else:
            clock_time_array = np.copy(time_array) * 1000000000.0
            clockUnit = 'nano seconds'
        if np.shape(planet_pos) == (nr_frames, 3):
            if self.debug:
                print('Dynamic planet position.')
            planet_pos = np.array(planet_pos, dtype=np.float64)
            planet_pos += origo_location
        elif np.shape(planet_pos) == (3,):
            if self.debug:
                print('Static planet position. Castig time axis.')
            planet_pos = np.array(planet_pos, dtype=np.float64)
            planet_pos = np.zeros(shape=(nr_frames, 3)) + planet_pos + origo_location
        else:
            raise IndexError("Parameter 'planet_pos' has shape %s, expected %s or %s" % (np.shape(planet_pos), (3,), (nr_frames, 3)))
        if planet_messages is None:
            planet_messages = [ '' for i in range(nr_frames) ]
        if camera_messages is None:
            camera_messages = [ '' for i in range(nr_frames) ]
        if np.shape(cam_pos) == (nr_frames, 3):
            stat_cam_pos = False
            if self.debug:
                print('Dynamic camera position.')
            cam_pos += origo_location
        elif np.shape(cam_pos) == (3,):
            stat_cam_pos = True
            if self.debug:
                print('Stationary camera position. Introducing time-axis manually.')
            cam_pos = np.zeros(shape=(nr_frames, 3)) + cam_pos + origo_location
        else:
            raise IndexError("Parameter 'cam_pos' has shape %s, expected %s or %s" % (np.shape(cam_pos), (3,), (nr_frames, 3)))
        dir_shape = np.shape(cam_dir)
        if dir_shape == (nr_frames, 2):
            stat_cam_uses_angles = False
            cam_uses_angles = True
            if self.debug:
                print('Dynamic camera angles.')
        elif dir_shape == (2,):
            stat_cam_uses_angles = True
            cam_uses_angles = True
            if self.debug:
                print('Static camera angles. Introducing time-axis manually.')
            cam_dir = np.zeros(shape=(nr_frames, 2)) + cam_dir
        elif dir_shape == (nr_frames, 3):
            stat_cam_dir_vec = False
            cam_uses_angles = False
            if self.debug:
                print('Dynamic camera vector.')
        elif dir_shape == (3,):
            stat_cam_dir_vec = True
            cam_uses_angles = False
            if self.debug:
                print('Static camera vector. Introducing time-axis manually.')
            cam_dir = np.zeros(shape=(nr_frames, 3)) + cam_dir
        else:
            raise IndexError("Parameter 'cam_dir' has shape %s. Expected %s or %s for angles or %s or %s for vector. Exiting." % (np.shape(cam_dir),
             (nr_frames, 2),
             (2,),
             (nr_frames, 3),
             (3,)))
        if cam_uses_angles:
            cam_dir_vec = np.zeros(shape=(nr_frames, 3))
            cam_dir_vec[:, 0] = np.cos(cam_dir[:, 1]) * np.sin(cam_dir[:, 0])
            cam_dir_vec[:, 1] = np.sin(cam_dir[:, 1]) * np.sin(cam_dir[:, 0])
            cam_dir_vec[:, 2] = np.cos(cam_dir[:, 0])
        else:
            cam_dir_vec = np.array(np.copy(cam_dir), dtype=np.float64)
            cam_dir_vec /= np.linalg.norm(cam_dir_vec, axis=1)[:, None]
        if str(up_vec) == 'auto':
            sun_vec = planet_pos
            up_vec = np.zeros(shape=(nr_frames, 3))
            for i in range(nr_frames):
                up_vec[i] = np.cross(sun_vec, cam_dir_vec[i])

            if up_vec[(0, 2)] < 0:
                up_vec = -up_vec
        if str(up_vec) == 'up':
            up_vec = cam_pos - planet_pos
        if np.shape(up_vec) == (nr_frames, 3):
            up_vec = np.array(up_vec, dtype=np.float64)
            up_vec /= np.linalg.norm(up_vec, axis=1)[:, None]
        elif np.shape(up_vec) == (3,):
            up_vec = np.array(up_vec, dtype=np.float64)
            up_vec /= np.linalg.norm(up_vec)
            up_vec = np.zeros(shape=(nr_frames, 3)) + up_vec
        else:
            raise IndexError('Dust')
        if ruler is not None:
            if len(ruler) != 4:
                raise IndexError('Ruler settings list has length %g, expected 4.' % len(ruler))
            rulerStart = ruler[0]
            rulerEnd = ruler[1]
            rulerTicks = ruler[2]
            rulerUnit = ruler[3]
        objects = etree.Element('Objects')
        star = etree.SubElement(objects, 'SerializedMCAstObject')
        if black_hole is True:
            etree.SubElement(star, 'category').text = str('black hole')
        else:
            etree.SubElement(star, 'category').text = str('star')
        etree.SubElement(star, 'pos_x').text = str(0)
        etree.SubElement(star, 'pos_z').text = str(0)
        etree.SubElement(star, 'pos_y').text = str(0)
        etree.SubElement(star, 'rot_y').text = str(0)
        etree.SubElement(star, 'radius').text = str(100000)
        etree.SubElement(star, 'temperature').text = str(4000)
        etree.SubElement(star, 'seed').text = str(int(self.seed * 1000 + 990))
        etree.SubElement(star, 'atmosphereDensity').text = str(10)
        etree.SubElement(star, 'atmosphereHeight').text = str(1.025)
        etree.SubElement(star, 'outerRadiusScale').text = str(1.0025)
        etree.SubElement(star, 'name').text = str('The star')
        star_objects = etree.SubElement(star, 'Objects')
        planet_pos_array = np.zeros(shape=(nr_frames, 3))
        planet_pos_array[:] = planet_pos
        if use_obj_scaling == 1:
            if self.debug:
                print('Scaling planet.')
            planet_scales_array = 1 / self.get_lorentz(time_array, cam_pos, planet_pos_array, c, object_name='planet')
        else:
            planet_scales_array = np.ones(nr_frames)
        planet = etree.SubElement(star_objects, 'SerializedMCAstObject')
        etree.SubElement(planet, 'pos_x').text = str(planet_pos[(0, 0)])
        etree.SubElement(planet, 'pos_z').text = str(planet_pos[(0, 1)])
        etree.SubElement(planet, 'pos_y').text = str(planet_pos[(0, 2)])
        etree.SubElement(planet, 'rot_y').text = str(0)
        etree.SubElement(planet, 'seed').text = str(int(self.seed * 1000 + chosen_planet))
        etree.SubElement(planet, 'radius').text = str(self.radius[chosen_planet])
        etree.SubElement(planet, 'temperature').text = str(self.temperature * np.sqrt(self.star_radius * km_to_AU / (2 * self.a[chosen_planet])))
        etree.SubElement(planet, 'atmosphereDensity').text = str(np.log(self.rho0[chosen_planet]) / np.log(25))
        etree.SubElement(planet, 'atmosphereHeight').text = str(1.025)
        etree.SubElement(planet, 'outerRadiusScale').text = str(1.0025)
        etree.SubElement(planet, 'category').text = str('planet')
        etree.SubElement(planet, 'name').text = str('planet ' + str(chosen_planet))
        frames = etree.SubElement(planet, 'Frames')
        for i in range(nr_frames):
            frame = etree.SubElement(frames, 'Frame')
            etree.SubElement(frame, 'id').text = str(i)
            etree.SubElement(frame, 'pos_x').text = str(planet_pos[i, 0])
            etree.SubElement(frame, 'pos_z').text = str(planet_pos[i, 1])
            etree.SubElement(frame, 'pos_y').text = str(planet_pos[i, 2])
            etree.SubElement(frame, 'displayMessage').text = str(planet_messages[i])
            etree.SubElement(frame, 'rot_y').text = str(0)
            etree.SubElement(frame, 'scale_z').text = str(planet_scales_array[i])

        for other_object in other_objects:
            obj_name = other_object[0]
            obj_string = other_object[1]
            obj_pos_array = other_object[2] + origo_location
            obj_scale = other_object[3]
            obj_color = other_object[4]
            obj_message_list = other_object[5]
            obj_sound_list = other_object[6]
            obj_visible = other_object[7]
            obj_orientation = other_object[8]
            if obj_message_list is None:
                obj_message_list = [ '' for i in range(nr_frames) ]
            if obj_sound_list is None:
                obj_sound_list = [ '' for i in range(nr_frames) ]
            if obj_visible is None or obj_visible is True:
                obj_visible = [ 1 for i in range(nr_frames) ]
            if np.shape(obj_pos_array) == (3,):
                obj_pos_array = np.zeros(shape=(nr_frames, 3)) + obj_pos_array
            if self.debug:
                print('Making', obj_string, 'object')
            if obj_string == 'Laser':
                obj_scales_array = laser_scale * np.ones(nr_frames)
            elif use_obj_scaling == 1:
                obj_scales_array = 1 / self.get_lorentz(time_array, cam_pos, obj_pos_array, c, object_name=obj_name)
            else:
                obj_scales_array = np.ones(nr_frames)
            obj = etree.SubElement(star_objects, 'SerializedMCAstObject')
            if obj_orientation is None:
                etree.SubElement(obj, 'autoOrient').text = str(1)
            elif np.shape(obj_orientation) == (3,):
                etree.SubElement(obj, 'autoOrient').text = str(0)
            else:
                raise ValueError('Expected shape (3,) or None for Orientation parameter on object %s. Got shape %s' % (obj_name, np.shape(obj_orientation)))
            etree.SubElement(obj, 'pos_x').text = str(obj_pos_array[(0, 0)])
            etree.SubElement(obj, 'pos_z').text = str(obj_pos_array[(0, 1)])
            etree.SubElement(obj, 'pos_y').text = str(obj_pos_array[(0, 2)])
            if obj_string == 'explosion':
                etree.SubElement(obj, 'category').text = str('explosion')
            else:
                etree.SubElement(obj, 'category').text = str('3dobject')
            etree.SubElement(obj, 'name').text = str(obj_name)
            if obj_string == 'explosion':
                pass
            elif obj_string == 'Laser':
                etree.SubElement(obj, 'objectMaterial').text = str('LaserMaterial')
                etree.SubElement(obj, 'objectString').text = str(obj_string)
            else:
                etree.SubElement(obj, 'objectMaterial').text = str('HullMaterial')
                etree.SubElement(obj, 'objectString').text = str(obj_string)
            etree.SubElement(obj, 'objectScale').text = str(obj_scale)
            etree.SubElement(obj, 'color_r').text = str(obj_color[0])
            etree.SubElement(obj, 'color_g').text = str(obj_color[1])
            etree.SubElement(obj, 'color_b').text = str(obj_color[2])
            frames = etree.SubElement(obj, 'Frames')
            for i in range(nr_frames):
                frame = etree.SubElement(frames, 'Frame')
                etree.SubElement(frame, 'id').text = str(i)
                etree.SubElement(frame, 'pos_x').text = str(obj_pos_array[i, 0])
                etree.SubElement(frame, 'pos_z').text = str(obj_pos_array[i, 1])
                etree.SubElement(frame, 'pos_y').text = str(obj_pos_array[i, 2])
                if obj_sound_list is not None:
                    etree.SubElement(frame, 'sound').text = str(obj_sound_list[i])
                if obj_message_list[i] != '':
                    etree.SubElement(frame, 'displayMessage').text = str(obj_message_list[i])
                etree.SubElement(frame, 'time').text = str(clock_time_array[i])
                if obj_string == 'explosion':
                    etree.SubElement(frame, 'color_r').text = str(obj_color[0])
                    etree.SubElement(frame, 'color_g').text = str(obj_color[1])
                    etree.SubElement(frame, 'color_b').text = str(obj_color[2])
                if obj_string == 'explosion':
                    etree.SubElement(frame, 'scale_x').text = str(obj_scales_array[i])
                else:
                    etree.SubElement(frame, 'scale_z').text = str(obj_scales_array[i])
                if obj_visible[i] == 0:
                    etree.SubElement(frame, 'visible').text = str(0)
                if obj_orientation is not None:
                    etree.SubElement(frame, 'rot_x').text = str(obj_orientation[0])
                    etree.SubElement(frame, 'rot_z').text = str(obj_orientation[1])
                    etree.SubElement(frame, 'rot_y').text = str(obj_orientation[2])

        cameras = etree.Element('Cameras')
        for i in range(nr_frames):
            camera = etree.SubElement(cameras, 'SerializedCamera')
            etree.SubElement(camera, 'cam_x').text = str(cam_pos[i, 0])
            etree.SubElement(camera, 'cam_z').text = str(cam_pos[i, 1])
            etree.SubElement(camera, 'cam_y').text = str(cam_pos[i, 2])
            etree.SubElement(camera, 'dir_x').text = str(cam_dir_vec[i, 0])
            etree.SubElement(camera, 'dir_z').text = str(cam_dir_vec[i, 1])
            etree.SubElement(camera, 'dir_y').text = str(cam_dir_vec[i, 2])
            etree.SubElement(camera, 'up_x').text = str(up_vec[i, 0])
            etree.SubElement(camera, 'up_z').text = str(up_vec[i, 1])
            etree.SubElement(camera, 'up_y').text = str(up_vec[i, 2])
            etree.SubElement(camera, 'fov').text = str(FOV)
            etree.SubElement(camera, 'scale_x').text = str(1)
            etree.SubElement(camera, 'scale_z').text = str(1)
            etree.SubElement(camera, 'scale_y').text = str(1)
            etree.SubElement(camera, 'displayMessage').text = camera_messages[i]
            etree.SubElement(camera, 'time').text = str(clock_time_array[i])
            etree.SubElement(camera, 'frame').text = str(i)

        with open(self.path + filename, 'w') as outfile:
            outfile.write('<?xml version="1.0" encoding="utf-8"?>\n')
            outfile.write('<SerializedWorld xmlns:xsi="http://www.w3.org/2001/')
            outfile.write('XMLSchema-instance"\n xmlns:xsd="http://www.w3.org/2001/XMLSchema">\n')
            outfile.write('<sun_intensity>0.100</sun_intensity>\n')
            outfile.write('<screenshot_width>900</screenshot_width>\n')
            outfile.write('<screenshot_height>900</screenshot_height>\n')
            outfile.write('<global_radius_scale>0.985</global_radius_scale>\n')
            outfile.write('<uuid>5acbd644-37c7-11e6-ac61-9e71128cae77</uuid>\n')
            outfile.write('<skybox>000</skybox>\n')
            if play_speed is None:
                outfile.write('<defaultPlaySpeed>0.6</defaultPlaySpeed>')
            else:
                outfile.write('<defaultPlaySpeed>%f</defaultPlaySpeed>' % play_speed)
            if ruler is not None:
                outfile.write('<rulerStart>%g</rulerStart>\n' % rulerStart)
                outfile.write('<rulerEnd>%g</rulerEnd>\n' % rulerEnd)
                outfile.write('<rulerTicks>%g</rulerTicks>\n' % rulerTicks)
                outfile.write('<rulerUnit>%s</rulerUnit>\n' % rulerUnit)
            outfile.write('<clockUnit>%s</clockUnit>\n' % clockUnit)
            outfile.write(str(etree.tostring(objects, pretty_print=True)))
            outfile.write(str(etree.tostring(cameras, pretty_print=True)))
            outfile.write('</SerializedWorld>')
        print('Video file written to %s. Open and view it in MCast!' % (self.path + filename))

    def get_lorentz(self, time_array, cam_pos_array, obj_pos_array, c = c_AU, object_name = 'unknown'):
        """ Compute lorentz length contraction for an object moving relative
        to the camera for each frame.
        Params:
        @ time_array     = Array with time points, shape (nr_frames)
        @ cam_pos_array  = Array with camera positions, shape (nr_frames, 3)
        @ obj_pos_array  = Array with object positions, shape (nr_frames, 3)
        @ c              = Light speed in the system
        Returns:
        @ lorentz        = Array with lorentz factor for each frame, shape (nr_frames)
        """
        obj_pos_arr = np.copy(obj_pos_array)
        obj_pos_arr = obj_pos_arr - cam_pos_array
        nr_frames = len(time_array)
        vel_array = np.zeros(shape=(nr_frames, 3))
        vel_array[1:, 0] = (obj_pos_arr[1:, 0] - obj_pos_arr[:-1, 0]) / (time_array[1:] - time_array[:-1])
        vel_array[1:, 1] = (obj_pos_arr[1:, 1] - obj_pos_arr[:-1, 1]) / (time_array[1:] - time_array[:-1])
        vel_array[1:, 2] = (obj_pos_arr[1:, 2] - obj_pos_arr[:-1, 2]) / (time_array[1:] - time_array[:-1])
        vel_array[(0, 0)] = (obj_pos_arr[(1, 0)] - obj_pos_arr[(0, 0)]) / (time_array[1] - time_array[0])
        vel_array[(0, 1)] = (obj_pos_arr[(1, 1)] - obj_pos_arr[(0, 1)]) / (time_array[1] - time_array[0])
        vel_array[(0, 2)] = (obj_pos_arr[(1, 2)] - obj_pos_arr[(0, 2)]) / (time_array[1] - time_array[0])
        absvels = np.linalg.norm(vel_array, axis=1)
        maxvel = np.max(absvels)
        if self.cheat_light_speed is False:
            if maxvel > c:
                raise ValueError('Detected maximum speed of %g in object %s, light speed is %g' % (maxvel, object_name, c))
        if self.cheat_light_speed is True:
            for i in range(len(absvels)):
                if absvels[i] > c:
                    absvels[i] = 0

        lorentz = 1 / np.sqrt(1 - absvels ** 2 / c ** 2)
        return lorentz

    def lorentz_transform(self, t, pos, v, c = c_AU):
        """
        Transformation of position/time 4-vector of event from rest-frame to moving frame.
        To transform an event from moving frame to rest frame, send in negative v.
        This function now only takes in 1-dimmentional movement along the axis of movement.
        All objects must be positioned, and move in the same one-dimmentional axis.
        INPUT
        @ pos  =  Position of events in rest frame along axis of movement. Shape = (nr_events)
        @ t    =  Time of event in rest frame. Shape = (nr_events).
        @ v    =  Speed of moving frame along movement axis. Shape = (nr_events) or scalar.
        RETURNS
        @ pos_marked   =  Poisitions of events in moving frame, along axis of movement. Shape = (nr_events)
        @ t_marked     =  Time of events in moving frame. Shape = (nr_events)
        """
        if np.sum(np.abs(v) >= c) > 0:
            print('Speed v = %f cannot exceed the speed of light, c = %f!' % (v, c))
            sys.exit()
        gamma = 1.0 / np.sqrt(1 - v ** 2 / c ** 2)
        t_marked = -v * gamma * pos / c ** 2 + gamma * t
        pos_marked = gamma * pos - v * gamma * t
        return (t_marked, pos_marked)

    def focus_tanhyper(self, NR_frames_used, start = None, time = None):
        """
        Hyperbolic function to change focus or position smoothly
        Returns V array from 0 to 1 of two posible shapes(NR_frames,) or (time,)
        @ NR_frames_used        Number of frames used to change focus
        @ time                  if provided, shape(V) = (len(time),) IF u want V
                                to go over the hole time, otherwise shape(V) = (NR_frames_used,)
        @ start                 if provided, start = indices of time, where change of focus
                                should start, and last over NR_frames_used, otherwise start at begining
        """
        x = np.linspace(-2.5, 2.5, NR_frames_used)
        if np.shape(time) == ():
            V = np.tanh(x)
            V += abs(np.min(V))
            V = V / np.max(V)
        else:
            V = np.zeros([len(time)])
            if start is not None:
                V[start:start + NR_frames_used] = np.tanh(x)
                V[start:start + NR_frames_used] += abs(np.min(V))
                V[start + NR_frames_used:] += V[start + NR_frames_used - 1]
            else:
                V[:NR_frames_used] = np.tanh(x)
                V[:NR_frames_used] += abs(np.min(V))
                V[NR_frames_used:] += V[NR_frames_used - 1]
            V = V / np.max(V)
        return V

    def spherical_to_cart(self, r, phi, theta = '2D'):
        """
        Transforms from spherical to cartesian coordinates. 2D or 3D.
        @ r      =  Radius in sphere
        @ phi    =  Degrees in the plane [0,2*pi]
        @ theta  =  Degrees from zenit   [0,pi].  Leave empty for 2D.
        (Note: phi = pi points forwards the star)
        """
        if theta == '2D':
            x = r * cos(phi)
            y = r * sin(phi)
            return np.array([x, y])
        else:
            x = r * cos(phi) * sin(theta)
            y = r * sin(theta) * sin(phi)
            z = r * cos(theta)
            return np.array([x, y, z])

    def velocity_transformation(self, v_observer, v_object):
        """
        WARNING: Only takes natural units (fractions of c)
        @ v_observer  =  velocity of new observer relative to rest frame(old observer).
        @ v_object  =  velocity of observed object relative to rest frame(old observer).
        RETURN  =  velocity of object relative to the new observer.
        """
        return (v_object - v_observer) / (1 - v_observer * v_object)

    def get_ruler_length(self, distance_to_object, FOV = 90):
        """
        @ distance_to_object  =  distance to object we wish to measure with our ruler.
        @ FOV  =  camera field of view in degrees.
        RETURN  =  length of ruler in the same units as distance_to_object
        """
        return 2 * distance_to_object * np.tan(FOV / 2 * (pi / 180))

    def ref_sys_interpolate_func(self, time_array, pos_array, v):
        """
        @ time_array  =  time of positions in original frame of refrence.
        @ pos_array  =  corresponding positions in original frame of refrence.
        @ v  =  velocity of new frame of refrence relative to original.
        Units in AU and seconds.
        RETURNS  =  Callable function which returns positions of event in new refrence frame given a timepoint.
        """
        new_time_array, new_pos_array = self.lorentz_transform(time_array, pos_array, v)
        return interpolate.interp1d(new_time_array, new_pos_array, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)

    def relativistic_doppler_shift(self, v, c = c_AU):
        """
        wl : Wavelength (lambda in the lecture notes)
        Formula from part 2b.
        delta wl / wl = ( sqrt ([1 + v]/[1 - v]) - 1 )
        @ v      = Speed of object that emits light relative to observer in AU/s.
        @ c      = Speed of light in your system
        RETURNS   = delta wl / wl. (Relative change of wavelength)
        """
        return np.sqrt((1 + v) / (1 - v)) - 1

    def relativistic_doppler_find_vel(self, dl, l):
        """
        Given wavelength and wavelength shift, find velocity of object.
        dl = Delta lambda
        l  = lambda
        RETURNS: Speed in light speeds (v/c)
        """
        u = (dl / l + 1) ** 2
        return (u - 1) / (1 + u)

    def object_list(self, obj_name, obj_type, pos, scale = 1, color_rgb = [1, 1, 1], msg_list = None, sound_list = None, visible = None, orient = None):
        """ Formats input about an object to a list which can be sent to the writer"""
        return [obj_name,
         obj_type,
         pos,
         scale,
         color_rgb,
         msg_list,
         sound_list,
         visible,
         orient]

    def A1(self, chosen_planet, N = 256, filename1 = 'A1_1.xml', filename2 = 'A1_2.xml', increase_height = False):
        """
        @ chosen_planet  =  Planet at which events are happening.
        @ N  =  Number of frames.
        M\xc3\x85 TESTES FOR SLUTT-H\xc3\x98YDE (end_radius). HVIS SOLEN SER RAR UT I SCENE 2 S\xc3\x85 ER DU INNI ATMOSF\xc3\x86REN.
        """
        random.seed(self.seed + 19)
        atmos_frame = randint(int(0.6 * N), int(0.7 * N))
        crash_frame = randint(int(0.9 * N), int(0.95 * N))
        cam_speed = random.uniform(0.98 * c, 0.99 * c) * km_to_AU
        planet_radius = self.radius[chosen_planet] * km_to_AU
        start_radius = planet_radius * 2
        if increase_height is False:
            end_radius = planet_radius * 1.001
        else:
            end_radius = planet_radius * 1.01
        cam_start_pos = [-start_radius, 0, 0]
        cam_end_pos = [-end_radius, 0, 0]
        cam_pos = np.zeros(shape=(N, 3))
        cam_pos[:crash_frame, 0] = np.linspace(cam_start_pos[0], cam_end_pos[0], crash_frame)
        cam_pos[crash_frame:] += cam_end_pos
        cam_dir = np.zeros(shape=(N, 3)) + np.array([1, 0, 0])
        cam_dir[crash_frame:] += np.array([-1, 1, 0])
        camera_messages = [ '' for i in range(N) ]
        distance_to_planet_in_sat_frame = (end_radius - cam_pos[atmos_frame, 0]) * sqrt(1 - (cam_speed / c_AU) ** 2) * AU_to_km
        print(distance_to_planet_in_sat_frame, (end_radius - cam_pos[atmos_frame, 0]) * AU_to_km)
        for i in range(atmos_frame, atmos_frame + int(N / 10.0)):
            camera_messages[i] = 'Spaceship has entered the atmosphere and sent message.\nConstant velocity %f m/s relative to planet.\nDistance to ground %f km in our frame of refrence.' % (cam_speed * AU_to_m, distance_to_planet_in_sat_frame)

        for i in range(crash_frame, N):
            camera_messages[i] = 'Spaceship has crashed.'

        expl_vis = np.zeros(N)
        expl_vis[crash_frame:] += 1
        expl_sound = [ '' for i in range(N) ]
        expl_sound[crash_frame] = 'explosion'
        expl_pos = cam_pos + np.array([0, 10 * km_to_AU, 0])
        other_objects = [['Exp0',
          'explosion',
          expl_pos,
          1000,
          [1, 0, 0],
          None,
          expl_sound,
          expl_vis,
          None]]
        crash_time = (cam_end_pos[0] - cam_start_pos[0]) / cam_speed
        end_time = N / float(crash_frame) * crash_time
        time_array = np.linspace(0, end_time, N)
        self.write_to_xml(time_array, cam_pos, cam_dir, other_objects=other_objects, camera_messages=camera_messages, chosen_planet=chosen_planet, filename=filename1)
        random.seed(self.seed + 20)
        atmos_frame = random.randint(int(0.6 * N), int(0.7 * N))
        crash_frame = random.randint(int(0.9 * N), int(0.95 * N))
        sat_speed = random.uniform(0.98 * c, 0.99 * c) * km_to_AU
        planet_radius = self.radius[chosen_planet] * km_to_AU
        start_radius = planet_radius * 2
        end_radius = planet_radius * 1.05
        cam_pos = np.zeros(shape=3) + np.array([-end_radius, 0, 0])
        sat_pos = np.zeros(shape=(N, 3))
        sat_start_pos = [-start_radius, 0, 0]
        sat_end_pos = [-end_radius, 0, 0]
        sat_pos[:crash_frame, 0] = np.linspace(sat_start_pos[0], sat_end_pos[0], crash_frame)
        sat_pos[crash_frame:] += sat_end_pos
        cam_dir = [-1, 0, 0]
        camera_messages = [ '' for i in range(N) ]
        distance_to_planet_in_planet_frame = (end_radius - sat_pos[atmos_frame, 0]) * AU_to_km
        for i in range(atmos_frame, atmos_frame + int(N / 10.0)):
            camera_messages[i] = 'Spaceship has entered the atmosphere and sent message.\nConstant velocity %f m/s relative to planet.\nDistance to ground %f km in our frame of refrence.' % (sat_speed * AU_to_m, distance_to_planet_in_planet_frame)

        for i in range(crash_frame, N):
            camera_messages[i] = 'Spaceship has crashed.'

        expl_vis = np.zeros(N)
        expl_vis[crash_frame:] += 1
        expl_sound = [ '' for i in range(N) ]
        expl_sound[crash_frame] = 'explosion'
        expl_pos = sat_end_pos - np.array([10 * km_to_AU, 0, 0])
        other_objects = [['Sat0',
          'Satellite',
          sat_pos,
          0.05,
          [0.5, 0.5, 0.5],
          None,
          None,
          None,
          None], ['Exp0',
          'explosion',
          expl_pos,
          1000,
          [1, 0, 0],
          None,
          expl_sound,
          expl_vis,
          None]]
        crash_time = (sat_end_pos[0] - sat_start_pos[0]) / sat_speed
        end_time = N / float(crash_frame) * crash_time
        time_array = np.linspace(0, end_time, N)
        self.write_to_xml(time_array, cam_pos, cam_dir, other_objects=other_objects, camera_messages=camera_messages, chosen_planet=chosen_planet, filename=filename2)

    def A2(self, friend_seed, chosen_planet, isYourSystem, N = 1024, filename1 = 'A2_1.xml', filename2 = 'A2_2.xml', increase_height = False):
        """
        @ friend_seed  =  The other students seed.
        @ chosen_planet  =  The chosen planet for the events.
        @ isYourSystem  =  True/False, your solar system or your friends.
        @ N  =  Number of frames.
        """
        if self.debug:
            print("Task A2 with own seed %d, friend seed %d, chosen planet %d, and 'Own system' = %s" % (self.seed,
             friend_seed,
             chosen_planet,
             isYourSystem))
        low_sat_speed = [0.8, 0.84]
        high_sat_speed = [0.86, 0.9]
        if isYourSystem:
            planet_radius = self.radius[chosen_planet] * km_to_AU
            distance_to_events = 400 * km_to_AU
            cam_dir = np.array([0, 0, -1])
        else:
            SS_friend = AST2000SolarSystem(friend_seed)
            planet_radius = SS_friend.radius[chosen_planet] * km_to_AU
            distance_to_events = -400 * km_to_AU
            cam_dir = np.array([0, 0, 1])
        if increase_height is False:
            event_radius = planet_radius * 1.005
        else:
            event_radius = planet_radius * 1.01
        if self.seed > friend_seed:
            seed_for_random = self.seed
            random.seed(self.seed + 4)
            v_cam_planframe = uniform(high_sat_speed[0], high_sat_speed[1]) * c_AU
            random.seed(friend_seed + 4)
            v_friend_planframe = uniform(low_sat_speed[0], low_sat_speed[1]) * c_AU
        else:
            seed_for_random = friend_seed
            random.seed(self.seed + 4)
            v_cam_planframe = uniform(low_sat_speed[0], low_sat_speed[1]) * c_AU
            random.seed(friend_seed + 4)
            v_friend_planframe = uniform(high_sat_speed[0], high_sat_speed[1]) * c_AU
        v_friend_camframe = self.velocity_transformation(v_cam_planframe, v_friend_planframe)
        if self.debug:
            print('Our velocity in planet frame = %gc\nFriend velocity in planet frame = %gc\nFriend velocity in our frame = %gc' % (v_cam_planframe, v_friend_planframe, v_friend_camframe))
        if not isYourSystem:
            self.seed = friend_seed
        sat_movement_length_planframe = 2000 * km_to_AU
        end_time_planframe = sat_movement_length_planframe / max([v_cam_planframe, v_friend_planframe])
        time_array_planframe = np.linspace(0, end_time_planframe, N)
        cam_pos_1D_camframe = np.zeros(N)
        friend_pos_1D_camframe = np.linspace(0, end_time_planframe * v_friend_camframe, N)
        cam_pos_1D_planframe = np.linspace(0, end_time_planframe * v_cam_planframe, N)
        friend_pos_1D_planframe = np.linspace(0, end_time_planframe * v_friend_planframe, N)
        end_time_camframe, _ = self.lorentz_transform(time_array_planframe[-1], cam_pos_1D_planframe[-1], v_cam_planframe)
        end_time_friendframe, _ = self.lorentz_transform(time_array_planframe[-1], cam_pos_1D_planframe[-1], v_friend_planframe)
        time_array_camframe = np.linspace(0, end_time_camframe, N)
        time_array_friendframe = np.linspace(0, end_time_friendframe, N)
        planet_pos_1D_camframe = np.linspace(0, -end_time_planframe * v_cam_planframe, N)
        cam_pos = np.zeros(shape=(N, 3))
        cam_pos[:, 2] += distance_to_events
        friend_pos_camframe = np.zeros(shape=(N, 3))
        friend_pos_camframe[:, 1] = friend_pos_1D_camframe
        friend_pos_camframe[:, 2] -= distance_to_events
        planet_pos_camframe = np.zeros(shape=(N, 3))
        planet_pos_camframe[:, 0] += event_radius
        planet_pos_camframe[:, 1] = planet_pos_1D_camframe
        ball1_pos_camframe = np.zeros(shape=(N, 3))
        ball1_visibility = np.zeros(shape=N)
        random.seed(seed_for_random)
        events_index_planframe = randint(N // 4, 3 * N // 4)
        print(events_index_planframe)
        events_time_planframe = time_array_planframe[events_index_planframe]
        avg_sat_pos_at_event = (cam_pos_1D_planframe[events_index_planframe] + friend_pos_1D_planframe[events_index_planframe]) / 2.0
        print(cam_pos_1D_planframe[events_index_planframe] * AU_to_km, friend_pos_1D_planframe[events_index_planframe] * AU_to_km)
        event1_pos_1D_planframe = avg_sat_pos_at_event + 110 * km_to_AU
        event2_pos_1D_planframe = avg_sat_pos_at_event - 110 * km_to_AU
        event1_time_camframe, event1_pos_1D_camframe = self.lorentz_transform(events_time_planframe, event1_pos_1D_planframe, v_cam_planframe)
        event2_time_camframe, event2_pos_1D_camframe = self.lorentz_transform(events_time_planframe, event2_pos_1D_planframe, v_cam_planframe)
        event1_time_friendframe, event1_pos_1D_friendframe = self.lorentz_transform(events_time_planframe, event1_pos_1D_planframe, v_friend_planframe)
        event2_time_friendframe, event2_pos_1D_friendframe = self.lorentz_transform(events_time_planframe, event2_pos_1D_planframe, v_friend_planframe)
        print(event1_pos_1D_camframe * AU_to_km, event1_time_camframe, event2_pos_1D_camframe * AU_to_km, event2_time_camframe)
        print(event1_pos_1D_friendframe * AU_to_km, event1_time_friendframe, event2_pos_1D_friendframe * AU_to_km, event2_time_friendframe)
        event1_index_camframe = np.abs(event1_time_camframe - time_array_planframe).argmin()
        event2_index_camframe = np.abs(event2_time_camframe - time_array_planframe).argmin()
        event1_visibility = np.zeros(shape=N)
        event2_visibility = np.zeros(shape=N)
        event1_visibility[event1_index_camframe:event1_index_camframe + N // 16] += 1
        event2_visibility[event2_index_camframe:event2_index_camframe + N // 16] += 1
        camera_messages = [ '' for i in range(N) ]
        event1_messages = [ '' for i in range(N) ]
        event2_messages = [ '' for i in range(N) ]
        event1_pos_camframe = np.zeros(shape=(N, 3))
        event2_pos_camframe = np.zeros(shape=(N, 3))
        event1_pos_camframe[:, 1] = event1_pos_1D_camframe
        event2_pos_camframe[:, 1] = event2_pos_1D_camframe
        event3_index_planframe = randint(N // 8, N // 4)
        event4_index_planframe = randint(3 * N // 4, 7 * N // 8)
        event3_time_planframe = time_array_planframe[event3_index_planframe]
        event4_time_planframe = time_array_planframe[event4_index_planframe]
        avg_sat_pos_at_event3 = (cam_pos_1D_planframe[event3_index_planframe] + friend_pos_1D_planframe[event3_index_planframe]) / 2.0
        avg_sat_pos_at_event4 = (cam_pos_1D_planframe[event4_index_planframe] + friend_pos_1D_planframe[event4_index_planframe]) / 2.0
        event3_pos_1D_planframe = avg_sat_pos_at_event3 - 70 * km_to_AU
        event4_pos_1D_planframe = avg_sat_pos_at_event4 + 90 * km_to_AU
        event3_time_camframe, event3_pos_1D_camframe = self.lorentz_transform(event3_time_planframe, event3_pos_1D_planframe, v_cam_planframe)
        event4_time_camframe, event4_pos_1D_camframe = self.lorentz_transform(event4_time_planframe, event4_pos_1D_planframe, v_cam_planframe)
        event3_index_camframe = np.abs(event3_time_camframe - time_array_planframe).argmin()
        event4_index_camframe = np.abs(event4_time_camframe - time_array_planframe).argmin()
        event3_visibility = np.zeros(shape=N)
        event4_visibility = np.zeros(shape=N)
        event3_visibility[event3_index_camframe:event3_index_camframe + N // 16] += 1
        event4_visibility[event4_index_camframe:event4_index_camframe + N // 16] += 1
        event3_messages = [ '' for i in range(N) ]
        event4_messages = [ '' for i in range(N) ]
        event3_pos_camframe = np.zeros(shape=(N, 3))
        event4_pos_camframe = np.zeros(shape=(N, 3))
        event3_pos_camframe[:, 1] = event3_pos_1D_camframe
        event4_pos_camframe[:, 1] = event4_pos_1D_camframe
        camera_messages = [ camera_messages[i] + 'Velocity of our satellite in planet frame = %fc' % (v_cam_planframe / c_AU) for i in range(N) ]
        camera_messages = [ camera_messages[i] + '\nVelocity of friendly satellite in our frame = %fc' % (v_friend_camframe / c_AU) for i in range(N) ]
        sat_messages = [ '\n\n\nFriend velocity relative to us = %fc' % (v_friend_camframe / c_AU) for i in range(N) ]
        ruler_length = self.get_ruler_length(distance_to_events * 2, 90) * AU_to_km
        ruler = [-ruler_length / 2.0,
         ruler_length / 2.0,
         20,
         'km']
        other_objects = [['Friend',
          'Satellite',
          friend_pos_camframe,
          0.08,
          [1, 1, 1],
          sat_messages,
          None,
          None,
          [0, -1, 0]],
         ['Ball1',
          'Sphere01',
          ball1_pos_camframe,
          10,
          [1, 1, 1],
          None,
          None,
          ball1_visibility,
          None],
         ['Event1',
          'explosion',
          event1_pos_camframe,
          400,
          [1, 0, 1],
          event1_messages,
          None,
          event1_visibility,
          None],
         ['Event2',
          'explosion',
          event2_pos_camframe,
          400,
          [0, 1, 1],
          event2_messages,
          None,
          event2_visibility,
          None],
         ['Event3',
          'explosion',
          event3_pos_camframe,
          400,
          [1, 1, 0],
          event3_messages,
          None,
          event3_visibility,
          None],
         ['Event4',
          'explosion',
          event4_pos_camframe,
          400,
          [0, 1, 0],
          event4_messages,
          None,
          event4_visibility,
          None]]
        if isYourSystem:
            SS = AST2000SolarSystem(self.seed, debug=self.debug, data_path=self.path)
        else:
            SS = AST2000SolarSystem(friend_seed, debug=self.debug, data_path=self.path)
        SS.write_to_xml(time_array_planframe, cam_pos, cam_dir, planet_pos_camframe, other_objects=other_objects, camera_messages=camera_messages, chosen_planet=chosen_planet, ruler=ruler, up_vec=[-1, 0, 0], filename=filename1, play_speed=0.55)
        for i in range(N):
            print(time_array_planframe[i], time_array_camframe[i])

        sat1_pos_planframe = np.zeros(shape=(N, 3))
        sat2_pos_planframe = np.zeros(shape=(N, 3))
        sat1_pos_planframe[:, 1] = cam_pos_1D_planframe
        sat2_pos_planframe[:, 1] = friend_pos_1D_planframe
        sat1_pos_planframe[:, 0] -= event_radius
        sat2_pos_planframe[:, 0] -= event_radius
        sat1_pos_planframe[:, 2] -= distance_to_events
        sat2_pos_planframe[:, 2] += distance_to_events
        event1_pos_planframe = np.zeros(shape=(N, 3))
        event2_pos_planframe = np.zeros(shape=(N, 3))
        event3_pos_planframe = np.zeros(shape=(N, 3))
        event4_pos_planframe = np.zeros(shape=(N, 3))
        event1_pos_planframe[:, 1] += event1_pos_1D_planframe
        event2_pos_planframe[:, 1] += event2_pos_1D_planframe
        event3_pos_planframe[:, 1] += event3_pos_1D_planframe
        event4_pos_planframe[:, 1] += event4_pos_1D_planframe
        event1_pos_planframe[:, 0] -= event_radius
        event2_pos_planframe[:, 0] -= event_radius
        event3_pos_planframe[:, 0] -= event_radius
        event4_pos_planframe[:, 0] -= event_radius
        event1_visibility2 = [0] * N
        event2_visibility2 = [0] * N
        event3_visibility2 = [0] * N
        event4_visibility2 = [0] * N
        event1_visibility2[events_index_planframe:events_index_planframe + N // 16] = [1] * (N // 16)
        event2_visibility2[events_index_planframe:events_index_planframe + N // 16] = [1] * (N // 16)
        event3_visibility2[event3_index_planframe:event3_index_planframe + N // 16] = [1] * (N // 16)
        event4_visibility2[event4_index_planframe:event4_index_planframe + N // 16] = [1] * (N // 16)
        cam_pos2 = np.zeros(shape=(N, 3))
        cam_pos2[:, 0] -= event_radius + 2500 * km_to_AU
        cam_dir2 = [1, 0, 0]
        planet_pos2 = [0, 0, 0]
        ruler_length = self.get_ruler_length(2500, FOV=90)
        ruler = [-ruler_length / 2,
         ruler_length / 2,
         20,
         'km']
        other_objects = [['Sat1',
          'Satellite',
          sat1_pos_planframe,
          0.16,
          [1, 1, 1],
          None,
          None,
          None,
          None],
         ['Sat2',
          'Satellite',
          sat2_pos_planframe,
          0.16,
          [1, 1, 1],
          None,
          None,
          None,
          None],
         ['Event1',
          'explosion',
          event1_pos_planframe,
          400,
          [1, 0, 1],
          event1_messages,
          None,
          event1_visibility2,
          None],
         ['Event2',
          'explosion',
          event2_pos_planframe,
          400,
          [0, 1, 1],
          event2_messages,
          None,
          event2_visibility2,
          None],
         ['Event3',
          'explosion',
          event3_pos_planframe,
          400,
          [1, 1, 0],
          event3_messages,
          None,
          event3_visibility2,
          None],
         ['Event4',
          'explosion',
          event4_pos_planframe,
          400,
          [0, 1, 0],
          event4_messages,
          None,
          event4_visibility2,
          None]]
        self.write_to_xml(time_array_planframe, cam_pos2, cam_dir2, planet_pos2, other_objects=other_objects, chosen_planet=chosen_planet, up_vec=[0, 0, 1], filename=filename2, ruler=ruler, play_speed=0.55)

    def A3(self, chosen_planet, N = 400, filename1 = 'A3_1.xml', filename2 = 'A3_2.xml'):
        """ Two spaceships moving relative to the planet. Both ships
        simultaneously fire lasers at each other, and explode simultaneously in
        their frame of reference.
        Scene 1: From the spaceships' frame of reference.
        Scene 2: From the planet's frame of reference."""
        random.seed(self.seed + 14)
        ref_frame_speed = random.uniform(0.87 * c_AU, 0.92 * c_AU)
        ref_frame_movement = np.linspace(-1000, 1000, N) * km_to_AU
        events_radius = self.radius[chosen_planet] * km_to_AU * 1.06
        spaceship_distance = 1200.0 * km_to_AU
        dist_from_star = np.hypot(self.x0[chosen_planet], self.y0[chosen_planet])
        planet_pos = np.array([0, 0, 0])
        rocket1_pos_1D = np.zeros(N) - spaceship_distance / 2.0
        rocket2_pos_1D = np.zeros(N) + spaceship_distance / 2.0
        cam_pos_1D = np.zeros(N)
        laser1_pos_1D = np.zeros(N)
        laser2_pos_1D = np.zeros(N)
        end_time = (ref_frame_movement[-1] - ref_frame_movement[0]) / ref_frame_speed
        time_array = np.linspace(0, end_time, N)
        dt = time_array[1] - time_array[0]
        dx = (time_array[-1] - time_array[-2]) * c * km_to_AU
        laser1_pos_1D[0] = rocket1_pos_1D[0]
        laser2_pos_1D[0] = rocket2_pos_1D[0]
        for i in range(N - 1):
            if laser1_pos_1D[i] > rocket2_pos_1D[i] or laser2_pos_1D[i] < rocket1_pos_1D[i]:
                laser1_pos_1D[i:] = laser1_pos_1D[i]
                laser2_pos_1D[i:] = laser2_pos_1D[i]
                expl_index = i + 1
                lasers_visible = np.ones(N)
                lasers_visible[expl_index:] = 0
                explosions_visible = np.copy(lasers_visible)
                explosions_visible = 1 - explosions_visible
                break
            else:
                laser1_pos_1D[i + 1] = laser1_pos_1D[i] + dx
                laser2_pos_1D[i + 1] = laser2_pos_1D[i] - dx

        cam_pos = np.zeros(shape=(N, 3))
        rocket1_pos = np.zeros(shape=(N, 3))
        rocket2_pos = np.zeros(shape=(N, 3))
        laser1_pos = np.zeros(shape=(N, 3))
        laser2_pos = np.zeros(shape=(N, 3))
        cam_pos[:, 2] = cam_pos_1D + ref_frame_movement
        rocket1_pos[:, 2] = rocket1_pos_1D + ref_frame_movement
        rocket2_pos[:, 2] = rocket2_pos_1D + ref_frame_movement
        laser1_pos[:, 2] = laser1_pos_1D + ref_frame_movement
        laser2_pos[:, 2] = laser2_pos_1D + ref_frame_movement
        laser1_pos[:, 2] -= 50 * km_to_AU
        laser2_pos[:, 2] -= 50 * km_to_AU
        cam_pos[:, 0] -= events_radius
        rocket1_pos[:, 0] -= events_radius
        rocket2_pos[:, 0] -= events_radius
        laser1_pos[:, 0] -= events_radius
        laser2_pos[:, 0] -= events_radius
        cam_pos[:, 1] -= 650 * km_to_AU
        cam_dir = [0, 1, 0]
        up_vec = [-1, 0, 0]
        other_objects = [['lsr1',
          'Laser',
          laser1_pos,
          4,
          [1, 1, 0],
          None,
          None,
          lasers_visible,
          None],
         ['lsr1',
          'Laser',
          laser2_pos,
          4,
          [1, 0, 0],
          None,
          None,
          lasers_visible,
          None],
         ['rocket1',
          'Satellite',
          rocket1_pos,
          0.1,
          [1, 1, 1],
          None,
          None,
          lasers_visible,
          None],
         ['rocket2',
          'Satellite',
          rocket2_pos,
          0.1,
          [1, 1, 1],
          None,
          None,
          lasers_visible,
          None],
         ['boom1',
          'explosion',
          rocket1_pos,
          1000,
          [1, 1, 1],
          None,
          None,
          explosions_visible,
          None],
         ['boom2',
          'explosion',
          rocket2_pos,
          1000,
          [1, 1, 1],
          None,
          None,
          explosions_visible,
          None]]
        ruler_length = self.get_ruler_length(650)
        ruler_start = 0
        ruler_end = ruler_length
        ticks = 11
        rulerUnit = 'km'
        ruler_s1 = [ruler_start,
         ruler_end,
         ticks,
         rulerUnit]
        self.write_to_xml(time_array, cam_pos, cam_dir, planet_pos, other_objects, ruler=ruler_s1, chosen_planet=chosen_planet, up_vec=up_vec, laser_scale=0.2, filename=filename1)
        rocket1_s2_func = self.ref_sys_interpolate_func(time_array, rocket1_pos_1D, -ref_frame_speed)
        rocket2_s2_func = self.ref_sys_interpolate_func(time_array, rocket2_pos_1D, -ref_frame_speed)
        laser1_s2_func = self.ref_sys_interpolate_func(time_array, laser1_pos_1D, -ref_frame_speed)
        laser2_s2_func = self.ref_sys_interpolate_func(time_array, laser2_pos_1D, -ref_frame_speed)
        t1 = self.lorentz_transform(time_array, rocket1_pos_1D, -ref_frame_speed)[0]
        gamma = 1 / np.sqrt(1 - ref_frame_speed ** 2 / c_AU ** 2)
        startTime = 0
        endTime = t1[-1] - t1[0]
        endTime *= 1.5
        times_s2 = np.linspace(startTime, endTime, N)
        t2 = self.lorentz_transform(time_array, rocket2_pos_1D, -ref_frame_speed)[0]
        rocket1_pos_1D_s2 = rocket1_s2_func(times_s2)
        rocket2_pos_1D_s2 = rocket2_s2_func(times_s2)
        laser1_pos_1D_s2 = laser1_s2_func(times_s2)
        laser2_pos_1D_s2 = laser2_s2_func(times_s2)
        time_diff_ticks = t2[0] - t1[0]
        dist_between_rockets_s2 = rocket2_pos_1D_s2[0] - rocket1_pos_1D_s2[0]
        dist_movement_s2 = ref_frame_speed * endTime
        ref_frame_movement_s2 = np.linspace(-dist_movement_s2 / 2, dist_movement_s2 / 2, N)
        rocket1_s2 = np.zeros(shape=(N, 3))
        rocket2_s2 = np.zeros(shape=(N, 3))
        cam_pos_s2 = np.zeros(shape=(N, 3))
        pos_of_ref_frame = np.zeros(shape=(N, 3))
        rocket1_s2[:, 0] -= events_radius
        rocket2_s2[:, 0] -= events_radius
        cam_pos_s2[:, 0] -= events_radius
        cam_pos_s2[:, 1] -= 3000 * km_to_AU
        rocket1_s2[:, 2] = ref_frame_movement_s2 - dist_between_rockets_s2 / 2
        rocket2_s2[:, 2] = ref_frame_movement_s2 + dist_between_rockets_s2 / 2
        pos_of_ref_frame[:] = 0.5 * rocket1_s2 + 0.5 * rocket2_s2
        cam_dir_s2 = pos_of_ref_frame - cam_pos_s2
        cam_dir_s2[:] = cam_dir_s2[int(N / 2), :]
        laser1_s2 = np.zeros(shape=(N, 3))
        laser2_s2 = np.zeros(shape=(N, 3))
        laser1_s2[0] = rocket1_s2[0]
        laser2_s2[0] = rocket2_s2[0]
        dt = times_s2[1] - times_s2[0]
        dx = c_AU * dt
        index_laser2_emission = np.argmin(np.abs(times_s2 - time_diff_ticks))
        rocket1_visible = np.ones(N)
        rocket2_visible = np.ones(N)
        laser2_visible = np.zeros(N)
        laser1_s2[:, 0] = rocket2_s2[:, 0]
        laser2_s2[:, 0] = rocket2_s2[:, 0]
        laser2_s2[:, 2] = rocket2_s2[:, 2]
        for i in range(N - 1):
            laser1_s2[i + 1, 2] = laser1_s2[i, 2] + dx
            if laser1_s2[i + 1, 2] > rocket2_s2[i + 1, 2]:
                rocket2_visible[i + 1] = 0
            if i >= index_laser2_emission:
                laser2_s2[i + 1, 2] = laser2_s2[i, 2] - dx
                laser2_visible[i + 1] = 1
            if laser2_s2[i + 1, 2] < rocket1_s2[i + 1, 2]:
                rocket1_visible[i + 1] = 0

        laser2_visible = laser2_visible - (1 - rocket1_visible)
        ship_scale_s2 = 0.5
        lasers_scale_s2 = 20
        other_objects_s2 = [self.object_list('rocket1', 'Satellite', rocket1_s2, ship_scale_s2, [1, 1, 1], visible=rocket1_visible),
         self.object_list('rocket2', 'Satellite', rocket2_s2, ship_scale_s2, [1, 1, 1], visible=rocket2_visible),
         self.object_list('laser1', 'Laser', laser1_s2, lasers_scale_s2, [1, 1, 0], visible=rocket2_visible),
         self.object_list('laser2', 'Laser', laser2_s2, lasers_scale_s2, [1, 0, 0], visible=laser2_visible),
         self.object_list('Booom1', 'explosion', rocket1_s2, 1000, [0, 1, 1], visible=1 - rocket1_visible),
         self.object_list('Boom2', 'explosion', rocket2_s2, 1000, [0, 1, 1], visible=1 - rocket2_visible)]
        self.write_to_xml(times_s2, cam_pos_s2, cam_dir_s2, planet_pos, other_objects_s2, ruler=ruler_s1, chosen_planet=chosen_planet, up_vec=up_vec, laser_scale=0.2, filename=filename2)

    def A4(self, chosen_planet, N = 1024, filename1 = 'A4_1.xml', filename2 = 'A4_2.xml'):
        distance_to_events = 300
        sat1_messages_satframe = [ '' for i in range(N) ]
        sat2_messages_satframe = [ '' for i in range(N) ]
        sat1_messages_planframe = [ '' for i in range(N) ]
        sat2_messages_planframe = [ '' for i in range(N) ]
        sat1_sounds_satframe = [ '' for i in range(N) ]
        sat2_sounds_satframe = [ '' for i in range(N) ]
        sat1_sounds_planframe = [ '' for i in range(N) ]
        sat2_sounds_planframe = [ '' for i in range(N) ]
        sat1_color = [1, 0, 0]
        sat2_color = [0, 0, 1]
        planet_radius = self.radius[chosen_planet] * km_to_AU
        random.seed(self.seed + 10)
        sat_frame_speed = uniform(0.65 * c, 0.65 * c) * km_to_AU
        sat_frame_movement = np.linspace(0, 1400, N) * km_to_AU
        events_radius = planet_radius * 1.5
        spaceship_distance = 400.0 * km_to_AU
        ball_speed_satframe = c_AU * 0.8
        sat1_pos_1D_satframe = np.zeros(N) - spaceship_distance / 2.0
        sat2_pos_1D_satframe = np.zeros(N) + spaceship_distance / 2.0
        cam_pos_1D_satframe = np.zeros(N)
        ball_pos_1D_satframe = np.zeros(shape=N)
        end_time = (sat_frame_movement[-1] - sat_frame_movement[0]) / sat_frame_speed
        time_array = np.linspace(0, end_time, N)
        dx = (time_array[-1] - time_array[-2]) * ball_speed_satframe
        ball_pos_1D_satframe[0] = sat1_pos_1D_satframe[0]
        ball_moving_forward = True
        for i in range(1, N):
            if ball_moving_forward:
                ball_pos_1D_satframe[i] = ball_pos_1D_satframe[i - 1] + dx
            else:
                ball_pos_1D_satframe[i] = ball_pos_1D_satframe[i - 1] - dx
            if ball_pos_1D_satframe[i] > sat2_pos_1D_satframe[i]:
                if self.debug:
                    print('ball bounced backwards at frame %d, [%es, %em] in satellite frame' % (i, time_array[i], ball_pos_1D_satframe[i]))
                ball_moving_forward = False
                ball_pos_1D_satframe[i] -= abs(ball_pos_1D_satframe[i] - sat2_pos_1D_satframe[i])
                sat2_messages_satframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat2_sounds_satframe[i] = 'bounce'
            if ball_pos_1D_satframe[i] < sat1_pos_1D_satframe[i]:
                if self.debug:
                    print('ball bounced forwards at frame %d, [%es, %em] in satellite frame' % (i, time_array[i], ball_pos_1D_satframe[i]))
                ball_moving_forward = True
                ball_pos_1D_satframe[i] += abs(ball_pos_1D_satframe[i] - sat1_pos_1D_satframe[i])
                sat1_messages_satframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat1_sounds_satframe[i] = 'bounce'

        cam_pos_satframe = np.zeros(shape=(N, 3))
        sat1_pos_satframe = np.zeros(shape=(N, 3))
        sat2_pos_satframe = np.zeros(shape=(N, 3))
        ball_pos_satframe = np.zeros(shape=(N, 3))
        cam_pos_satframe[:, 1] = cam_pos_1D_satframe + sat_frame_movement
        sat1_pos_satframe[:, 1] = sat1_pos_1D_satframe + sat_frame_movement
        sat2_pos_satframe[:, 1] = sat2_pos_1D_satframe + sat_frame_movement
        ball_pos_satframe[:, 1] = ball_pos_1D_satframe + sat_frame_movement
        cam_pos_satframe[:, 0] -= events_radius
        sat1_pos_satframe[:, 0] -= events_radius
        sat2_pos_satframe[:, 0] -= events_radius
        ball_pos_satframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_pos_satframe[:, 0] -= distance_to_events * km_to_AU
        cam_dir = [1, 0, 0]
        up_vec = [0, 0, 1]
        other_objects = [['Ball',
          'Sphere01',
          ball_pos_satframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None], ['sat1',
          'Satellite',
          sat1_pos_satframe,
          0.1,
          sat1_color,
          sat1_messages_satframe,
          sat1_sounds_satframe,
          None,
          [0, -1, 0]], ['sat2',
          'Satellite',
          sat2_pos_satframe,
          0.1,
          sat2_color,
          sat2_messages_satframe,
          sat2_sounds_satframe,
          None,
          [0, 1, 0]]]
        ruler_length = self.get_ruler_length(distance_to_events)
        ruler = [0,
         ruler_length,
         15,
         'km']
        self.write_to_xml(time_array, cam_pos_satframe, cam_dir, other_objects=other_objects, chosen_planet=chosen_planet, up_vec=up_vec, filename=filename1, ruler=ruler)
        ball_speed_planframe_forward = self.velocity_transformation(-sat_frame_speed / c_AU, ball_speed_satframe / c_AU) * c_AU
        ball_speed_planframe_backward = self.velocity_transformation(-sat_frame_speed / c_AU, -ball_speed_satframe / c_AU) * c_AU
        time_array_planframe1, sat1_pos_1D_planframe = self.lorentz_transform(time_array, sat1_pos_1D_satframe, -sat_frame_speed)
        time_array_planframe2, sat2_pos_1D_planframe = self.lorentz_transform(time_array, sat2_pos_1D_satframe, -sat_frame_speed)
        time_array_planframe = np.linspace(0, time_array_planframe1[-1] - time_array_planframe1[0], N)
        sat1_pos_1D_planframe = interpolate.interp1d(time_array_planframe1, sat1_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_array_planframe)
        sat2_pos_1D_planframe = interpolate.interp1d(time_array_planframe2, sat2_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_array_planframe)
        ball_pos_1D_planframe = np.zeros(shape=N)
        dx_forward = (time_array_planframe[-1] - time_array_planframe[-2]) * ball_speed_planframe_forward
        dx_backward = (time_array_planframe[-1] - time_array_planframe[-2]) * ball_speed_planframe_backward
        ball_pos_1D_planframe[0] = sat1_pos_1D_planframe[0]
        ball_moving_forward = True
        for i in range(1, N):
            if ball_moving_forward:
                ball_pos_1D_planframe[i] = ball_pos_1D_planframe[i - 1] + dx_forward
            else:
                ball_pos_1D_planframe[i] = ball_pos_1D_planframe[i - 1] + dx_backward
            if ball_pos_1D_planframe[i] > sat2_pos_1D_planframe[i]:
                if self.debug:
                    print('ball bounced backwards at frame %d, [%es, %em] in planet frame' % (i, time_array_planframe[i], ball_pos_1D_planframe[i]))
                ball_moving_forward = False
                ball_pos_1D_planframe[i] -= abs(ball_pos_1D_planframe[i] - sat2_pos_1D_planframe[i])
                sat2_messages_planframe[i:i + int(N / 16)] = [ 'Bounce! x = %.4fkm' % (ball_pos_1D_planframe[i] * AU_to_km) for j in range(int(N / 16)) ]
                sat2_sounds_planframe[i] = 'bounce'
            if ball_pos_1D_planframe[i] < sat1_pos_1D_planframe[i]:
                if self.debug:
                    print('ball bounced forwards at frame %d, [%es, %em] in planet frame' % (i, time_array_planframe[i], ball_pos_1D_planframe[i]))
                ball_moving_forward = True
                ball_pos_1D_planframe[i] += abs(ball_pos_1D_planframe[i] - sat1_pos_1D_planframe[i])
                sat1_messages_planframe[i:i + int(N / 16)] = [ 'Bounce! x = %.4fkm' % (ball_pos_1D_planframe[i] * AU_to_km) for j in range(int(N / 16)) ]
                sat1_sounds_planframe[i] = 'bounce'

        cam_pos_planframe = np.zeros(shape=(N, 3))
        sat1_pos_planframe = np.zeros(shape=(N, 3))
        sat2_pos_planframe = np.zeros(shape=(N, 3))
        ball_pos_planframe = np.zeros(shape=(N, 3))
        sat1_pos_planframe[:, 1] = sat1_pos_1D_planframe
        sat2_pos_planframe[:, 1] = sat2_pos_1D_planframe
        ball_pos_planframe[:, 1] = ball_pos_1D_planframe
        cam_pos_planframe[:, 0] -= events_radius + 1200 * km_to_AU
        cam_pos_planframe[:, 2] += uniform(-500, 500) * km_to_AU
        sat1_pos_planframe[:, 0] -= events_radius
        sat2_pos_planframe[:, 0] -= events_radius
        ball_pos_planframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_dir_planframe = sat1_pos_planframe - cam_pos_planframe
        up_vec_planframe = [0, 0, 1]
        other_objects_planframe = [['lsr',
          'Sphere01',
          ball_pos_planframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None], ['sat1',
          'Satellite',
          sat1_pos_planframe,
          0.1,
          sat1_color,
          sat1_messages_planframe,
          sat1_sounds_planframe,
          None,
          [0, -1, 0]], ['sat2',
          'Satellite',
          sat2_pos_planframe,
          0.1,
          sat2_color,
          sat2_messages_planframe,
          sat2_sounds_planframe,
          None,
          [0, 1, 0]]]
        self.write_to_xml(time_array_planframe, cam_pos_planframe, cam_dir_planframe, other_objects=other_objects_planframe, chosen_planet=chosen_planet, up_vec=up_vec_planframe, laser_scale=0.1, filename=filename2)

    def A5(self, chosen_planet, N = 1024, filename1 = 'A5_1.xml', filename2 = 'A5_2.xml'):
        sat1_messages_satframe = [ '' for i in range(N) ]
        sat2_messages_satframe = [ '' for i in range(N) ]
        sat1_messages_planframe = [ '' for i in range(N) ]
        sat2_messages_planframe = [ '' for i in range(N) ]
        sat1_sounds_satframe = [ '' for i in range(N) ]
        sat2_sounds_satframe = [ '' for i in range(N) ]
        sat1_sounds_planframe = [ '' for i in range(N) ]
        sat2_sounds_planframe = [ '' for i in range(N) ]
        sat1_color = [1, 0, 0]
        sat2_color = [0, 0, 1]
        planet_radius = self.radius[chosen_planet] * km_to_AU
        random.seed(self.seed + 10)
        sat_frame_speed = uniform(0.65 * c, 0.65 * c) * km_to_AU
        sat_frame_movement = np.linspace(0, 1400, N) * km_to_AU
        events_radius = planet_radius * 2
        spaceship_distance = 400.0 * km_to_AU
        sat1_pos_1D_satframe = np.zeros(N) - spaceship_distance / 2.0
        sat2_pos_1D_satframe = np.zeros(N) + spaceship_distance / 2.0
        cam_pos_1D_satframe = np.zeros(N)
        laser_pos_1D_satframe = np.zeros(shape=N)
        end_time = (sat_frame_movement[-1] - sat_frame_movement[0]) / sat_frame_speed
        time_array = np.linspace(0, end_time, N)
        dx = (time_array[-1] - time_array[-2]) * c * km_to_AU
        laser_pos_1D_satframe[0] = sat1_pos_1D_satframe[0]
        laser_moving_forward = True
        for i in range(1, N):
            if laser_moving_forward:
                laser_pos_1D_satframe[i] = laser_pos_1D_satframe[i - 1] + dx
            else:
                laser_pos_1D_satframe[i] = laser_pos_1D_satframe[i - 1] - dx
            if laser_pos_1D_satframe[i] > sat2_pos_1D_satframe[i]:
                if self.debug:
                    print('Laser bounced backwards at frame %d, [%es, %em] in satellite frame' % (i, time_array[i], laser_pos_1D_satframe[i]))
                laser_moving_forward = False
                laser_pos_1D_satframe[i] -= abs(laser_pos_1D_satframe[i] - sat2_pos_1D_satframe[i])
                sat2_messages_satframe[i:i + int(N / 16)] = [ 'Bounce!' for j in range(int(N / 16)) ]
                sat2_sounds_satframe[i] = 'laser'
            if laser_pos_1D_satframe[i] < sat1_pos_1D_satframe[i]:
                if self.debug:
                    print('Laser bounced forwards at frame %d, [%es, %em] in satellite frame' % (i, time_array[i], laser_pos_1D_satframe[i]))
                laser_moving_forward = True
                laser_pos_1D_satframe[i] += abs(laser_pos_1D_satframe[i] - sat1_pos_1D_satframe[i])
                sat2_messages_satframe[i:i + int(N / 16)] = [ 'Bounce!' for j in range(int(N / 16)) ]
                sat1_sounds_satframe[i] = 'laser'

        cam_pos_satframe = np.zeros(shape=(N, 3))
        sat1_pos_satframe = np.zeros(shape=(N, 3))
        sat2_pos_satframe = np.zeros(shape=(N, 3))
        laser_pos_satframe = np.zeros(shape=(N, 3))
        cam_pos_satframe[:, 1] = cam_pos_1D_satframe + sat_frame_movement
        sat1_pos_satframe[:, 1] = sat1_pos_1D_satframe + sat_frame_movement
        sat2_pos_satframe[:, 1] = sat2_pos_1D_satframe + sat_frame_movement
        laser_pos_satframe[:, 1] = laser_pos_1D_satframe + sat_frame_movement
        cam_pos_satframe[:, 0] -= events_radius
        sat1_pos_satframe[:, 0] -= events_radius
        sat2_pos_satframe[:, 0] -= events_radius
        laser_pos_satframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_pos_satframe[:, 0] -= 300 * km_to_AU
        cam_dir = [1, 0, 0]
        up_vec = [0, 0, 1]
        other_objects = [['lsr',
          'Laser',
          laser_pos_satframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None], ['sat1',
          'Satellite',
          sat1_pos_satframe,
          0.1,
          sat1_color,
          sat1_messages_satframe,
          sat1_sounds_satframe,
          None,
          [0, -1, 0]], ['sat2',
          'Satellite',
          sat2_pos_satframe,
          0.1,
          sat2_color,
          sat2_messages_satframe,
          sat2_sounds_satframe,
          None,
          [0, 1, 0]]]
        self.write_to_xml(time_array, cam_pos_satframe, cam_dir, other_objects=other_objects, chosen_planet=chosen_planet, up_vec=up_vec, laser_scale=0.1, filename=filename1, play_speed=0.6)
        time_array_planframe1, sat1_pos_1D_planframe = self.lorentz_transform(time_array, sat1_pos_1D_satframe, -sat_frame_speed)
        time_array_planframe2, sat2_pos_1D_planframe = self.lorentz_transform(time_array, sat2_pos_1D_satframe, -sat_frame_speed)
        time_array_planframe = np.linspace(0, time_array_planframe1[-1] - time_array_planframe1[0], N)
        sat1_pos_1D_planframe = interpolate.interp1d(time_array_planframe1, sat1_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_array_planframe)
        sat2_pos_1D_planframe = interpolate.interp1d(time_array_planframe2, sat2_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_array_planframe)
        laser_pos_1D_planframe = np.zeros(shape=N)
        dx = (time_array_planframe[-1] - time_array_planframe[-2]) * c * km_to_AU
        laser_pos_1D_planframe[0] = sat1_pos_1D_planframe[0]
        laser_moving_forward = True
        for i in range(1, N):
            if laser_moving_forward:
                laser_pos_1D_planframe[i] = laser_pos_1D_planframe[i - 1] + dx
            else:
                laser_pos_1D_planframe[i] = laser_pos_1D_planframe[i - 1] - dx
            if laser_pos_1D_planframe[i] > sat2_pos_1D_planframe[i]:
                if self.debug:
                    print('Laser bounced backwards at frame %d, [%es, %em] in planet frame' % (i, time_array_planframe[i], laser_pos_1D_planframe[i]))
                laser_moving_forward = False
                laser_pos_1D_planframe[i] -= abs(laser_pos_1D_planframe[i] - sat2_pos_1D_planframe[i])
                sat2_messages_planframe[i:i + int(N / 16)] = [ 'Bounce! x = %.4fkm' % (laser_pos_1D_planframe[i] * AU_to_km) for j in range(int(N / 16)) ]
                sat2_sounds_planframe[i] = 'laser'
            if laser_pos_1D_planframe[i] < sat1_pos_1D_planframe[i]:
                if self.debug:
                    print('Laser bounced forwards at frame %d, [%es, %em] in planet frame' % (i, time_array_planframe[i], laser_pos_1D_planframe[i]))
                laser_moving_forward = True
                laser_pos_1D_planframe[i] += abs(laser_pos_1D_planframe[i] - sat1_pos_1D_planframe[i])
                sat1_messages_planframe[i:i + int(N / 16)] = [ 'Bounce! x = %.4fkm' % (laser_pos_1D_planframe[i] * AU_to_km) for j in range(int(N / 16)) ]
                sat1_sounds_planframe[i] = 'laser'

        cam_pos_planframe = np.zeros(shape=(N, 3))
        sat1_pos_planframe = np.zeros(shape=(N, 3))
        sat2_pos_planframe = np.zeros(shape=(N, 3))
        laser_pos_planframe = np.zeros(shape=(N, 3))
        sat1_pos_planframe[:, 1] = sat1_pos_1D_planframe
        sat2_pos_planframe[:, 1] = sat2_pos_1D_planframe
        laser_pos_planframe[:, 1] = laser_pos_1D_planframe
        cam_pos_planframe[:, 0] -= events_radius + 1200 * km_to_AU
        cam_pos_planframe[:, 2] += uniform(-500, 500) * km_to_AU
        sat1_pos_planframe[:, 0] -= events_radius
        sat2_pos_planframe[:, 0] -= events_radius
        laser_pos_planframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_dir_planframe = sat1_pos_planframe - cam_pos_planframe
        up_vec_planframe = [0, 0, 1]
        other_objects_planframe = [['lsr',
          'Laser',
          laser_pos_planframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None], ['sat1',
          'Satellite',
          sat1_pos_planframe,
          0.1,
          sat1_color,
          sat1_messages_planframe,
          sat1_sounds_planframe,
          None,
          [0, -1, 0]], ['sat2',
          'Satellite',
          sat2_pos_planframe,
          0.1,
          sat2_color,
          sat2_messages_planframe,
          sat2_sounds_planframe,
          None,
          [0, 1, 0]]]
        self.write_to_xml(time_array_planframe, cam_pos_planframe, cam_dir_planframe, other_objects=other_objects_planframe, chosen_planet=chosen_planet, up_vec=up_vec_planframe, laser_scale=0.1, filename=filename2, play_speed=0.6)

    def A6(self, chosen_planet = 0, N = 500, height = 100):
        """
        Scene 1:
        Sat frame: set satelittes to be still, and planet moving to get right orientation for laser. This makes sat frame the rest frame!
        Scene 2:
        Planet frame: Moving frame.
        """
        random.seed(self.seed + 13)
        planet_radius = self.radius[chosen_planet] * km_to_AU
        plan_frame_speed = uniform(0.8 * c, 0.9 * c) * km_to_AU
        ball_speed = uniform(0.92 * c, 0.98 * c) * km_to_AU
        plan_frame_movement = np.linspace(-1000, 1000, N) * km_to_AU
        events_radius = planet_radius * 1.05 + height * km_to_AU
        spaceship_distance = 200.0 * km_to_AU
        sat1_messages_satframe = [ '' for i in range(N) ]
        sat2_messages_satframe = [ '' for i in range(N) ]
        sat1_messages_planframe = [ '' for i in range(N) ]
        sat2_messages_planframe = [ '' for i in range(N) ]
        sat1_sounds_satframe = [ '' for i in range(N) ]
        sat2_sounds_satframe = [ '' for i in range(N) ]
        sat1_sounds_planframe = [ '' for i in range(N) ]
        sat2_sounds_planframe = [ '' for i in range(N) ]
        endtime = (plan_frame_movement[-1] - plan_frame_movement[0]) / plan_frame_speed
        time_sat_frame = np.linspace(0, endtime, N)
        sat1_pos_satframe = np.zeros((N, 3))
        sat2_pos_satframe = np.zeros((N, 3))
        cam_pos_satframe = np.zeros((N, 3))
        ball_pos_satframe = np.zeros((N, 3))
        plan_pos_satframe = np.zeros((N, 3))
        sat2_pos_satframe[:, 2] += spaceship_distance / 2
        sat1_pos_satframe[:, 2] -= spaceship_distance / 2
        dz = (time_sat_frame[-1] - time_sat_frame[-2]) * ball_speed
        ball_pos_satframe[0, :] = sat1_pos_satframe[0, :]
        ball_moving_upward = True
        for i in range(1, N):
            if ball_moving_upward:
                ball_pos_satframe[i, 2] = ball_pos_satframe[i - 1, 2] + dz
            else:
                ball_pos_satframe[i, 2] = ball_pos_satframe[i - 1, 2] - dz
            if ball_pos_satframe[i, 2] > sat2_pos_satframe[i, 2]:
                ball_moving_upward = False
                ball_pos_satframe[i, 2] -= abs(ball_pos_satframe[i, 2] - sat2_pos_satframe[i, 2])
                sat2_messages_satframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat2_sounds_satframe[i] = 'laser'
            if ball_pos_satframe[i, 2] < sat1_pos_satframe[i, 2]:
                ball_moving_upward = True
                ball_pos_satframe[i, 2] += abs(ball_pos_satframe[i, 2] - sat1_pos_satframe[i, 2])
                sat1_messages_satframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat1_sounds_satframe[i] = 'laser'

        plan_pos_satframe[:, 1] += plan_frame_movement
        sat1_pos_satframe[:, 0] -= events_radius
        sat2_pos_satframe[:, 0] -= events_radius
        ball_pos_satframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_pos_satframe[:, 0] -= events_radius + 500 * km_to_AU
        sat1_color = [1, 0.5, 0]
        sat2_color = [0, 0.5, 1]
        obj_list_satframe = [['sat1',
          'Satellite',
          sat1_pos_satframe,
          0.1,
          sat1_color,
          sat1_messages_satframe,
          sat1_sounds_satframe,
          None,
          None,
          None], ['sat2',
          'Satellite',
          sat2_pos_satframe,
          0.1,
          sat2_color,
          sat2_messages_satframe,
          sat2_sounds_satframe,
          None,
          None,
          None], ['lsr',
          'Sphere01',
          ball_pos_satframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None]]
        camdir_satframe = np.array([1, 0, 0])
        upvec_satframe = np.array([0, 0, 1])
        self.write_to_xml(time_sat_frame, cam_pos_satframe, camdir_satframe, planet_pos=plan_pos_satframe, other_objects=obj_list_satframe, chosen_planet=chosen_planet, filename='A6_1.xml')
        time_planframe1, sat1_pos_1D_planframe = self.lorentz_transform(time_sat_frame, sat1_pos_satframe[:, 1], plan_frame_speed)
        time_planframe2, sat2_pos_1D_planframe = self.lorentz_transform(time_sat_frame, sat2_pos_satframe[:, 1], plan_frame_speed)
        time_planframe = np.linspace(0, time_planframe1[-1] - time_planframe1[0], N)
        sat1_pos_1D_planframe = interpolate.interp1d(time_planframe1, sat1_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_planframe)
        sat2_pos_1D_planframe = interpolate.interp1d(time_planframe2, sat2_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_planframe)
        sat1_pos_planframe = np.zeros((N, 3))
        sat2_pos_planframe = np.zeros((N, 3))
        ball_pos_planframe = np.zeros((N, 3))
        cam_pos_planframe = np.zeros((N, 3))
        sat1_pos_planframe[:, 1] = sat1_pos_1D_planframe
        sat2_pos_planframe[:, 1] = sat2_pos_1D_planframe
        sat1_pos_planframe[:, 2] -= spaceship_distance / 2
        sat2_pos_planframe[:, 2] += spaceship_distance / 2
        v_z = sqrt(ball_speed ** 2 - plan_frame_speed ** 2)
        dz = (time_planframe[-1] - time_planframe[-2]) * v_z
        dy = (time_planframe[-1] - time_planframe[-2]) * plan_frame_speed
        ball_pos_planframe[0, :] = sat1_pos_planframe[0, :]
        ball_moving_upward = True
        for i in range(1, N):
            ball_pos_planframe[i, 1] = ball_pos_planframe[i - 1, 1] - dy
            if ball_moving_upward:
                ball_pos_planframe[i, 2] = ball_pos_planframe[i - 1, 2] + dz
            else:
                ball_pos_planframe[i, 2] = ball_pos_planframe[i - 1, 2] - dz
            if ball_pos_planframe[i, 2] > sat2_pos_planframe[i, 2]:
                if self.debug:
                    print('ball bounced downwards at frame %d, [%es, %em] in satellite frame' % (i, time_sat_frame[i], ball_pos_planframe[i, 2]))
                ball_moving_upward = False
                ball_pos_planframe[i, 2] -= abs(ball_pos_planframe[i, 2] - sat2_pos_planframe[i, 2])
                sat2_messages_planframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat2_sounds_planframe[i] = 'laser'
            if ball_pos_planframe[i, 2] < sat1_pos_planframe[i, 2]:
                if self.debug:
                    print('ball bounced upwards at frame %d, [%es, %em] in satellite frame' % (i, time_sat_frame[i], ball_pos_planframe[i, 2]))
                ball_moving_upward = True
                ball_pos_planframe[i, 2] += abs(ball_pos_planframe[i, 2] - sat1_pos_planframe[i, 2])
                sat1_messages_planframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat1_sounds_planframe[i] = 'laser'

        sat1_pos_planframe[:, 0] -= events_radius
        sat2_pos_planframe[:, 0] -= events_radius
        ball_pos_planframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_pos_planframe[:, 0] -= planet_radius + height * km_to_AU
        sat1_color = [1, 0.5, 0]
        sat2_color = [0, 0.5, 1]
        obj_list_planframe = [['sat1',
          'Satellite',
          sat1_pos_planframe,
          0.1,
          sat1_color,
          sat1_messages_planframe,
          sat1_sounds_planframe,
          None,
          None], ['sat2',
          'Satellite',
          sat2_pos_planframe,
          0.1,
          sat2_color,
          sat2_messages_planframe,
          sat2_sounds_planframe,
          None,
          None], ['ball',
          'Sphere01',
          ball_pos_planframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None]]
        camdir_planframe = np.copy(sat1_pos_planframe - cam_pos_planframe)
        camdir_planframe[:, 2] += spaceship_distance / 2
        upvec_planframe = np.array([0, 0, 1])
        self.write_to_xml(time_planframe, cam_pos_planframe, camdir_planframe, other_objects=obj_list_planframe, chosen_planet=chosen_planet, filename='A6_2.xml')

    def A7(self, chosen_planet = 0, N = 500, height = 100):
        """
        Scene 1:
        Sat frame: set satelittes to be still, and planet moving to get right orientation for laser. This makes sat frame the rest frame!
        Scene 2:
        Planet frame: Moving frame.
        """
        random.seed(self.seed + 12)
        planet_radius = self.radius[chosen_planet] * km_to_AU
        plan_frame_speed = uniform(0.85 * c, 0.92 * c) * km_to_AU
        plan_frame_movement = np.linspace(-1000, 1000, N) * km_to_AU
        events_radius = planet_radius * 1.05 + height * km_to_AU
        spaceship_distance = 200.0 * km_to_AU
        sat1_messages_satframe = [ '' for i in range(N) ]
        sat2_messages_satframe = [ '' for i in range(N) ]
        sat1_messages_planframe = [ '' for i in range(N) ]
        sat2_messages_planframe = [ '' for i in range(N) ]
        sat1_sounds_satframe = [ '' for i in range(N) ]
        sat2_sounds_satframe = [ '' for i in range(N) ]
        sat1_sounds_planframe = [ '' for i in range(N) ]
        sat2_sounds_planframe = [ '' for i in range(N) ]
        endtime = (plan_frame_movement[-1] - plan_frame_movement[0]) / plan_frame_speed
        time_sat_frame = np.linspace(0, endtime, N)
        sat1_pos_satframe = np.zeros((N, 3))
        sat2_pos_satframe = np.zeros((N, 3))
        cam_pos_satframe = np.zeros((N, 3))
        laser_pos_satframe = np.zeros((N, 3))
        plan_pos_satframe = np.zeros((N, 3))
        sat2_pos_satframe[:, 2] += spaceship_distance / 2
        sat1_pos_satframe[:, 2] -= spaceship_distance / 2
        dz = (time_sat_frame[-1] - time_sat_frame[-2]) * c * km_to_AU
        laser_pos_satframe[0, :] = sat1_pos_satframe[0, :]
        laser_moving_upward = True
        for i in range(1, N):
            if laser_moving_upward:
                laser_pos_satframe[i, 2] = laser_pos_satframe[i - 1, 2] + dz
            else:
                laser_pos_satframe[i, 2] = laser_pos_satframe[i - 1, 2] - dz
            if laser_pos_satframe[i, 2] > sat2_pos_satframe[i, 2]:
                laser_moving_upward = False
                laser_pos_satframe[i, 2] -= abs(laser_pos_satframe[i, 2] - sat2_pos_satframe[i, 2])
                sat2_messages_satframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat2_sounds_satframe[i] = 'laser'
            if laser_pos_satframe[i, 2] < sat1_pos_satframe[i, 2]:
                laser_moving_upward = True
                laser_pos_satframe[i, 2] += abs(laser_pos_satframe[i, 2] - sat1_pos_satframe[i, 2])
                sat1_messages_satframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat1_sounds_satframe[i] = 'laser'

        plan_pos_satframe[:, 1] += plan_frame_movement
        sat1_pos_satframe[:, 0] -= events_radius
        sat2_pos_satframe[:, 0] -= events_radius
        laser_pos_satframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_pos_satframe[:, 0] -= events_radius + 500 * km_to_AU
        sat1_color = [1, 0.5, 0]
        sat2_color = [0, 0.5, 1]
        obj_list_satframe = [['sat1',
          'Satellite',
          sat1_pos_satframe,
          0.1,
          sat1_color,
          sat1_messages_satframe,
          sat1_sounds_satframe,
          None,
          None], ['sat2',
          'Satellite',
          sat2_pos_satframe,
          0.1,
          sat2_color,
          sat2_messages_satframe,
          sat2_sounds_satframe,
          None,
          None], ['lsr',
          'Laser',
          laser_pos_satframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None]]
        camdir_satframe = np.array([1, 0, 0])
        upvec_satframe = np.array([0, 0, 1])
        self.write_to_xml(time_sat_frame, cam_pos_satframe, camdir_satframe, planet_pos=plan_pos_satframe, other_objects=obj_list_satframe, chosen_planet=chosen_planet, filename='A7_1.xml', laser_scale=0.1)
        time_planframe1, sat1_pos_1D_planframe = self.lorentz_transform(time_sat_frame, sat1_pos_satframe[:, 1], plan_frame_speed)
        time_planframe2, sat2_pos_1D_planframe = self.lorentz_transform(time_sat_frame, sat2_pos_satframe[:, 1], plan_frame_speed)
        time_planframe = np.linspace(0, time_planframe1[-1] - time_planframe1[0], N)
        sat1_pos_1D_planframe = interpolate.interp1d(time_planframe1, sat1_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_planframe)
        sat2_pos_1D_planframe = interpolate.interp1d(time_planframe2, sat2_pos_1D_planframe, kind='linear', bounds_error=False, fill_value='extrapolate', assume_sorted=True)(time_planframe)
        sat1_pos_planframe = np.zeros((N, 3))
        sat2_pos_planframe = np.zeros((N, 3))
        laser_pos_planframe = np.zeros((N, 3))
        cam_pos_planframe = np.zeros((N, 3))
        sat1_pos_planframe[:, 1] = sat1_pos_1D_planframe
        sat2_pos_planframe[:, 1] = sat2_pos_1D_planframe
        sat1_pos_planframe[:, 2] -= spaceship_distance / 2
        sat2_pos_planframe[:, 2] += spaceship_distance / 2
        v_z = sqrt(c_AU ** 2 - plan_frame_speed ** 2)
        dz = (time_planframe[-1] - time_planframe[-2]) * v_z
        dy = (time_planframe[-1] - time_planframe[-2]) * plan_frame_speed
        laser_pos_planframe[0, :] = sat1_pos_planframe[0, :]
        laser_moving_upward = True
        for i in range(1, N):
            laser_pos_planframe[i, 1] = laser_pos_planframe[i - 1, 1] - dy
            if laser_moving_upward:
                laser_pos_planframe[i, 2] = laser_pos_planframe[i - 1, 2] + dz
            else:
                laser_pos_planframe[i, 2] = laser_pos_planframe[i - 1, 2] - dz
            if laser_pos_planframe[i, 2] > sat2_pos_planframe[i, 2]:
                if self.debug:
                    print('Laser bounced downwards at frame %d, [%es, %em] in satellite frame' % (i, time_sat_frame[i], laser_pos_planframe[i, 2]))
                laser_moving_upward = False
                laser_pos_planframe[i, 2] -= abs(laser_pos_planframe[i, 2] - sat2_pos_planframe[i, 2])
                sat2_messages_planframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat2_sounds_planframe[i] = 'laser'
            if laser_pos_planframe[i, 2] < sat1_pos_planframe[i, 2]:
                if self.debug:
                    print('Laser bounced upwards at frame %d, [%es, %em] in satellite frame' % (i, time_sat_frame[i], laser_pos_planframe[i, 2]))
                laser_moving_upward = True
                laser_pos_planframe[i, 2] += abs(laser_pos_planframe[i, 2] - sat1_pos_planframe[i, 2])
                sat1_messages_planframe[i:i + int(N / 64)] = [ 'Bounce!' for j in range(int(N / 64)) ]
                sat1_sounds_planframe[i] = 'laser'

        sat1_pos_planframe[:, 0] -= events_radius
        sat2_pos_planframe[:, 0] -= events_radius
        laser_pos_planframe[:, 0] -= events_radius + 10 * km_to_AU
        cam_pos_planframe[:, 0] -= planet_radius + height * km_to_AU
        sat1_color = [1, 0.5, 0]
        sat2_color = [0, 0.5, 1]
        obj_list_planframe = [['sat1',
          'Satellite',
          sat1_pos_planframe,
          0.1,
          sat1_color,
          sat1_messages_planframe,
          sat1_sounds_planframe,
          None,
          None], ['sat2',
          'Satellite',
          sat2_pos_planframe,
          0.1,
          sat2_color,
          sat2_messages_planframe,
          sat2_sounds_planframe,
          None,
          None], ['lsr',
          'Laser',
          laser_pos_planframe,
          8,
          [1, 0, 0],
          None,
          None,
          None,
          None]]
        camdir_planframe = np.copy(sat1_pos_planframe - cam_pos_planframe)
        camdir_planframe[:, 2] += spaceship_distance / 2
        upvec_planframe = np.array([0, 0, 1])
        self.write_to_xml(time_planframe, cam_pos_planframe, camdir_planframe, other_objects=obj_list_planframe, chosen_planet=chosen_planet, filename='A7_2.xml', laser_scale=0.1)

    def B1_1(self, chosen_planet, filename = 'B1_2.xml', N = 400):
        """ Part 2B tasak 1 scene 1. A spaceship passes by a planet and emits a
        light beam. Seen from the planet frame of reference"""
        random.seed(self.seed - 14)
        ref_frame_speed = random.uniform(0.97 * c, 0.995 * c)
        ref_frame_movement = np.linspace(-1000, 1000, N) * km_to_AU
        events_radius = self.radius[chosen_planet] * km_to_AU * 1.06
        dist_from_star = np.hypot(self.x0[chosen_planet], self.y0[chosen_planet])
        planet_pos = np.array([0, 0, 0])
        rocket_pos_1D = np.zeros(N)
        laser_pos_1D = np.zeros(N)
        end_time = AU_to_km * (ref_frame_movement[-1] - ref_frame_movement[0]) / ref_frame_speed
        time_array = np.linspace(0, end_time, N)
        dt = time_array[1] - time_array[0]
        dx = (time_array[-1] - time_array[-2]) * c * km_to_AU
        laser_pos_1D[0] = rocket_pos_1D[0]
        for i in range(N - 1):
            laser_pos_1D[i + 1] = laser_pos_1D[i] + dx

        cam_pos = np.zeros(shape=(3,))
        rocket_pos = np.zeros(shape=(N, 3))
        laser_pos = np.zeros(shape=(N, 3))
        rocket_pos[:, 2] = rocket_pos_1D + ref_frame_movement
        laser_pos[:, 2] = laser_pos_1D + rocket_pos[(0, 2)]
        cam_pos[2] = np.average(rocket_pos[:, 2])
        cam_pos[1] -= 700 * km_to_AU
        laser_pos[:, 2] += 25 * km_to_AU
        cam_pos[0] -= events_radius
        rocket_pos[:, 0] -= events_radius
        laser_pos[:, 0] -= events_radius
        cam_dir = [0, 1, 0]
        cam_upvec = [-1, 0, 0]
        other_objects = [['lsr1',
          'Laser',
          laser_pos,
          8,
          [1, 1, 0],
          None,
          None,
          None,
          None], ['rocket1',
          'Satellite',
          rocket_pos,
          0.1,
          [1, 1, 1],
          None,
          None,
          None,
          None]]
        self.write_to_xml(time_array, cam_pos, cam_dir, [1, 0, 0], other_objects, chosen_planet=chosen_planet, up_vec=cam_upvec, laser_scale=0.2, filename=filename)

    def B1_2(self, chosen_planet, filename = 'B1_2.xml', N = 400):
        """ Part 2B task 1 scene 2. A spaceship passes by a planet and emits a
        light beam. Seen from the spaceship frame of reference. """
        random.seed(self.seed - 14)
        ref_frame_speed = random.uniform(0.97 * c, 0.995 * c)
        ref_frame_movement = np.linspace(-1000, 1000, N) * km_to_AU
        events_radius = self.radius[chosen_planet] * km_to_AU * 1.06
        dist_from_star = np.hypot(self.x0[chosen_planet], self.y0[chosen_planet])
        planet_pos = np.array([0, 0, 0])
        rocket_pos_1D = np.zeros(N) - 200 * km_to_AU
        cam_pos_1D = np.zeros(N)
        laser_pos_1D = np.zeros(N)
        laser2_pos_1D = np.zeros(N)
        end_time = AU_to_km * (ref_frame_movement[-1] - ref_frame_movement[0]) / ref_frame_speed
        time_array = np.linspace(0, end_time, N)
        dt = time_array[1] - time_array[0]
        dx = (time_array[-1] - time_array[-2]) * c * km_to_AU
        laser_pos_1D[0] = rocket_pos_1D[0]
        laser2_pos_1D[0] = rocket_pos_1D[0]
        laser2_visibility = np.zeros(N)
        for i in range(N - 1):
            laser_pos_1D[i + 1] = laser_pos_1D[i] + dx
            if i > int(N / 2.0):
                laser2_pos_1D[i + 1] = laser_pos_1D[i] + dx
                laser2_visibility[i + 1] = 1

        cam_pos = np.zeros(shape=(N, 3))
        rocket_pos = np.zeros(shape=(N, 3))
        laser_pos = np.zeros(shape=(N, 3))
        laser2_pos = np.zeros(shape=(N, 3))
        cam_pos[:, 2] = cam_pos_1D + ref_frame_movement
        rocket_pos[:, 2] = rocket_pos_1D + ref_frame_movement
        laser_pos[:, 2] = laser_pos_1D + ref_frame_movement
        laser2_pos[:, 2] = laser2_pos_1D + ref_frame_movement
        laser_pos[:, 2] += 25 * km_to_AU
        laser2_pos[:, 2] += 25 * km_to_AU
        cam_pos[:, 0] -= events_radius
        rocket_pos[:, 0] -= events_radius
        laser_pos[:, 0] -= events_radius
        laser2_pos[:, 0] -= events_radius
        cam_pos[:, 1] -= 400 * km_to_AU
        cam_dir = [0, 1, 0]
        up_vec = [-1, 0, 0]
        other_objects = [['lsr1',
          'Laser',
          laser_pos,
          4,
          [1, 1, 0],
          None,
          None,
          None,
          None], ['rocket1',
          'Satellite',
          rocket_pos,
          0.1,
          [1, 1, 1],
          None,
          None,
          None,
          None]]
        self.write_to_xml(time_array, cam_pos, cam_dir, [1, 0, 0], other_objects, chosen_planet=chosen_planet, up_vec=up_vec, laser_scale=0.2, filename=filename)

    def lorentz_vel(self, v, w, c):
        return (v + w) / (1.0 + v * w / c ** 2)

    def B2(self, chosen_planet, chosen_seed, reference_system):
        system = AST2000SolarSystem(chosen_seed)
        radius = system.radius[chosen_planet]
        N = 300
        dt = 5e-05
        ships = np.zeros((4, N, 3))
        planet_pos0 = np.array((0, 1 * AU_to_km, 0))
        planet_pos = np.zeros((N, 3))
        vel = np.zeros((4, 3))
        vel[(0, 0)] = 0.9 * c
        vel[(1, 0)] = 0.7 * c
        vel[(2, 0)] = 0.6 * c
        vel[(3, 0)] = 0.95 * c
        cam_pos = np.zeros((N, 3))
        cam_dir = np.zeros((N, 3))
        for i in range(4):
            if i == 0:
                ships[0, 0, :] = np.array((2000, -1.5 * radius, 0))
            if i == 1:
                ships[1, 0, :] = np.array((8000, -1.5 * radius, 0))
            if i == 2:
                ships[2, 0, :] = np.array((-4000, -1.5 * radius, 0))
            else:
                ships[3, 0, :] = np.array((-8000, -1.5 * radius, 0))

        cam_dist = 12000
        ships_center = np.array((-1000, ships[(0, 0, 1)], ships[(0, 0, 2)]))
        planet_vel = np.array((system.vx0[chosen_planet], system.vy0[chosen_planet], 0)) * AU_pr_year_to_m_pr_s / 1000
        if reference_system == 0:
            cam_pos[0, :] = ships_center + np.array((0, -cam_dist, 0))
            cam_dir[0, :] = ships_center - cam_pos[0, :]
            v = vel[0, :]
            planet_rel_vel = -self.lorentz_vel(v, planet_vel, c)
            for i in range(N - 1):
                planet_pos[i + 1, :] = planet_pos[i, :] + dt * planet_rel_vel
                ships[0, i + 1, :] = ships[0, 0, :]
                cam_pos[i + 1, :] = ships_center + np.array((0, -cam_dist, 0))
                cam_dir[i + 1, :] = ships_center - cam_pos[i + 1, :]
                for ship in range(1, 4):
                    if ship == 1:
                        v = -v
                    w = vel[ship, :]
                    rel_vel = self.lorentz_vel(v, w, c)
                    ships[ship, i + 1, :] = ships[ship, i, :] + rel_vel * dt

        elif reference_system == 1:
            ship_count = [0, 2, 3]
            cam_pos[0, :] = ships_center + np.array((0, -cam_dist, 0))
            cam_dir[0, :] = ships_center - cam_pos[0, :]
            v = vel[1, :]
            planet_rel_vel = -self.lorentz_vel(v, planet_vel, c)
            for i in range(N - 1):
                planet_pos[i + 1, :] = planet_pos[i, :] + dt * planet_rel_vel
                ships[1, i + 1, :] = ships[1, 0, :]
                cam_pos[i + 1, :] = ships_center + np.array((0, -cam_dist, 0))
                cam_dir[i + 1, :] = ships_center - cam_pos[i + 1, :]
                for ship in ship_count:
                    w = vel[ship, :]
                    rel_vel = self.lorentz_vel(v, w, c)
                    ships[ship, i + 1, :] = ships[ship, i, :] + rel_vel * dt

        else:
            print('You can only choose the reference system as 0 (student nr. 1) or 1 (student nr. 2)')
            sys.exit()
        time_array = np.linspace(0, dt * N, N)
        object_name = ['student 1',
         'student 2',
         'ship 3',
         'ship 4']
        object_string = ['satellite',
         'satellite',
         'satellite',
         'satellite']
        planet_pos *= km_to_AU
        ships *= km_to_AU
        cam_pos *= km_to_AU
        cam_dir *= km_to_AU
        ships1 = np.zeros((4, 2, 3))
        cam_pos1 = np.zeros((2, 3))
        cam_dir1 = np.zeros((2, 3))
        planet_pos1 = np.zeros((2, 3))
        time_array = np.array((time_array[0], time_array[-1]))
        planet_pos1[0, :] = planet_pos[0, :]
        planet_pos1[1, :] = planet_pos[-1, :]
        ships1[:, 0, :] = ships[:, 0, :]
        ships1[:, 1, :] = ships[:, -1, :]
        cam_pos1[0, :] = cam_pos[0, :]
        cam_pos1[1, :] = cam_pos[-1, :]
        cam_dir1[0, :] = cam_dir[0, :]
        cam_dir1[1, :] = cam_dir[-1, :]
        other_objects = []
        message1 = []
        message2 = []
        message3 = []
        message4 = []
        for i in range(N):
            message1.append('sat 1')
            message2.append('sat 2')
            message3.append('sat 3')
            message4.append('sat 4')

        ruler_length = self.get_ruler_length(cam_dist * 2)
        ruler = [-ruler_length / 2.0,
         ruler_length / 2.0,
         10,
         'km']
        other_objects.append([object_name[0],
         object_string[0],
         ships1[0, :, :],
         1.5,
         [1, 0, 0],
         message1,
         None,
         None,
         None])
        other_objects.append([object_name[1],
         object_string[1],
         ships1[1, :, :],
         1.5,
         [1, 0, 0],
         message2,
         None,
         None,
         None])
        other_objects.append([object_name[2],
         object_string[2],
         ships1[2, :, :],
         1.5,
         [1, 0, 0],
         message3,
         None,
         None,
         None])
        other_objects.append([object_name[3],
         object_string[3],
         ships1[3, :, :],
         1.5,
         [1, 0, 0],
         message4,
         None,
         None,
         None])
        system.write_to_xml(time_array, cam_pos1, cam_dir1, planet_pos1, other_objects, chosen_planet=chosen_planet, ruler=ruler)

    def B4(self, chosen_planet, chosen_seed, reference_system):
        """
        @ chosen planet         background planet
        @ chosen seed           the chosen system in which to place the satellites
        @ reference system      can only chose 1 or 2 for the 2 different possibilities
        """
        system = AST2000SolarSystem(chosen_seed)
        radius = system.radius[chosen_planet]
        N = 300
        dt = 0.0005
        ships = np.zeros((2, N, 3))
        planet_pos0 = np.array((system.x0[chosen_planet], system.y0[chosen_planet], 0))
        planet_pos = np.zeros((N, 3))
        planet_pos[:, :] = planet_pos0
        vel = np.zeros((2, 3))
        vel[(0, 0)] = 0
        vel[(1, 0)] = -0.95 * c
        cam_pos = np.zeros((N, 3))
        cam_dir = np.zeros((N, 3))
        mass1 = 1500.0
        mass2 = 1500.0
        gamma2 = 1.0 / np.sqrt(1.0 - vel[(1, 0)] ** 2 / c ** 2)
        momentum1 = 0
        momentum2 = mass2 * vel[(1, 0)] * gamma2
        ships[0, :, :] = np.array((0, -1.5 * radius, 0))
        ships[1, 0, :] = np.array((8000, -1.5 * radius, 0))
        cam_dist = 12000
        ships_center = ships[0, 0, :]
        planet_vel = 0
        cam_pos[:, :] = ships_center + np.array((0, -cam_dist, 0))
        cam_dir[:, :] = ships_center - cam_pos[0, :]
        collision_index = 52
        for i in range(collision_index):
            ships[1, i + 1, :] = ships[1, i, :] + dt * vel[1, :]

        collision_vel = gamma2 * mass2 * vel[(1, 0)] / (mass1 + gamma2 * mass2)
        gamma3 = 1.0 / np.sqrt(1.0 - (collision_vel / c) ** 2)
        m3 = gamma2 * mass2 * vel[(1, 0)] / (gamma3 * collision_vel)
        for i in range(collision_index, N - 1):
            ships[1, i + 1, :] = ships[1, i, :] + dt * np.array((collision_vel, 0, 0))

        time_array = np.linspace(0, dt * N, N)
        object_name = ['ship 1', 'ship 2', 'invisible ship']
        object_string = ['satellite', 'satellite', 'satellite']
        for i in range(N):
            ships[:, i, :] -= planet_pos[i, :]
            cam_pos[i, :] -= planet_pos[i, :]

        if reference_system == 0:
            print('Reference system of the planet chosen')
        elif reference_system == 1:
            print('Reference system of the mass after collision chosen')
            cam_pos[0:collision_index, :] = ships[0, 0:collision_index, :] + np.array((0, -cam_dist, 0))
            cam_pos[collision_index:, :] = ships[1, collision_index:, :] + np.array((0, -cam_dist, 0))
        else:
            print('You can only choose the reference system as 0 (student nr. 1) or 1 (student nr. 2)')
            sys.exit()
        planet_pos *= km_to_AU
        ships *= km_to_AU
        cam_pos *= km_to_AU
        cam_dir *= km_to_AU
        time_array = np.array((time_array[0], time_array[collision_index], time_array[-1]))
        planet_pos[0, :] = planet_pos[0, :]
        planet_pos[1, :] = planet_pos[collision_index, :]
        planet_pos[2, :] = planet_pos[-1, :]
        planet_pos = planet_pos[0:3, :]
        ships[:, 0, :] = ships[:, 0, :]
        ships[:, 1, :] = ships[:, collision_index, :]
        ships[:, 2, :] = ships[:, -1, :]
        ships = ships[:, 0:3, :]
        cam_pos[0, :] = cam_pos[0, :]
        cam_pos[1, :] = cam_pos[collision_index, :]
        cam_pos[2, :] = cam_pos[-1, :]
        cam_pos = cam_pos[0:3, :]
        cam_dir[0, :] = cam_dir[0, :]
        cam_dir[1, :] = cam_dir[collision_index, :]
        cam_dir[2, :] = cam_dir[-1, :]
        cam_dir = cam_dir[0:3, :]
        invisible_ship_pos = np.zeros((3, 3))
        invisible_ship_pos[:, 0] = -600 * km_to_AU
        invisible_ship_pos[:, :] += ships[1, :, :]
        visibility1 = [1, 0, 0]
        visibility2 = [0, 1, 1]
        other_objects = []
        message1 = []
        message2 = []
        message3 = []
        for i in range(N):
            message1.append('sat 1')
            message2.append('sat 2')
            message3.append('invisible ship')

        ruler_length = self.get_ruler_length(cam_dist * 2)
        ruler = [-ruler_length / 2.0,
         ruler_length / 2.0,
         10,
         'km']
        other_objects.append([object_name[0],
         object_string[0],
         ships[0, :, :],
         1.5,
         [1, 0, 0],
         None,
         None,
         visibility1,
         None])
        other_objects.append([object_name[1],
         object_string[1],
         ships[1, :, :],
         1.5,
         [1, 0, 0],
         None,
         None,
         None,
         None])
        other_objects.append([object_name[2],
         object_string[2],
         invisible_ship_pos,
         1.5,
         [1, 0, 0],
         None,
         None,
         visibility2,
         None])
        system.write_to_xml(time_array, cam_pos, cam_dir, planet_pos, other_objects, chosen_planet=chosen_planet, cheat_light_speed=True, ruler=ruler)

    def A9(self, pl_nr, height, c = 1):
        """
        @ pl_nr         planet number for observer 3 to stand on
        @ height        add to radial distance from planetsurface to avoid mountains
        @ c             set speed of light for twin paradox to have markable effect. Give in fraction of real speed of light in AU/yr
        OBS these scenes use the XML writer without scaling any objects or planets to more freely set speeds
        """
        N = 500
        pos3 = np.zeros((N, 3))
        pos2 = np.zeros((N, 3))
        pos1 = np.zeros((N, 3))
        times3 = np.linspace(0, 1000, N)
        times2 = np.zeros(N)
        times1 = np.zeros(N)
        other_objects = []
        sc1_pos = np.zeros(shape=(N, 3))
        sc2_pos = np.zeros(shape=(N, 3))
        cam_dir2 = np.zeros(shape=(N, 3))
        expl_pos = np.zeros(shape=(N, 3))
        scale1 = 0.8
        obj_list = []
        P_orbit = self.a[pl_nr]
        P_pos = np.array([0, 0, 0])
        P_radius = self.radius[pl_nr] * km_to_AU
        print('orbit', P_orbit)
        sc1_pos += self.spherical_to_cart(1.05 * P_radius + height * km_to_AU, pi / 2, 1.5)
        sc1_pos[:, 2] += np.linspace(0, 0.5, N)
        sc1_pos[:, 0] -= 1000 * km_to_AU
        colour1 = np.array([0, 0.8, 0.6])
        message_list = [ '' for i in range(N) ]
        sound_list = None
        visible = None
        sc1_list = ['Kerbal',
         'Satellite',
         sc1_pos,
         scale1,
         colour1,
         message_list,
         sound_list,
         visible,
         None]
        obj_list.append(sc1_list)
        sc2_pos = np.flipud(np.copy(sc1_pos))
        sc2_pos[:, 0] += 1000 * km_to_AU
        colour2 = np.array([1, 0, 0])
        sc2_list = ['Herbal',
         'Satellite',
         sc2_pos,
         scale1,
         colour2,
         message_list,
         sound_list,
         visible,
         None]
        obj_list.append(sc2_list)
        pos1 += np.copy(sc2_pos)
        pos1[:, 2] += 2000 * km_to_AU
        pos1[:, 1] += 1000 * km_to_AU
        cam_dir1 = np.copy(sc2_pos) - np.copy(pos1)
        pos3 += self.spherical_to_cart(P_radius + height * km_to_AU, pi / 2, pi / 2)
        cam_dir3 = np.copy(sc2_pos - pos3)
        cam_dir3 = cam_dir3 / np.linalg.norm(cam_dir3, axis=1)[:, None]
        up_vec3 = np.copy(pos3[0, :])
        expl_pos[:, :] = pos3[1, :]
        expl_pos[:, 2] += 0.0002
        colour3 = np.array([1, 1, 1])
        visible_expl = np.zeros(N)
        visible_expl[250:265] = 1
        expl_list = ['Flash',
         'explosion',
         expl_pos,
         50000,
         colour3,
         message_list,
         sound_list,
         visible_expl,
         None]
        obj_list.append(expl_list)
        print('expl_pos', expl_pos[0, :])
        V = self.focus_tanhyper(50, start=200, time=times3)
        for i in range(3):
            pos2[:, i] = (1 - V) * sc1_pos[:, i] + V * sc2_pos[:, i]
            if i == 2:
                pos2[:, i] = (1 - V) * (sc1_pos[:, i] - 4000 * km_to_AU) + V * (sc2_pos[:, i] + 2000 * km_to_AU)

        pos2[:, 1] += 2000 * km_to_AU
        cam_dir2 = np.copy(sc2_pos - pos2)
        cam_dir2 = cam_dir2 / np.linalg.norm(cam_dir2, axis=1)[:, None]
        up_vec2 = np.copy(pos3[0, :])
        cam2_messages = [ '' for i in range(N) ]
        for i in range(3):
            for j in range(140 + 20 * i, 160 + 20 * i):
                cam2_messages[j] = 'Launch in %i...' % (3 - i)

        for i in range(10):
            cam2_messages[200 + i] = 'Launch!'
            cam2_messages[245 + i] = 'Securly attached to spacecraft 2'

        dt_3 = times3[1] - times3[0]
        lightspeed = c * c_AU
        self.cheat_light_speed = False
        gamma = self.get_lorentz(times3, pos3, pos1, c=lightspeed)
        dt_marked = dt_3 / gamma[0]
        for i in range(N - 1):
            times2[i + 1] = times2[i] + dt_marked

        times1 = np.copy(times2)
        print('gamma', gamma[0])
        print('dt3', dt_3)
        print('dtmarked', dt_marked)
        print('endtime 3', times3[-1])
        print('endtime 2', times2[-1])
        self.write_to_xml(times3, pos3, cam_dir3, P_pos, obj_list, up_vec=up_vec3, filename='A9_3.xml', chosen_planet=pl_nr, c=lightspeed, use_obj_scaling=0)
        self.write_to_xml(times2, pos2, cam_dir2, P_pos, obj_list, up_vec=up_vec2, filename='A9_2.xml', chosen_planet=pl_nr, camera_messages=cam2_messages, c=lightspeed, use_obj_scaling=0)
        self.write_to_xml(times1, pos1, cam_dir1, P_pos, obj_list, up_vec=up_vec3, filename='A9_1.xml', chosen_planet=pl_nr, c=lightspeed, use_obj_scaling=0)

    def B3(self, chosen_planet, N = 400):
        """ Two spaceships crash, all the particles annihilate and the energy
        is converted into only photons, all with a frequency nu. Frames of reference:
        1: Planet frame of reference
        2: One of the ships' frame of reference
        3: A frame of reference moving with some speed towards the explosion
        ARGS:
        chosen_planet = Index of planet the events will happen on
        path          = Path to where you want the xmls to end up
        N             = Number of frames/time steps in video
        """
        filename_s1 = 'B3_1'
        filename_s2 = 'B3_2'
        filename_s3 = 'B3_3'
        random.seed(self.seed - 3)
        ship_speeds = random.uniform(0.97 * c_AU, 0.995 * c_AU)
        events_radius = self.radius[chosen_planet] * km_to_AU * 1.03
        dist_from_star = np.hypot(self.x0[chosen_planet], self.y0[chosen_planet])
        planet_pos = np.array([0, 0, 0])
        ship_travel_dist = 1000
        rocket_pos_1D = np.linspace(-ship_travel_dist, ship_travel_dist, N) * km_to_AU
        rocket2_pos_1D = np.linspace(ship_travel_dist, -ship_travel_dist, N) * km_to_AU
        end_time = np.abs(rocket_pos_1D[-1] - rocket_pos_1D[0]) / ship_speeds
        time_array = np.linspace(0, end_time, N)
        dt = time_array[1] - time_array[0]
        i_crash = np.argmin(np.abs(rocket_pos_1D - rocket2_pos_1D))
        ships_visible = np.ones(N)
        ships_visible[i_crash:] = 0
        explosions_visible = 1 - ships_visible
        cam_pos = np.zeros(shape=(3,))
        rocket_pos = np.zeros(shape=(N, 3))
        rocket2_pos = np.zeros(shape=(N, 3))
        expl_pos = np.zeros(shape=(N, 3))
        rocket_pos[:, 2] = rocket_pos_1D
        rocket2_pos[:, 2] = rocket2_pos_1D
        cam_pos[2] = 0
        cam_pos[1] -= 750 * km_to_AU
        cam_pos[0] -= events_radius
        rocket_pos[:, 0] -= events_radius
        rocket2_pos[:, 0] -= events_radius
        expl_pos[:, 0] -= events_radius
        cam_dir = [0, 1, 0]
        cam_upvec = [-1, 0, 0]
        colors = ['green', 'yellow']
        colors_wl = [510, 580]
        colors_rgbs = [(0, 1, 0), (1, 1, 0)]
        explosion_color_i = random.randint(0, 1)
        explosion_color_wl = colors_wl[explosion_color_i]
        explosion_rgb = colors_rgbs[explosion_color_i]
        m_rest = 10000
        vel = c_m * ship_speeds / c_AU
        gamma = 1 / np.sqrt(1 - vel ** 2 / c_m ** 2)
        if self.debug:
            print('Ship velocities = %g' % (vel / c_m))
        if self.debug:
            print('Gamma = %g' % gamma)
        p_rel = m_rest * vel * gamma
        spaceship_total_energy = np.sqrt(p_rel ** 2 * c_m ** 2 + (m_rest * c_m ** 2) ** 2)
        E_sys = 2 * spaceship_total_energy
        h = 6.626069934e-34
        E_photon = h * c_m / explosion_color_wl
        n_photons = E_sys / E_photon
        print('The spaceships collide and all their atoms annihilate! Information you will need to solve the tasks:')
        print('Spaceship speed        = %g' % (vel / c_m))
        print('Spaceship rest mass    = %g kg' % m_rest)
        print('Number of photons      = %g' % n_photons)
        print('(The mass and speed is the same for both spaceships.)')
        other_objects = [self.object_list('rocket2', 'Satellite', rocket2_pos, 0.15, [1, 1, 0], visible=ships_visible), self.object_list('rocket1', 'Satellite', rocket_pos, 0.15, [1, 1, 1], visible=ships_visible), self.object_list('Expl', 'explosion', expl_pos, 1000, explosion_rgb, visible=1 - ships_visible)]
        self.write_to_xml(time_array, cam_pos, cam_dir, planet_pos, other_objects, chosen_planet=chosen_planet, up_vec=cam_upvec, laser_scale=0.2, filename=filename_s1, ruler=[0,
         2 * ship_travel_dist,
         12,
         'km'])
        rocket1_s2_1D_func = self.ref_sys_interpolate_func(time_array, rocket_pos_1D, ship_speeds)
        rocket2_s2_1D_func = self.ref_sys_interpolate_func(time_array, rocket2_pos_1D, ship_speeds)
        time_array1_s2 = self.lorentz_transform(time_array, rocket_pos_1D, ship_speeds)[0]
        time_array2_s2 = self.lorentz_transform(time_array, rocket2_pos_1D, ship_speeds)[0]
        ts = time_array1_s2
        rocket1_s2_1D = rocket1_s2_1D_func(ts)
        rocket2_s2_1D = rocket2_s2_1D_func(ts)
        time_in_s2 = ts[-1] - ts[0]
        dist_ref_frame2 = ship_speeds * time_in_s2
        start = -dist_ref_frame2 / 2
        stop = -start
        ref_sys_movement = np.linspace(start, stop, N)
        rocket1_s2 = np.zeros(shape=(N, 3))
        rocket2_s2 = np.zeros(shape=(N, 3))
        rocket1_s2[:, 2] = rocket1_s2_1D - rocket1_s2_1D[0]
        rocket2_s2[:, 2] = rocket2_s2_1D - rocket1_s2_1D[0]
        rocket2_rel2_rocket1 = np.copy(rocket2_s2)
        cam_pos_s2 = np.copy(rocket1_s2)
        rocket1_s2[:, 2] += ref_sys_movement
        rocket2_s2[:, 2] += ref_sys_movement
        cam_pos_s2[:, 2] += ref_sys_movement
        cam_pos_s2[:, 0] -= events_radius
        rocket1_s2[:, 0] -= events_radius
        rocket2_s2[:, 0] -= events_radius
        rocket2_s2[:, 2] = rocket1_s2[:, 2] + rocket2_rel2_rocket1[:, 2]
        cam_pos_s2[:, 0] -= 100 * km_to_AU
        cam_pos_s2[:, 2] -= 25 * km_to_AU
        cam_pos_s2[:, 1] += 120 * km_to_AU
        i_crash_s2 = np.argmin(np.abs(rocket1_s2_1D - rocket2_s2_1D))
        crash_site_s2 = rocket1_s2[i_crash_s2]
        cam_dir_s2 = crash_site_s2 - cam_pos_s2
        cam_upvec_s2 = cam_pos_s2 / np.linalg.norm(cam_pos_s2)
        expl_rgb_s2 = (random.uniform(0.9, 1), 0, 0)
        other_objects_s2 = [['rocket2',
          'Satellite',
          rocket2_s2,
          0.05,
          [1, 1, 1],
          None,
          None,
          ships_visible,
          None], ['rocket1',
          'Satellite',
          rocket1_s2,
          0.05,
          [1, 1, 1],
          None,
          None,
          ships_visible,
          None]]
        self.write_to_xml(ts - ts[0], cam_pos_s2, cam_dir_s2, planet_pos, other_objects_s2, chosen_planet=chosen_planet, up_vec=cam_upvec_s2, filename=filename_s2)
        wl_blue = 400 + random.uniform(-20, 20)
        original_wl = 580
        delta_wl = original_wl - wl_blue
        rgb_s3 = [0, 0, 1]
        speed_s3 = self.relativistic_doppler_find_vel(original_wl - wl_blue, original_wl)
        if self.debug:
            print('Camera speed in scene 3:', speed_s3)
        print('In scene 3, the camera is moving towards the explosion with velocity %g' % speed_s3)
        speed_s3 *= c_AU
        cam_movement_dist = 10000 * km_to_AU
        explosion_loc_1D = np.zeros(N)
        cam_movement_1D = np.linspace(-2 * cam_movement_dist, -cam_movement_dist, N)
        endtime = cam_movement_dist / speed_s3
        times = np.linspace(0, endtime, N)
        cam_pos_s3 = np.zeros(shape=(N, 3))
        cam_pos_s3[:, 0] -= events_radius
        cam_pos_s3[:, 0] += cam_movement_1D
        cam_dir_s3 = cam_pos_s3[-1, :] - cam_pos_s3[0, :]
        cam_upvec_s3 = [0, 0, 1]
        explosion_loc = np.zeros(shape=(N, 3))
        other_objects_s3 = []
        other_objects_s3.append(self.object_list('expl', 'explosion', explosion_loc, 10000, visible=1 - ships_visible))
        self.write_to_xml(times, cam_pos_s3, cam_dir_s3, planet_pos, other_objects_s3, chosen_planet=chosen_planet, up_vec=cam_upvec_s3, filename=filename_s3)

    def C4(self, chosen_planet = 0, N = 400, filename = 'C4.xml'):
        """ Global Positioning System. The video is the same no matter what, but the
        positions shown are as if it is somewhere else on the planet.
        Lag posisjoner og tidspunkter for satelittene. Lag nye arrays med posisjoner til
        satelittene sett fra personen p\xc3\xa5 planeten (tidspunkter f\xc3\xb8rst for n\xc3\xa5r korresponderende
        informasjon treffer, og posisjoner til disse tidene. Interpoler s\xc3\xa5 posisjonene). """
        G = 6.67408e-11
        c_m = c * 1000.0
        M_kg_to_m = G / c_m ** 2
        M_m_to_kg = 1 / M_kg_to_m

        def time_diff_shells(M, r1, r2, v2):
            """ Returns t1/t2, the factor at which time runs differently, at the
            shell a distance r1, relative
            to somebody on the shell at a distance r2 moving with tangential velocity v2."""
            M_m = M * M_kg_to_m
            teller = 1 - 2 * M_m / r1
            nevner = 1 - 2 * M_m / r2 - v2 ** 2 / c_m ** 2
            print(2 * M_m / r2, v2 ** 2 / c_m ** 2)
            return np.sqrt(teller / nevner)

        def sat_speed(M, r):
            """ Speed of orbiting satellite at a radius r in kilometers, orbiting a planet
            with mass M in kilograms"""
            print(np.sqrt(G * M / r))
            return np.sqrt(G * M / r)

        random.seed(int(self.seed / 2 + 716))
        M = self.mass[chosen_planet] * sun_mass_to_kg
        R = self.radius[chosen_planet]
        sat_orbit_radius = 2 * R
        v_sat = sat_speed(M, sat_orbit_radius * 1000.0) / 1000.0
        omega_sat = v_sat / sat_orbit_radius
        T_sat = 2 * np.pi / omega_sat
        T_sat_hr = T_sat / 3600
        theta = random.uniform(0, 2 * np.pi)
        person_pos0 = R * np.array([cos(theta), sin(theta), 0])
        sat_sys_angular_dist = pi
        ang_dist_between_sats = pi / 7
        tot_time = sat_sys_angular_dist / omega_sat
        times_earth = np.linspace(0, tot_time, N)
        times_sat = np.copy(times_earth) / time_diff_shells(M, R, sat_orbit_radius, v_sat)
        cam_pos = np.zeros(shape=(N, 3))
        sat1_pos = np.zeros(shape=(N, 3))
        sat2_pos = np.zeros(shape=(N, 3))
        sat_sys_pos = np.zeros(shape=(N, 3))
        cam_pos[:, 0] -= 1.01 * R * km_to_AU
        sat_sys_pos[:, 0] = sat_orbit_radius * np.cos(pi + sat_sys_angular_dist / 2 - omega_sat * times_earth) * km_to_AU
        sat_sys_pos[:, 1] = sat_orbit_radius * np.sin(pi + sat_sys_angular_dist / 2 - omega_sat * times_earth) * km_to_AU
        sat1_pos[:, 0] = sat_orbit_radius * np.cos(pi + sat_sys_angular_dist / 2 - omega_sat * times_earth - ang_dist_between_sats) * km_to_AU
        sat1_pos[:, 1] = sat_orbit_radius * np.sin(pi + sat_sys_angular_dist / 2 - omega_sat * times_earth - ang_dist_between_sats) * km_to_AU
        sat2_pos[:, 0] = sat_orbit_radius * np.cos(pi + sat_sys_angular_dist / 2 - omega_sat * times_earth + ang_dist_between_sats) * km_to_AU
        sat2_pos[:, 1] = sat_orbit_radius * np.sin(pi + sat_sys_angular_dist / 2 - omega_sat * times_earth + ang_dist_between_sats) * km_to_AU
        cam_dir = np.copy(sat_sys_pos) - cam_pos[:]
        cam_up = np.zeros(shape=(N, 3))
        sat_normalized = sat_sys_pos / np.linalg.norm(sat_sys_pos, axis=1)[:, None]
        cam_up[:int(N / 2), 0] = np.linspace(sat_normalized[(0, 1)], 0, int(N / 2))
        cam_up[int(N / 2):, 0] = np.linspace(0, sat_normalized[(0, 1)], int(N / 2))
        cam_up[:int(N / 2), 1] = np.linspace(-sat_normalized[(0, 0)], 0, int(N / 2))
        cam_up[int(N / 2):, 1] = np.linspace(0, -sat_normalized[(0, 0)], int(N / 2))
        z_upvec = np.zeros(N)
        z_upvec[:int(N / 2)] = self.focus_tanhyper(int(N / 2))
        z_upvec[int(N / 2):] = np.flipud(z_upvec[:int(N / 2)])
        cam_up[:, 2] = z_upvec[:]
        person_pos = np.zeros(shape=(N, 3))
        person_pos[:, :] = person_pos0
        t_received1_arr = np.zeros(N)
        t_received2_arr = np.zeros(N)
        for i in xrange(N):
            d1 = np.linalg.norm(person_pos[i] - sat1_pos[i])
            d2 = np.linalg.norm(person_pos[i] - sat2_pos[i])
            t1 = d1 / c_AU
            t2 = d2 / c_AU
            t_received1_arr[i] = times_sat[i] + t1
            t_received2_arr[i] = times_sat[i] + t2

        pos_sat1_received = interpolate.interp1d(t_received1_arr, sat1_pos, axis=0, fill_value='extrapolate')
        pos_sat2_received = interpolate.interp1d(t_received2_arr, sat2_pos, axis=0, fill_value='extrapolate')
        pos_sat1_seen_from_plan = pos_sat1_received(times_earth)
        pos_sat2_seen_from_plan = pos_sat2_received(times_earth)
        sat1_messages = [ str(pos_sat1_seen_from_plan[i]) for i in range(N) ]
        sat2_messages = [ str(pos_sat2_seen_from_plan[i]) for i in range(N) ]
        other_objects = [self.object_list('sat1', 'Satellite', sat1_pos, 1.5, [1, 1, 1]), self.object_list('sat2', 'Satellite', sat2_pos, 1.5, [1, 1, 1])]
        self.write_to_xml(times_earth, cam_pos, cam_dir, [0, 0, 0], other_objects, up_vec=cam_up, chosen_planet=chosen_planet, filename=filename)

    def schwarzschild_line_element(self, dt, dr, r, M, dtheta = 0):
        """
        dt  =  Time interval on far away observer clock
        dr  =  Change in distance from black hole for shell observer
        r   =  Radius of shell observer (assumes small dr, and relatively constant r)
        M   =  Mass of black hole
        dtheta  =  Angular movement around black hole
        RETURNS  =  Line element of shell observer
        """
        return (1 - 2 * M / r) * dt ** 2 - dr / (1 - 2 * M / r) - r ** 2 * dtheta

    def C2_1(self, chosen_planet):
        M = sun_mass_to_kg * kg_to_m
        M *= 40000000.0
        r_AU = 1
        planet_pos = [r_AU, 0, 0]
        r_m = r_AU * AU_to_m
        sat_start_speed = 0.1
        gamma = 1 / sqrt(1 - sat_start_speed ** 2)
        print('r_AU = ', r_AU)
        print('r_m = ', r_m)
        print('M_SM = ', M * m_to_AU)
        print('M_m = ', M)
        obs_start_radius = r_m - self.radius[chosen_planet] * 1000 * 1.04
        sat_start_radius = r_m - self.radius[chosen_planet] * 1000 * 3.0
        dt = 1000000000.0
        time_array_sat = [0]
        time_array_obs = [0]
        sat_pos_1D = [sat_start_radius]
        while sat_pos_1D[-1] > 2 * M * 1.01:
            time_array_obs.append(time_array_obs[-1] + dt)
            dtau = (1 - 2 * M / sat_pos_1D[-1]) / gamma * dt
            time_array_sat.append(time_array_sat[-1] + dtau)
            sat_pos_1D.append(sat_pos_1D[-1] - sqrt(gamma ** 2 - (1 - 2 * M / sat_pos_1D[-1])) * dtau)

        sat_pos_1D = np.array(sat_pos_1D[:-1])
        time_array_sat = np.array(time_array_sat[:-1])
        time_array_obs = np.array(time_array_obs[:-1])
        if self.debug:
            print('Task C2\nSatellite start pos: %e\nend pos: %e\nschwarzschild radius: %e' % (sat_pos_1D[0] * m_to_AU, sat_pos_1D[-1] * m_to_AU, 2 * M * m_to_AU))
        N3 = len(sat_pos_1D)
        N2 = N3 // 4
        N1 = N3 // 2
        N = N1 + N2 + N3
        if self.debug:
            print('Frames = %d' % N)
        obs_pos = np.zeros(shape=(N, 3))
        obs_pos[:, 0] += obs_start_radius
        light_interval = (time_array_sat[-1] - time_array_sat[0]) / 18
        send_light_times_sat = np.array([ i * light_interval for i in range(2, 19) ])
        light_indexes = []
        for i in range(len(send_light_times_sat)):
            light_indexes.append(np.abs(send_light_times_sat[i] - time_array_sat).argmin())

        print('light indexes ', np.array(light_indexes[1:]) - np.array(light_indexes[:-1]))
        light_obj_pos = np.zeros(shape=(N, 3)) + np.array([r_AU - 0.001, 0, 0])
        light_obj_visibility = np.zeros(N)
        for i in light_indexes:
            light_obj_visibility[N1 + N2 + i:N1 + N2 + i + N // 128] += 1
            print(N1 + N2 + i, N1 + N2 + i + N // 128)

        time_array_XML = np.zeros(N)
        time_array_XML[:N1] = np.linspace(0, dt * N1, N1)
        time_array_XML[N1 - 1:N1 + N2 + 1] = np.linspace(time_array_XML[N1 - 1], time_array_XML[N1 - 1] * 100, N2 + 2)
        time_array_XML[N1 + N2:] = time_array_obs + time_array_XML[N1 + N2]
        for i in range(1, N):
            print('%e' % (time_array_XML[i] / c_m))

        sat_pos_fake = np.zeros(shape=(N, 3))
        sat_pos_fake[:, 0] = np.logspace(log(sat_start_radius * m_to_AU), log(sat_start_radius * m_to_AU - 1e-05), N)
        print(sat_pos_fake)
        other_objects = [['rocket1',
          'Satellite',
          sat_pos_fake,
          3.0,
          [1, 1, 1],
          None,
          None,
          None,
          None], ['light',
          'Sphere01',
          light_obj_pos,
          100000,
          [1, 1, 1],
          None,
          None,
          light_obj_visibility,
          None]]
        cam_dir_s1 = [-1, 0, 0]
        self.write_to_xml(time_array_XML / c_m, obs_pos * m_to_AU, cam_dir_s1, planet_pos=planet_pos, other_objects=other_objects, chosen_planet=chosen_planet, filename='C2_1.xml', origo_location=np.array([0, 0, 0]))


if __name__ == '__main__':
    from sys import argv
    if len(argv) > 1:
        seed = int(argv[1])
    else:
        print('This module is intended to be imported into your python code.')
        print('Please provide a seed if you want to print info about a solar system.')
        print('Now I am just generating a random seed and printing it for fun...')
        seed = randint(0, 99999)
        print('Randomly generated seed', seed, '. THIS IS NOT YOUR SEED, testing only.')
        input('Press enter to agree...')
    system = AST2000SolarSystem(seed)
    system.print_info()
