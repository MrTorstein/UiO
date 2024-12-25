#coding=utf-8
from numpy import sqrt, sin, pi, cos
from ast2000solarsystem_27_v2 import AST2000SolarSystem as AST2000SolarSystem
from Part1 import part1
from Part4 import lsf, comp_v, PosSat
seed = 12827
MittSolSystem = AST2000SolarSystem(seed)

tilfelle = part1()

R = MittSolSystem.radius[0]*1e3
M = MittSolSystem.mass[0]*1.989*10**30
G = 6.67428*10**(-11)
v_esc = sqrt((2*G*M)/R)
anpart = 7.7e13


v_end, n_box, pos, sec, fuel = tilfelle.launch_sim(no_boxes = 3e+13, v_end = v_esc, A = 4.4e-9, no_particles_sec = anpart)	


init_sat_pos = (MittSolSystem.x0[0] + MittSolSystem.radius[0]/149597870.700, MittSolSystem.y0[0])
MittSolSystem.engine_settings(4.4e-9, n_box, anpart, 10000, sec, init_sat_pos, 0)

calculated_y =  1.61407848e-04
calculated_x = 3.50934789e+00
pos_after_launch = calculated_x, calculated_y
a = MittSolSystem.mass_needed_launch(pos_after_launch, test = True)

MittSolSystem.send_satellite("sat_commands.txt")
















tilfelle = part1()
p, anpart, pos = tilfelle.simulate_box(n = 10**5, N = 1000, dt = 10**(-9)/1000., box_dim = 10**(-6), T = 10**4)
print "p = ", p		#3.03841102114e-20
print "hullteller =", anpart	#3290

R = MittSolSystem.radius[0]*1e3
M = MittSolSystem.mass[0]*1.989*10**30
G = 6.67428*10**(-11)
v_esc = sqrt((2*G*M)/R)
anpart = 7.7e13



v_end, n_box, pos, sec, fuel = tilfelle.launch_sim(no_boxes = 3e+13, v_end = v_esc, A = 4.4e-9, no_particles_sec = anpart)	

print 10577060.0089*6.68458712*10**-12

print tilfelle.boost_sim(no_boxes = n_box, v_start = 0, v_end = 1000, A = 4.4e-9, no_particles_sec = anpart, fuel_mass = 500)

"""
drivs_trengs = 91.06026455567189
tid_brukt = 11.399999999999801
"""

init_sat_pos = (MittSolSystem.x0[0] + MittSolSystem.radius[0]/149597870.700, MittSolSystem.y0[0])
MittSolSystem.engine_settings(4.4e-9, n_box, anpart, 10000, sec, init_sat_pos, 0)

