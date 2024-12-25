#coding=utf-8
from numpy import sqrt
from ast2000solarsystem_27_v2 import AST2000SolarSystem as AST2000SolarSystem1
from AST2000SolarSystem import AST2000SolarSystem
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
MittSolSystem.send_satellite("sat_commands.txt")