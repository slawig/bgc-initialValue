#!/usr/bin/env python
# -*- coding: utf8 -*-

from Random_Tracer import Random_Tracer as rt
import numpy as np
import argparse
import os


#Command line parameter
parser = argparse.ArgumentParser()
parser.add_argument("start", type=int, help="Startindex")
parser.add_argument("last", type=int, help="Lastindex")
parser.add_argument("tracer_count", type=int, help="Tracer count")
parser.add_argument("distribution", type=str, help="Distribution")
parser.add_argument("-set_mass_distribution", "--set_mass_distribution", action="store_true")
args = parser.parse_args()

start = args.start
last = args.last
tracer_count = args.tracer_count
distribution = args.distribution
assert distribution == 'lognormal' or distribution == 'normal' or distribution == 'uniform' or distribution == 'OneBox'

#Path to save initial tracer
path_metos3d = '/sfs/fs5/home-sh/sunip350/.metos3d'
path = '/sfs/fs2/work-sh1/sunip350/metos3d/InitialTracer'
distribution_path = {'lognormal' : 'LognormalDistribution', 'normal' : 'NormalDistribution', 'uniform' : 'UniformDistribution', 'OneBox' : 'OneBoxDistribution' }
set_mass_distribution_path = {True : 'set_mass_distribution', False : 'random_mass_distribution'}
tracer_count_path = {1 : 'One_Tracer', 2 : 'Two_Tracer', 3 : 'Three_Tracer', 4 : 'Four_Tracer', 5 : 'Five_Tracer'}

path_tracer = os.path.join(path, distribution_path[distribution], set_mass_distribution_path[args.set_mass_distribution], tracer_count_path[tracer_count])

mass = 2.17 + 10**(-4) * (tracer_count - 1)

tracer = rt(tracer_count, path_metos3d, distribution=distribution, mass=mass)
tracername = ['InitialValue_Tracer_0_', 'InitialValue_Tracer_1_', 'InitialValue_Tracer_2_', 'InitialValue_Tracer_3_', 'InitialValue_Tracer_4_']
end = '.petsc'

#Distribution of the single tracer
mass_distribution = np.zeros(tracer_count, dtype='>f8')
mass_distribution[0]=2.17 / mass
for i in range(1, tracer_count):
    mass_distribution[i]=10**(-4) / mass

for i in range(start, last):
    if args.set_mass_distribution:
        tracer.set_mass_distribution(mass_distribution)
        tracer.calculate_mean_sigma()
    else:
        if i%10 == 0:
            tracer.random_mass_distribution()

    tracer.calculate_mu_sigma()
    tracer.create_random_tracer()
    
    #Adapt mass of the randomly created tracer
    for j in range(tracer_count):
        tracer.masscorrection(j, iteration=1e4)
    
    filename = [tracername[0] + str(i).zfill(3) + end, tracername[1] + str(i).zfill(3) + end, tracername[2] + str(i).zfill(3) + end, tracername[3] + str(i).zfill(3) + end, tracername[4] + str(i).zfill(3) + end]

    tracer.write_all_tracer_to_file(path_tracer, filename)
    print('Created tracer i={}\n'.format(i))

