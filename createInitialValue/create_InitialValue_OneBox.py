#!/usr/bin/env python
# -*- coding: utf8 -*-

from Random_Tracer import Random_Tracer as rt
import numpy as np

#Berechnung mehrer zufälliger initialer Tracerkonzentrationen
start=0
last=100
tracer_count = 5
distribution = 'OneBox'
path_metos3d = '/sfs/fs5/home-sh/sunip350/.metos3d'
path_tracer = '/sfs/fs2/work-sh1/sunip350/metos3d/InitialTracer/OneBoxDistribution/Five_Tracer'

tracer = rt(tracer_count, path_metos3d, distribution=distribution)
tracername = ['InitialValue_Tracer_1_', 'InitialValue_Tracer_2_', 'InitialValue_Tracer_3_', 'InitialValue_Tracer_4_', 'InitialValue_Tracer_5_']
end = '.petsc'

#Verteilung der Masse auf die einzelnen Tracer
mass_distribution = np.zeros(tracer_count, dtype='>f8')
mass_distribution[0]=1.0
#mass_distribution[1]=0.0
#mass_distribution[2]=0.0
#mass_distribution[3]=0.0
#mass_distribution[4]=0.0

#Erwartungswerte und Standardabweichungen für die einzelnen Tracer
mean = np.ones(tracer_count, dtype='>f8')
mean[0] = 3.5
#mean[1] = 2.0
sigma = np.ones(tracer_count, dtype='>f8')
sigma[0] = 1.0
#sigma[1] = 0.5

for i in range(start, last):
    if i%10 == 0:
        tracer.random_mass_distribution()

    #tracer.set_mass_distribution(mass_distribution)
    #tracer.calculate_mean_sigma()
    #tracer.set_mean(mean)
    #tracer.set_sigma(sigma)
    tracer.calculate_mu_sigma()
    tracer.create_random_tracer()
    
    #Masse der zufällig erstellten Tracer anpassen
    for j in range(tracer_count):
        tracer.masscorrection(j)
    
    filename = [tracername[0] + str(i).zfill(3) + end, tracername[1] + str(i).zfill(3) + end, tracername[2] + str(i).zfill(3) + end, tracername[3] + str(i).zfill(3) + end, tracername[4] + str(i).zfill(3) + end]

    tracer.write_all_tracer_to_file(path_tracer, filename)
    print('Tracer fuer i={} erstellt\n'.format(i))
