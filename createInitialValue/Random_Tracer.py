#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np

class Random_Tracer():

    def __init__(self, tracer_count, path_metos3d, distribution='uniform', length=52749, mass=2.17):
        self.__tracer_count = tracer_count
        self.__path_metos3d = path_metos3d
        self.__length = length
        self.__tracer = np.zeros(shape=(self.__tracer_count, self.__length), dtype='>f8')
        self.__distribution = distribution

        self.__read_boxvolumen_from_file()

        self.random_mass_distribution()        
        self.__calculate_total_mass(mass)
        
        self.calculate_mean_sigma()


    def write_tracer_to_file(self, number, path, filename):
        tracer_filename = path + '/' + filename
        self.__petsc_vector['data'] = self.__tracer[number]
        self.__petsc_vector.tofile(tracer_filename)
        print('Write tracer number {} as petsc vector to file {}.'.format(number, tracer_filename))


    def write_all_tracer_to_file(self, path, filename=['N_Initial_Tracer.petsc', 'DOP_Initial_Tracer.petsc', 'P_Initial_Tracer.petsc', 'Z_Initial_Tracer.petsc', 'D_Initial_Tracer.petsc']):
        for i in range(self.__tracer_count):
            self.write_tracer_to_file(i, path, filename[i])


    def __create_random_tracer(self, number, mean, sigma):
        if self.__mass_distribution[number] == 0.0:
            self.__tracer[number, :] = np.zeros(self.__length, dtype='>f8')
        elif self.__distribution == 'uniform':
            self.__tracer[number, :] = np.random.uniform(max(0, mean - sigma), max(0, mean + sigma), self.__length)
        elif self.__distribution == 'lognormal':
            self.__tracer[number, :] = np.random.lognormal(mean, sigma, self.__length)
        elif self.__distribution == 'normal':
            self.__tracer[number, :] = np.random.normal(mean, sigma, self.__length)
            # eliminate negative entrys
            for i in range(0,self.__length):
                if self.__tracer[number, i] < 0.0:
                    self.__tracer[number, i] = 0.0
        elif self.__distribution == 'OneBox':
            self.__tracer[number, :] = np.zeros(self.__length)
            self.__tracer[number, np.random.randint(self.__length)] = 1.0
        print('Create random tracer for tracer number {} with distribution {} and mean={}, sigma={}.'.format(number, self.__distribution, mean, sigma))


    def create_random_tracer(self):
        for i in range(self.__tracer_count):
            self.__create_random_tracer(i, self.__mean[i], self.__sigma_value[i])


    def set_tracer(self, number, boxnumber):
        """
        Initialize tracer with whole mass in the given box
        """
        self.__tracer[number,:] = np.zeros(self.__length)
        self.__tracer[number, boxnumber] = 1.0


    def __read_boxvolumen_from_file(self):
        petsc_vector_type = np.dtype([('header',('>i4',2)),('data','>f8',self.__length)])
        self.__petsc_vector = np.empty(shape=(1), dtype=petsc_vector_type)
        
        filename = self.__path_metos3d + '/data/data/TMM/2.8/Geometry/volumes.petsc'
        f = open(filename, 'rb')
        self.__petsc_vector['header'] = np.fromfile(f, dtype='>i4', count=2)
        self.__boxvolumes = np.empty(shape=(self.__length), dtype='>f8')
        self.__boxvolumes = np.fromfile(f, dtype='>f8', count=self.__length)
        f.close()


    def random_mass_distribution(self):
        self.__mass_distribution = np.ones(shape=(self.__tracer_count), dtype='>f8')
        if self.__tracer_count != 1:
            self.__mass_distribution[0] = np.random.rand()
            rest = 1 - self.__mass_distribution[0]
            for i in range(1, self.__tracer_count-1):
                self.__mass_distribution[i] = np.random.rand() * rest
                rest -= self.__mass_distribution[i]
            self.__mass_distribution[self.__tracer_count-1] = rest
            self.__mass_distribution = np.random.permutation(self.__mass_distribution)
        while np.sum(self.__mass_distribution) != 1.0:
            self.__mass_distribution[np.random.randint(self.__tracer_count)] += max(1.0 - np.sum(self.__mass_distribution), 0.0)
        print('Random distribution of mass for the tracer: {}'.format(self.__mass_distribution))


    def set_mass_distribution(self, mass_distribution):
        if ((self.__tracer_count == mass_distribution.size) and (np.sum(mass_distribution) == 1.0)):
            self.__mass_distribution = mass_distribution
        else:
            print("Sum of all elements in the distribution must one!")
        print('Set random distribution of mass for the tracer: {}'.format(self.__mass_distribution))


    def __calculate_total_mass(self, mass):
        self.__total_mass = np.sum(mass) * np.sum(self.__boxvolumes)
        print('Total mass: %1.20e' % (self.__total_mass))


    def set_total_mass(self,mass):
        self.__total_mass = mass
        print('Set total mass: %1.20e' % (self.__total_mass))


    def calculate_mean_sigma(self):
        self.__mean = np.zeros(shape=(self.__tracer_count), dtype='>f8')
        self.__sigma_value = np.ones(shape=(self.__tracer_count), dtype='>f8')
        for i in range(self.__tracer_count):
            self.__mean[i] = self.__total_mass * self.__mass_distribution[i] / np.sum(self.__boxvolumes)
            if self.__distribution == 'uniform':
                self.__sigma_value[i] = min(1.0, self.__mean[i])
            elif self.__distribution == 'lognormal':
                self.__sigma_value[i] = min(self.__sigma_value[i], self.__mean[i] / 2.0)
            elif self.__distribution == 'normal':
                self.__sigma_value[i] = self.__mean[i] / 3.0


    def calculate_mu_sigma(self):
        if self.__distribution == 'lognormal':
            for i in range(self.__tracer_count):
                mu = self.__mean[i]
                sigma = self.__sigma_value[i]
                self.__mean[i] = np.log(np.square(mu) / np.sqrt(sigma + np.square(mu)))
                self.__sigma_value[i] = np.sqrt(np.log(sigma / np.square(mu) + 1))
            

    def set_mean(self, mean):
        if self.__tracer_count == mean.size:
            self.__mean = mean


    def set_sigma(self, sigma):
        if self.__tracer_count == sigma.size:
            self.__sigma_value = sigma


    def __calculate_tracer_mass(self, number):
        tracer_mass = 0.0
        for i in range(self.__length):
            tracer_mass += self.__tracer[number, i] * self.__boxvolumes[i]
        return tracer_mass


    def masscorrection(self, number, eps=1e-20, iteration=1e3):
        if self.__mass_distribution[number] == 0:
            for i in range(self.__length):
                self.__tracer[number, i] = 0.0
            print('Change mass of tracer number {}: relativ error {}, absolut error {}'.format(number, 0,0))
        elif self.__distribution == 'OneBox':
            i = np.argmax(self.__tracer[number,:])
            self.__tracer[number,i] = self.__total_mass * self.__mass_distribution[number] / self.__boxvolumes[i]
            print('Change mass of tracer number {}: relativ error {}, absolut error {}'.format(number, 0,0))
        else:
            tracer_mass = self.__calculate_tracer_mass(number)
            mass = self.__total_mass * self.__mass_distribution[number]

            l = 0
            Change_Volumes_tracer = 0
            while ((abs((tracer_mass - mass) / mass) > eps) and (l < iteration)):
                i = np.random.randint(1, self.__length)
                changeEntrys = np.random.randint(self.__length, size=(i))
                diff = tracer_mass - mass
                part = diff / (i * 1.0)
                for j in range(i):
                    part_j = part * np.random.rand()
                    c = part_j / self.__boxvolumes[changeEntrys[j]]
                    if ((self.__tracer[number, changeEntrys[j]] - c >= 0.0) and (abs(tracer_mass - mass - part_j) < abs(diff))):
                        self.__tracer[number, changeEntrys[j]] -= c
                        tracer_mass -= part_j;
                        diff -= part_j
                        Change_Volumes_tracer += 1
                l = l+1

            l = 0
            while ((abs((tracer_mass - mass) / mass) > eps) and (l < iteration)):
                i = np.random.randint(self.__length)
                diff = tracer_mass - mass
                diff_c = diff * np.random.rand()
                c = diff_c / self.__boxvolumes[i]
                if (self.__tracer[number,i] - c >= 0):
                    self.__tracer[number,i] -= c
                    tracer_mass -= diff_c
                    Change_Volumes_tracer += 1
                l += 1

            while ((abs((tracer_mass - mass) / mass) > eps) and (l < 2*iteration)):
                i = np.random.randint(self.__length)
                diff = tracer_mass - mass
                c = diff / self.__boxvolumes[i]
                if (self.__tracer[number,i] - c >= 0):
                    self.__tracer[number,i] -= c
                    tracer_mass -= diff
                    Change_Volumes_tracer += 1
            print('Change mass of tracer number {}: Count of changes {}, relativ error {}, absolut error {}'.format(number, Change_Volumes_tracer, abs((tracer_mass - mass) / mass), tracer_mass - mass))

