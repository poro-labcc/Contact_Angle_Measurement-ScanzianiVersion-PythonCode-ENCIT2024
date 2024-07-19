#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:42:20 2024

@author: usuario
"""

#K_sub     -> Length of the cubic sub-volume for extracking 2d image
#k_line    -> Length of the fluid-fluid interface
#k_ROI     -> Size of the region for the moving average procedure
#N         -> Number of points for the moving average procedure
#Parallel  -> Number of Cores for Parallelization of the code


import numpy as np
from argparse import ArgumentParser
import os.path
import Library_Scanziani_Par as S

Radiuss = [  7,  14,  28,  56]
Size   = [ 38,  52,  80, 136]
Types  = [  0,   1,   2,   3]
Corner = [  0]
Angle  = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
Directory = '/home/usuario/ENCIT_2024/Benchmark'
Directory_out = '/home/usuario/ENCIT_2024/Medições/Scanziani_New2'

for i in range(len(Radiuss)):
    radius = Radiuss[i]
    size = Size[i]
    for types in Types:
        for angle in Angle:
            input_file = f'{Directory}/T{types}R{radius}x{size}/T{types}_R{radius}_A{angle}_x{size}.raw'
            if not os.path.exists(f'{Directory_out}/T{types}R{radius}x{size}'):
                os.mkdir(f'{Directory_out}/T{types}R{radius}x{size}')
            output_file_name = f'Scanziani_{types}_R{radius}_A{angle}_x{size}'
            npimg = np.fromfile(input_file, dtype=np.uint8)
            imageSize = (size, size, size)
            npimg = npimg.reshape(imageSize)
            #Result in radians
            
            try:
                theta, Radius, angles_position = S.Scanziani(npimg, k_sub = size, k_line = 50, k_ROI = 15, N = 4, Parallel = 6, directory = Directory, alpha = 1)
            except IndexError():
                print('err')

            np.save(f"{Directory_out}/T{types}R{radius}x{size}/{output_file_name}_Angle", theta)
            np.save(f"{Directory_out}/T{types}R{radius}x{size}/{output_file_name}_Position", angles_position)


            if not os.path.exists(f'{Directory_out}/T{types}R{radius}x{size}_inv'):
                os.mkdir(f'{Directory_out}/T{types}R{radius}x{size}_inv')
            npimg = np.where(npimg == 1, 2, np.where(npimg == 2, 1, npimg))


            try:
                theta, Radius, angles_position = S.Scanziani(npimg, k_sub = size, k_line = 50, k_ROI = 15, N = 4, Parallel = 6, directory = Directory, alpha = 1)
            except IndexError():
                print('err')             

            np.save(f"{Directory_out}/T{types}R{radius}x{size}_inv/{output_file_name}_Angle", theta)
            np.save(f"{Directory_out}/T{types}R{radius}x{size}_inv/{output_file_name}_Position", angles_position)



