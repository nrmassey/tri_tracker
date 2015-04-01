#!/usr/bin/env python
#******************************************************************************
#** Program : tri_area_levels.py
#** Author  : Neil Massey
#** Date    : 20/01/15
#** Purpose : calculate the area of a triangle for a number of different levels
#**           so they can be directly compared with grid box areas in gridpoint
#**           models
#******************************************************************************

from indexed_tri_3D import indexed_tri_3D
from vector_3D import vector_3D
from point_cloud import point_cloud
from tri_grid import *
from Geodetic import hD
from geo_convert import cart_to_model

import numpy
import os, getopt, sys

#******************************************************************************

def split_triangle(tri, pc):
    # spit the triangle into 4 new triangles and insert the new points into
    # the point cloud
    pc_idx = numpy.zeros(3, 'i')
    for i in range(0,3):
        if i < 2:
            j = i + 1
        else:
            j = 0
        # create a new point as the midpoint
        pn = (tri[i] + tri[j]) * 0.5
        # normalise the midpoint so that it is projected to a sphere
        pn = pn * float(1.0 / pn.mag())
        pc.point_cloud_store.append(pn)
        pc_idx[i] = len(pc.point_cloud_store) - 1
    # create four triangles from the vectors
    nt0 = indexed_tri_3D(pc)
    nt0.pc_idx = numpy.array([tri.pc_idx[0], pc_idx[0], pc_idx[2]],'i')

    nt1 = indexed_tri_3D(pc)
    nt1.pc_idx = numpy.array([pc_idx[0], tri.pc_idx[1], pc_idx[1]],'i')

    nt2 = indexed_tri_3D(pc)
    nt2.pc_idx = numpy.array([pc_idx[1], tri.pc_idx[2], pc_idx[2]],'i')

    nt3 = indexed_tri_3D(pc)
    nt3.pc_idx = numpy.array([pc_idx[0],     pc_idx[1], pc_idx[2]],'i')
    
    return nt0, nt3     # only need one triangle to carry on

#******************************************************************************

def calc_area(T, d, Ar):
    Tr = numpy.radians(T)
    dr = numpy.radians(d)
    EARTH_R = 6371 * 1000
#    print d, d*Ar
    area = abs(2*dr * numpy.cos(Tr) * numpy.sin(0.5*(dr*Ar)) * EARTH_R**2)
    return area

#******************************************************************************

def find_equivalent_resolution_area(area, T, Ar):
    # find the equivalent grid point resolution for a triangle level
    # Ar = aspect ratio
    d = 1.0 # initial guess, 1x1*Ar degree

    diff = 0.001
    guess = True
     
    while guess:
        A = calc_area(T, d, Ar)
        if A > area:
            d = d - diff * area/A
        if A < area:
            d = d + diff * area/A
        if abs(1.0 - A/area) < 1e-3:
            guess = False
    return d

#******************************************************************************

def calc_distance(tri0, tri1):
    # calculate the distance between the 
    C0 = tri0.get_centroid()
    C1 = tri1.get_centroid()
    # use the haversine distance
    M0 = cart_to_model(C0)
    M1 = cart_to_model(C1)
    # convert to radians
    d = hD(M0[1], M0[0], M1[1], M1[0])
    return d * 1e-3

#******************************************************************************

def find_equivalent_resolution_distance(T_dist, T, Ar):
    d = 1.0 # initial guess, 1x1*Ar degree
    diff = 0.001
    guess = True
    G_dist = hD(T, 0.0, T, d) * 1e-3
    
    while guess:
        # get the distance between the grid boxes
        G_dist = hD(T, 0.0, T, d) * 1e-3
        if G_dist > T_dist:
            d = d - diff * T_dist / G_dist
        if G_dist < T_dist:
            d = d + diff * T_dist / G_dist
        if abs(1.0 - T_dist/G_dist) < 1e-3:
            guess = False
    return d

#******************************************************************************

if __name__ == "__main__":
    # initialise the point cloud with 3 points.  The coordinates and reasoning
    # can be found in Massey 2012 (Comp and Geo).
    b   = numpy.cos(numpy.pi/5.0);
    c   = numpy.sqrt(1+4*b*b);
    phi = (2 * numpy.cos(numpy.pi/5.0)) / c;
    the = 1.0 / c;

    # create the 3 vectors
    p0 = vector_3D( 0.0, -phi,  the)        # idx 1
    p1 = vector_3D( 0.0, -phi, -the)        # idx 3
    p2 = vector_3D( phi, -the,  0.0)        # idx 10
    p3 = vector_3D( the,  0.0,  phi)        # idx 4

    # create the points in the point cloud
    pc = point_cloud()
    pc.point_cloud_store.extend([p0,p1,p2,p3])

    # now create the triangle
    tri0_O = indexed_tri_3D(pc)
    tri0_O.pc_idx = numpy.array([0,1,2], 'i')
    tri1_O = indexed_tri_3D(pc)
    tri1_O.pc_idx = numpy.array([3,0,2], 'i')
    
    tri0 = tri0_O
    tri1 = tri1_O

    # split the triangle at a number of levels to output the surface area
    # at each triangle level
    d = 1.0
    Ar = 1.0
    Area = True
    if Area:
        for L in range(0,11):
            T0 = tri0.surface_area()
            d0 = find_equivalent_resolution_area(T0, 80, Ar)
            d1 = find_equivalent_resolution_area(T0, 60, Ar)
            d2 = find_equivalent_resolution_area(T0, 40, Ar)
            d3 = find_equivalent_resolution_area(T0, 20, Ar)
            d4 = find_equivalent_resolution_area(T0, 0, Ar)
        
            A0 = calc_area(80, d0, Ar) / (1000**2)
            A1 = calc_area(60, d1, Ar) / (1000**2)
            A2 = calc_area(40, d2, Ar) / (1000**2)
            A3 = calc_area(20, d3, Ar) / (1000**2)
            A4 = calc_area(0,  d4, Ar) / (1000**2)
        
            print str(L) +" & %i" % int(0.5+T0 / (1000**2)) + \
                          " & %.2f" % d0 +\
                          " & %.2f" % d1 +\
                          " & %.2f" % d2 +\
                          " & %.2f" % d3 +\
                          " & %.2f" % d4 +\
                          " \\\\"
            tri0, tri1 = split_triangle(tri0, pc)

    # reset the triangles
    tri0 = tri0_O
    tri1 = tri1_O

    for L in range(0,11):
        D0 = calc_distance(tri0, tri1)
        if L<=0:
            d0 = 360.0
        else:
            d0 = find_equivalent_resolution_distance(D0, 80, Ar)
        d1 = find_equivalent_resolution_distance(D0, 60, Ar)
        d2 = find_equivalent_resolution_distance(D0, 40, Ar)
        d3 = find_equivalent_resolution_distance(D0, 20, Ar)
        d4 = find_equivalent_resolution_distance(D0, 0, Ar)
        print str(L) +" & %i" % int(D0) + \
                      " & %.2f" % d0 +\
                      " & %.2f" % d1 +\
                      " & %.2f" % d2 +\
                      " & %.2f" % d3 +\
                      " & %.2f" % d4 +\
                      " \\\\"
        tri0, tri1 = split_triangle(tri0, pc)
