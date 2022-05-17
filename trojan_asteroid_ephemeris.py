from astroquery.jplhorizons import Horizons
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.io import ascii
import imageio
import pandas as pd



telescope = '500@0'

# You can change these at will, but be warned that Horizons won't give you an ephemeris that has too many entries
# t0: animation start time
t0 = '2000-01-01'
# t1: animation stop time
t1 = '2013-01-01'
# dt: timestep between days
dt = '10d'
# these are the Lucy targets, but you can change these to whatever asteroids you want
# they don't even have to be Trojans!
asteroids = ['Patroclus', 'Leucus', 'Eurybates', 'Orus', 'Polymele']


# NOTE: Horizons requires an internet connection to work
def grab_data(ob, telescope, t0, t1, dt, convert_dms = True, full_eph = False):
    '''This grabs asteroid data specifically, converts it to plottable x & y coords'''
    obj = Horizons(id=ob, location=telescope, epochs={'start':t0, 'stop':t1, 'step':dt})
    eph = obj.ephemerides()
    radius = eph['r']
    long = eph['EclLon']
    x = radius * np.cos(long*2*np.pi/360)
    y = radius * np.sin(long*2*np.pi/360)
    return(x, y)

# This portion grabs Jupiter's coordinates. 
obj = Horizons(id='599', location=telescope, epochs={'start':t0, 'stop':t1, 'step':dt}, id_type='majorbody')
eph = obj.ephemerides()
radius = eph['r']
long = eph['EclLon']
j_x = radius * np.cos(long*2*np.pi/360)
j_y = radius * np.sin(long*2*np.pi/360)

xs = [j_x]
ys = [j_y]

names = ['Jupiter']

# This portion grabs ephemeris data for the asteroids
for i in range(len(asteroids)):
    aster = asteroids[i]
    name = asteroids[i]
    x, y = grab_data(aster, telescope, t0, t1, dt, True, True)
    xs.append(x)
    ys.append(y)
    names.append(name)
    # you can comment out the print statement, but I like to keep track of progress. 
    print(i)

# This saves the data to an offline file so you don't have to query Horizons every time
all_data = {
    'xs' : xs,
    'ys' : ys,
    'names' : names
}

import pickle
pickle.dump( all_data, open( "save.p", "wb" ) )

# Everything past this point can be run offline if you have the 'save.p' file

all_data = pd.read_pickle('save.p')

xs = np.transpose(all_data['xs'])
ys = np.transpose(all_data['ys'])


len(xs)
for i in range(len(xs)):
    plt.figure(figsize=(10, 10))
    plt.fill([-6, -6, 6, 6], [-6, 6, 6, -6], color='k')

    plt.ylim(-6, 6)
    plt.xlim(-6, 6)
    # This plots all the asteroids (and Jupiter)
    plt.plot(xs[i], ys[i],'.w')
    # This plots the sun at 0,0
    plt.plot(0, 0, 'o' , markersize = 12, color = '#ffff00')
    # This plots Jupiter's orbit
    plt.plot(all_data['xs'][0], all_data['ys'][0], 'r')
    # This plots Jupiter, but large and red so it's more obvious
    plt.plot(xs[i][0], ys[i][0], 'or')
    
    # This is an example of plotting the orbit of an individual asteroid
    # For this example, I'm using Patroclus, which has an index of 1
    # you can change which asteroid this is by changing the index, using the names variable for reference
    plt.plot(all_data['xs'][1], all_data['ys'][1], 'b')
    plt.plot(xs[i][1], ys[i][1], 'ob')
    
    
    plt.axis('off')
    axes=plt.gca()
    axes.set_aspect(1)
    plt.savefig("asteroid" + str(i)+".png", bbox_inches='tight')
    plt.close()
    
# to make gifs
images = []
filenames = np.asarray(range(len(xs)))
for filename in filenames:
    f = 'asteroid'+str(filename)+'.png'
    images.append(imageio.imread(f))
imageio.mimsave('asteroid.gif', images,  duration=0.05)
