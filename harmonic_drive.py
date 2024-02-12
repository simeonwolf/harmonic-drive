import numpy as np
from scipy import optimize
import scipy.special

import shapely as sp

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpl_patches

from tqdm import tqdm

def ROT(phi):
    '''2x2 rotation matrix with angle phi'''
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s],
                     [s,  c]])

def ROTATE(vec_arr, phi):
    '''
    rotates a curve with angle phi around the origin/z-axis
    curve --> ([x,y],)-vectors in an (:,2)-array'''
    return np.array([np.dot(ROT(phi), xy) for xy in vec_arr])

def MIRROR(vec_arr, axis = 'x'):
    '''
    mirrors a curve on the x or y-axis
    curve --> ([x,y],)-vectors in an (:,2)-array'''
    if axis == 'x':
        return np.flip(np.array([np.array([x,-y]) for x,y in vec_arr]), axis = 0)
    if axis == 'y':
        return np.flip(np.array([np.array([-x,y]) for x,y in vec_arr]), axis = 0)

def length_vec_array(vec_arr):
    '''
    length of a curve (approximation)
    curve --> ([x,y],)-vectors in an (:,2)-array'''
    tangent = np.gradient(vec_arr)[0]
    return np.sum( np.array( [np.linalg.norm(t) for t in tangent] ) )

def gear_from_flank(flank, z):
    '''
    Generates the contour of the gear from a tooth flank and the number of teeth by rotating and mirroring
    tooth flank --> curve in the first quadrant on the x-axis'''
    gear = flank
    mirrored = MIRROR(flank)

    for idx in range(1,int(z)):
        gear = np.append(gear, ROTATE(mirrored, 2*np.pi/z*idx), axis = 0)
        gear = np.append(gear, ROTATE(flank, 2*np.pi/z*idx), axis = 0)

    gear = np.append(gear, mirrored, axis = 0)

    return gear

def plot_polygon(ax, poly, **kwargs):
    '''
    plots a shapely polygon to pyplot axis, supports polygons with holes'''
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])
    
    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)
    ax.add_collection(collection, autolim=True)
    
    return collection

def plot_va(ax, vec_arr, **kwargs):
    '''
    plots a curve to pyplot axis
    curve --> ([x,y],)-vectors in an (:,2)-array'''
    ax.plot(vec_arr[:,0], vec_arr[:,1], **kwargs)

def plot_grid(ax):
    '''
    adds a grid to pyplot axis'''
    ax.axis('equal')
    ax.grid(which='major', color='#DDDDDD', linewidth=1)
    ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
    ax.minorticks_on()
    ax.set_axisbelow(True)

class Flexspline():
    '''
    represents the Flexspline (FS) of a Harmonic Drive gearbox
    --> generates the initial tooth flank geometry of the Flexspline'''
    def __init__(self):
        #parameters that trigger a parameter update when changed:
        self.update_parameter = {'z', 'd_i', 's_st', 'd', 
                                 'd_h', 'alpha', 'c', 'r_fh', 'r_ff', 
                                 'r_hr', 'r_fr'}
        
        #---general parameters---:
        self.z      = 20    #number of teeth
        self.d_i    = 75    #Inner diameter
        self.s_st   = 2.5   #Sprocket Thickness
        self.d_nf   = []    #Neutral Fiber Diameter
        self.u_nf   = []    #circumference neutral fiber
        self.d_f    = []    #Foot circle diameter
        self.d      = 81    #Pitch circle diameter
        self.d_h    = 82    #Head circle diameter

        #---tooth parameters---:
        self.alpha  = 20*np.pi/180  #Tooth flank angle
        self.c      = 1             #Tooth space ratio
        self.r_fh   = 5             #Head flank radius
        self.r_ff   = 5             #Foot flank radius
        self.r_hr   = 0.3           #Head rounding
        self.r_fr   = 0.3           #Foot rounding
        self.h_t    = []            #Tooth Height

        self.flank = [] #curve of tooth-flank
        self.tooth = [] #curve of tooth
    
        self.r_z_i         = []
        self.CS2_i         = []
        self.inside_fibre  = []

        self.update()

    def __setattr__(self, name, value):
        if name == 'update_parameter':
            self.__dict__[name] = value
            return
        
        if name in self.update_parameter:
            self.__dict__[name] = value
            try:
                self.update()
            except:
                pass
        else:
            self.__dict__[name] = value
        
    def update(self):
        self.d_nf = self.d_i + self.s_st
        self.d_f  = self.d_i + self.s_st * 2 
        self.u_nf = self.d_nf * np.pi
        self.h_t  = 0.5*(self.d_h-self.d_f)
        
        self.flank = self.calculate_flank(num_of_discretization = 100)
        self.tooth = self.calculate_tooth()
        
    def calculate_flank(self, num_of_discretization = 100):
        '''
        calculates the tooth flank from 6 tangential arc sections with the radii r_fl_i
        Part of the tooth flank are the head and foot circle (r_fl_1, r_fl_6)'''
        r_fl_1 = self.d_h/2
        r_fl_2 = self.r_hr

        r_fl_3 = self.r_fh
        r_fl_4 = self.r_ff
        r_fl_5 = self.r_fr
        r_fl_6 = self.d_f/2
        
        beta_6 = np.pi/self.z
        beta_3 = beta_6/2*self.c
        
        s_34 = self.d/2 * np.array([np.cos(beta_3),np.sin(beta_3)])
        m_3 = s_34 - r_fl_3 * np.array([np.sin(self.alpha),np.cos(self.alpha)])
        m_4 = s_34 + r_fl_4 * np.array([np.sin(self.alpha),np.cos(self.alpha)])
          
        #intersection m_2
        gamma = np.arccos( ((m_3[0]**2+m_3[1]**2) + (r_fl_1-r_fl_2)**2 - (r_fl_3-r_fl_2)**2 )/(2*np.linalg.norm(m_3)*(r_fl_1-r_fl_2)))
        delta = np.abs(np.arctan(m_3[1]/m_3[0]))
        
        m_2 = (r_fl_1-r_fl_2) * np.array([np.cos(gamma-delta),np.sin(gamma-delta)])
        beta_1 = np.abs(np.arctan(m_2[1]/m_2[0]))
        
        #intersection s_23
        s_23 = m_3 + (m_2-m_3)/np.linalg.norm(m_2-m_3) * r_fl_3
        beta_2 = np.abs(np.arctan(s_23[1]/s_23[0]))
        
        #intersection m_5
        gamma = np.arccos( ((m_4[0]**2+m_4[1]**2) + (r_fl_6+r_fl_5)**2 - (r_fl_4-r_fl_5)**2 )/(2*np.linalg.norm(m_4)*(r_fl_6+r_fl_5)))
        delta = np.abs(np.arctan(m_4[1]/m_4[0]))
        
        m_5 = (r_fl_6+r_fl_5) * np.array([np.cos(delta-gamma),np.sin(delta-gamma)])
        beta_5 = np.abs(np.arctan(m_5[1]/m_5[0]))
        
        #intersection s_23
        s_45 = m_5 + (m_5-m_4)/np.linalg.norm(m_5-m_4) * r_fl_5
        beta_4 = np.abs(np.arctan(s_45[1]/s_45[0]))
        
        circular_arcs = [
            {"start": 0,     "end": beta_1, "m": np.array([0,0]), "r": r_fl_1},
            {"start": beta_1,"end": beta_2, "m": m_2            , "r": r_fl_2, "gamma": np.arctan(m_2[1]/m_2[0])},
            {"start": beta_2,"end": beta_3, "m": m_3            , "r": r_fl_3, "gamma": np.arctan(m_3[1]/m_3[0])},
            {"start": beta_3,"end": beta_4, "m": m_4            , "r": -r_fl_4, "gamma": np.arctan(m_4[1]/m_4[0])},
            {"start": beta_4,"end": beta_5, "m": m_5            , "r": -r_fl_5, "gamma": np.arctan(m_5[1]/m_5[0])},
            {"start": beta_5,"end": beta_6, "m": np.array([0,0]), "r": r_fl_6}
        ]
        
        beta_flank = np.linspace(0, beta_6, num_of_discretization)
        flank = np.zeros((len(beta_flank),2))
        
        for idx in range(len(flank)):
            if circular_arcs[0]["start"] <= beta_flank[idx] and beta_flank[idx] < circular_arcs[0]["end"]:
                flank[idx] = circular_arcs[0]["m"] + circular_arcs[0]["r"] * np.array([np.cos(beta_flank[idx]),np.sin(beta_flank[idx])])
                continue
            if circular_arcs[1]["start"] <= beta_flank[idx] and beta_flank[idx] < circular_arcs[1]["end"]:        
                epsilon = np.arcsin(np.sin(beta_flank[idx] - circular_arcs[1]["gamma"]) * np.linalg.norm(circular_arcs[1]["m"]) / circular_arcs[1]["r"])
                delta = beta_flank[idx] + epsilon
                flank[idx] = circular_arcs[1]["m"] + circular_arcs[1]["r"] * np.array([np.cos(delta), np.sin(delta)])
                continue
            if circular_arcs[2]["start"] <= beta_flank[idx] and beta_flank[idx] < circular_arcs[2]["end"]:        
                epsilon = np.arcsin(np.sin(beta_flank[idx] - circular_arcs[2]["gamma"]) * np.linalg.norm(circular_arcs[2]["m"]) / circular_arcs[2]["r"])
                delta = beta_flank[idx] + epsilon
                flank[idx] = circular_arcs[2]["m"] + circular_arcs[2]["r"] * np.array([np.cos(delta), np.sin(delta)])
                continue
            if circular_arcs[3]["start"] <= beta_flank[idx] and beta_flank[idx] < circular_arcs[3]["end"]:
                epsilon = np.arcsin(np.sin(beta_flank[idx] - circular_arcs[3]["gamma"]) * np.linalg.norm(circular_arcs[3]["m"]) / circular_arcs[3]["r"])
                delta = beta_flank[idx] + epsilon
                flank[idx] = circular_arcs[3]["m"] + circular_arcs[3]["r"] * np.array([np.cos(delta), np.sin(delta)])
                continue
            if circular_arcs[4]["start"] <= beta_flank[idx] and beta_flank[idx] < circular_arcs[4]["end"]:
                epsilon = np.arcsin(np.sin(beta_flank[idx] - circular_arcs[4]["gamma"]) * np.linalg.norm(circular_arcs[4]["m"]) / circular_arcs[4]["r"])
                delta = beta_flank[idx] + epsilon
                flank[idx] = circular_arcs[4]["m"] + circular_arcs[4]["r"] * np.array([np.cos(delta), np.sin(delta)])
                continue
            if circular_arcs[5]["start"] <= beta_flank[idx] and beta_flank[idx] <= circular_arcs[5]["end"]:
                flank[idx] = circular_arcs[5]["m"] + circular_arcs[5]["r"] * np.array([np.cos(beta_flank[idx]),np.sin(beta_flank[idx])])
                continue
        
        self.circular_arcs = circular_arcs
        
        self.r_fl_2 = r_fl_2    
        self.r_fl_3 = r_fl_3
        self.r_fl_4 = r_fl_4
        self.r_fl_5 = r_fl_5
        
        self.m_2 = m_2
        self.m_3 = m_3
        self.m_4 = m_4
        self.m_5 = m_5
        
        self.s_34 = s_34
        
        return flank
    
    def calculate_flank_normal(self, num_of_discretization = 100):
        '''
        returns the tooth flank and the corresponding normal vectors (normalized)
        in the tooth coordinate system
        ! recalculation of the tooth flank according to num_of_discretization !'''
        flank = self.calculate_flank(num_of_discretization = num_of_discretization)
        flank = np.array([np.array([-y,x-self.d_nf/2]) for x,y in flank])
        
        #TODO More precise way to calculate the normal vector from the 
        #definition of the tooth flank:
        tangent = np.gradient(flank)[0]
        normal = np.array([[y,-x]for x,y in tangent])
        normal_flank = np.array([1/np.linalg.norm(n)*n for n in normal])

        return flank, normal_flank
    
    def calculate_tooth(self):
        '''
        returns the tooth (two flanks) in the tooth coordinate system'''
        tooth_flank_mirrored = MIRROR(self.flank)
        tooth = np.append(tooth_flank_mirrored, self.flank, axis = 0)
                        
        return np.array([np.array([-y,x-self.d_nf/2]) for x,y in tooth])
    
    def gear(self):
        '''
        contour of the gear'''
        return gear_from_flank(self.flank, self.z)
    
    def polygon(self):
        '''
        displays the gear contour as a plygon (e.g. to display the gear wheel in color)'''
        gear = self.gear()
        gear_polygon = sp.Polygon(gear).buffer(0)
        
        ALPHA = np.linspace(0,2*np.pi,int(len(gear)/10))
        inner_circle_polygon = sp.Polygon([[self.d_i/2*np.cos(alpha), self.d_i/2*np.sin(alpha)] for alpha in ALPHA])
        
        return gear_polygon.difference(inner_circle_polygon)
    
    def gear_deformed(self):
        '''
        deformed contour of the gear
        ! only works after r_z_i (vectors to toot i) and CS2_i (coordinate system of toot i)
        have been set'''
        gear_deformed = self.r_z_i[0] + np.array([np.dot(self.CS2_i[0], xy) for xy in self.tooth])
        for i in range(1,int(self.z)):
            tooth_i = np.array([self.r_z_i[i] + np.dot(self.CS2_i[i], xy) for xy in self.tooth])
            gear_deformed = np.append(gear_deformed, tooth_i, axis = 0)
        
        return gear_deformed
    
    def polygon_deformed(self):
        '''
        displays the deformed contour of the gear as a plygon (e.g. to display the gear wheel in color)
        ! only works after r_z_i, CS2_i and inside_fibre have been set'''
        gear_deformed = self.gear_deformed()
        polygon_gear_deformed = sp.Polygon(gear_deformed).buffer(0)
        
        inside_fibre_polygon = sp.Polygon(self.inside_fibre)
        
        return polygon_gear_deformed.difference(inside_fibre_polygon)
        
class Wavegenerator():
    '''
    represents the Wavegenerator (WG) of a Harmonic Drive gearbox
    --> generates the initial profil of the wave generator
    currently 2 available profil shapes:
        1) Elliptical       (shape = 1)
        2) 345-Polynomial   (shape = 2)''' 
    def __init__(self, shape = 1):
        #parameters that trigger a parameter update when changed:
        self.update_parameter = {'a', 'b', 'arc', 'shape'}
        
        self.a = 1      #Large semi-axis
        self.b = 0.9    #Small semi-axis
        self.arc = 2 * np.pi/180 #angle of the arc, only used for 345polynomial profil
        
        if shape == 1:
            self.shape = 'elliptical'
        elif shape == 2:
            self.shape = '345polynomial'
        
        self.profil = []
        
    def __setattr__(self, name, value):
        if name == 'update_parameter':
            self.__dict__[name] = value
            return
        
        if name in self.update_parameter:
            self.__dict__[name] = value
            try:
                self.update()
            except:
                pass
        else:
            self.__dict__[name] = value
        

    def update(self):
        self.profil = self.calculate_profil(num_of_discretization = 1000)
    
    def calculate_profil(self, num_of_discretization = 1000):
        '''
        calculates the oval profile of the wave generator
        1) Elliptical or
        2) 345-Polynomial
        --> returns curve in parameter representation
        (curve --> ([x,y],)-vectors in an (:,2)-array)'''
        
        if self.shape == 'elliptical':
            t = np.linspace(0,2*np.pi,num_of_discretization)
            return np.array([np.array([self.a * np.cos(t),
                                       self.b * np.sin(t)]) for t in t])
            
        if self.shape == '345polynomial':
            #composed of 345 polynomial and circular arc
            T = np.linspace(0,2*np.pi,num_of_discretization)
            quarter = int(num_of_discretization/4)
            rho = np.zeros(int(quarter))
            for i in range(int(quarter)):
                if T[i] < self.arc:
                    rho[i] = self.a
                else:
                    z = 1-(T[i]-self.arc)/(np.pi/2-self.arc)
                    rho[i] = self.b+(self.a-self.b)*( 10*z**3 - 15*z**4 + 6*z**5)
            rho = np.append(rho, [np.flip(rho), rho, np.flip(rho)])
            
            return np.array([[r*np.cos(t), r*np.sin(t)] for r, t in zip(rho, T)])
        
    def calculate_profil_length(self, a, b):
        '''
        calculates the length of the profile of the wave generator
        1) Elliptical       --> numerical approximation via elliptic integral (scipy.special.ellipe))
        2) 345-Polynomial   --> numerical approximation via Newton-Cotes(mid-point rule)'''
        
        if self.shape == 'elliptical':
            m = 1 - b**2/a**2
            return 4*a*scipy.special.ellipe(m)
        
        if self.shape == '345polynomial':
            l_1 = a * self.arc
                        
            num_of_discretization = 10000
            dt = (np.pi/2-self.arc)/num_of_discretization
            
            l_2 = 0
            for i in range(num_of_discretization):
                t_i = self.arc + dt/2 + dt * i
                
                z = 1-(t_i-self.arc)/(np.pi/2-self.arc)
                z_dash = -1/(np.pi/2-self.arc)
                rho = b+(a-b)*( 10*z**3 - 15*z**4 + 6*z**5)
                rho_dash = (a-b)*(30*z**2*z_dash-60*z**3*z_dash+30*z**4*z_dash)
                
                l_2 = l_2 + np.sqrt( rho_dash**2 + rho**2 ) * dt
            
            return 4*(l_1+l_2)
            
    def aequidistant(self, q, num_of_discretization = 1000):
        '''
        returns the curve equidistant to the profile of the wave generator (distance q)
        (curve --> ([x,y],)-vectors in an (:,2)-array)'''
        aequidistant = np.zeros((num_of_discretization, 2))
        
        if self.shape == 'elliptical':
            T = np.linspace(0,2*np.pi,num_of_discretization)
            for i in range(num_of_discretization):
                x = self.a * np.cos(T[i])
                y = self.b * np.sin(T[i])
                n = np.array([x/self.a**2,y/self.b**2])
                n = 1/np.linalg.norm(n) * n
                aequidistant[i] = np.array([x, y] + n * q)
            return aequidistant
        
        if self.shape == '345polynomial':
            T = np.linspace(0,2*np.pi,num_of_discretization)            
            quarter = int(num_of_discretization/4)
            aequidistant_quarter = np.zeros((quarter, 2))
            for i in range(quarter):
                if T[i] < self.arc:
                    aequidistant_quarter[i] = np.array([np.cos(T[i]), np.sin(T[i])])*(self.a + q)
                else:
                    z = 1-(T[i]-self.arc)/(np.pi/2-self.arc)
                    z_dash = -1/(np.pi/2-self.arc)
                    rho = self.b+(self.a-self.b)*( 10*z**3 - 15*z**4 + 6*z**5)
                    rho_dash = (self.a-self.b)*(30*z**2*z_dash-60*z**3*z_dash+30*z**4*z_dash)
                    c = np.cos(T[i])
                    s = np.sin(T[i])
                    n = rho_dash*np.array([s,-c])+rho*np.array([c,s])
                    n = 1/np.linalg.norm(n) * n
                    aequidistant_quarter[i] = rho * np.array([c,s]) + n * q
            aequidistant[0:quarter]         = aequidistant_quarter
            aequidistant[quarter:2*quarter] = np.flip([np.array([-x,y]) for x,y in aequidistant_quarter], axis = 0)
            aequidistant[2*quarter:3*quarter] = [np.array([-x,-y]) for x,y in aequidistant_quarter]
            aequidistant[3*quarter:] = np.flip([np.array([x,-y]) for x,y in aequidistant_quarter], axis = 0)
            
            return aequidistant
    
    def point_aequidistant(self, q, t):
        '''
        returns a point to the angle t of a curve equidistant to the profile of
        the wave generator (distance q)'''
        
        if self.shape == 'elliptical':
            x = self.a * np.cos(t)
            y = self.b * np.sin(t)
            n = np.array([x/self.a**2,y/self.b**2])
            n = 1/np.linalg.norm(n) * n
            
            return np.array([x, y] + n * q)
        
        if self.shape == '345polynomial':
            if t <= np.pi/2:
                if t < self.arc:
                    return np.array([np.cos(t), np.sin(t)])*(self.a + q)
                else:
                    z = 1-(t-self.arc)/(np.pi/2-self.arc)
                    z_dash = -1/(np.pi/2-self.arc)
                    rho = self.b+(self.a-self.b)*( 10*z**3 - 15*z**4 + 6*z**5)
                    rho_dash = (self.a-self.b)*(30*z**2*z_dash-60*z**3*z_dash+30*z**4*z_dash)
                    c = np.cos(t)
                    s = np.sin(t)
                    n = rho_dash*np.array([s,-c])+rho*np.array([c,s])
                    n = 1/np.linalg.norm(n) * n
                    return rho * np.array([c,s]) + n * q
            if t < np.pi:
                x,y = self.point_aequidistant(q, np.pi -t)
                return np.array([-x,y])
            if t < 3/2*np.pi:
                x,y = self.point_aequidistant(q, t - np.pi)
                return np.array([-x,-y])
            if t < 2*np.pi:
                x,y = self.point_aequidistant(q, 2*np.pi-t)
                return np.array([x,-y])
            
    def tangent_aequidistant(self, q, t):
        '''
        returns the tangent to a point with the angle t of a curve equidistant to the profile of
        the wave generator (distance q)'''
        
        if self.shape == 'elliptical':
            x = self.a * np.cos(t)
            y = self.b * np.sin(t)
            tangent = np.array([-y/self.b**2, x/self.a**2])
            return 1/np.linalg.norm(tangent) * tangent
        
        if self.shape == '345polynomial':
            if t <= np.pi/2:
                if t < self.arc:
                    tangent = np.array([-np.sin(t), np.cos(t)])*(self.a + q)
                    return 1/np.linalg.norm(tangent) * tangent
                else:
                    z = 1-(t-self.arc)/(np.pi/2-self.arc)
                    z_dash = -1/(np.pi/2-self.arc)
                    rho = self.b+(self.a-self.b)*( 10*z**3 - 15*z**4 + 6*z**5)
                    rho_dash = (self.a-self.b)*(30*z**2*z_dash-60*z**3*z_dash+30*z**4*z_dash)
                    c = np.cos(t)
                    s = np.sin(t)
                    tangent = rho_dash*np.array([c, s])+rho*np.array([-s, c])
                    return 1/np.linalg.norm(tangent) * tangent
            if t < np.pi:
                x,y = self.tangent_aequidistant(q, np.pi -t)
                return np.array([x,-y])
            if t < 3/2*np.pi:
                x,y = self.tangent_aequidistant(q, t - np.pi)
                return np.array([-x,-y])
            if t < 2*np.pi:
                x,y = self.tangent_aequidistant(q, 2*np.pi-t)
                return np.array([-x,y])

    def profil_rotated(self, phi):
        '''
        returns the rotated profil of the Wavegenerator (angle phi)'''
        return ROTATE(self.profil, phi)
    
    def polygon(self, phi = 0, index = True):
        '''displays rotated profil of the Wavegenerator (angle phi)
        index: Possibility to add a hole to indicate the orientation'''
        if index == True:
            ALPHA = np.linspace(0,2*np.pi,100)
            m = np.dot(ROT(phi), self.profil[0]) * 0.95
            index_circle_polygon = sp.Polygon([[m[0] + np.cos(alpha), m[1] + np.sin(alpha)] for alpha in ALPHA])
            return sp.Polygon(ROTATE(self.profil, phi)).buffer(0).difference(index_circle_polygon)
        else:
            return sp.Polygon(ROTATE(self.profil, phi)).buffer(0)
        

class CircularSpline():
    '''
    represents the Circular Spline (CS) of a Harmonic Drive gearbox'''
    def __init__(self):
        self.z      = []
        self.flank  = []
                
    def gear(self):
        '''
        contour of the gear'''
        flank = np.array([xy for xy in self.flank if np.arctan(xy[1]/xy[0]) < np.pi/self.z*2 ])
        return gear_from_flank(flank, int(self.z))
        #return gear_from_flank(self.flank, int(self.z))
        
    def polygon(self):
        '''
        displays the gear contour as a plygon (e.g. to display the gear wheel in color)'''
        gear = self.gear()
        gear_polygon = sp.Polygon(gear).buffer(0)
        
        ALPHA = np.linspace(0,2*np.pi,int(len(gear)/100))
        d_o = (gear[0,0]*2) * 1.1
        outer_circle_polygon = sp.Polygon([[d_o/2*np.cos(alpha), d_o/2*np.sin(alpha)] for alpha in ALPHA])
        
        return outer_circle_polygon.difference(gear_polygon)
    
class DynamicSpline():
    '''
    represents the Dynamic Spline (DS) of a Harmonic Drive gearbox'''
    def __init__(self):
        self.z      = []
        self.flank  = []
                
    def gear(self):
        '''
        contour of the gear'''
        flank = np.array([xy for xy in MIRROR(self.flank) if np.arctan(xy[1]/xy[0]) < np.pi/self.z ])
        return gear_from_flank(flank, self.z)
        
    def polygon(self, phi = 0, index = True):
        '''
        displays the gear contour as a plygon (e.g. to display the gear wheel in color)
        index: Possibility to add a hole to indicate the orientation'''
        gear = self.gear()
        if index == True:
            ALPHA = np.linspace(0,2*np.pi,int(len(gear)/100))
            d_o = (gear[0,0]*2) * 1.05 #outer diameter
            outer_circle_polygon = sp.Polygon([[d_o/2*np.cos(alpha), d_o/2*np.sin(alpha)] for alpha in ALPHA])
            m = np.dot(ROT(phi), gear[0]) * 1.05
            index_circle_polygon = sp.Polygon([[m[0] + np.cos(alpha), m[1] + np.sin(alpha)] for alpha in ALPHA])
        
            gear_polygon = sp.Polygon(ROTATE(gear, phi)).buffer(0)
        
            return outer_circle_polygon.difference(gear_polygon).difference(index_circle_polygon)
        else:
            return sp.Polygon(ROTATE(gear, phi)).buffer(0)

class HarmonicDrive():
    def __init__(self, flexspline_inst, wavegenerator_inst, circularspline_inst, dynamicspline_inst):
        #parameters that trigger a parameter update when changed:
        self.update_parameter = {'i', 'd_br'}
        
        #objects
        self.fs = flexspline_inst
        self.wg = wavegenerator_inst
        self.cs = circularspline_inst
        self.ds = dynamicspline_inst
        
        #parameter
        self.i = -20        #Transmission Ratio (i_wg->ds)
        self.phi_wg = 0     #Drive Angle Wavegenerator
        self.phi_ds = 0     #Output Angle Dynamic Spline
        
        self.d_br = 2       #Diameter Wavegenertor bearing Rolling element
        self.q_nf = []      #Equidistant distance neutral fiber
                
    def __setattr__(self, name, value):
        if name == 'update_parameter':
            self.__dict__[name] = value
            return
        
        if name in self.update_parameter:
            self.__dict__[name] = value
            try:               
                self.update()
            except:
                pass
        elif name == 'phi_wg':
            self.__dict__[name] = value
            self.phi_ds = self.phi_wg/self.i
        else:
            self.__dict__[name] = value

    def update(self):
        print("hd --> update")
        #Number of teeth
        self.fs.z = 2*(-self.i)
        self.cs.z = self.fs.z + 2
        self.ds.z = self.fs.z
        
        self.q_nf = self.d_br + self.fs.s_st/2
        
        #calculate a and b according to the kinematic boundary conditions
        self.wg.a = self.fs.u_nf * (-self.i + 1)/-self.i * 1/(2*np.pi) - self.q_nf
        self.wg.b = optimize.newton(func = self.optimize_u_nf__b, x0 = self.wg.a/2, args=(self.wg.a, self.q_nf, self.fs.u_nf),)
            
    def optimize_u_nf__b(self, b, *args):
        '''
        optimization function for calculating a and b'''
        (a, q_nf, u_nf) = args
        lenght_neutral_fibre = self.wg.calculate_profil_length(a,b) + 2*np.pi*q_nf
        return lenght_neutral_fibre - u_nf
    
    def neutral_fibre(self, num_of_discretization = 1000):
        '''
        returns the neutral fiber (in the middle of the flexspline sprocket),
        which corresponds to the aequidistant to the profile of the wave generator'''
        return self.wg.aequidistant(self.q_nf, num_of_discretization = num_of_discretization)
    
    def neutral_fibre_rotated(self, phi):
        '''
        returns the rotated neutral fibre (angle phi)'''
        neutral_fibre = self.neutral_fibre()
        return ROTATE(neutral_fibre, phi)
    
    def inside_fibre_rotated(self, phi, num_of_discretization = 1000):
        '''
        returns the rotated inside fibre of the Flexspline (angle phi)'''
        q_nf = self.q_nf
        s_st = self.fs.s_st
        inside_fibre = self.wg.aequidistant(q_nf-s_st/2, num_of_discretization = num_of_discretization)
        return ROTATE(inside_fibre, phi)
    
    def kincematics(self, num_of_discretization = 100000, phi_wg_lim  = 2*np.pi):
        '''
        calculates the kinematics implicitly using the euler-forward-method, 
        the result is the trace of the point Z on the neutral fiber r_z
        --> phi_wg_lim: up to which drive angle r_z is calculated'''
        a = self.wg.a
        q_nf = self.q_nf
        u_nf = self.fs.u_nf
        i = self.i

        #lookup-table
        neutral_fibre = self.wg.aequidistant(self.q_nf, num_of_discretization = num_of_discretization)
        gradient = np.gradient(neutral_fibre)[0]
        
        l_nf = np.cumsum( np.array([np.linalg.norm(xy) for xy in gradient]) ) #length of neutral fibre
        l_nf_t = np.linspace(0,2*np.pi,len(neutral_fibre)) #t to length of neutral fibre
        
        #euler forward
        d_phi_wg = phi_wg_lim/num_of_discretization
        
        phi_wg_i = 0
        s_i = 0
        r_z = np.zeros((num_of_discretization,2))
        r_z[0] = np.array([a + q_nf,0])
        e_2_1 = np.zeros((num_of_discretization,2))
        e_2_2 = np.zeros((num_of_discretization,2))
        
        for idx in range(num_of_discretization-1):
            while s_i >= u_nf:
                s_i = s_i - u_nf
            t = np.interp(s_i,l_nf,l_nf_t)
            
            tangent = self.wg.tangent_aequidistant(q_nf, t)
            
            v_r = tangent * (a + q_nf)       
            v_r = np.dot(ROT(phi_wg_i), v_r)
            
            v_t = np.array([r_z[idx,1], -r_z[idx,0]])
            
            v_abs = v_r + v_t
            
            phi_wg_i = phi_wg_i - d_phi_wg
            s_i = s_i + (a + q_nf) * d_phi_wg
            
            r_z[idx+1] = r_z[idx] + v_abs * d_phi_wg
            
            e_2_1[idx] = -1/np.linalg.norm(v_r)*v_r
            e_2_2[idx] = np.array([-e_2_1[idx,1],e_2_1[idx,0]])
        
        phi_wg_r_z = np.linspace(0,phi_wg_lim,num_of_discretization)
        phi_ds_r_z = -phi_wg_r_z/i
                
        _3_r_z = np.array([ np.dot(ROT(-phi),xy) for phi, xy in zip(phi_ds_r_z, r_z)])
        _3_e_2_1 = np.array([ np.dot(ROT(-phi),xy) for phi, xy in zip(phi_ds_r_z, e_2_1)])
        _3_e_2_2 = np.array([ np.dot(ROT(-phi),xy) for phi, xy in zip(phi_ds_r_z, e_2_2)])
        
        return r_z, e_2_1, e_2_2, _3_r_z, _3_e_2_1, _3_e_2_2, l_nf, l_nf_t, phi_wg_r_z, phi_ds_r_z
    
    def calculate_kinematics(self, num_of_discretization = 10000):
        print("calculate kinematics...")

        (r_z, e_2_1, e_2_2, _3_r_z, _3_e_2_1, _3_e_2_2, l_nf, l_nf_t, phi_wg_r_z, phi_ds_r_z) =  self.kincematics(num_of_discretization, phi_wg_lim = 2*np.pi)

        self.r_z = r_z
        self.e_2_1 = e_2_1
        self.e_2_2 = e_2_2
        
        self._3_r_z = _3_r_z
        self._3_e_2_1 = _3_e_2_1
        self._3_e_2_2 = _3_e_2_2
        
        self.l_nf = l_nf
        self.l_nf_t = l_nf_t
        
        self.phi_wg_r_z = phi_wg_r_z
        self.phi_ds_r_z = phi_ds_r_z
        
    def calculate_CS2_i(self):
        u_nf = self.fs.u_nf
        z_fs = self.fs.z
        l_nf = self.l_nf
        l_nf_t = self.l_nf_t
        a = self.wg.a
        b = self.wg.b
        q_nf = self.q_nf
        phi_wg = self.phi_wg
        r_z = self.r_z
        phi_wg_r_z = self.phi_wg_r_z
        e_2_1 = self.e_2_1
        e_2_2 = self.e_2_2
        i = self.i
        
        self.phi_ds = phi_wg/i
        r_z_i = np.zeros((int(z_fs),2))
        CS2_i = np.zeros((int(z_fs),2,2))
        
        #Tooth 1:
        s_0 = -phi_wg*(a+q_nf) #way tooth 1
                
        r_z_i[0] = np.array([np.interp(-phi_wg,phi_wg_r_z,r_z[:,0]),
                             np.interp(-phi_wg,phi_wg_r_z,r_z[:,1])])
        r_z_i[0 + int(z_fs/2)] = -r_z_i[0]
        
        e_2_1_i = np.array([np.interp(-phi_wg,phi_wg_r_z,e_2_1[:,0]),
                            np.interp(-phi_wg,phi_wg_r_z,e_2_1[:,1])])
        e_2_2_i = np.array([np.interp(-phi_wg,phi_wg_r_z,e_2_2[:,0]),
                            np.interp(-phi_wg,phi_wg_r_z,e_2_2[:,1])])
        
        CS2_i[0] = np.array([[e_2_1_i[0],e_2_2_i[0]],
                             [e_2_1_i[1],e_2_2_i[1]]])
        
        CS2_i[0 + int(z_fs/2)] = -CS2_i[0]
        
        for i_tooth in range(1,int(z_fs/2)):
            s_i = u_nf/z_fs * i_tooth + s_0
            
            if s_i > u_nf/2:
                s_i = s_i - u_nf/2
                s_0 = s_0 - u_nf/2
                phi_wg = phi_wg + np.pi
                            
            t_i = np.interp(s_i, l_nf, l_nf_t) 
            p_i = self.wg.point_aequidistant(q_nf, t_i)
            
            r_z_i[i_tooth] = np.dot(ROT(phi_wg), p_i)
            r_z_i[i_tooth + int(z_fs/2)] = -r_z_i[i_tooth]
            
            tangent = self.wg.tangent_aequidistant(q_nf, t_i)
            
            e_2_1_i = np.dot(ROT(phi_wg), -tangent)
            e_2_2_i = np.array([-e_2_1_i[1],e_2_1_i[0]])
            
            CS2_i[i_tooth] = np.array([[e_2_1_i[0],e_2_2_i[0]],
                                       [e_2_1_i[1],e_2_2_i[1]]])
            
            CS2_i[i_tooth + int(z_fs/2)] = -CS2_i[i_tooth]
        
        self.r_z_i = r_z_i
        self.CS2_i = CS2_i
            
        self.fs.r_z_i = r_z_i
        self.fs.CS2_i = CS2_i
        
        self.fs.inside_fibre = self.inside_fibre_rotated(phi_wg, num_of_discretization = 100)
        
    def calculate_flank_cs(self):
        print('HD     >calculate_flank_cs()')
        flank, normal = self.fs.calculate_flank_normal(num_of_discretization=50)  
    
        kinematics =  self.kincematics(num_of_discretization=10000, phi_wg_lim = np.pi/2)
        (r_z, e_2_1, e_2_2) = kinematics[:3]
        
        r_flank = np.zeros((len(flank), len(r_z), 2))
        n_flank = np.zeros((len(flank), len(r_z), 2))
    
        sc_pr_r_n = np.zeros((len(flank), len(r_z)))
    
        flank_cs = []
    
        for idx_flank in tqdm(range(len(flank))):
        #for idx_flank in range(len(flank)):  
            r_flank[idx_flank] = [xy_3 + np.dot(flank[idx_flank], np.vstack((xy_1,xy_2))) for xy_1, xy_2, xy_3  in zip(e_2_1, e_2_2, r_z)]
            n_flank[idx_flank] = [np.dot(normal[idx_flank], np.vstack((xy_1,xy_2))) for xy_1, xy_2  in zip(e_2_1, e_2_2)]
            
            v_r_flank = np.gradient(r_flank[idx_flank])[0]    

            sc_pr_r_n[idx_flank] = np.array([np.dot(xy_1, xy_2) for xy_1, xy_2 in zip(v_r_flank, n_flank[idx_flank])])
            
            sc = sc_pr_r_n[idx_flank][0:-2]
                
            idx = np.where(sc[:-1] * sc[1:] < 0)[0] + 1
            
            if len(idx) == 3:
                flank_cs.append(r_flank[idx_flank,idx[1]])
            if len(idx) == 2:
                flank_cs.append(r_flank[idx_flank,idx[0]])
                
        flank_cs = np.array(flank_cs)
        self.cs.flank = flank_cs
        return flank_cs

    def calculate_flank_ds(self):
        print('HD     >calculate_flank_ds()')
        flank, normal = self.fs.calculate_flank_normal(num_of_discretization=50)  
        
        flank = MIRROR(flank, axis = 'y')
        normal = MIRROR(normal, axis = 'y')
            
        kinematics =  self.kincematics(num_of_discretization=10000, phi_wg_lim = np.pi/2)
        (_3_r_z, _3_e_2_1, _3_e_2_2) = kinematics[3:6]
    
        r_flank = np.zeros((len(flank), len(_3_r_z), 2))
        n_flank = np.zeros((len(flank), len(_3_r_z), 2))
       
        flank_ds = []
        flank_ds_idx = []
    
        idx_0 = []
        for idx_flank in tqdm(range(len(flank))):
        #for idx_flank in range(len(flank)): 
            r_flank[idx_flank] = np.array([xy_3 + np.dot(flank[idx_flank], np.vstack((xy_1,xy_2))) for xy_1, xy_2, xy_3  in zip(_3_e_2_1, _3_e_2_2, _3_r_z)])
            n_flank[idx_flank] = np.array([np.dot(normal[idx_flank], np.vstack((xy_1,xy_2))) for xy_1, xy_2  in zip(_3_e_2_1, _3_e_2_2)])
            
            v_r_flank = np.gradient(r_flank[idx_flank])[0]    
            
            sc_pr_r_n = np.array([np.dot(xy_1, xy_2) for xy_1, xy_2 in zip(v_r_flank, n_flank[idx_flank])])[0:-2]
            
            idx = np.where(sc_pr_r_n[:-1] * sc_pr_r_n[1:] < 0)[0] + 1
                        
            if idx_flank == 0:
                idx_0 = idx[1]
            else:
                if len(idx) == 1:
                    idx_0 = idx[0]
                else:
                    diff = [abs(idx-idx_0) for idx in idx]
                    idx_0 = idx[np.argmin(diff)]
      
            flank_ds.append(r_flank[idx_flank,idx_0]) 
            flank_ds_idx.append(idx_0)
                        
        flank_ds = np.array(flank_ds)
        self.ds.flank = flank_ds
        return flank_ds