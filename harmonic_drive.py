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
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s],
                     [s,  c]])

def ROTATE(vec_array, phi):      
     return np.array([np.dot(ROT(phi), xy) for xy in vec_array])

def MIRROR(vec_arr, axis = 'x'):
    if axis == 'x':
        return np.flip(np.array([np.array([x,-y]) for x,y in vec_arr]), axis = 0)
    if axis == 'y':
        return np.flip(np.array([np.array([-x,y]) for x,y in vec_arr]), axis = 0)

def aequidistant(profil, q):
    tangent = np.gradient(profil)[0]
    normal = np.array([np.array([y,-x])/np.linalg.norm(np.array([x,y])) for x,y in tangent])
    return profil + normal * q

def length_vec_array(vec_array):
    tangent = np.gradient(vec_array)[0]
    return np.sum( np.array( [np.linalg.norm(t) for t in tangent] ) )

def ellipsis(a, b, num_of_discretization = 1000):
    t = np.linspace(0,2*np.pi,num_of_discretization)
    
    return np.array([np.array([a * np.cos(t),
                               b * np.sin(t)]) for t in t])

def point_aequidistant_ellipsis(a, b, q, t):
    x = a * np.cos(t)
    y = b * np.sin(t)
    
    n = np.array([x/a**2,y/b**2])
    n = 1/np.linalg.norm(n) * n
    
    return np.array([x, y] + n * q)

def tangent_aequidistant_ellipsis(a, b, q, t):
    x = a * np.cos(t)
    y = b * np.sin(t)
    
    n = np.array([x/a**2,y/b**2])
    n = 1/np.linalg.norm(n) * n
    
    return np.array([-n[1], n[0]])

def gear_from_flank(flank, z):
    gear = flank
    mirrored = MIRROR(flank)

    for idx in range(1,int(z)):
        gear = np.append(gear, ROTATE(mirrored, 2*np.pi/z*idx), axis = 0)
        gear = np.append(gear, ROTATE(flank, 2*np.pi/z*idx), axis = 0)

    gear = np.append(gear, mirrored, axis = 0)

    return gear

def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)
    
    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection

class Flexspline():
    
    def __init__(self):
        self.__is_init = True
        
        self.parameter = {'z', 'd_i', 's_st', 'd_nf', 'd_f', 'd', 'd_h', 'alpha', 'c', 'r_fh', 'r_ff', 'r_hr', 'r_fr', 'h_t'}
        
        #general:
        self.z      = 20    #number of teeth
        self.d_i    = 75     #Inner diameter
        self.s_st   = 2.5    #Sprocket Thickness
        self.d_nf   = []    #Neutral Fiber Diameter
        self.u_nf   = []      #circumference neutral fiber
        self.d_f    = []    #Foot circle diameter
        self.d      = 81    #Pitch circle diameter
        self.d_h    = 82    #Head circle diameter

        #tooth:
        self.alpha  = 20*np.pi/180 #Tooth flank angle
        self.c      = 1     #Tooth space ratio
        self.r_fh   = 5     #Head flank radius
        self.r_ff   = 5     #Foot flank radius
        self.r_hr   = 0.3   #Head rounding
        self.r_fr   = 0.3   #Foot rounding
        self.h_t    = []    #Tooth Height

        self.flank = []
        self.tooth = []
    
        self.r_z_i         = []
        self.CS2_i         = []
        self.inside_fibre  = []
    
        self.__is_init = False
        
        self.update_parameter()

    def __setattr__(self, name, value):
        #object.__setattr__(self, name, value)
        if name == '_Flexspline__is_init':
            self.__dict__[name] = value
        if self.__is_init:
            self.__dict__[name] = value   
        else:
            if name in self.parameter:
                self.__dict__[name] = value
                self.update_parameter()
            else:
                self.__dict__[name] = value
                
    def update_parameter(self):        
        #object.__setattr__(self, name, value)
        self.__dict__['d_nf']   = self.d_i + self.s_st
        self.__dict__['d_f']    = self.d_i + self.s_st * 2 
        self.__dict__['u_nf']   = self.d_nf * np.pi
        self.__dict__['h_t']    = 0.5*(self.d_h-self.d_f)
        
        self.flank              = self.calculate_flank(num_of_discretization = 50)
        self.tooth              = self.calculate_tooth()
        
    def calculate_flank(self, num_of_discretization = 1000):
        
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
        flank = self.calculate_flank(num_of_discretization = 100)
        flank = np.array([np.array([-y,x-self.d_nf/2]) for x,y in flank])
        
        tangent = np.gradient(flank)[0]
        
        normal = np.array([[y,-x]for x,y in tangent])
        normal_flank = np.array([1/np.linalg.norm(n)*n for n in normal])

        return flank, normal_flank
    
    def calculate_tooth(self):
        flank = self.flank
        tooth_flank_mirrored = MIRROR(flank)
        tooth = np.append(tooth_flank_mirrored, flank, axis = 0)
                
        tooth = np.array([np.array([-y,x-self.d_nf/2]) for x,y in tooth])
        
        return tooth
    
    def gear(self):
        return gear_from_flank(self.flank, self.z)
    
    def polygon(self):
        gear = self.gear()
        gear_polygon = sp.Polygon(gear).buffer(0)
        
        ALPHA = np.linspace(0,2*np.pi,int(len(gear)/10))
        inner_circle_polygon = sp.Polygon([[self.d_i/2*np.cos(alpha), self.d_i/2*np.sin(alpha)] for alpha in ALPHA])
        
        return gear_polygon.difference(inner_circle_polygon)
    
    def gear_deformed(self):
        gear_deformed = self.r_z_i[0] + np.array([np.dot(self.CS2_i[0], xy) for xy in self.tooth])
        for i in range(1,int(self.z)):
            tooth_i = np.array([self.r_z_i[i] + np.dot(self.CS2_i[i], xy) for xy in self.tooth])
            gear_deformed = np.append(gear_deformed, tooth_i, axis = 0)
        
        return gear_deformed
    
    def polygon_deformed(self):
        gear_deformed = self.gear_deformed()
        polygon_gear_deformed = sp.Polygon(gear_deformed).buffer(0)
        
        inside_fibre_polygon = sp.Polygon(self.inside_fibre)
        
        return polygon_gear_deformed.difference(inside_fibre_polygon)
        
class Wavegenerator():
    def __init__(self, shape = 1):
        self.__dict__['a'] = 1 #Large elliptical semi-axis
        self.__dict__['b'] = 0.9 #Small elliptical semi-axis
        
        if shape == 1:
            self.__dict__['shape'] = 'elliptical'
        elif shape == 2:
            self.__dict__['shape'] = '345polynomial'
            self.__dict__['arc'] = np.pi/10
        
        self.profil = []
        
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        self.update_parameter()

    def update_parameter(self):
        self.__dict__['profil'] = self.calculate_profil()
    
    def calculate_profil(self, num_of_discretization = 1000):
        if self.shape == 'elliptical':
            return ellipsis(self.a, self.b, num_of_discretization = num_of_discretization)
        
        if self.shape == '345polynomial':
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
        return ROTATE(self.profil, phi)
    
    def polygon(self, phi = 0):
        ALPHA = np.linspace(0,2*np.pi,100)
        m = np.dot(ROT(phi), self.profil[0]) * 0.95
        index_circle_polygon = sp.Polygon([[m[0] + np.cos(alpha), m[1] + np.sin(alpha)] for alpha in ALPHA])
        
        return sp.Polygon(ROTATE(self.profil, phi)).buffer(0).difference(index_circle_polygon)

class CircularSpline():
    def __init__(self):
        self.z      = []
        self.flank  = []
                
    def gear(self):
        #flank = np.array([xy for xy in self.flank if np.arctan(xy[1]/xy[0]) < np.pi/self.z*2 ])
        #return gear_from_flank(flank, int(self.z))
        return gear_from_flank(self.flank, int(self.z))
        
    def polygon(self):
        gear = self.gear()
        gear_polygon = sp.Polygon(gear).buffer(0)
        
        ALPHA = np.linspace(0,2*np.pi,int(len(gear)/100))
        d_o = (gear[0,0]*2) * 1.1
        outer_circle_polygon = sp.Polygon([[d_o/2*np.cos(alpha), d_o/2*np.sin(alpha)] for alpha in ALPHA])
        
        return outer_circle_polygon.difference(gear_polygon)
    
class DynamicSpline():
    def __init__(self):
        self.z      = []
        self.flank  = []
                
    def gear(self):
        flank = np.array([xy for xy in MIRROR(self.flank) if np.arctan(xy[1]/xy[0]) < np.pi/self.z ])
        return gear_from_flank(flank, self.z)
        
    def polygon(self, phi = 0):        
        gear = self.gear()
        gear_polygon = sp.Polygon(ROTATE(gear, phi)).buffer(0)
        
        ALPHA = np.linspace(0,2*np.pi,int(len(gear)/100))
        d_o = (gear[0,0]*2) * 1.05
        outer_circle_polygon = sp.Polygon([[d_o/2*np.cos(alpha), d_o/2*np.sin(alpha)] for alpha in ALPHA])
        
        m = np.dot(ROT(phi), gear[0]) * 1.05
        index_circle_polygon = sp.Polygon([[m[0] + np.cos(alpha), m[1] + np.sin(alpha)] for alpha in ALPHA])
        
        return outer_circle_polygon.difference(gear_polygon).difference(index_circle_polygon)

class HarmonicDrive():
    def __init__(self, flexspline_inst, wavegenerator_inst, circularspline_inst, dynamicspline_inst):
        
        #objects
        self.fs = flexspline_inst
        self.wg = wavegenerator_inst
        self.cs = circularspline_inst
        self.ds = dynamicspline_inst
        
        #parameter
        self.i = -20 #Transmission Ratio (i_wg->ds) | i
        self.d_br = 2 #Diameter Wavegenertor bearing Rolling element
        self.phi_wg = 0 #Drive Angle Wavegenerator | phi_wg
        
        self.q_nf = [] #Equidistant distance neutral fiber
        
        self.l_nf = []
        self.l_nf_t = []
        self.r_z = []
        self.phi_wg_r_z = []
        self.e_2_1 = []
        self.e_2_2 = []
        
        self.r_z_i = []
        self.e_2_1_i = []
        self.e_2_2_i = []
        
        self._3_r_z = []
        self._3_e_2_1 = []
        self._3_e_2_2 = []
        
        self.phi_ds = []
        
        self.circular_spline = []
        self.dynamic_spline = []
        
        self.update_parameter()
        self.fs.update_parameter()

    def update_parameter(self):
        self.fs.z = 2*(-self.i)
        self.cs.z = self.fs.z + 2
        self.ds.z = self.fs.z
        
        self.q_nf = self.d_br + self.fs.s_st/2
        
        self.wg.a = self.fs.u_nf * (-self.i + 1)/-self.i * 1/(2*np.pi) - self.q_nf
        self.wg.b = optimize.newton(func = self.optimize_u_nf__b, x0 = self.wg.a/2, args=(self.wg.a, self.q_nf, self.fs.u_nf),)
            
    def optimize_u_nf__b(self, b, *args):
        (a, q_nf, u_nf) = args
        lenght_neutral_fibre = self.wg.calculate_profil_length(a,b) + 2*np.pi*q_nf
        return lenght_neutral_fibre - u_nf
    
    def neutral_fibre(self, num_of_discretization = 1000):
        return self.wg.aequidistant(self.q_nf, num_of_discretization = num_of_discretization)
    
    def neutral_fibre_rotated(self, phi):
        neutral_fibre = self.neutral_fibre()
        return ROTATE(neutral_fibre, phi)
    
    def inside_fibre_rotated(self, phi, num_of_discretization = 1000):
        q_nf = self.q_nf
        s_st = self.fs.s_st
        inside_fibre = self.wg.aequidistant(q_nf-s_st/2, num_of_discretization = num_of_discretization)
        return ROTATE(inside_fibre, phi)
    
    def kincematics(self, num_of_discretization = 100000, phi_wg_lim  = 2*np.pi):
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
            
            tangent = self.wg.tangent_aequidistant(q_nf, t) #TODO
            
            tangent = np.array([np.interp(t,l_nf_t,gradient[:,0]),
                                np.interp(t,l_nf_t,gradient[:,1])])
            
            tangent = 1/np.linalg.norm(tangent) * tangent
            
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
        flank, normal = self.fs.calculate_flank_normal(num_of_discretization=10)  
    
        kinematics =  self.kincematics(num_of_discretization=1000, phi_wg_lim = np.pi/2)
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
        flank, normal = self.fs.calculate_flank_normal(num_of_discretization=10)  
        
        flank = MIRROR(flank, axis = 'y')
        normal = MIRROR(normal, axis = 'y')
            
        kinematics =  self.kincematics(num_of_discretization=1000, phi_wg_lim = np.pi/2)
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
    
    
    
#%% DEBUG
'''
def präsi():
    hd.i = -20
    
    fs.d_i = 75
    fs.s_st = 2.5
    fs.d_nf = 77.5
    fs.d_f = 80
    fs.d = 83
    fs.d_h = 86
    
    fs.c = 1
    fs.r_fh = 8
    fs.r_ff = 8
    fs.r_hr = 1
    fs.r_fr = 1
    
    hd.d_br = 0.3
    
fs = Flexspline()
wg = Wavegenerator(shape = 2)
wg.arc = np.pi/180 * 10
cs = CircularSpline()
ds = DynamicSpline()

hd = HarmonicDrive(fs, wg, cs, ds)

präsi()

hd.update_parameter()
hd.calculate_kinematics(num_of_discretization=10000)
#hd.update_parameter()

def pl_a(array, **kwargs):
    ax.plot(array[:,0], array[:,1], **kwargs)

plt.close('all')
fig = plt.gcf()
ax = plt.gca()


# b = np.linspace(10, 50, 100)
# res = np.zeros(len(b))

# for i in range(len(b)):
#     res[i] = hd.optimize_u_nf__b(b[i], wg.a, hd.q_nf, fs.u_nf)
#ax.plot(b,res)


hd.calculate_CS2_i()

patch_fs = plot_polygon(ax, fs.polygon_deformed(),
                        ec = "black",
                        lw = 1,
                        fc = "#e30066")

pl_a(hd.r_z)
pl_a(hd.neutral_fibre())



ax.axis('equal')
ax.grid(which='major', color='#DDDDDD', linewidth=1)
ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
ax.minorticks_on()

fig.tight_layout()
'''
