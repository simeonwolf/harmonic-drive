import numpy as np

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.animation as animation
import matplotlib.patches as mpl_patches

import tkinter as tk
from tkinter import ttk

from harmonic_drive import *

def create_entry_frame(parent_frame, label_text, default_value, row):
    ttk.Label(parent_frame, text=label_text).grid(row=row, column=0, sticky=tk.W)
    entry = ttk.Entry(parent_frame)
    entry.grid(row=row, column=1, sticky=tk.E)
    entry.insert(0, default_value)
    entry.bind("<Return>", (lambda event: update_plot()))
    return entry

def plot_settings(ax):
    ax.axis('equal')
    ax.grid(which='major', color='#DDDDDD', linewidth=1)
    ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    
    ax.set_axis_off()
    

def plot_circle(ax, r, center = [0,0], **kwargs):
    ALPHA = np.linspace(0,2*np.pi,1000)
    ax.plot(center[0] + r*np.cos(ALPHA), center[1] + r*np.sin(ALPHA), **kwargs)
    
def plot_coordinate_system(ax, length, phi, x_o = 0, y_o = 0, **kwargs):
    CS = np.array([[np.cos(phi), -np.sin(phi)],
                   [np.sin(phi),  np.cos(phi)]])
    e_1 = np.dot(CS,[1,0])*length
    e_2 = np.dot(CS,[0,1])*length
    
    ax.arrow(x_o ,y_o, e_1[0], e_1[1], lw=0.5, head_length=2, head_width=0.8, length_includes_head=True, zorder = 2, edgecolor = 'black', **kwargs)
    ax.arrow(x_o, y_o, e_2[0], e_2[1], lw=0.5, head_length=2, head_width=0.8, length_includes_head=True, zorder = 2, edgecolor = 'black', **kwargs)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def close():
   root.destroy()
   root.quit()

def update_checkboxes(current_var):
    if current_var.get() == 1:
        if current_var == check_var_tooth:
            check_var_gear.set(0)
            check_var_harmonic_drive.set(0)
            check_var_array_flexspline_tooth.set(0)
            check_var_3_array_flexspline_tooth.set(0)
        elif current_var == check_var_gear:
            check_var_tooth.set(0)
            check_var_harmonic_drive.set(0)
            check_var_array_flexspline_tooth.set(0)
            check_var_3_array_flexspline_tooth.set(0)
        elif current_var == check_var_harmonic_drive:
            check_var_tooth.set(0)
            check_var_gear.set(0)
            check_var_array_flexspline_tooth.set(0)
            check_var_3_array_flexspline_tooth.set(0)
        elif current_var == check_var_array_flexspline_tooth:
            check_var_tooth.set(0)
            check_var_gear.set(0)
            check_var_harmonic_drive.set(0)
            check_var_3_array_flexspline_tooth.set(0)
        elif current_var == check_var_3_array_flexspline_tooth:
            check_var_tooth.set(0)
            check_var_gear.set(0)
            check_var_harmonic_drive.set(0)
            check_var_array_flexspline_tooth.set(0)
            
def update_checkboxes_wg(current_var):
    if current_var.get() == 1:
        if current_var == check_var_elliptical:
            check_var_345polynomial.set(0)
            entry_arc.configure(state="readonly")
        elif current_var == check_var_345polynomial:
            check_var_elliptical.set(0)
            entry_arc.configure(state="normal")

def update_plot(calculate = True):
    print("update...")
    
    #wg:
    if check_var_elliptical.get() == 1:
        wg.shape = 'elliptical'
    elif check_var_345polynomial.get() == 1:
        wg.shape = '345polynomial'

    wg.arc = float(entry_arc.get())*np.pi/180
    
    #fs
    fs.__dict__['d_i']    = float(entry_d_i.get())
    fs.__dict__['s_st']   = float(entry_s_st.get())
    fs.__dict__['d']      = float(entry_d.get())
    fs.__dict__['d_h']    = float(entry_d_h.get())
    
    fs.__dict__['alpha']  = float(entry_alpha.get())*np.pi/180
    fs.__dict__['c']      = float(entry_c.get())
    fs.__dict__['r_fh']   = float(entry_r_fh.get())
    fs.__dict__['r_ff']   = float(entry_r_ff.get())
    fs.__dict__['r_hr']   = float(entry_r_hr.get())
    fs.__dict__['r_fr']   = float(entry_r_fr.get())
    fs.update()

    hd.__dict__['i']      = -float(entry_i.get())
    
    br.d_br               = float(entry_d_br.get())
    br.n                  = float(entry_n.get())
    
    hd.phi_wg             = -float(entry_phi_wg.get())*np.pi/180    
    hd.update()
    
    #update entry values
    update_entrys()
    
    #update selected plot
    ax.clear()
    
    if check_var_tooth.get() == 1:
        plot_view_tooth()
    elif check_var_gear.get() == 1:
        plot_view_fs_gear()
    elif check_var_harmonic_drive.get() == 1:
        plot_view_harmonic_drive(calculate = calculate)
    elif check_var_array_flexspline_tooth.get() == 1:
        plot_view_array_flexspline_tooth()
    elif check_var_3_array_flexspline_tooth.get() == 1:
        plot_view_3_array_flexspline_tooth()

    canvas.draw()
    print("...ready")
    
def update_entrys():
    update_entry(entry_d_nf , fs.d_nf)
    update_entry(entry_d_f  , fs.d_f)
    update_entry(entry_z    , fs.z)
    
    update_entry(entry_a, wg.a)
    update_entry(entry_b, wg.b)
    
    update_entry(entry_z_fs , fs.z)
    update_entry(entry_z_cs , cs.z)
    update_entry(entry_z_ds , ds.z)
    update_entry(entry_q_nf , hd.q_nf)

def update_entry(entry, value):
    entry.configure(state="normal")
    entry.delete(0, "end")
    entry.insert(0, value)
    entry.configure(state="readonly")

def plot_view_tooth():
    fs_flank = fs.flank
    ax.plot(fs_flank[:,0], fs_flank[:,1], c = "darkblue", zorder = 3, lw = 2)
    
    ALPHA = np.linspace(0,np.pi/fs.z,1000)
    ax.plot(fs.d_i/2*np.cos(ALPHA), fs.d_i/2*np.sin(ALPHA), lw = 2, c = "darkblue")
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    gear = fs.gear()
    ax.plot(gear[:,0], gear[:,1], lw = 1, c = "darkblue")
    
    #Circles 
    plot_circle(ax, fs.d_i/2, lw = 1, c = "darkblue")
    plot_circle(ax, fs.d_nf/2, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d_f/2, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d/2, ls = "dashdot", lw = 1, c = "black")
    plot_circle(ax, fs.d_h/2, ls = "dashed", lw = 1, c = "black")
    
    for i in range(2, 6):
        plot_circle(ax, getattr(fs, f'r_fl_{i}'), center = getattr(fs, f'm_{i}'), ls = "dashed", lw = 1, c = "lightcoral")
    
    fa_1 = fs.s_34 - (fs.d_h-fs.d_i)/2 * np.array([np.cos(fs.alpha), -np.sin(fs.alpha)])
    fa_2 = fs.s_34 + (fs.d_h-fs.d_i)/2 * np.array([np.cos(fs.alpha), -np.sin(fs.alpha)])
    ax.plot([fa_1[0],fa_2[0]],[fa_1[1],fa_2[1]], ls = "dashed", lw = 2, c = "lightcoral")

    ax.plot([fs.d_i/2-fs.d_i/100,fs.d_h/2+fs.d_h/100],[0,0], ls = "dashdot", lw = 1, c = "k")

    plot_settings(ax)
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def plot_view_fs_gear():
    plot_polygon(ax, fs.polygon(), 
                      ec = "black",
                      lw = 1,
                      fc = "lightblue")
    
    plot_circle(ax, fs.d_nf/2, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d_f/2, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d/2, ls = "dashdot", lw = 1, c = "black")
    plot_circle(ax, fs.d_h/2, ls = "dashed", lw = 1, c = "black")
    
    plot_settings(ax)

def plot_view_harmonic_drive(calculate = True):
    if calculate:
        hd.calculate_kinematics()
        hd.calculate_flank_cs()
        hd.calculate_flank_ds()
        hd.calculate_bearing_kinematics()
    
    phi_wg = hd.phi_wg

    update_rotation(phi_wg)
            
    R = (wg.a + hd.q_nf + fs.s_st/2 + fs.h_t)*1.11
    plot_settings(ax) 
    ax.set_xlim([-R,R])
    ax.set_ylim([-R,R])
    ax.set_aspect('equal', 'box')

def update_rotation(phi_wg, plot_circular_spline = True):
    
    hd.calculate_CS2_i()
    hd.calculate_bearing()

    #---Wavegenerator---
    color = 'deepskyblue'
    plot_coordinate_system(ax, length=wg.b/2, phi=phi_wg, facecolor = color) #CS1
    #profil_wg_phi_wg = wg.profil_rotated(phi_wg)
    #ax.plot(profil_wg_phi_wg[:,0], profil_wg_wg[:,1])
    #ax.plot([0,np.cos(phi_wg)*a_el],[0,np.sin(phi_wg)*a_el], ls = 'dashed',lw = 1, c = "darkviolet" , alpha = 0.5)
    try:
        plot_polygon(ax, wg.polygon(phi=phi_wg),
                            ec = "black",
                            lw = 1,
                            fc = (color, 0.8))
    except:
        print('Plot Error Wavegenerator')

    #---Bearing---
    color = 'grey'
    try:
        [plot_polygon(ax, polygon,
                        ec = "black",
                        lw = 1,
                        fc = (color, 0.8)) for polygon in br.polygons()]
    except:
        print('Plot Error Bearing')

    #---Flexspline---
    color = 'springgreen'
    #neutral_fibre_phi_wg = hd.neutral_fibre_rotated(phi_wg)
    #ax.plot(neutral_fibre_phi_wg[:,0], neutral_fibre_phi_wg[:,1], ls = "-", lw = 1, c = "black")
    try:
        plot_polygon(ax, fs.polygon_deformed(),
                            ec = "black",
                            lw = 1,
                            fc = (color, 0.8))
    except:
        print('Plot Error Flexspline')

    #---Flexspline first tooth---
    r_z = hd.r_z
    r_z_i = hd.r_z_i
    phi_wg_r_z = hd.phi_wg_r_z
    e_2_1 = hd.e_2_1
    e_2_2 = hd.e_2_2
    
    e_2_1_phi_wg = np.array([np.interp(-phi_wg,phi_wg_r_z,e_2_1[:,0]),
                             np.interp(-phi_wg,phi_wg_r_z,e_2_1[:,1])])
    
    e_2_2_phi_wg = np.array([np.interp(-phi_wg,phi_wg_r_z,e_2_2[:,0]),
                             np.interp(-phi_wg,phi_wg_r_z,e_2_2[:,1])])
    
    #trace r_z:
    ax.plot(r_z[:,0], r_z[:,1], '--', lw = 1, c = "r")
    ax.arrow(0,0,r_z_i[0,0],r_z_i[0,1], color = 'red', lw=1, head_length=2, head_width=0.8, length_includes_head=True)
    
    cs_2_phi_wg = np.array([[e_2_1_phi_wg[0], e_2_2_phi_wg[0]],
                            [e_2_1_phi_wg[1], e_2_2_phi_wg[1]]])
    
    tooth_phi_wg = r_z_i[0] + np.array([np.dot(cs_2_phi_wg, xy) for xy in fs.tooth])
    ax.plot(tooth_phi_wg[:,0], tooth_phi_wg[:,1], lw = 2, c = "darkblue")
    
    #coordinate system:
    ax.arrow(r_z_i[0,0],r_z_i[0,1],e_2_1_phi_wg[0]*4,e_2_1_phi_wg[1]*4, facecolor = color, lw=0.5, head_length=1, head_width=0.4, length_includes_head=True, edgecolor = 'black')
    ax.arrow(r_z_i[0,0],r_z_i[0,1],e_2_2_phi_wg[0]*4,e_2_2_phi_wg[1]*4, facecolor = color, lw=0.5, head_length=1, head_width=0.4, length_includes_head=True, edgecolor = 'black')
    
    #---Circular Spline---
    color = 'darkgrey'
    plot_coordinate_system(ax, length=wg.b/2, phi=0, facecolor = color) #CS0
    try:
        plot_polygon(ax, cs.polygon(),
                          ec = "black",
                          lw = 1,
                          fc = (color, 0.8))
    except:
        print('Plot Error Circular Spline')
            
    #---Dynamic Spline---
    color = 'fuchsia'
    plot_coordinate_system(ax, length=wg.b/2, phi=hd.phi_ds, facecolor = color) #CS3
    try:
        plot_polygon(ax, ds.polygon(phi = hd.phi_ds),
                            ec = "black",
                            lw = 1,
                            fc = (color, 0.5))
    except:
        print('Plot Error Dynamic Spline')
    
       
def plot_view_array_flexspline_tooth():
    tooth = fs.tooth
    
    flank_cs = hd.calculate_flank_cs()       
    
    kinematics =  hd.kincematics(num_of_discretization=10000, phi_wg_lim = np.pi)
    (r_z, e_2_1, e_2_2) = kinematics[:3]
    
    ax.plot(r_z[:,0], r_z[:,1], '--', lw = 1, c = "r")
    
    for i in range(0,len(r_z), int(len(r_z)/100)):
        CS_2_i = np.array([[e_2_1[i,0],e_2_2[i,0]],
                           [e_2_1[i,1],e_2_2[i,1]]])
        tooth_i = np.array([ r_z[i] + np.dot(CS_2_i, xy) for xy in tooth])
        ax.plot(tooth_i[:,0], tooth_i[:,1], lw = 1, c = "darkblue")
    
    ax.plot(flank_cs[:,0], flank_cs[:,1], '+r')
    
    plot_settings(ax)
    
def plot_view_3_array_flexspline_tooth():
    fs.calculate_tooth()
    tooth = fs.tooth
    
    flank_ds = hd.calculate_flank_ds()
    
    kinematics =  hd.kincematics(num_of_discretization=10000, phi_wg_lim = np.pi)
    (_3_r_z, _3_e_2_1, _3_e_2_2) = kinematics[3:6]
    
    ax.plot(_3_r_z[:,0], _3_r_z[:,1], '--', lw = 1, c = "r")
    
    for i in range(0,len(_3_r_z), int(len(_3_r_z)/100)):
        _3_CS_2_i = np.array([[_3_e_2_1[i,0],_3_e_2_2[i,0]],
                              [_3_e_2_1[i,1],_3_e_2_2[i,1]]])
        _3_tooth_i = np.array([ _3_r_z[i] + np.dot(_3_CS_2_i, xy) for xy in tooth])
        ax.plot(_3_tooth_i[:,0], _3_tooth_i[:,1], lw = 1, c = "darkblue")
        
    ax.plot(flank_ds[:,0], flank_ds[:,1], '+r')
        
    plot_settings(ax)

def update_animation(frame):
    print("frame", frame)

    for artist in ax.get_children():
        if isinstance(artist, mpl_patches.FancyArrow) or isinstance(artist, plt.Line2D) or isinstance(artist, PatchCollection):
            artist.remove()
    
    phi_wg = -frame/2*np.pi/180
    hd.phi_wg = phi_wg
    update_rotation(phi_wg)
    entry_phi_wg.delete(0, "end")
    entry_phi_wg.insert(0, int(-phi_wg*180/np.pi))

    return ax.get_children()

def start_stop_animation():
    global animation_running, animation_instance
    
    if animation_running == False:
        animation_running = True
        if not animation_instance:
            
            animation_instance = animation.FuncAnimation(fig=fig, func=update_animation,
                                                              frames = 360*2, interval=5, repeat=True)#, blit=True)
            root.update()
            canvas.draw()
        
            if check_var_safe_animation.get() == 1:
                print("Calculate Animation")
                writervideo = animation.FFMpegWriter(fps=30, bitrate=30000)#, metadata=dict(artist='Me'), bitrate=1800, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
                animation_instance.save('animation.mp4', writer=writervideo) 
                print("Save Animation")
        else:
            animation_instance.resume()
    else:
        animation_running = False
        if animation_instance:
            animation_instance.pause()
            
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def initial():
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
    
    br.d_br = 5
    br.n = 20

#objects
fs = Flexspline()
wg = Wavegenerator(shape = 2)
cs = CircularSpline()
ds = DynamicSpline()
br = Bearing()
hd = HarmonicDrive(fs, wg, cs, ds, br)

#change inital paramter
initial()

#root window
root = tk.Tk()
root.iconbitmap("icon.ico")
root.title("Harmonic Drive")

#figure and axes for Matplotlib
fig = plt.Figure(figsize=(10, 10), dpi=100)
ax = fig.add_subplot()

#animation
animation_running = False
animation_instance = []
update_list = []

#FigureCanvasTkAgg to integrate the Matplotlib graphic into Tkinter
frame_plot = ttk.Frame(root)
frame_plot.pack(side=tk.LEFT, fill = 'y')

canvas = FigureCanvasTkAgg(fig, frame_plot)
canvas_widget = canvas.get_tk_widget()

toolbar = NavigationToolbar2Tk(canvas, frame_plot)
toolbar.update()

canvas_widget.pack()
toolbar.pack()

# Entrys
# frame_flexspline
frame_flexspline = ttk.Frame(root)
frame_flexspline.pack(side=tk.TOP, fill='x', padx=10, pady=2)
ttk.Label(frame_flexspline, text="Flexspline Parameter").grid(row=0, column=0, columnspan=2)

entry_d_i = create_entry_frame(frame_flexspline, "Inner Diameter | d_i", fs.d_i, 1)
entry_s_st = create_entry_frame(frame_flexspline, "Sprocket Thickness | s_st", fs.s_st, 2)
entry_d_nf = create_entry_frame(frame_flexspline, "Neutral Fiber Diameter | d_nf", fs.d_nf, 3)
entry_d_nf.configure(state="readonly")
entry_d_f = create_entry_frame(frame_flexspline, "Foot Circle Diameter | d_f", fs.d_f, 4)
entry_d_f.configure(state="readonly")
entry_d = create_entry_frame(frame_flexspline, "Pitch Circle Diameter | d", fs.d, 5)
entry_d_h = create_entry_frame(frame_flexspline, "Head Circle Diameter | d_h", fs.d_h, 6)
entry_z = create_entry_frame(frame_flexspline, "Number of Teeth | z", fs.z, 7)
entry_z.configure(state="readonly")

# frame_tooth
frame_tooth = ttk.Frame(root)
frame_tooth.pack(side=tk.TOP, fill='x', padx=10, pady=2)
ttk.Label(frame_tooth, text="Flexspline Tooth Parameter").grid(row=0, column=0, columnspan=2)

entry_alpha = create_entry_frame(frame_tooth, "Tooth Flank Angle | alpha", fs.alpha * 180/np.pi, 1)
entry_c = create_entry_frame(frame_tooth, "Tooth Space Ratio | c", fs.c, 2)
entry_r_fh = create_entry_frame(frame_tooth, "Head Flank Radius | r_fh", fs.r_fh, 3)
entry_r_ff = create_entry_frame(frame_tooth, "Foot Flank Radius | r_ff", fs.r_ff, 4)
entry_r_hr = create_entry_frame(frame_tooth, "Head Rounding | r_hr", fs.r_hr, 5)
entry_r_fr = create_entry_frame(frame_tooth, "Foot Rounding | r_fr", fs.r_fr, 6)

# frame_wavegenerator
frame_wavegenerator = ttk.Frame(root)
frame_wavegenerator.pack(side=tk.TOP, fill='x', padx=10, pady=2)
ttk.Label(frame_wavegenerator, text="Wavegenerator Parameter").grid(row=0, column=0, columnspan=2)

entry_a = create_entry_frame(frame_wavegenerator, "Large Semi-Axis | a", wg.a, 1)
entry_a.configure(state="readonly")
entry_b = create_entry_frame(frame_wavegenerator, "Small Semi-Axis | b", wg.b, 2)
entry_b.configure(state="readonly")
entry_arc = create_entry_frame(frame_wavegenerator, "Arc | gamma", wg.arc * 180/np.pi, 3)

check_var_elliptical = tk.BooleanVar(root)
checkbox_elliptical = ttk.Checkbutton(root, text="Elliptical", variable=check_var_elliptical,
                                      command=lambda: update_checkboxes_wg(check_var_elliptical))
checkbox_elliptical.pack(padx=10, pady=2)

check_var_345polynomial = tk.BooleanVar(root)
checkbox_345polynomial = ttk.Checkbutton(root, text="345polynomial", variable=check_var_345polynomial,
                                      command=lambda: update_checkboxes_wg(check_var_345polynomial))
checkbox_345polynomial.pack(padx=10, pady=2)

if wg.shape == 'elliptical':
    check_var_345polynomial.set(0)
    check_var_elliptical.set(1)
    entry_arc.configure(state="readonly")
elif wg.shape == '345polynomial':
    check_var_345polynomial.set(1)
    check_var_elliptical.set(0)
    entry_arc.configure(state="normal")

# frame_harmonic_drive
frame_harmonic_drive = ttk.Frame(root)
frame_harmonic_drive.pack(side=tk.TOP, fill='x', padx=10, pady=2)
ttk.Label(frame_harmonic_drive, text="Harmonic Drive Parameter").grid(row=0, column=0, columnspan=2)

entry_i = create_entry_frame(frame_harmonic_drive, "Transmission Ratio (|i_wg->ds|) | i", -hd.i, 1)
entry_z_fs = create_entry_frame(frame_harmonic_drive, "Number of teeth Flexspline | z_fs", fs.z, 2)
entry_z_fs.configure(state="readonly")
entry_z_cs = create_entry_frame(frame_harmonic_drive, "Number of teeth Circular Spline | z_cs", cs.z, 3)
entry_z_cs.configure(state="readonly")
entry_z_ds = create_entry_frame(frame_harmonic_drive, "Number of teeth Dynamic Spline | z_ds", ds.z, 4)
entry_z_ds.configure(state="readonly")
entry_q_nf = create_entry_frame(frame_harmonic_drive, "Equidistant Distance Neutral Fiber | q_nf", hd.q_nf, 5)
entry_q_nf.configure(state="readonly")

# frame_bearing
frame_bearing = ttk.Frame(root)
frame_bearing.pack(side=tk.TOP, fill='x', padx=10, pady=2)
ttk.Label(frame_bearing, text="Bearing Parameter").grid(row=0, column=0, columnspan=2)

entry_d_br = create_entry_frame(frame_bearing, "Diameter Wavegenertor Bearing Rolling element | d_br", br.d_br, 1)
entry_n    = create_entry_frame(frame_bearing, "Number of bearing rollers | n", br.n, 2)

#CheckBoxes
check_var_tooth = tk.BooleanVar(root)
checkbox_tooth = ttk.Checkbutton(root, text="View Flexspline Tooth", variable=check_var_tooth,
                                      command=lambda: update_checkboxes(check_var_tooth))
checkbox_tooth.pack(padx=10, pady=2)

check_var_gear = tk.BooleanVar(root)
checkbox_gear = ttk.Checkbutton(root, text="View Flexspline", variable=check_var_gear,
                                     command=lambda: update_checkboxes(check_var_gear))
checkbox_gear.pack(padx=10, pady=2)

check_var_harmonic_drive = tk.BooleanVar(root)
checkbox_harmonic_drive = ttk.Checkbutton(root, text="View Harmonic Drive", variable=check_var_harmonic_drive,
                                     command=lambda: update_checkboxes(check_var_harmonic_drive))
checkbox_harmonic_drive.pack(padx=10, pady=2)

#Entry Plot Harmonic Drive
frame_plot_h_d = ttk.Frame(root)
frame_plot_h_d.pack(side=tk.TOP, fill='x', padx=10, pady=2)
ttk.Label(frame_plot_h_d, text="Plot Harmonic Drive Parameter").grid(row=0, column=0, columnspan=2)

entry_phi_wg = create_entry_frame(frame_plot_h_d, "Drive Angle Wavegenerator | phi_wg", hd.phi_wg * 180/np.pi, 1)
entry_phi_wg.bind("<Return>", (lambda event: update_plot(calculate = False)))

# Button for Animation
button_start_stop = ttk.Button(root, text="Start/Stop", command=start_stop_animation)
button_start_stop.pack(padx=10, pady=2)

check_var_safe_animation = tk.BooleanVar(root)
checkbox_safe_animation = ttk.Checkbutton(root, text="Safe Animation", variable=check_var_safe_animation)
checkbox_safe_animation.pack(padx=10, pady=2)
        
#---

check_var_array_flexspline_tooth = tk.BooleanVar(root)
checkbox_array_flexspline_tooth = ttk.Checkbutton(root, text="Array of Flexspline Tooth", variable=check_var_array_flexspline_tooth,
                                     command=lambda: update_checkboxes(check_var_array_flexspline_tooth))
checkbox_array_flexspline_tooth.pack(padx=10, pady=2)

check_var_3_array_flexspline_tooth = tk.BooleanVar(root)
checkbox_3_array_flexspline_tooth = ttk.Checkbutton(root, text="Array of Flexspline Tooth (CS3)", variable=check_var_3_array_flexspline_tooth,
                                     command=lambda: update_checkboxes(check_var_3_array_flexspline_tooth))
checkbox_3_array_flexspline_tooth.pack(padx=10, pady=2)

#Buttons
button_close = ttk.Button(root, text="close", command=close)
button_close.pack(padx=10, pady=2)

button_update = ttk.Button(root, text="update", command=update_plot)
button_update.pack(padx=10, pady=2)

root.mainloop()