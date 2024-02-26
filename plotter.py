import numpy as np

import os, fnmatch

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.animation as animation
import matplotlib.patches as mpl_patches
import matplotlib.colors as mcolors

import tkinter as tk
import tkinter.ttk as ttk
from ttkthemes import ThemedTk

from harmonic_drive import *

def plot_settings(ax):
    ax.axis('equal')
    ax.grid(which='major', color='#DDDDDD', linewidth=1)
    ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    
    ax.set_axis_off()
          
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

def create_entry_frame(parent_frame, description, parameter, default_value, row, readonly = False):
    frame = ttk.Frame(master=parent_frame)
    frame.grid(row=row, column=0, sticky="w")
    ttk.Label(master=frame, text=description).pack()
    
    frame = ttk.Frame(master=parent_frame)
    frame.grid(row=row, column=1, sticky="e")
    ttk.Label(master=frame, text=parameter).pack()
    
    frame = ttk.Frame(master=parent_frame)
    frame.grid(row=row, column=2, sticky="w")
    ent = ttk.Entry(master=frame)
    ent.insert(0, default_value)
    ent.bind("<Return>", (lambda event: update_plot()))
    ent.pack()
    
    if readonly: ent.configure(state="readonly")
    
    return ent

def update_plot(calc = True):
    print("update...")
    
    #hd
    hd.i      = -float(entry_i.get())
    
    #fs
    fs.d_i    = float(parameter_fs['d_i']['widget'].get())
    fs.s_st   = float(parameter_fs['s_st']['widget'].get())
    fs.d      = float(parameter_fs['d']['widget'].get())
    fs.d_h    = float(parameter_fs['d_h']['widget'].get())
    
    fs.alpha  = float(parameter_fs['alpha']['widget'].get())*np.pi/180
    fs.c      = float(parameter_fs['c']['widget'].get())
    fs.r_fh   = float(parameter_fs['r_fh']['widget'].get())
    fs.r_ff   = float(parameter_fs['r_ff']['widget'].get())
    fs.r_hr   = float(parameter_fs['r_hr']['widget'].get())
    fs.r_fr   = float(parameter_fs['r_fr']['widget'].get())
    
    #br
    br.d_br   = float(entry_d_br.get())
    br.n      = float(entry_n.get())
    
    #wg: 
    wg.a     = float(parameter_wg['a']['widget'].get())
    wg.b     = float(parameter_wg['b']['widget'].get())
    wg.arc   = float(parameter_wg['arc']['widget'].get())*np.pi/180    
    wg.shape = shape_var.get()            
                         
    hd.update()
    update_entrys()
    
    #update selected plot
    ax.clear()
    
    if view_var.get() == views[0]:
        plot_view_tooth()
    elif view_var.get() == views[1]:
        plot_view_fs_gear()
    elif view_var.get() == views[2]:
        hd.calc_flexspline_kinematics()
        hd.calc_circular_spline_flank()
        plot_view_array_flexspline_tooth()
    elif view_var.get() == views[3]:
        hd.calc_flexspline_kinematics()
        hd.calc_dynamic_spline_flank()
        plot_view_3_array_flexspline_tooth()
    elif view_var.get() == views[4]:
        if calc: hd.calc()
        hd.phi_wg = -float(ent_phi_wg.get())*np.pi/180
        plot_view_harmonic_drive(calc = calc)
    
    canvas.draw()
    print("...ready")
    
def update_entrys():
    for p in parameter_fs.values():
        if p['abbr'] == 'alpha':
            update_entry(p['widget'], fs.__dict__[p['abbr']]*180/np.pi)
        else:
            update_entry(p['widget'], fs.__dict__[p['abbr']])
    
    for p in parameter_wg.values():
        if p['abbr'] == 'arc':
            update_entry(p['widget'], wg.__dict__[p['abbr']]*180/np.pi)
        else:
            update_entry(p['widget'], wg.__dict__[p['abbr']])
    
    
    update_entry(entry_i, -hd.i)
    update_entry(entry_z_fs, fs.z)
    update_entry(entry_z_cs, cs.z)
    update_entry(entry_z_ds, ds.z)
    update_entry(entry_q_nf, hd.q_nf)

    update_entry(entry_d_br, br.d_br)
    update_entry(entry_n, br.n)

def update_entry(entry, value):
    if entry.cget("state") == 'readonly':
        entry.configure(state="normal")
        entry.delete(0, tk.END)
        entry.insert(0, value)
        entry.configure(state="readonly")
    else:
        entry.delete(0, tk.END)
        entry.insert(0, value)
        
def plot_view_tooth():
    fs_flank = fs.flank
    plot_va(ax, fs_flank, c = "darkblue", zorder = 3, lw = 2)
    
    ALPHA = np.linspace(0, np.pi/fs.z, 100)
    plot_va(ax, fs.d_i/2 * np.array([np.cos(ALPHA), np.sin(ALPHA)]).T, c = "darkblue")
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    gear = fs.gear_undeformed()
    plot_va(ax, fs.gear_undeformed(), lw = 1, c = "darkblue")
    
    #Circles
    plot_circle(ax, fs.d_i/2, num_of_discretization = 1000, lw = 1, c = "darkblue")
    plot_circle(ax, fs.d_nf/2, num_of_discretization = 1000, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d_f/2, num_of_discretization = 1000, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d/2, num_of_discretization = 1000, ls = "dashdot", lw = 1, c = "black")
    plot_circle(ax, fs.d_h/2, num_of_discretization = 1000, ls = "dashed", lw = 1, c = "black")
    
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
    plot_polygon(ax, fs.polygon_undeformed(), 
                      ec = "black",
                      lw = 1,
                      fc = "lightblue")
    
    plot_circle(ax, fs.d_nf/2, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d_f/2, ls = "dashed", lw = 1, c = "black")
    plot_circle(ax, fs.d/2, ls = "dashdot", lw = 1, c = "black")
    plot_circle(ax, fs.d_h/2, ls = "dashed", lw = 1, c = "black")
    
    plot_settings(ax)

def plot_view_harmonic_drive(calc = True):
    global patch_cs
    if not animation_running:
        #---Circular Spline---
        plot_coordinate_system(ax, length=wg.b/2, phi=0, facecolor = colors['cs']) #CS0
        try:
            patch_cs = plot_polygon(ax, cs.polygon(),
                                      ec = "black",
                                      lw = 1,
                                      fc = (colors['cs'], 0.8))
        except:
            print('Plot Error Circular Spline')
    
    update_rotation()
    
    R = (wg.a + hd.q_nf + fs.s_st/2 + fs.h_t)*1.11
    plot_settings(ax) 
    ax.set_xlim([-R,R])
    ax.set_ylim([-R,R])
    ax.set_aspect('equal', 'box')

def update_rotation():
    
    #---Wavegenerator---
    plot_coordinate_system(ax, length=wg.b/2, phi=hd.phi_wg, facecolor = colors['wg']) #CS1
    #profil_wg_phi_wg = wg.profil_rotated(phi_wg)
    #ax.plot(profil_wg_phi_wg[:,0], profil_wg_wg[:,1])
    #ax.plot([0,np.cos(phi_wg)*a_el],[0,np.sin(phi_wg)*a_el], ls = 'dashed',lw = 1, c = "darkviolet" , alpha = 0.5)
    try:
        plot_polygon(ax, wg.polygon(),
                            ec = "black",
                            lw = 1,
                            fc = (colors['wg'], 0.8))
    except Exception as error:
        print('Plot Error Wavegenerator' + error)

    #---Bearing---
    try:
        [plot_polygon(ax, polygon,
                        ec = "black",
                        lw = 1,
                        fc = (colors['br'], 0.8)) for polygon in br.polygons()]
    except Exception as error:
        print('Plot Error Bearing: ' + error)

    #---Flexspline---
    #neutral_fibre_phi_wg = hd.neutral_fibre_rotated(phi_wg)
    #ax.plot(neutral_fibre_phi_wg[:,0], neutral_fibre_phi_wg[:,1], ls = "-", lw = 1, c = "black")
    try:
        plot_polygon(ax, fs.polygon(),
                            ec = "black",
                            lw = 1,
                            fc = (colors['fs'], 0.8))
    except:
        print('Plot Error Flexspline')

    #---Flexspline first tooth---

    #trace r_z:
    r_z = hd.r_z
    ax.plot(r_z[:,0], r_z[:,1], '--', lw = 1, c = "r")
    ax.arrow(0,0,fs.r_z_i[0,0],fs.r_z_i[0,1], color = 'red', lw=1, head_length=2, head_width=0.8, length_includes_head=True)
    
    #TODO
    #TRANSFORM(fs.toooth, )
    tooth_phi_wg = fs.r_z_i[0] + np.array([np.dot(fs.CS_2_i[0], xy) for xy in fs.tooth])
    ax.plot(tooth_phi_wg[:,0], tooth_phi_wg[:,1], lw = 2, c = "darkblue")
    
    #coordinate system:
    [e_2_1_phi_wg, e_2_2_phi_wg] = np.hsplit(fs.CS_2_i[0], 2)
    e_2_1_phi_wg = np.squeeze(e_2_1_phi_wg)
    e_2_2_phi_wg = np.squeeze(e_2_2_phi_wg)
    
    ax.arrow(fs.r_z_i[0,0], fs.r_z_i[0,1], e_2_1_phi_wg[0]*4, e_2_1_phi_wg[1]*4, facecolor = colors['fs'], lw=0.5, head_length=1, head_width=0.4, length_includes_head=True, edgecolor = 'black')
    ax.arrow(fs.r_z_i[0,0], fs.r_z_i[0,1], e_2_2_phi_wg[0]*4, e_2_2_phi_wg[1]*4, facecolor = colors['fs'], lw=0.5, head_length=1, head_width=0.4, length_includes_head=True, edgecolor = 'black')
    
    #---Dynamic Spline---
    plot_coordinate_system(ax, length=wg.b/2, phi=hd.phi_ds, facecolor = colors['ds']) #CS3
    ds.polygon()
    try:
        plot_polygon(ax, ds.polygon(),
                            ec = "black",
                            lw = 1,
                            fc = (colors['ds'], 0.5))
    except:
        print('Plot Error Dynamic Spline')
       
def plot_view_array_flexspline_tooth():    
    plot_va(ax, hd.r_z, ls = 'dashed', lw = 1, c = 'r')        
    
    for i in range(0,len(hd.r_z), int(len(hd.r_z)/100)):
        tooth_i = TRANSFORM(fs.tooth, hd.CS_2[i],  hd.r_z[i])
        plot_va(ax, tooth_i, lw = 1, c = 'darkblue')
    
    plot_va(ax, cs.flank, marker = '+', color = 'red')

    plot_settings(ax)
    
def plot_view_3_array_flexspline_tooth():
    plot_va(ax, hd._3_r_z, ls = 'dashed', lw = 1, c = 'r')
    
    for i in range(0,len(hd._3_r_z), int(len(hd._3_r_z)/100)):
        _3_tooth_i = TRANSFORM(fs.tooth, hd._3_CS_2[i], hd._3_r_z[i])
        plot_va(ax, _3_tooth_i, lw = 1, c = 'darkblue')
    
    plot_va(ax, ds.flank, marker = '+', color = 'red')
        
    plot_settings(ax)

def update_animation(frame):

    for artist in ax.get_children():
        if not patch_cs == artist:
            if isinstance(artist, mpl_patches.FancyArrow) or isinstance(artist, plt.Line2D) or isinstance(artist, PatchCollection):
                artist.remove()
    
    phi_wg = -frame/4*np.pi/180
    hd.phi_wg = phi_wg
    update_rotation()
    ent_phi_wg.delete(0, tk.END)
    ent_phi_wg.insert(0, int(-phi_wg*180/np.pi))
    window.update()
    
    return ax.get_children()

def start_stop_animation():
    global animation_running, animation_instance
    
    if animation_running == False:
        animation_running = True
        if not animation_instance:
            
            animation_instance = animation.FuncAnimation(fig=fig, func=update_animation,
                                                              frames = 360*4, interval=5, repeat=True)#, blit=True)
            window.update()
            canvas.draw()

            if safe_animation_var.get() == 1:
                print("Calculate Animation")
                writervideo = animation.FFMpegWriter(fps=30, bitrate=1800)#, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
                animation_instance.save('animation/animation.gif', writer=writervideo) 
                print("Save Animation")
        else:
            animation_instance.resume()
    else:
        animation_running = False
        if animation_instance:
            animation_instance.pause()

#%%initialization            
#------------------------------------------------------------------------------
#-----------------------------INITIALIZATION-----------------------------------
#------------------------------------------------------------------------------

#objects
wg = Wavegenerator(shape = 1)
fs = Flexspline()
br = Bearing()
cs = CircularSpline()
ds = DynamicSpline()

hd = HarmonicDrive(fs, wg, cs, ds, br)

#inital paramter
hd.initial()

#colors
colors = {'wg': 'deepskyblue',
          'fs': 'springgreen',
          'br': 'grey',
          'cs': 'darkgrey',
          'ds': 'fuchsia'}

#window
window = ThemedTk(theme="arc")
window.iconbitmap("icon.ico")
window.title("Harmonic Drive")

#figure and axes
fig = plt.Figure(figsize=(8, 8), dpi=100)
fig.tight_layout()
ax = fig.add_subplot()

#animation
animation_running = False
animation_instance = []
patch_cs = []

#%%GUI
#------------------------------------------------------------------------------
#----------------------------------GUI-----------------------------------------
#------------------------------------------------------------------------------

#%%%FigureCanvasTkAgg to integrate the Matplotlib graphic into Tkinter
frame_plot = ttk.Frame(window)
canvas = FigureCanvasTkAgg(fig, frame_plot)

toolbar = NavigationToolbar2Tk(canvas, frame_plot)
toolbar.update()
toolbar.pack()

frame_plot.pack(side='left', fill='both', expand=True)
canvas.get_tk_widget().pack(side='left', fill='both', expand=True)

#%%%Notebook for the part-specific parameters
lf_parameter = ttk.LabelFrame(window, text='Parameter')
lf_parameter.pack(side='top', padx=3, pady=3)

notebook = ttk.Notebook(lf_parameter)
notebook.pack(padx=1, pady=1, side='top')

#flexspline
frame_fs = ttk.Frame(notebook)
frame_fs.pack()
notebook.add(frame_fs, text='Flexspline')

parameter_fs = dict(d_i   = dict(txt='Inner Diameter'        , abbr='d_i'  , readonly=False, widget=[]),
                    s_st  = dict(txt='Sprocket Thickness'    , abbr='s_st' , readonly=False, widget=[]),
                    d_nf  = dict(txt='Neutral Fiber Diameter', abbr='d_nf' , readonly=True , widget=[]),
                    d_f   = dict(txt='Foot Circle Diameter'  , abbr='d_f'  , readonly=True , widget=[]),
                    d     = dict(txt='Pitch Circle Diameter' , abbr='d'    , readonly=False, widget=[]),
                    d_h   = dict(txt='Head Circle Diameter'  , abbr='d_h'  , readonly=False, widget=[]),
                    z     = dict(txt='Number of Teeth'       , abbr='z'    , readonly=True , widget=[]),
                    alpha = dict(txt='Tooth Flank Angle'     , abbr='alpha', readonly=False, widget=[]),
                    c     = dict(txt='Tooth Space Ratio'     , abbr='c'    , readonly=False, widget=[]),
                    r_fh  = dict(txt='Head Flank Radius'     , abbr='r_fh' , readonly=False, widget=[]),
                    r_ff  = dict(txt='Foot Flank Radius'     , abbr='r_ff' , readonly=False, widget=[]),
                    r_hr  = dict(txt='Head Rounding'         , abbr='r_hr' , readonly=False, widget=[]),
                    r_fr  = dict(txt='Foot Rounding'         , abbr='r_fr' , readonly=False, widget=[]))

for row, p in enumerate(parameter_fs.values()):
    #text
    ttk.Label(frame_fs, text=p['txt']).grid(row=row, column=0, sticky="w")
    
    #abbreviation
    ttk.Label(frame_fs, text=p['abbr']).grid(row=row, column=1, sticky="e")
    
    #entry
    ent = ttk.Entry(frame_fs)
    if p['abbr'] == 'alpha':
        default_value = fs.__dict__[p['abbr']]*180/np.pi
    else:
        default_value = fs.__dict__[p['abbr']]
    ent.insert(0, default_value)
    ent.grid(row=row, column=2, sticky="w")
    if p['readonly']: 
        ent.configure(state="readonly")
    else:
        ent.bind("<Return>", (lambda event: update_plot()))
    p['widget'] = ent

#wavegenerator
frame_wg = ttk.Frame(notebook)
frame_wg.pack()
notebook.add(frame_wg, text='Wavegenerator')

parameter_wg = dict(a     = dict(txt='Large Semi-Axis'     , abbr='a'  , readonly=True , widget=[]),
                    b     = dict(txt='Small Semi-Axis'     , abbr='b'  , readonly=True , widget=[]),
                    arc   = dict(txt='Circle angle at apex', abbr='arc', readonly=False, widget=[]))

for row, p in enumerate(parameter_wg.values()):
    #text
    ttk.Label(frame_wg, text=p['txt']).grid(row=row, column=0, sticky="w")
    
    #abbreviation
    ttk.Label(frame_wg, text=p['abbr']).grid(row=row, column=1, sticky="e")
    
    #entry
    ent = ttk.Entry(frame_wg)
    if p['abbr'] == 'arc':
        default_value = wg.__dict__[p['abbr']]*180/np.pi
    else:
        default_value = wg.__dict__[p['abbr']]
    ent.insert(0, default_value)
    ent.grid(row=row, column=2, sticky="w")
    if p['readonly']: 
        ent.configure(state="readonly")
    else:
        ent.bind("<Return>", (lambda event: update_plot()))
    p['widget'] = ent

shapes = ('elliptical',
          '345polynomial')

shape_var = tk.StringVar(frame_wg)

for row, v in enumerate(shapes):
    r = ttk.Radiobutton(frame_wg,
                        text=v,
                        value=v,
                        variable=shape_var)
    r.grid(row = row+len(parameter_wg), column = 0, columnspan = 3)
    
    if wg.shape == v:
        r.invoke()
  
#TODO  
#harmonic drive
frame_hd = ttk.Frame(notebook)
frame_hd.pack()
notebook.add(frame_hd, text='Harmonic Drive')

entry_i    = create_entry_frame(frame_hd, "Transmission Ratio (|i_wg->ds|)", "i", -hd.i, 0)
entry_z_fs = create_entry_frame(frame_hd, "Number of teeth Flexspline", "z_fs", fs.z, 1, readonly = True)
entry_z_cs = create_entry_frame(frame_hd, "Number of teeth Circular Spline", "z_cs", cs.z, 2, readonly = True)
entry_z_ds = create_entry_frame(frame_hd, "Number of teeth Dynamic Spline", "z_ds", ds.z, 3, readonly = True)
entry_q_nf = create_entry_frame(frame_hd, "Equidistant Distance Neutral Fiber", "q_nf", hd.q_nf, 4, readonly = True)

#TODO
# frame_bearing
frame_br = ttk.Frame(notebook)
frame_br.pack()
notebook.add(frame_br, text='Bearing')

entry_d_br = create_entry_frame(frame_br, "Diameter Bearing Rolling element", "d_br", br.d_br, 0)
entry_n    = create_entry_frame(frame_br, "Number of bearing rollers", "n", br.n, 1)

#%%%View selection
lf_view_select = ttk.LabelFrame(window, text='View Selection')
lf_view_select.pack(side='top', padx=3, pady=3)

views = ('Flexspline Tooth',
         'Flexspline (undeformed)',
         'Family of Flexspline Tooth (for Circular Spline Tooth)',
         'Family of Flexspline Tooth (for Dynamic Spline Tooth)',
         'HARMONIC DRIVE')

view_var = tk.StringVar(lf_view_select)

for v in views:
    r = ttk.Radiobutton(lf_view_select,
                        text=v,
                        value=v,
                        variable=view_var)
    r.pack()

button_update = ttk.Button(lf_view_select, text="update", command=update_plot)
button_update.pack()

#%%%Plot Settings View Harmonic Drive
lf_plot_settings = ttk.LabelFrame(window, text='Plot Settings View Harmonic Drive')
lf_plot_settings.pack(side='top', padx=3, pady=3)

#colors
ttk.Label(lf_plot_settings, text='Change Colors').grid(row=0, column=0, columnspan = 3, sticky="w")

def drop_down_color(part, text, row):
    def set_c(*args):
        print(dd_color.get())
        colors[part] = dd_color.get()
    
    dd_color = ttk.Combobox(lf_plot_settings)
    dd_color.bind("<<ComboboxSelected>>", set_c)
    dd_color['values'] = list(mcolors.TABLEAU_COLORS.keys())
    dd_color.set(colors[part])
    
    ttk.Label(lf_plot_settings, text = text).grid(row = row, column = 0, columnspan=2, sticky="e")  
    dd_color.grid(row=row, column=2)

drop_down_color('wg', "Wavegenerator", row = 1) 
drop_down_color('fs', "Flexspline", row = 2)    
drop_down_color('br', "Bearing", row = 3) 
drop_down_color('cs', "Circular Spline", row = 4) 
drop_down_color('ds', "Dynamic Spline", row = 5) 

#phi_wg
#text
ttk.Label(lf_plot_settings, text='Drive Angle Wavegenerator').grid(row=6, column=0, sticky="w")

#abbreviation
ttk.Label(lf_plot_settings, text='|phi_wg|').grid(row=6, column=1, sticky="e")

#entry
ent_phi_wg = ttk.Entry(lf_plot_settings)
default_value = -hd.phi_wg * 180/np.pi
ent_phi_wg.insert(0, default_value)
ent_phi_wg.grid(row=6, column=2, sticky="w")
ent_phi_wg.bind("<Return>", (lambda event: update_plot(calc = False)))

#animation start/stop
button_start_stop = ttk.Button(lf_plot_settings, text="Start/Stop Animation", command=start_stop_animation)
button_start_stop.grid(row = 7, column = 0, columnspan=3)

#checkbox safe animation
safe_animation_var = tk.BooleanVar()
ttk.Checkbutton(lf_plot_settings,
                text='Play and Save Animation',
                variable=safe_animation_var).grid(row=8, column = 0, columnspan = 3)


#%%%Frame Load/Save Configuration
lf_config = ttk.LabelFrame(window, text='Load/Save Configuration')
lf_config.pack(side='top', padx=3, pady=3)

def save_config():
    hd.save_config(file_name = entry_save_config.get())
    drop_down['values'] = fnmatch.filter(os.listdir('config/'), '*.csv')
    
def load_config(*args):
    hd.load_config(file_name = drop_down.get())
    update_entrys()
    update_plot()

#Save
ttk.Label(lf_config, text='Save Config').grid(row=0, column=0, sticky="e")

entry_save_config = ttk.Entry(lf_config)
entry_save_config.insert(0, "---Name + Enter---")
entry_save_config.bind("<Return>", (lambda event: save_config()))
entry_save_config.grid(row=0, column=1)

#Load
ttk.Label(lf_config, text = "Load Config").grid(row = 1, column = 0, sticky="e")

drop_down_variable = tk.StringVar()
drop_down = ttk.Combobox(lf_config, textvariable=drop_down_variable)
drop_down.bind("<<ComboboxSelected>>", load_config)
drop_down['values'] = fnmatch.filter(os.listdir('config/'), '*.csv')
drop_down.set("---Select---")

drop_down.grid(row=1, column=1)

#%%% Close Button
def close():
   window.destroy()
   window.quit()

button_close = ttk.Button(window, text="close", command=close)
button_close.pack(padx=1, pady=1)

window.mainloop()