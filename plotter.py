import numpy as np

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.animation as animation
import matplotlib.patches as mpl_patches

import tkinter as tk
import tkinter.ttk as ttk

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
    frame = tk.Frame(master=parent_frame)
    frame.grid(row=row, column=0, sticky="w")
    tk.Label(master=frame, text=description).pack()
    
    frame = tk.Frame(master=parent_frame)
    frame.grid(row=row, column=1, sticky="e")
    tk.Label(master=frame, text=parameter).pack()
    
    frame = tk.Frame(master=parent_frame)
    frame.grid(row=row, column=2, sticky="w")
    ent = tk.Entry(master=frame)
    ent.insert(0, default_value)
    ent.bind("<Return>", (lambda event: update_plot()))
    ent.pack()
    
    if readonly: ent.configure(state="readonly")
    
    return ent

def create_checkbox(parent_frame, description, command, row):
    check_var = tk.BooleanVar(window)
    
    frame = tk.Frame(master=parent_frame)
    frame.grid(row=row, column=0, sticky="e", columnspan = 2)
    tk.Label(master=frame, text=description).pack()
    
    if command == "update_checkbox_wg":
        checkbox = ttk.Checkbutton(master = parent_frame,
                                   variable=check_var,
                                   command=lambda:update_checkbox_wg(check_var))
    elif command == "update_checkbox_view":
        checkbox = ttk.Checkbutton(master = parent_frame,
                                   variable=check_var,
                                   command=lambda:update_checkbox_view(check_var))
    else:
         checkbox = ttk.Checkbutton(master = parent_frame,
                                    variable=check_var)
        
    checkbox.grid(row=row, column=2, sticky = "w")

    return check_var

def update_checkbox_view(current_var):
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
            
def update_checkbox_wg(current_var):
    if current_var.get() == 1:
        if current_var == check_var_elliptical:
            check_var_345polynomial.set(0)
            entry_arc.configure(state="readonly")
        elif current_var == check_var_345polynomial:
            check_var_elliptical.set(0)
            entry_arc.configure(state="normal")

def update_plot(calc = True):
    print("update...")
    
    #hd
    hd.i      = -float(entry_i.get())
    
    #fs
    fs.d_i    = float(entry_d_i.get())
    fs.s_st   = float(entry_s_st.get())
    fs.d      = float(entry_d.get())
    fs.d_h    = float(entry_d_h.get())
    
    fs.alpha  = float(entry_alpha.get())*np.pi/180
    fs.c      = float(entry_c.get())
    fs.r_fh   = float(entry_r_fh.get())
    fs.r_ff   = float(entry_r_ff.get())
    fs.r_hr   = float(entry_r_hr.get())
    fs.r_fr   = float(entry_r_fr.get())

    #br
    br.d_br   = float(entry_d_br.get())
    br.n      = float(entry_n.get())
    
    #wg:
    if check_var_elliptical.get() == 1:
        wg.shape = 'elliptical'
    elif check_var_345polynomial.get() == 1:
        wg.shape = '345polynomial'

    wg.arc    = float(entry_arc.get())*np.pi/180
    
    hd.update()
    update_entrys()
    
    #update selected plot
    ax.clear()
    
    if check_var_tooth.get() == 1:
        plot_view_tooth()
    elif check_var_gear.get() == 1:
        plot_view_fs_gear()
    elif check_var_harmonic_drive.get() == 1:
        if calc: hd.calc()
        hd.phi_wg = -float(entry_phi_wg.get())*np.pi/180
        plot_view_harmonic_drive(calc = calc)
    elif check_var_array_flexspline_tooth.get() == 1:
        hd.calc_flexspline_kinematics()
        hd.calc_circular_spline_flank()
        plot_view_array_flexspline_tooth()
    elif check_var_3_array_flexspline_tooth.get() == 1:
        hd.calc_flexspline_kinematics()
        hd.calc_dynamic_spline_flank()
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
    entry.delete(0, tk.END)
    entry.insert(0, value)
    entry.configure(state="readonly")

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
        color = 'darkgrey'
        plot_coordinate_system(ax, length=wg.b/2, phi=0, facecolor = color) #CS0
        try:
            patch_cs = plot_polygon(ax, cs.polygon(),
                                      ec = "black",
                                      lw = 1,
                                      fc = (color, 0.8))
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
    color = 'deepskyblue'
    plot_coordinate_system(ax, length=wg.b/2, phi=hd.phi_wg, facecolor = color) #CS1
    #profil_wg_phi_wg = wg.profil_rotated(phi_wg)
    #ax.plot(profil_wg_phi_wg[:,0], profil_wg_wg[:,1])
    #ax.plot([0,np.cos(phi_wg)*a_el],[0,np.sin(phi_wg)*a_el], ls = 'dashed',lw = 1, c = "darkviolet" , alpha = 0.5)
    try:
        plot_polygon(ax, wg.polygon(),
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
        plot_polygon(ax, fs.polygon(),
                            ec = "black",
                            lw = 1,
                            fc = (color, 0.8))
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
    
    ax.arrow(fs.r_z_i[0,0], fs.r_z_i[0,1], e_2_1_phi_wg[0]*4, e_2_1_phi_wg[1]*4, facecolor = color, lw=0.5, head_length=1, head_width=0.4, length_includes_head=True, edgecolor = 'black')
    ax.arrow(fs.r_z_i[0,0], fs.r_z_i[0,1], e_2_2_phi_wg[0]*4, e_2_2_phi_wg[1]*4, facecolor = color, lw=0.5, head_length=1, head_width=0.4, length_includes_head=True, edgecolor = 'black')
    
    #---Dynamic Spline---
    color = 'fuchsia'
    plot_coordinate_system(ax, length=wg.b/2, phi=hd.phi_ds, facecolor = color) #CS3
    ds.polygon()
    try:
        plot_polygon(ax, ds.polygon(),
                            ec = "black",
                            lw = 1,
                            fc = (color, 0.5))
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
    
    phi_wg = -frame/2*np.pi/180
    hd.phi_wg = phi_wg
    update_rotation()
    entry_phi_wg.delete(0, tk.END)
    entry_phi_wg.insert(0, int(-phi_wg*180/np.pi))
    window.update()
    
    return ax.get_children()

def start_stop_animation():
    global animation_running, animation_instance
    
    if animation_running == False:
        animation_running = True
        if not animation_instance:
            
            animation_instance = animation.FuncAnimation(fig=fig, func=update_animation,
                                                              frames = 360*2, interval=5, repeat=True)#, blit=True)
            window.update()
            canvas.draw()

            if check_var_safe_animation.get() == 1:
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
            
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#objects
fs = Flexspline()
wg = Wavegenerator(shape = 1)
cs = CircularSpline()
ds = DynamicSpline()
br = Bearing()
hd = HarmonicDrive(fs, wg, cs, ds, br)

#change inital paramter
hd.initial()

#window
window = tk.Tk()
window.iconbitmap("icon.ico")
window.title("Harmonic Drive")

#figure and axes for Matplotlib
fig = plt.Figure(figsize=(10, 10), dpi=100)
ax = fig.add_subplot()

#animation
animation_running = False
animation_instance = []
patch_cs = []

#FigureCanvasTkAgg to integrate the Matplotlib graphic into Tkinter
frame_plot = ttk.Frame(window)
frame_plot.pack(side=tk.LEFT, fill = 'y')

canvas = FigureCanvasTkAgg(fig, frame_plot)
canvas_widget = canvas.get_tk_widget()

toolbar = NavigationToolbar2Tk(canvas, frame_plot)
toolbar.update()

canvas_widget.pack()
toolbar.pack()

# Entrys
frame_entrys = tk.Frame(window)
frame_entrys.pack(side=tk.TOP, padx=1, pady=1)

# frame_flexspline
frame_flexspline = tk.Frame(frame_entrys,
                            highlightbackground= 'springgreen',
                            highlightthickness=2)
frame_flexspline.grid(row=0, column=0, sticky="w")
tk.Label(frame_flexspline, text="Flexspline Parameter").grid(row=0, column=0, columnspan=3)

entry_d_i  = create_entry_frame(frame_flexspline, "Inner Diameter", "d_i", fs.d_i, 1)
entry_s_st = create_entry_frame(frame_flexspline, "Sprocket Thickness", "s_st", fs.s_st, 2)
entry_d_nf = create_entry_frame(frame_flexspline, "Neutral Fiber Diameter", "d_nf", fs.d_nf, 3, readonly = True)
entry_d_f  = create_entry_frame(frame_flexspline, "Foot Circle Diameter", "d_f", fs.d_f, 4, readonly = True)
entry_d    = create_entry_frame(frame_flexspline, "Pitch Circle Diameter", "d", fs.d, 5)
entry_d_h  = create_entry_frame(frame_flexspline, "Head Circle Diameter","d_h", fs.d_h, 6)
entry_z    = create_entry_frame(frame_flexspline, "Number of Teeth", "z",fs.z, 7, readonly = True)

# frame_tooth
frame_tooth = tk.Frame(frame_entrys,
                       highlightbackground= 'springgreen',
                       highlightthickness=2)
frame_tooth.grid(row=0, column=1, sticky="w")
tk.Label(frame_tooth, text="Flexspline Tooth Parameter").grid(row=0, column=0, columnspan=3)

entry_alpha = create_entry_frame(frame_tooth, "Tooth Flank Angle", "alpha", fs.alpha * 180/np.pi, 1)
entry_c     = create_entry_frame(frame_tooth, "Tooth Space Ratio", "c", fs.c, 2)
entry_r_fh  = create_entry_frame(frame_tooth, "Head Flank Radius", "r_fh", fs.r_fh, 3)
entry_r_ff  = create_entry_frame(frame_tooth, "Foot Flank Radius", "r_ff", fs.r_ff, 4)
entry_r_hr  = create_entry_frame(frame_tooth, "Head Rounding", "r_hr", fs.r_hr, 5)
entry_r_fr  = create_entry_frame(frame_tooth, "Foot Rounding", "r_fr", fs.r_fr, 6)

# frame_wavegenerator
frame_wavegenerator = tk.Frame(window,
                               highlightbackground= 'deepskyblue',
                               highlightthickness=2)
frame_wavegenerator.pack(side=tk.TOP, padx=1, pady=1)
tk.Label(frame_wavegenerator, text="Wavegenerator Parameter").grid(row=0, column=0, columnspan=3)

entry_a   = create_entry_frame(frame_wavegenerator, "Large Semi-Axis", "a", wg.a, 1, readonly = True)
entry_b   = create_entry_frame(frame_wavegenerator, "Small Semi-Axis", "b", wg.b, 2, readonly = True)
entry_arc = create_entry_frame(frame_wavegenerator, "Arc", "gamma", wg.arc * 180/np.pi, 3)
check_var_elliptical    = create_checkbox(frame_wavegenerator, "Elliptical", "update_checkbox_wg", 4)
check_var_345polynomial = create_checkbox(frame_wavegenerator, "345polynomial", "update_checkbox_wg", 5)

if wg.shape == 'elliptical':
    check_var_345polynomial.set(0)
    check_var_elliptical.set(1)
    entry_arc.configure(state="readonly")
elif wg.shape == '345polynomial':
    check_var_345polynomial.set(1)
    check_var_elliptical.set(0)
    entry_arc.configure(state="normal")

# frame_harmonic_drive
frame_harmonic_drive = tk.Frame(window,
                                highlightbackground= 'black',
                                highlightthickness=2)
frame_harmonic_drive.pack(side=tk.TOP, padx=1, pady=1)
tk.Label(frame_harmonic_drive, text="Harmonic Drive Parameter").grid(row=0, column=0, columnspan=3)

entry_i    = create_entry_frame(frame_harmonic_drive, "Transmission Ratio (|i_wg->ds|)", "i", -hd.i, 1)
entry_z_fs = create_entry_frame(frame_harmonic_drive, "Number of teeth Flexspline", "z_fs", fs.z, 2, readonly = True)
entry_z_cs = create_entry_frame(frame_harmonic_drive, "Number of teeth Circular Spline", "z_cs", cs.z, 3, readonly = True)
entry_z_ds = create_entry_frame(frame_harmonic_drive, "Number of teeth Dynamic Spline", "z_ds", ds.z, 4, readonly = True)
entry_q_nf = create_entry_frame(frame_harmonic_drive, "Equidistant Distance Neutral Fiber", "q_nf", hd.q_nf, 5, readonly = True)

# frame_bearing
frame_bearing = tk.Frame(window,
                         highlightbackground= 'grey',
                         highlightthickness=2)
frame_bearing.pack(side=tk.TOP, padx=1, pady=1)
ttk.Label(frame_bearing, text="Bearing Parameter").grid(row=0, column=0, columnspan=3)

entry_d_br = create_entry_frame(frame_bearing, "Diameter Bearing Rolling element", "d_br", br.d_br, 1)
entry_n    = create_entry_frame(frame_bearing, "Number of bearing rollers", "n", br.n, 2)

# frame_selection
frame_selection = tk.Frame(window,
                           highlightbackground= 'black',
                           highlightthickness=2)
frame_selection.pack(side=tk.TOP, padx=1, pady=1)
tk.Label(frame_selection, text="Select View").grid(row=0, column=0, columnspan=3)

check_var_gear                     = create_checkbox(frame_selection, "Flexspline", "update_checkbox_view", 1)
check_var_tooth                    = create_checkbox(frame_selection, "Flexspline Tooth", "update_checkbox_view", 2)
check_var_array_flexspline_tooth   = create_checkbox(frame_selection, "Array of Flexspline Tooth", "update_checkbox_view", 3)
check_var_3_array_flexspline_tooth = create_checkbox(frame_selection, "Array of Flexspline Tooth", "update_checkbox_view", 4)
check_var_harmonic_drive           = create_checkbox(frame_selection, "Harmonic Drive Gear", "update_checkbox_view", 5)

button_update = ttk.Button(frame_selection, text="update", command=update_plot)
button_update.grid(row = 6, column = 0, columnspan=3)

# frame_plot_settings
frame_plot_settings = tk.Frame(window,
                               highlightbackground= 'black',
                               highlightthickness=2)
frame_plot_settings.pack(side=tk.TOP, padx=1, pady=1)
ttk.Label(frame_plot_settings, text="Plot Settings").grid(row=0, column=0, columnspan=3)

entry_phi_wg = create_entry_frame(frame_plot_settings, "Drive Angle Wavegenerator", "|phi_wg|", -hd.phi_wg * 180/np.pi, 1)
entry_phi_wg.bind("<Return>", (lambda event: update_plot(calc = False)))

button_start_stop = ttk.Button(frame_plot_settings, text="Start/Stop Animation", command=start_stop_animation)
button_start_stop.grid(row = 2, column = 0, columnspan=3)

check_var_safe_animation = create_checkbox(frame_plot_settings, "Safe animation.mp4", "-y", 3)
        
#---

# close button
def close():
   window.destroy()
   window.quit()

button_close = ttk.Button(window, text="close", command=close)
button_close.pack(padx=1, pady=1)


window.mainloop()