from tkinter import *
from tkinter import messagebox
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def calculate_page(button_number):

    new_window = Toplevel(root)
    Label(new_window, text=f"Enter values for eqation {button_number}").pack()
    result_label = Label(new_window, text="")
    result_label.pack()

    if button_number == 1:
        entries = [Entry(new_window, width=8) for _ in range(4)]
        labels = ["Enter K (md)", "Enter Ks (md)", "Enter rs (ft):", "Enter rw (ft):"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_1(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 2:
        entries = [Entry(new_window, width=8) for _ in range(6)]
        labels = ["Enter K (md)", "Enter Kp (md) _ permeability in the particle-invaded zone", "Enter Kf (md) _ permeability in the filtrate invaded zone", "Enter rw (ft)", "Enter rp (ft)", "Enter rf (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_2(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 3:
        entries = [Entry(new_window, width=8) for _ in range(4)]
        labels = ["Enter h (ft)", "Enter hw (ft) _ completion thickness", "Enter re (ft):", "Enter rw (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_3(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 4:
        entries = [Entry(new_window, width=8) for _ in range(4)]
        labels = ["Enter h (ft)", "Enter hw (ft) _ completion thickness", "Enter re (ft)", "Enter rw (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_4(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 5:
        entries = [Entry(new_window, width=8) for _ in range(6)]
        labels = ["Enter h (ft)", "Enter hp (ft)", "kh (md)", "Enter kv (md)", "Enter y (ft) _ distance between the top of the open interval to the top of the reservoir", "Enter rw (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_5(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 6:
        entries = [Entry(new_window, width=8) for _ in range(6)]
        labels = ["Enter h (ft)", "Enter hw (ft) _ completion thickness", "Enter h1 (ft) _ distance between the top of the open interval to the top of the reservoir)", "Enter Kv (md)", "Enter Kh (md)", "Enter rw (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_6(entries, result_label)).pack()
        Button(new_window, text="Sensitivity analysis on h1D", command=lambda: sens6_h1_input(entries)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()
        
    if button_number == 7:
        entries = [Entry(new_window, width=8) for _ in range(5)]
        labels = ["Enter θ (degree)", "Enter h (ft)", "Enter Kv (md)", "Enter Kh (md)", "Enter rw (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_7(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 8:
        entries = [Entry(new_window, width=8) for _ in range(3)]
        labels = ["Enter rw (ft)", "Enter θ (degree)", "Enter h (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_8(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 9:
        entries = [Entry(new_window, width=8) for _ in range(4)]
        labels = ["Enter rw (ft)", "Enter I_ani", "Enter h (ft)", "Enter θ (degree)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_9(entries, result_label)).pack()
        Button(new_window, text="Sensitivity analysis on θ", command=lambda: sens9_I_ani_input(entries)).pack()        
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 10:
        entries = [Entry(new_window, width=8) for _ in range(5)]
        labels = ["Enter k_H (md)", "Enter k_sH (md)", "Enter I_ani", "Enter r_sH (ft)", "Enter rw (ft)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_10(entries, result_label)).pack()
        Button(new_window, text="Sensitivity analysis on K_H/Ks_H", command=lambda: sens10_KH_KsH_input(entries)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 11:
        entries = [Entry(new_window, width=8) for _ in range(7)]
        labels = ["Enter θ (degree)", "Enter rw (ft)", "Enter l_perf (ft)", "Enter h_perf (ft) or 1/SPF", "Enter r_perf (ft)", "Enter Kh (md)", "Enter Kv (md)"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_11(entries, result_label)).pack()
        Button(new_window, text="Sensitivity analysis on SPF", command=lambda: sens11_SPF_input(entries)).pack()
        Button(new_window, text="Sensitivity analysis on θ", command=lambda: sens11_theta_input(entries)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()

    if button_number == 12:
        entries = [Entry(new_window, width=8) for _ in range(6)]
        labels = ["Enter K (md)", "Enter Kg (md) _ gravel pack permeability", "Enter Ks (md) _ damaged zone permeability", "Enter rw (ft)", "Enter rs (ft)", "Enter rg (ft) _ radius of screen"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_2(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()
    if button_number == 13:
        entries = [Entry(new_window, width=8) for _ in range(9)]
        labels = ["Enter rw (ft)", "Enter rt (ft) _ radius of the tunnel", "Enter rp (ft)", 
                "Enter hp (ft) _ perforation spacing or 1/SPF", 
                "Enter lt (ft) _ length of the tunnel through the casing and cement", 
                "Enter K (md)", "Enter Kg (md) _ gravel pack permeability", 
                "Enter θ (degree) _ (enter 0 instead of 360)", 
                "Enter sp _ perforation skin factor"]
        for i, entry in enumerate(entries):
            Label(new_window, text=labels[i]).pack()
            entry.pack()
        Button(new_window, text="Calculate", command=lambda: formula_13(entries, result_label)).pack()
        Button(new_window, text="Exit", command=new_window.destroy).pack()


def formula_1(entries, result_label):
    try:
        k = float(entries[0].get())
        k_s = float(entries[1].get())
        r_s = float(entries[2].get())
        r_w = float(entries[3].get())
        s = ((k / k_s) - 1) * math.log(r_s / r_w)
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_2(entries, result_label):
    try:
        k = float(entries[0].get())
        k_p = float(entries[1].get())
        k_f = float(entries[2].get())
        r_w = float(entries[3].get())
        r_p = float(entries[4].get())
        r_f = float(entries[5].get())
        s = ((k / k_p) * math.log(r_p / r_w) + (k / k_f) * math.log(r_f / r_p) - math.log(r_f / r_w))
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_3(entries, result_label):
    try:
        h = float(entries[0].get())
        h_w = float(entries[1].get())
        r_e = float(entries[2].get())
        r_w = float(entries[3].get())
        s = ((h / h_w) / (1 + 7*math.sqrt((r_w / (2 * h_w)) * math.cos((math.pi * h_w) / (2 * h)))) - 1) * math.log(r_e / r_w)
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_4(entries, result_label):
    try:
        h = float(entries[0].get())
        h_w = float(entries[1].get())
        r_e = float(entries[2].get())
        r_w = float(entries[3].get())
        s = ((h / h_w) / (1 + 7*math.sqrt((r_w / h_w) * math.cos((math.pi * h_w) / (2 * h)))) - 1) * math.log(r_e / r_w)
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_5(entries, result_label):
    try:
        h = float(entries[0].get())
        h_p = float(entries[1].get())
        k_h = float(entries[2].get())
        K_v = float(entries[3].get())
        y = float(entries[4].get())
        r_w = float(entries[5].get())
        z_m = y + h_p / 2
        r_wc = r_w * math.exp(0.2126 * (z_m / h + 2.753))

        s = 1.35 * (h / h_p - 1) ** 0.835 * (math.log((h * math.sqrt(k_h / K_v)) + 7) - 1.95 - math.log(r_wc) * (0.49 + 0.1 * math.log(h * math.sqrt(k_h / K_v))))
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_6(entries, result_label):
    try:
        h = float(entries[0].get())
        h_w = float(entries[1].get())
        h_1 = float(entries[2].get())
        k_v = float(entries[3].get())
        k_H = float(entries[4].get())
        r_w = float(entries[5].get())
        h_wD = (h_w / h)
        r_D = (r_w / h) * math.sqrt (k_v / k_H)
        h_1D = (h_1 / h)
        A = 1 / (h_1D + h_wD / 4)
        B = 1 / (h_1D + (3 * h_wD / 4))

        s = (1 / h_wD - 1) * math.log(math.pi / (2 * r_D)) + (1 / h_wD) * math.log(h_wD / (2 + h_wD) * math.sqrt((A-1)/(B-1)))
    
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def sens6_h1_input(entries):
    new_window = Toplevel(root)
    new_window.title("Enter values for Sensitivity Analysis")

    Label(new_window, text="Enter h1/h values (comma-separated):").pack()
    h1_h_entry = Entry(new_window, width=20)
    h1_h_entry.pack()

    Label(new_window, text="Enter I_ani range:").pack()
    i_ani_min = Entry(new_window, width=8)
    i_ani_min.pack()
    Label(new_window, text="to").pack()
    i_ani_max = Entry(new_window, width=8)
    i_ani_max.pack()

    Button(new_window, text="Plot", command=lambda: sens6_anisotropy(entries, h1_h_entry, i_ani_min, i_ani_max)).pack()
    Button(new_window, text="Exit", command=new_window.destroy).pack()

def sens6_anisotropy(entries, h1_h_entry, i_ani_min, i_ani_max):
    try:
        h = float(entries[0].get())
        h_w = float(entries[1].get())
        r_w = float(entries[5].get())

        h1_h_values = [float(x.strip()) for x in h1_h_entry.get().split(',')]
        i_ani_min = float(i_ani_min.get())
        i_ani_max = float(i_ani_max.get())

        new_window = Toplevel(root)
        new_window.title("Sensitivity Analysis on Anisotropy")

        fig, ax = plt.subplots(figsize=(10, 8))

        i_ani_values = np.linspace(i_ani_min, i_ani_max, 100)

        for h1_h in h1_h_values:
            skin_values = []
            for i_ani in i_ani_values:
                
                k_v = 1
                k_H = i_ani ** 2
                h_wD = h_w / h
                r_D = (r_w / h) * math.sqrt(k_v / k_H)
                h_1D = h1_h
                A = 1 / (h_1D + h_wD / 4)
                B = 1 / (h_1D + (3 * h_wD / 4))
                s = (1 / h_wD - 1) * math.log(math.pi / (2 * r_D)) + (1 / h_wD) * math.log(h_wD / (2 + h_wD) * math.sqrt((A-1)/(B-1)))
                skin_values.append(s)

            ax.plot(i_ani_values, skin_values, label=f'h1/h = {h1_h}')

        ax.set_xlabel('I_ani')
        ax.set_ylabel('Partial Completion Skin Factor, s')
        ax.set_title('The effect of anisotropy on the partial completion skin factor')
        ax.legend()
        ax.grid(True)

        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()

        Button(new_window, text="Exit", command=new_window.destroy).pack()

    except ValueError:
        messagebox.showerror("Error", "Please fill in all the boxes with valid numbers")

def formula_7(entries, result_label):
    try:
        theta = float(entries[0].get())
        h = float(entries[1].get())
        k_v = float(entries[2].get())
        k_h = float(entries[3].get())
        r_w = float(entries[4].get())
        theta_rad = math.radians(theta)
        theta_prime_rad = math.atan(math.sqrt(k_v / k_h) * math.tan(theta_rad))
        theta_prime_deg = math.degrees(theta_prime_rad)

        s = (-(theta_prime_deg / 41) ** 2.06) - ((theta_prime_deg / 56) ** 1.865) * math.log10(h / (100 * r_w) * math.sqrt(k_h / k_v))
    
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_8(entries, result_label):
    try:
        r_w = float(entries[0].get())
        theta = float(entries[1].get())
        h = float(entries[2].get())
        theta_rad = math.radians(theta)

        s = math.log((4 * r_w * math.cos(theta_rad)) / h) + math.cos(theta_rad) * math.log(h / (4 * r_w * (math.cos(theta_rad))**0.5))

        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_9(entries, result_label):
    try:
        r_w = float(entries[0].get())
        I_ani = float(entries[1].get())
        h = float(entries[2].get())
        theta_deg = float(entries[3].get())
         
        theta_rad = math.radians(theta_deg)
        y = math.sqrt((1/I_ani**2) + (math.cos(theta_rad))**2 * (1 - (1/ I_ani**2)))

        s = math.log(1 / (I_ani * y)*(4 * r_w * math.cos(theta_rad) / h)) + \
            (math.cos(theta_rad) / y) * math.log(((2 * I_ani * math.sqrt(y)) \
            / (1 + (1 / y))) * (h / (4 * r_w * math.sqrt(math.cos(theta_rad)))))        
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def sens9_I_ani_input(entries):
    new_window = Toplevel(root)
    new_window.title("Enter values for Sensitivity Analysis")

    Label(new_window, text="Enter I_ani range:").pack()
    i_ani_min = Entry(new_window, width=8)
    i_ani_min.pack()
    Label(new_window, text="to").pack()
    i_ani_max = Entry(new_window, width=8)
    i_ani_max.pack()

    Label(new_window, text="Enter theta values (comma-separated):").pack()
    theta_entry = Entry(new_window, width=20)
    theta_entry.pack()

    Button(new_window, text="Plot", command=lambda: sens9_I_ani(entries, i_ani_min, i_ani_max, theta_entry)).pack()
    Button(new_window, text="Exit", command=new_window.destroy).pack()

def sens9_I_ani(entries, i_ani_min, i_ani_max, theta_entry):
    try:
        r_w = float(entries[0].get())
        h = float(entries[2].get())

        i_ani_min = float(i_ani_min.get())
        i_ani_max = float(i_ani_max.get())
        theta_values = [float(x.strip()) for x in theta_entry.get().split(',')]

        new_window = Toplevel(root)
        new_window.title("Sensitivity Analysis on I_ani")

        fig, ax = plt.subplots(figsize=(10, 8))

        i_ani_values = np.linspace(i_ani_min, i_ani_max, 100)

        for theta_deg in theta_values:
            skin_values = []
            theta_rad = math.radians(theta_deg)
            for I_ani in i_ani_values:
                y = math.sqrt((1/I_ani**2) + (math.cos(theta_rad))**2 * (1 - (1/ I_ani**2)))
                s = math.log(1 / (I_ani * y)*(4 * r_w * math.cos(theta_rad) / h)) + \
                    (math.cos(theta_rad) / y) * math.log(((2 * I_ani * math.sqrt(y)) \
                    / (1 + (1 / y))) * (h / (4 * r_w * math.sqrt(math.cos(theta_rad)))))
                skin_values.append(s)

            ax.plot(i_ani_values, skin_values, label=f'θ = {theta_deg}°')

        ax.set_xlabel('Anisotropy Index, I_ani')
        ax.set_ylabel('Slant Skin, s_θ')
        ax.set_title('The effect of anisotropy on the slant skin factor')
        ax.legend()
        ax.grid(True)

        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()

        Button(new_window, text="Exit", command=new_window.destroy).pack()

    except ValueError:
        messagebox.showerror("Error", "Please fill in all the boxes with valid numbers")

def formula_10(entries, result_label):
    try:
        k_H = float(entries[0].get())
        k_sH = float(entries[1].get())
        I_ani = float(entries[2].get())
        r_sH = float(entries[3].get())
        r_w = float(entries[4].get())
        s = ((k_H / k_sH) - 1) * math.log(1 / (I_ani + 1) * (r_sH / r_w + math.sqrt((r_sH / r_w)**2 + I_ani ** 2 -1)))
    
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def sens10_KH_KsH_input(entries):
    new_window = Toplevel(root)
    new_window.title("Enter values for Sensitivity Analysis")

    Label(new_window, text="Enter K_H/Ks_H values (comma-separated):").pack()
    kh_ksh_entry = Entry(new_window, width=20)
    kh_ksh_entry.pack()

    Label(new_window, text="Enter r_sH range:").pack()
    rsh_min = Entry(new_window, width=8)
    rsh_min.pack()
    Label(new_window, text="to").pack()
    rsh_max = Entry(new_window, width=8)
    rsh_max.pack()

    Button(new_window, text="Plot", command=lambda: sens10_KH_KsH(entries, kh_ksh_entry, rsh_min, rsh_max)).pack()
    Button(new_window, text="Exit", command=new_window.destroy).pack()

def sens10_KH_KsH(entries, kh_ksh_entry, rsh_min, rsh_max):
    try:
        I_ani = float(entries[2].get())
        r_w = float(entries[4].get())

        kh_ksh_values = [float(x.strip()) for x in kh_ksh_entry.get().split(',')]
        rsh_min = float(rsh_min.get())
        rsh_max = float(rsh_max.get())

        new_window = Toplevel(root)
        new_window.title("Sensitivity Analysis on K_H/Ks_H")

        fig, ax = plt.subplots(figsize=(10, 8))

        rsh_values = np.linspace(rsh_min, rsh_max, 100)

        for kh_ksh in kh_ksh_values:
            skin_values = []
            for r_sH in rsh_values:
                s = ((kh_ksh) - 1) * math.log(1 / (I_ani + 1) * (r_sH / r_w + math.sqrt((r_sH / r_w)**2 + I_ani ** 2 -1)))
                skin_values.append(s)

            ax.plot(rsh_values, skin_values, label=f'K_H/Ks_H = {kh_ksh}')

        ax.set_xlabel('Penetration of Damage (r_sH)')
        ax.set_ylabel('Skin Factor (s)')
        ax.set_title('The effect of K_H/Ks_H on the Skin Factor')
        ax.legend()
        ax.grid(True)

        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()

        Button(new_window, text="Exit", command=new_window.destroy).pack()

    except ValueError:
        messagebox.showerror("Error", "Please fill in all the boxes with valid numbers")

def formula_11(entries, result_label):
    try:
        theta = float(entries[0].get())
        r_w = float(entries[1].get())
        l_perf = float(entries[2].get())
        h_perf = float(entries[3].get())
        r_perf = float(entries[4].get())
        k_h = float(entries[5].get())
        k_v = float(entries[6].get())

        angle_params = {
            0: {'a0': 0.250, 'a1': -2.091, 'a2': 0.0453, 'b1': 5.1313, 'b2': 1.8672, 'c1': 0.16, 'c2': 2.675},
            180: {'a0': 0.5, 'a1': -2.025, 'a2': 0.0943, 'b1': 3.0373, 'b2': 1.8115, 'c1': 0.026, 'c2': 4.532},
            120: {'a0': 0.648, 'a1': -2.018, 'a2': 0.0634, 'b1': 1.6136, 'b2': 1.7770, 'c1': 0.0066, 'c2': 5.320},
            90: {'a0': 0.726, 'a1': -1.905, 'a2': 0.1038, 'b1': 1.5674, 'b2': 1.6935, 'c1': 0.0019, 'c2': 6.155},
            60: {'a0': 0.813, 'a1': -1.898, 'a2': 0.1023, 'b1': 1.3654, 'b2': 1.6490, 'c1': 0.0003, 'c2': 7.509},
            45: {'a0': 0.860, 'a1': -1.788, 'a2': 0.2398, 'b1': 1.1915, 'b2': 1.6392, 'c1': 0.000046, 'c2': 8.791}
        }

        if theta not in angle_params:
            raise ValueError("Invalid theta. Please enter one of: 0, 45, 60, 90, 120, 180")

        params = angle_params[theta]

        # Calculate horizontal skin factor (S_H)
        rw_prime = params['a0'] * (r_w + l_perf) if theta != 0 else l_perf / 4
        S_H = math.log(r_w / rw_prime)

        # Calculate vertical skin factor (S_V)
        h_D = (h_perf / l_perf) * math.sqrt(k_h / k_v)
        r_D = (r_perf / (2 * h_perf)) * (1 + math.sqrt(k_v / k_h))
        a = params['a1'] * math.log10(r_D) + params['a2']
        b = params['b1'] * r_D + params['b2']
        S_V = (10**a) * (h_D**(b-1)) * (r_D**b)

        # Calculate wellbore skin factor (S_wb)
        r_wD = r_w / (l_perf + r_w)
        S_wb = params['c1'] * math.exp(params['c2'] * r_wD)

        # Calculate total perforation skin
        S_p = S_H + S_V + S_wb

        result_label.config(text=f"S_H: {S_H:.4f}\nS_V: {S_V:.4f}\nS_wb: {S_wb:.4f}\nTotal perforation skin: {S_p:.4f}", bg="greenyellow")
    except ValueError as e:
        result_label.config(text=f"Error: {str(e)}", bg="orangered")


def sens11_SPF_input(entries):
    new_window = Toplevel(root)
    new_window.title("Enter values for Sensitivity Analysis")

    Label(new_window, text="Enter SPF values (comma-separated):").pack()
    spf_entry = Entry(new_window, width=20)
    spf_entry.pack()

    Label(new_window, text="Enter kH/kv range:").pack()
    kh_kv_min = Entry(new_window, width=8)
    kh_kv_min.pack()
    Label(new_window, text="to").pack()
    kh_kv_max = Entry(new_window, width=8)
    kh_kv_max.pack()

    Button(new_window, text="Plot", command=lambda: sens11_SPF(entries, spf_entry, kh_kv_min, kh_kv_max)).pack()
    Button(new_window, text="Exit", command=new_window.destroy).pack()

def sens11_SPF(entries, spf_entry, kh_kv_min, kh_kv_max):
    try:
        theta = float(entries[0].get())
        r_w = float(entries[1].get())
        l_perf = float(entries[2].get())
        r_perf = float(entries[4].get())

        spf_values = [float(x.strip()) for x in spf_entry.get().split(',')]
        kh_kv_min = float(kh_kv_min.get())
        kh_kv_max = float(kh_kv_max.get())

        new_window = Toplevel(root)
        new_window.title("Sensitivity Analysis on SPF")

        fig, ax = plt.subplots(figsize=(12, 8))

        kh_kv_values = np.linspace(kh_kv_min, kh_kv_max, 100)

        angle_params = {
            0: {'a0': 0.250, 'a1': -2.091, 'a2': 0.0453, 'b1': 5.1313, 'b2': 1.8672, 'c1': 0.16, 'c2': 2.675},
            180: {'a0': 0.5, 'a1': -2.025, 'a2': 0.0943, 'b1': 3.0373, 'b2': 1.8115, 'c1': 0.026, 'c2': 4.532},
            120: {'a0': 0.648, 'a1': -2.018, 'a2': 0.0634, 'b1': 1.6136, 'b2': 1.7770, 'c1': 0.0066, 'c2': 5.320},
            90: {'a0': 0.726, 'a1': -1.905, 'a2': 0.1038, 'b1': 1.5674, 'b2': 1.6935, 'c1': 0.0019, 'c2': 6.155},
            60: {'a0': 0.813, 'a1': -1.898, 'a2': 0.1023, 'b1': 1.3654, 'b2': 1.6490, 'c1': 0.0003, 'c2': 7.509},
            45: {'a0': 0.860, 'a1': -1.788, 'a2': 0.2398, 'b1': 1.1915, 'b2': 1.6392, 'c1': 0.000046, 'c2': 8.791}
        }

        if theta not in angle_params:
            raise ValueError("Invalid theta. Please enter one of: 0, 45, 60, 90, 120, 180")

        params = angle_params[theta]
        
        for spf in spf_values:
            s_v_values = []
            h_perf = 1 / spf  # SPF = 1 / h_perf
            for kh_kv in kh_kv_values:
                h_D = (h_perf / l_perf) * math.sqrt(kh_kv)
                r_D = (r_perf / (2 * h_perf)) * (1 + math.sqrt(1 / kh_kv))

                a = params['a1'] * math.log10(r_D) + params['a2']
                b = params['b1'] * r_D + params['b2']
                S_V = (10**a) * (h_D**(b-1)) * (r_D**b)
                s_v_values.append(S_V)

            ax.plot(kh_kv_values, s_v_values, label=f'SPF = {spf}')

        ax.set_xlabel('kH/kv')
        ax.set_ylabel('Vertical Skin Factor, S_V')
        ax.set_title(f'The effect of kH/kv on the Vertical Skin Factor for θ = {theta}°')
        ax.legend()
        ax.grid(True)

        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()

        Button(new_window, text="Exit", command=new_window.destroy).pack()

    except ValueError as e:
        messagebox.showerror("Error", str(e))

def sens11_theta_input(entries):
    new_window = Toplevel(root)
    new_window.title("Enter values for Sensitivity Analysis")

    Label(new_window, text="Enter θ values (comma-separated):").pack()
    theta_entry = Entry(new_window, width=20)
    theta_entry.pack()

    Label(new_window, text="Enter kH/kv range:").pack()
    kh_kv_min = Entry(new_window, width=8)
    kh_kv_min.pack()
    Label(new_window, text="to").pack()
    kh_kv_max = Entry(new_window, width=8)
    kh_kv_max.pack()

    Button(new_window, text="Plot", command=lambda: sens11_theta(entries, theta_entry, kh_kv_min, kh_kv_max)).pack()
    Button(new_window, text="Exit", command=new_window.destroy).pack()

def sens11_theta(entries, theta_entry, kh_kv_min, kh_kv_max):
    try:
        r_w = float(entries[1].get())
        l_perf = float(entries[2].get())
        h_perf = float(entries[3].get())
        r_perf = float(entries[4].get())

        theta_values = [float(x.strip()) for x in theta_entry.get().split(',')]
        kh_kv_min = float(kh_kv_min.get())
        kh_kv_max = float(kh_kv_max.get())

        new_window = Toplevel(root)
        new_window.title("Sensitivity Analysis on θ")

        fig, axs = plt.subplots(2, 2, figsize=(16, 12))
        axs = axs.ravel()

        kh_kv_values = np.linspace(kh_kv_min, kh_kv_max, 100)

        angle_params = {
            0: {'a0': 0.250, 'a1': -2.091, 'a2': 0.0453, 'b1': 5.1313, 'b2': 1.8672, 'c1': 0.16, 'c2': 2.675},
            180: {'a0': 0.5, 'a1': -2.025, 'a2': 0.0943, 'b1': 3.0373, 'b2': 1.8115, 'c1': 0.026, 'c2': 4.532},
            120: {'a0': 0.648, 'a1': -2.018, 'a2': 0.0634, 'b1': 1.6136, 'b2': 1.7770, 'c1': 0.0066, 'c2': 5.320},
            90: {'a0': 0.726, 'a1': -1.905, 'a2': 0.1038, 'b1': 1.5674, 'b2': 1.6935, 'c1': 0.0019, 'c2': 6.155},
            60: {'a0': 0.813, 'a1': -1.898, 'a2': 0.1023, 'b1': 1.3654, 'b2': 1.6490, 'c1': 0.0003, 'c2': 7.509},
            45: {'a0': 0.860, 'a1': -1.788, 'a2': 0.2398, 'b1': 1.1915, 'b2': 1.6392, 'c1': 0.000046, 'c2': 8.791}
        }

        for theta in theta_values:
            if theta not in angle_params:
                raise ValueError(f"Invalid theta: {theta}. Please enter one of: 0, 45, 60, 90, 120, 180")

            params = angle_params[theta]
            s_v_values, s_h_values, s_wb_values, s_p_values = [], [], [], []

            for kh_kv in kh_kv_values:
                k_h = math.sqrt(kh_kv)
                k_v = 1 / math.sqrt(kh_kv)

                # Calculate S_H
                rw_prime = params['a0'] * (r_w + l_perf) if theta != 0 else l_perf / 4
                S_H = math.log(r_w / rw_prime)

                # Calculate S_V
                h_D = (h_perf / l_perf) * math.sqrt(k_h / k_v)
                r_D = (r_perf / (2 * h_perf)) * (1 + math.sqrt(k_v / k_h))
                a = params['a1'] * math.log10(r_D) + params['a2']
                b = params['b1'] * r_D + params['b2']
                S_V = (10**a) * (h_D**(b-1)) * (r_D**b)

                # Calculate S_wb
                r_wD = r_w / (l_perf + r_w)
                S_wb = params['c1'] * math.exp(params['c2'] * r_wD)

                # Calculate total perforation skin
                S_p = S_H + S_V + S_wb

                s_v_values.append(S_V)
                s_h_values.append(S_H)
                s_wb_values.append(S_wb)
                s_p_values.append(S_p)

            axs[0].plot(kh_kv_values, s_v_values, label=f'θ = {theta}°')
            axs[1].plot(kh_kv_values, s_h_values, label=f'θ = {theta}°')
            axs[2].plot(kh_kv_values, s_wb_values, label=f'θ = {theta}°')
            axs[3].plot(kh_kv_values, s_p_values, label=f'θ = {theta}°')

        titles = ['Vertical Skin Factor, S_V', 'Horizontal Skin Factor, S_H', 
                  'Wellbore Skin Factor, S_wb', 'Total Perforation Skin, S_p']
        
        for i, ax in enumerate(axs):
            ax.set_xlabel('kH/kv')
            ax.set_ylabel('Skin Factor')
            ax.set_title(titles[i])
            ax.legend()
            ax.grid(True)

        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()

        Button(new_window, text="Exit", command=new_window.destroy).pack()

    except ValueError as e:
        messagebox.showerror("Error", str(e))

def formula_12(entries, result_label):
    try:
        k = float(entries[0].get())
        k_g = float(entries[1].get())
        k_s = float(entries[2].get())
        r_w = float(entries[3].get())
        r_s = float(entries[4].get())
        r_g = float(entries[4].get())

        s = k / k_g * math.log(r_w / r_g) + (k / k_s - 1) * math.log(r_s / r_w)    
        result_label.config(text=f"The calculated skin value is {s}", bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in the boxes with numbers", bg="orangered")

def formula_13(entries, result_label):
    try:
        r_w = float(entries[0].get())
        r_t = float(entries[1].get())
        r_p = float(entries[2].get())
        h_p = float(entries[3].get())
        l_t = float(entries[4].get())
        k = float(entries[5].get())
        k_g = float(entries[6].get())
        theta = float(entries[7].get())
        s_p = float(entries[8].get())

        r_tD = r_t / r_w
        h_pD = h_p / r_w
        l_tD = l_t / r_w
        k_gD = k_g / k
        r_pD = r_p / r_w

        # The skin component caused by flow through the perforation tunnel through the cement and casing (S_CG,ic)
        S_CGic = (2 * h_pD / r_tD**2) * (l_tD / k_gD)

        # Perforated liner skin factor
        if theta == 0:
            S_pl = ((3 * h_pD) / (2 * r_pD)) + math.log(1.5**2 / (h_pD**2 * (1 + 1.5))) - 0.61
        else:
            v = math.sin(math.pi * theta / 360)
            S_pl = ((3 * h_pD) / (2 * r_pD)) + math.log(v**2 / (h_pD**2 * (1 + v))) - 0.61

        # The skin factor caused by flow to and within the perforation cavity extending into the formation (S_CG,oc)
        S_CGoc = (1 - k_gD**(-0.5)) * s_p + (k_gD**(-0.5)) * S_pl

        # Cased-hole gravel pack skin factor
        S_CG = S_CGic + S_CGoc

        result = f"S_CG,ic: {S_CGic:.4f}\nS_CG,oc: {S_CGoc:.4f}\nCased-hole gravel pack skin factor (S_CG): {S_CG:.4f}"
        result_label.config(text=result, bg="greenyellow")
    except ValueError:
        result_label.config(text="Please fill in all the boxes with valid numbers", bg="orangered")

from tkinter import *

def on_enter(event, button):
    button.config(borderwidth=1, relief="solid", highlightthickness=0, bg="#f0e6d2") 

def on_leave(event, button):
    button.config(borderwidth=1, relief="raised", highlightthickness=0, bg="ivory") 

root = Tk()
root.title("SkinCalc")
root.iconbitmap("F:\\software_icon.ico")

mylabel = Label(root, text="Choose the equation", bg="lightgrey")
mylabel.pack()

buttons_info = [
    ("1. Hawkins' formula", 1),
    ("2. Hawkins formula for Concentric Radial Damage Zones", 2),
    ("3. Muskat formula for calculating the partially completion skin factor", 3),
    ("4. Muskat (hw in the middle of a reservoir)", 4),
    ("5. Odeh equation", 5),
    ("6. Papatzacos formula", 6),
    ("7. Well Deviation skin (Cinco)", 7),
    ("8. Well Deviation skin (Besson isotropic)", 8),
    ("9. Well Deviatian skin (Besson anisotropic)", 9),
    ("10. Horizontal well damage skin", 10),
    ("11. Perforation skin", 11),
    ("12. Open hole gravel pack skin factor", 12),
    ("13. Cased hole gravel pack skin factor", 13)
]

for text, page in buttons_info:
    button = Button(root, text=text, bg="ivory", width=100, height=2, border=1, highlightthickness=0, command=lambda page=page: calculate_page(page))
    button.pack()
    button.bind("<Enter>", lambda event, b=button: on_enter(event, b))
    button.bind("<Leave>", lambda event, b=button: on_leave(event, b))

root.mainloop()