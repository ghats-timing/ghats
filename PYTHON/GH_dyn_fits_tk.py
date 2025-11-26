#!/usr/bin/env python3
import tkinter as tk
from tkinter import ttk
import argparse
from astropy.io import fits
from numpy import array, arange, linspace, interp, log10, quantile, array2string, multiply
from numpy import diff, min, where, append, split, ones, zeros, int32, zeros
from numpy import abs as npabs, nan as npnan
import matplotlib
matplotlib.use('TkAgg')  # Set the backend to TkAgg
from matplotlib.colors import LogNorm, Normalize as Norm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

# GH_dyn_fits_tk.py
# Makes dynPDS and lightcurves from GHATS FITS files exported using GH_dyn_fits.pro
# Version 2025/07/29 updated by Federico Garcia
#   Solved a bug in the ytick formatting when --savePDF/PNG or --noshow are used
# Version 2025/07/19 updated by Federico Garcia (+M. Méndez)
#   Now captures some arguments from GH_dyn_fits.pro and GH_dyn_reim_fits.pro
#   Better number representation and Leahy Units added to colorbar
# Version 2025/07/18 updated by Federico Garcia
#   Added the hability to show the GTI gaps with vertical lines using --showGaps
#   Fixed slight bug on the highlighted freq and times when using tmin and tmax
# Version 2025/07/17 updated by Federico Garcia
#   Added Crosshair now also prints de Power, together with PDS, time and Freq.
#   Added Interpolator option with arguments -interp or --interpolator
#   Fixed slight bug on the use of tmin tmax for the lightcurve.
# Version 2025/07/15 updated by Federico Garcia
#   Added Colorbar option with arguments -cb or --colorbar
# Version 2025/06/30 updated by Federico Garcia
#   Added interactive Mouse Coordinates and Freq and Time Highlights
# Version 2024/02/24 updated by Federico Garcia
#   Added TZERO to include a reference time for the X axis
#   Added the possibility to fill GTI gaps and plot the full timespan --fillGaps
# Version 2023/07/28 updated by Federico Garcia
#   Updated interactivity using tkinter Scales, Textboxes and Buttons
# Version 2023/07/27 created by Federico Garcia
title = 'GH_DYN_FITS v2025/07/29'

#Prepare the plot scheme
plt.rcParams['axes.linewidth'] = 1.3
plt.rcParams['axes.axisbelow'] = False
plt.rcParams['grid.linewidth'] = 1.2
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.3
plt.rcParams['ytick.major.width'] = 1.3
plt.rcParams['xtick.minor.width'] = 1.3
plt.rcParams['ytick.minor.width'] = 1.3
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.minor.size'] = 5

plt.rcParams['figure.figsize'] = (9,7.5)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
#plt.rcParams['axes.labelsize'] = 7
#plt.rcParams['xtick.labelsize'] = 7
#plt.rcParams['ytick.labelsize'] = 7


# Parse arguments
parser = argparse.ArgumentParser(prog=title,
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          description='Make dynPDS and lightcurve from GHATS FITS file exported using GH_dyn_fits.')

parser.add_argument("filename", help="Input filename. The program assumes GH_DYN_FITS format.", type=str, default='dynPDS.fits')
parser.add_argument('-title','--title', help="Add a title to the tkinter window", type=str, default=None, required=None)
parser.add_argument('-t0','--TZERO', help="MJD0 for dynPDS", type=float, default=None, required=None)
parser.add_argument('-cmap','--cmap', help="Matplotlib colorbar (i.e. gray, jet, rainbow, heat, ...)", type=str, default='viridis', required=None)
parser.add_argument('-interp','--interpolator', help="Matplotlib imshow interpolator (i.e. none, linear, gaussian, ...)", type=str, default='gaussian', required=None)
parser.add_argument('-pmin','--pmin', help="Minimum power in the colorbar", type=float, default=None, required=None)
parser.add_argument('-pmax','--pmax', help="Maximum power in the colorbar", type=float, default=None, required=None)
parser.add_argument('-fmin','--fmin', help="Minimum frequency in the dynPDS", type=float, default=None, required=None)
parser.add_argument('-fmax','--fmax', help="Maximum frequency in the dynPDS", type=float, default=None, required=None)
parser.add_argument('-ft','--tick-freqs', help="Ticked frequencies in the dynPDS", nargs='+', default=None, required=None)
parser.add_argument('-fh','--highlight-freqs', help="Highlighted frequencies in the dynPDS", nargs='+', default=None, required=None)
parser.add_argument('-hc','--highlight-color', help="Color used to highlight times and frequencies", type=str, default='C3', required=None)
parser.add_argument('-tmin','--tmin', help="Minimum time in the dynPDS and lightcurve", type=float, default=None, required=None)
parser.add_argument('-tmax','--tmax', help="Maximum time in the dynPDS and lightcurve", type=float, default=None, required=None)
parser.add_argument('-tt','--tick-times', help="Ticked times in the lightcurve", nargs='+', default=None, required=None)
parser.add_argument('-th','--highlight-times', help="Highlighted times in the dynPDS and lightcurve", nargs='+', default=None, required=None)
parser.add_argument('-trebin','--trebin', help="Rebin factor in time, for the xaxis label", type=int, default=None, required=None)
parser.add_argument('--savepdf', help="Save figure as PDF", action='store_true', default=False, required=None)
parser.add_argument('--savepng', help="Save figure as PDF", action='store_true', default=False, required=None)
parser.add_argument('--noshow', help="Do not display figure", action='store_true', default=False, required=None)
parser.add_argument('--dynPDSonly', help="Hide the Lightcurve", action='store_true', default=False, required=None)
parser.add_argument('--nupnu', help="Plot nuPnu instead of Pnu", action='store_true', default=False, required=None)
parser.add_argument('--linear', help="Plot linear cmap", action='store_true', default=False, required=None)
parser.add_argument("-cb", "--colorbar", action="store_true", help="Show colorbar at startup")
parser.add_argument('--fillGaps', help="Plot full timespan with gaps due to GTIs", action='store_true', default=False, required=None)
parser.add_argument('--showGaps', help="Plot vertical lines at the GTI gaps", action='store_true', default=False, required=None)

args = parser.parse_args()

# Add a title to the titlebar of the tkinter window
if args.title:
    title += ' | '
    title += args.title

# Auxiliary formatting functions
def format_significant_3(x):
    s = f"{x:.3g}"  # hasta 3 dígitos significativos
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s

def format_significant_2(x):
    s = f"{x:.2g}"  # hasta 2 dígitos significativos
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s

def format_significant_1(x):
    s = f"{x:.1g}"  # hasta 1 dígitos significativos
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s

def format_significant_0(x):
    s = f"{x:.0f}"  # hasta 0 dígitos significativos
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s

# Read FITS file
f = fits.open(args.filename)
image = f[0].data
freqs = f[1].data
times = f[2].data #- f[2].data[0]
rates = f[3].data

if args.nupnu:
    image = (freqs*image.T).T

# Highlight the gaps due to GTIs if showGaps = True
def show_gaps(times, rates, image):
    time_diff = diff(times)
    timestep = min(time_diff)
    threshold = timestep*1.5
    jump_indices = where(time_diff > threshold)[0]
    return jump_indices

if args.showGaps:
    jump_indices = show_gaps(times, rates, image)


# Fill the gaps due to GTIs if fillGaps = True
def fill_gaps(times, rates, image):
    time_diff = diff(times)
    timestep = min(time_diff)
    threshold = timestep*1.5
    jump_indices = where(time_diff > threshold)[0]
    jumps = int32(time_diff[jump_indices]//timestep)

    time_intervals = split(times, jump_indices + 1)
    image_intervals = split(image, jump_indices + 1, axis=1)
    rates_intervals = split(rates, jump_indices + 1)

    FS = image.shape[0]

    nimage = image_intervals[0]
    nrates = rates_intervals[0]
    for i, time_interval in enumerate(time_intervals[1:]):
        nimage = append(nimage, ones((FS,jumps[i])), axis=1)
        nimage = append(nimage, image_intervals[i+1], axis=1)
        nrates = append(nrates, zeros(jumps[i]))
        nrates = append(nrates, rates_intervals[i+1])

    ntimes = arange(times[0],times[-1],times[-1]/nimage.shape[1])
    return ntimes, nrates, nimage

if args.fillGaps:
    times, rates, image = fill_gaps(times, rates, image)


# Prepare interpolators to handle frequency and compacted-time ticks
def freq_to_idx(yfreqs, freqs):
    return interp(yfreqs, freqs, arange(len(freqs)))

def idx_to_freq(idxs, freqs):
    return interp(idxs, arange(len(freqs)), freqs)

def time_to_idx(ytimes, times):
    return interp(ytimes, times, arange(len(times)))

def idx_to_time(idxs, times):
    return interp(idxs, arange(len(times)), times)

# Read arguments related to plot specifications
if args.pmin:
    pmin = args.pmin
else:
    pmin = quantile(image, 0.001)
    print('Using 0.1% as minimum power = ', pmin)

if args.pmax:
    pmax = args.pmax
else:
    pmax = quantile(image, 0.999)
    print('Using 99.9% as maximum power = ', pmax)

if args.tmin:
    tmin = args.tmin
else:
    tmin = times[0]

if args.tmax:
    tmax = args.tmax
else:
    tmax = times[-1]

if args.fmin:
    fmin = args.fmin
else:
    fmin = freqs[0]

if args.fmax:
    fmax = args.fmax
else:
    fmax = freqs[-1]


#Make subplots
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2.5, 1]}, constrained_layout=True)
fig.canvas.manager.set_window_title(title)
fig.align_ylabels()

# Plot dynamical PDS
if args.linear:
    pcm = ax1.imshow(image, cmap=args.cmap, norm=Norm(vmin=pmin, vmax=pmax), 
           origin='lower', interpolation=args.interpolator, aspect='auto')
else:
    pcm = ax1.imshow(image, cmap=args.cmap, norm=LogNorm(vmin=pmin, vmax=pmax), 
           origin='lower', interpolation=args.interpolator, aspect='auto')

# Create the colorbar
if args.colorbar:
    cbar = fig.colorbar(pcm, ax=ax1, orientation='vertical', pad=0.01, aspect=20)
    if args.nupnu:
        cbar.set_label("Frequency × Power (Leahy units)", fontsize=10)
    else:
        cbar.set_label("Power (Leahy units)", fontsize=10)    
    cbar.ax.tick_params(labelsize=9)


if args.tick_freqs:
#    yfreqs = array(args.tick_freqs, dtype=float)
    yfreqs = array([float(y) for y in args.tick_freqs])
    ax1.set_yticks(freq_to_idx(yfreqs, freqs), args.tick_freqs)
else:
    yfreqs = idx_to_freq(linspace(freq_to_idx(fmin,freqs),freq_to_idx(fmax,freqs),5), freqs)
    ax1.set_yticks(freq_to_idx(yfreqs, freqs),['{format_significant_3(float(yfreq))}' for yfreq in yfreqs])

if args.tick_times:
    ytimes = array(args.tick_times, dtype=float)
else:
    ytimes = idx_to_time(linspace(time_to_idx(tmin, times), time_to_idx(tmax, times), 5), times)

if args.trebin:
    trebin = args.trebin
else:
    trebin = 1

label_xaxis = 'Segment #'
if trebin > 1 :
    label_xaxis = label_xaxis + ' (rebinned by a factor of ' + str(int(trebin)) + ')'

ax1.set_xlabel(label_xaxis)
ax1.set_xlim(time_to_idx(tmin, times),time_to_idx(tmax, times))
ax1.xaxis.set_tick_params(labeltop=False)
ax1.xaxis.set_tick_params(labelbottom=True)
ax1.xaxis.set_label_position('bottom')

ax1.set_ylabel('Frequency (Hz)')
ax1.set_ylim(freq_to_idx(fmin,freqs),freq_to_idx(fmax,freqs))
#ax1.set_yticks(freq_to_idx(yfreqs, freqs),['{format_significant_3(float(yfreq))}' for yfreq in yfreqs])

# Plot lightcurve
ax2.plot(arange(len(rates)), 1e-3*rates, color='k', lw=2)

ax2.grid(zorder=-100)
ax2.set_ylabel('Rate (10$^3$ cts/s)')
nzrates = rates[rates > 0]
ax2.set_ylim(1e-3*(nzrates.min()-0.1*nzrates.mean()),1e-3*(nzrates.max()+0.1*nzrates.mean()))

if args.fillGaps:
    if args.TZERO:
        ax2.set_xlabel(f'time (s) since MJD {args.TZERO:.5f}')
    else:
        ax2.set_xlabel('time (s)')
else:
    ax2.set_xlabel('Compactified time (s)')    
    
ax2.set_xlim(time_to_idx(tmin, times),time_to_idx(tmax, times))
ax2.set_xticks(time_to_idx(ytimes, times),['{:.0f}'.format(ytime).rstrip('0').rstrip('.') for ytime in ytimes])
ax2.xaxis.set_label_position('bottom')

# Show GAPS
if args.showGaps:
    for ji in jump_indices:
        ax1.axvline(ji+0.5, c=args.highlight_color, ls='--', alpha=0.75)
        ax2.axvline(ji+0.5, c=args.highlight_color, ls='--', alpha=0.75)

# Hide the Lightcurve and Set the ticks to bottom in the dynPDS
if args.dynPDSonly:
    ax2.set_visible(False)
    ax1.xaxis.set_tick_params(labeltop=False)
    ax1.xaxis.set_tick_params(labelbottom=True)
    ax1.xaxis.set_label_position('bottom')


# Pre-format the yticks when the GUI is hidden
if args.savepdf or args.savepng:
    if args.tick_freqs:
        yfreqs = array([float(y) for y in args.tick_freqs])
        ax1.set_yticks(freq_to_idx(yfreqs, freqs), [format_significant_3(float(y)) for y in yfreqs])
    else:
        yfreqs = idx_to_freq(linspace(freq_to_idx(fmin,freqs),freq_to_idx(fmax,freqs),5), freqs)
        ax1.set_yticks(freq_to_idx(yfreqs, freqs), [format_significant_3(float(y)) for y in yfreqs])

    if args.tick_times:
        ytimes = array([float(y) for y in args.tick_times])
        ax2.set_xticks(time_to_idx(ytimes, times), [format_significant_3(float(t)) for t in ytimes])


# Funtion to save the figure as a PDF
def save_as_pdf():
    print('Exporting figure to {}.pdf'.format(args.filename.strip('.fits')))
    plt.savefig('{}.pdf'.format(args.filename.strip('.fits')), dpi=200)

# Funtion to save the figure as a PDF
def save_as_png():
    print('Exporting figure to {}.png'.format(args.filename.strip('.fits')))
    plt.savefig('{}.png'.format(args.filename.strip('.fits')), dpi=200)

# Export figure
if args.savepdf:
    save_as_pdf()
if args.savepng:
    save_as_png()

# Avoid displaying the figure on the screen
if args.noshow:
    exit(0)

# Create the main Tkinter window
root = tk.Tk()
root.title(title)

# Build the tKinter frames
ticks_frame = ttk.Frame(root)
scale_frame = ttk.Frame(root)
coord_frame = ttk.Frame(root)
button_frame = ttk.Frame(root)
ticks_frame.grid(row=0, column=0)
coord_frame.grid(row=2, column=0)
scale_frame.grid(row=3, column=0)
button_frame.grid(row=4, column=0)

def update_colormap_from_sliders(val=None):
    if args.linear:
        pmin = scale_pmin.get()
        pmax = scale_pmax.get()
    else:
        pmin = 10**scale_pmin.get()
        pmax = 10**scale_pmax.get()
    pcm.set_clim(vmin=pmin, vmax=pmax)      
    if 'cbar' in globals():
        cbar.update_normal(pcm)  # this is needed to sync the cbar
    fig.canvas.draw_idle()

# Create two Tkinter scales for the sliders
if args.linear:
    scale_pmin = tk.Scale(scale_frame, label='pmin', from_=quantile(image, 0.001), to=quantile(image, 0.999), resolution=0.01, orient=tk.HORIZONTAL, length=350)
    scale_pmax = tk.Scale(scale_frame, label='pmax', from_=quantile(image, 0.001), to=quantile(image, 0.999), resolution=0.01, orient=tk.HORIZONTAL, length=350)
else:
    scale_pmin = tk.Scale(scale_frame, label='log10(pmin)', from_=log10(quantile(image, 0.001)), to=log10(quantile(image, 0.999)), resolution=0.01, orient=tk.HORIZONTAL, length=350)
    scale_pmax = tk.Scale(scale_frame, label='log10(pmax)', from_=log10(quantile(image, 0.001)), to=log10(quantile(image, 0.999)), resolution=0.01, orient=tk.HORIZONTAL, length=350)

scale_pmin.config(command=update_colormap_from_sliders)
scale_pmax.config(command=update_colormap_from_sliders)

# Set the Tkinter scales to initial values
if args.linear:
    scale_pmin.set(pmin)
    scale_pmax.set(pmax)
else:
    scale_pmin.set(log10(pmin))
    scale_pmax.set(log10(pmax))

# Function to update the color scheme when scales are adjusted
def update_scale_pmin(event):
    if args.linear:
        pmin = scale_pmin.get()
        pmax = scale_pmax.get()
        if pmin < pmax:
            pcm.set_clim(vmin=pmin)
        else:
            pmin = pmax-0.01
            pcm.set_clim(vmin=pmin)
            scale_pmin.set(pmin)
        fig.canvas.draw_idle()

    else:
        pmin = 10**scale_pmin.get()
        pmax = 10**scale_pmax.get()
        if pmin < pmax:
            pcm.set_clim(vmin=pmin)
        else:
            pmin = pmax/1.01
            pcm.set_clim(vmin=pmin)
            scale_pmin.set(log10(pmin))
        fig.canvas.draw_idle()

def update_scale_pmax(event):
    if args.linear:
        pmin = scale_pmin.get()
        pmax = scale_pmax.get()
        if pmax > pmin:
            pcm.set_clim(vmax=pmax)
        else:
            pmax = pmin+0.01
            pcm.set_clim(vmax=pmax)
            scale_pmax.set(pmax)
        fig.canvas.draw_idle()

    else:
        pmin = 10**scale_pmin.get()
        pmax = 10**scale_pmax.get()
        if pmax > pmin:
            pcm.set_clim(vmax=pmax)
        else:
            pmax = pmin*1.01
            pcm.set_clim(vmax=pmax)
            scale_pmax.set(log10(pmax))
        fig.canvas.draw_idle()

# Bind the update functions to the scales' adjustments
scale_pmin.bind("<ButtonRelease-1>", update_scale_pmin)
scale_pmax.bind("<ButtonRelease-1>", update_scale_pmax)

# Pack the scales to the bottom of the window
scale_pmin.grid(row=0, column=1, padx=25, pady=5)
scale_pmax.grid(row=0, column=2, padx=25, pady=5)

# Create the Save PDF and PNG buttons
save_pdf_button = tk.Button(button_frame, text="Save as PDF", command=save_as_pdf)
save_pdf_button.grid(row=0, column=3, padx=10, pady=2)
save_png_button = tk.Button(button_frame, text="Save as PNG", command=save_as_png)
save_png_button.grid(row=0, column=4, padx=10, pady=2)

# Update ticks and highlights
def update_ticks():
    global highlight_freq_lines, highlight_time_lines_ax1, highlight_time_lines_ax2, highlight_texts

    # Frequency ticks
    tbtf = text_box_tf.get("1.0", tk.END).strip()
    if tbtf:
        try:      
            tick_labels_tf = tbtf.split()
            yfreqs = array([float(y) for y in tick_labels_tf])
            ax1.set_yticks(freq_to_idx(yfreqs, freqs), tick_labels_tf)
        except Exception as e:
            print("Error parsing frequency ticks:", e)

    # Frequency range
    tbts = text_box_fs.get("1.0", tk.END).strip()
    try:
        fs = array(tbts.split(), dtype=float)
        fmin = fs[0]
        fmax = fs[-1]
        ax1.set_ylim(freq_to_idx(fs[0],freqs),freq_to_idx(fs[-1],freqs))
    except Exception as e:
        print("Error parsing frequency range:", e)

    # Time ticks
    tbtt = text_box_tt.get("1.0", tk.END).strip()
    if tbtt:
        try:
            tick_labels_tt = tbtt.split()
            ytimes = array([float(y) for y in tick_labels_tt])
            ax2.set_xticks(time_to_idx(ytimes, times), tick_labels_tt)
        except Exception as e:
            print("Error parsing time ticks:", e)

    # === Clear previous highlight lines and texts ===
    for item in highlight_freq_lines + highlight_time_lines_ax1 + highlight_time_lines_ax2 + highlight_texts:
        item.remove()
    highlight_freq_lines.clear()
    highlight_time_lines_ax1.clear()
    highlight_time_lines_ax2.clear()
    highlight_texts.clear()

    # Highlight frequencies
    hf_text = text_box_hf.get("1.0", tk.END).strip()
    if hf_text:
        try:
            highlight_freqs = array(hf_text.split(), dtype=float)
            for f in highlight_freqs:
                y = freq_to_idx(f, freqs)
                line = ax1.axhline(y, c=args.highlight_color, ls='--', alpha=0.75)
                text = ax1.text(len(times), y, s=f'{format_significant_3(f)} Hz', c=args.highlight_color, alpha=1,
                                ha='left', va='center', fontsize=9)
                highlight_freq_lines.append(line)
                highlight_texts.append(text)
        except Exception as e:
            print("Error parsing highlight freqs:", e)

    # Highlight times
    ht_text = text_box_ht.get("1.0", tk.END).strip()
    if ht_text:
        try:
            highlight_times = array(ht_text.split(), dtype=float)
            for t in highlight_times:
                x = time_to_idx(t, times)
                line1 = ax1.axvline(x, c=args.highlight_color, ls='--', alpha=0.75)
                line2 = ax2.axvline(x, c=args.highlight_color, ls='--', alpha=0.75)
                text = ax1.text(x, freq_to_idx(fmax, freqs), s=f' {format_significant_3(t/1000)} ks',
                                c=args.highlight_color, alpha=1, rotation=90,
                                ha='center', va='bottom', fontsize=9)
                highlight_time_lines_ax1.append(line1)
                highlight_time_lines_ax2.append(line2)
                highlight_texts.append(text)
        except Exception as e:
            print("Error parsing highlight times:", e)

    fig.canvas.draw_idle()


def clear_all_highlights():
    global highlight_freq_lines, highlight_time_lines_ax1, highlight_time_lines_ax2, highlight_texts, consolidated_lines

    # Remove lines and dynamic texts
    for item in highlight_freq_lines + highlight_time_lines_ax1 + highlight_time_lines_ax2 + highlight_texts:
        item.remove()
    highlight_freq_lines.clear()
    highlight_time_lines_ax1.clear()
    highlight_time_lines_ax2.clear()
    highlight_texts.clear()

    # Remove crosshairs by click
    for entry in consolidated_lines:
        for item in entry[2:]:  # v_ax1, h_ax1, v_ax2, text_f, text_t
            item.remove()
    consolidated_lines.clear()

    # Clean up the textboxes
    text_box_hf.delete('1.0', tk.END)
    text_box_ht.delete('1.0', tk.END)

    fig.canvas.draw_idle()


# Create the textboxes to introduce frequency ticks
text_box_tf_label = tk.Label(ticks_frame, text="Freq. ticks (Hz):")
text_box_tf_label.grid(row=0, column=0, padx=5, pady=5)
text_box_tf = tk.Text(ticks_frame, height=1, width=40)
#text_box_tf.insert(tk.END," ".join(array2string(yfreqs, precision=2, formatter={'float_kind':lambda x: "%.2f" % x}).lstrip("[").rstrip("]").split()))
text_box_tf.insert(tk.END," ".join(array2string(yfreqs, precision=2, formatter={'float_kind':lambda x: "%.2g" % x}).lstrip("[").rstrip("]").split()))
text_box_tf.grid(row=0, column=1, padx=5, pady=5)

# Create the textboxes to introduce time ticks
text_box_tt_label = tk.Label(ticks_frame, text="Time ticks (s):")
text_box_tt_label.grid(row=1, column=0, padx=5, pady=5)
text_box_tt = tk.Text(ticks_frame, height=1, width=40)
text_box_tt.insert(tk.END,array2string(ytimes, formatter={'float_kind':lambda x: "%.0f" % x}).lstrip("[").rstrip("]"))
text_box_tt.grid(row=1, column=1, padx=5, pady=5)

# Highlighted Frequencies
text_box_hf_label = tk.Label(ticks_frame, text="Highlight Freqs. (Hz):")
text_box_hf_label.grid(row=0, column=4, padx=5, pady=5)
text_box_hf = tk.Text(ticks_frame, height=1, width=25)
text_box_hf.grid(row=0, column=5, padx=5, pady=5)

# Highlighted Times
text_box_ht_label = tk.Label(ticks_frame, text="Highlight Times (s):")
text_box_ht_label.grid(row=1, column=4, padx=5, pady=5)
text_box_ht = tk.Text(ticks_frame, height=1, width=25)
text_box_ht.grid(row=1, column=5, padx=5, pady=5)

# Fill-up the highlighted textboxes with the passed arguments
if args.highlight_freqs:
    hf_str = " ".join(f"{format_significant_3(float(f))}" for f in args.highlight_freqs)
    text_box_hf.insert(tk.END, hf_str)

if args.highlight_times:
    ht_str = " ".join(f"{format_significant_3(float(t))}" for t in args.highlight_times)
    text_box_ht.insert(tk.END, ht_str)

# Create the textboxes to introduce frequency range
text_box_fs_label = tk.Label(button_frame, text="Freq. Min Max (Hz):")
text_box_fs_label.grid(row=0, column=0, padx=5, pady=5)
text_box_fs = tk.Text(button_frame, height=1, width=30)
text_box_fs.insert(tk.END," ".join(array2string(array([fmin,fmax]), precision=2, formatter={'float_kind':lambda x: "%.2f" % x}).lstrip("[").rstrip("]").split()))
text_box_fs.grid(row=0, column=1, padx=5, pady=5)

# Bind <Return> and <FocusOut> to update_ticks for all text boxes
for textbox in [text_box_tf, text_box_tt, text_box_fs, text_box_hf, text_box_ht]:
    textbox.bind("<Return>", lambda event: update_ticks())
    textbox.bind("<FocusOut>", lambda event: update_ticks())

# Create the UPDATE button to update the plot
read_button = tk.Button(button_frame, text='UPDATE', font='serif 11 bold', command=update_ticks)
read_button.grid(row=0, column=5, padx=5, pady=5)

# Create the CLEAR button to clean the highlights of the plot
clear_button = tk.Button(button_frame, text='CLEAR', font='serif 11 bold', command=clear_all_highlights)
clear_button.grid(row=0, column=6, padx=5, pady=5)

# Create the canvas for the plot
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=1,column=0)

# Function to close the plot and release resources
def on_closing():
    exit(0)

# Bind the function to the window's close event
root.protocol("WM_DELETE_WINDOW", on_closing)

# Print Mouse Coordinates
coord_label = tk.Label(coord_frame, text="PDS# = ---   ν = ---   P = ---", font='bold')
coord_label.grid(row=0, column=0, sticky='w', padx=10, pady=5)

# Mouse Crosshair
crosshair_vline = ax1.axvline(0, color=args.highlight_color, linestyle='--', lw=2, alpha=0.6, visible=False)
crosshair_hline = ax1.axhline(0, color=args.highlight_color, linestyle='--', lw=2, alpha=0.6, visible=False)
crosshair_vline_ax2 = ax2.axvline(0, color='C3', linestyle='--', lw=2, alpha=0.6, visible=False)

# Mouse motion control: update crosshair lines
def on_mouse_move(event):
    if event.inaxes == ax1:
        x = event.xdata
        y = event.ydata
        if x is not None and y is not None:
            t = idx_to_time(x, times)
            f = idx_to_freq(y, freqs)
            
            ix = int(x)
            iy = int(y)
            if 0 <= ix < image.shape[1] and 0 <= iy < image.shape[0]:
                value = image[iy, ix]
                if args.nupnu:
                    value_str = f"νPν = {value:.2f}"
                else:
                    value_str = f"P = {value:.2f}"
            else:
                value_str = "P = ---"

            coord_label.config(text=f"PDS# = {x:.2f} ({t:.0f} s)    ν = {f:.2f} Hz    {value_str}", font='bold')

#FG         coord_label.config(text="PDS# = {:.2f} ({:.0f} s)    Freq = {:.2f} Hz".format(x, t, f), font='bold')
            crosshair_vline.set_xdata([x])
            crosshair_hline.set_ydata([y])
            crosshair_vline_ax2.set_xdata([x])
            crosshair_vline.set_visible(True)
            crosshair_hline.set_visible(True)
            crosshair_vline_ax2.set_visible(True)
            fig.canvas.draw_idle()
    else:
        crosshair_vline.set_visible(False)
        crosshair_hline.set_visible(False)
        crosshair_vline_ax2.set_visible(False)
        coord_label.config(text="")
        fig.canvas.draw_idle()

# Consolidated crosshair lines
consolidated_lines = []  # lista de tuplas: (x, y, vline_ax1, hline_ax1, vline_ax2)

highlight_freq_lines = []
highlight_time_lines_ax1 = []
highlight_time_lines_ax2 = []
highlight_texts = []

# Mouse click control: consolidate or remove crosshair highlights
def on_click(event):
    if event.inaxes == ax1:
        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            return

        # Check if there is a crosshair closeby
        to_remove = None
        for entry in consolidated_lines:
            x0, y0 = entry[0], entry[1]
            if abs(x - x0) < 2 and abs(y - y0) < 0.3:
                to_remove = entry
                break

        if to_remove:
            # Remove lines and texts
            to_remove[2].remove()  # v_ax1
            to_remove[3].remove()  # h_ax1
            to_remove[4].remove()  # v_ax2
            to_remove[5].remove()  # text_f
            to_remove[6].remove()  # text_t
            consolidated_lines.remove(to_remove)
        else:
            # Create lines
            v_ax1 = ax1.axvline(x, color=args.highlight_color, linestyle='--')
            h_ax1 = ax1.axhline(y, color=args.highlight_color, linestyle='--')
            v_ax2 = ax2.axvline(x, color=args.highlight_color, linestyle='--')

            fmax = idx_to_freq(ax1.get_ylim()[1], freqs)
            tmax = idx_to_time(ax1.get_xlim()[1], times)

            # Add Frequency label (right)                        
            f_label = idx_to_freq(y, freqs)
            text_f = ax1.text(time_to_idx(tmax, times), y, s=f' {format_significant_3(float(f_label))} Hz',
                                c=args.highlight_color, alpha=1, ha='left', va='center', fontsize=9)
                                
            # Add Time label (top)           
            t_label = idx_to_time(x, times)
            text_t = ax1.text(x, freq_to_idx(fmax, freqs), s=f'  {format_significant_3(t_label/1000)} ks',
                  c=args.highlight_color, alpha=1, rotation=90, ha='center', va='bottom', fontsize=9)

            consolidated_lines.append((x, y, v_ax1, h_ax1, v_ax2, text_f, text_t))

        fig.canvas.draw_idle()

    # Updated textbox of highlighted frequencies
    fh_list = text_box_hf.get("1.0", tk.END).strip().split()
    f_str = format_significant_3(f_label)
    if f_str not in fh_list:
        fh_list.append(f_str)
        fh_list_sorted = sorted(set(float(f) for f in fh_list))
        text_box_hf.delete("1.0", tk.END)
        text_box_hf.insert(tk.END, " ".join(format_significant_3(f) for f in fh_list_sorted))

    # Updated textbox of highlighted times
    th_list = text_box_ht.get("1.0", tk.END).strip().split()
    t_str = format_significant_3(t_label)
    if t_str not in th_list:
        th_list.append(t_str)
        th_list_sorted = sorted(set(float(t) for t in th_list))
        text_box_ht.delete("1.0", tk.END)
        text_box_ht.insert(tk.END, " ".join(format_significant_3(t) for t in th_list_sorted))


# Connect Mouse events to the canvas
fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)
fig.canvas.mpl_connect('button_press_event', on_click)

# Start the Tkinter event loop
update_ticks()
root.mainloop()
