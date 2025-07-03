import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from rtlsdr import RtlSdr
from numpy.fft import fft, fftshift
from datetime import datetime, timedelta
import csv
import os
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, Galactic
from astropy.time import Time
import astropy.units as u

# Constants from low-cost 21 cm paper
R0 = 8.5            # kpc (Sun to Galactic Center)
V0 = 220            # km/s (Solar orbital velocity)
F0 = 1420.405751    # MHz (21cm rest frequency)
C = 299792.458      # km/s (speed of light)
A = 14.8            # km/s/kpc (Oort A)
B = A - (V0 / R0)   # Oort B
def capture_rtl_data(center_freq, num_samples, chunk_size=1024000):
    sdr = RtlSdr()
    sdr.sample_rate = 2.4e6
    sdr.center_freq = center_freq
    sdr.gain = 'auto'
    samples = np.array([], dtype=np.complex64)
    try:
        for _ in range(num_samples // chunk_size):
            samples = np.concatenate((samples, sdr.read_samples(chunk_size)))
    finally:
        sdr.close()
    return samples

def return_averaged_spectras(ipdata, chNo, nAvgerages, nSets, npt):
    ipdata = np.asarray(ipdata)
    row = len(ipdata)
    totalSpectrasAsked = nAvgerages * nSets
    NoOfAvailableSpectras = (row // npt) - 1

    if totalSpectrasAsked <= NoOfAvailableSpectras:
        aspecA = np.zeros((npt, nSets))
        for set_idx in range(nSets):
            for I in range(nAvgerages):
                startNo = (I * npt) + (set_idx * nAvgerages * npt) + 1
                endNo = startNo + npt
                if endNo > len(ipdata): 
                    break
                segment = ipdata[startNo:endNo] if ipdata.ndim == 1 else ipdata[startNo:endNo, chNo]
                spc = np.abs(fft(segment)) ** 2
                aspecA[:, set_idx] += spc[:npt]
            aspecA[:, set_idx] = fftshift(aspecA[:, set_idx] / nAvgerages)
        return aspecA
    else:
        st.error("Error: Not enough data for averaging.")
        return np.zeros((2, 1))

def process_data(aa, bb, colNo=1, nfft=128):
    navg = len(aa) // nfft - 3
    nsets = 1
    avgps = return_averaged_spectras(aa, colNo, navg, nsets, nfft)
    avgps2 = return_averaged_spectras(bb, colNo, navg, nsets, nfft)
    return avgps, avgps2

def calculate_velocity_and_distance(obs_freq_mhz, L):
    L_rad = np.radians(L)
    Vr = ((obs_freq_mhz / F0) * (C)) - C
    denominator = A * np.sin(2 * L_rad)
    d = Vr / denominator if abs(denominator) > 1e-6 else np.nan
    Vt = d * (A * np.cos(2 * L_rad) + B)
    Ur = Vr + V0 * np.sin(L_rad)
    Ut = Vt + V0 * np.cos(L_rad)
    V = np.sqrt(Ur**2 + Ut**2)
    R = np.sqrt(R0**2 + d**2 - 2*R0*d*np.cos(L_rad))
    return {
        'Vr': round(Vr, 1),
        'd': round(d, 1),
        'R': round(R, 1),
        'Ur': round(Ur, 1),
        'Ut': round(Ut, 1),
        'V': round(V, 1)
    }

def create_plots(avgps, avgps2, num, L, nfft=128):
    Fcenter = 1420.405751
    freq_hz = np.linspace(Fcenter*1e6 - 1e6, Fcenter*1e6 + 1e6, nfft)
    diff_spec = avgps2 - avgps
    peak = np.argmax(diff_spec)
    obs_freq_mhz = freq_hz[peak] / 1e6

    results = calculate_velocity_and_distance(obs_freq_mhz, L)

    plt.figure(figsize=(10, 10))
    plt.subplot(3, 1, 1)
    plt.plot(freq_hz/1e6, avgps, 'b-', label='1422 MHz')
    plt.plot(freq_hz/1e6, avgps2, 'r--', label='1420 MHz')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power in counts')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.subplot(3, 1, 2)
    plt.plot(freq_hz/1e6, diff_spec, 'g-')
    plt.axvline(Fcenter, color='k', linestyle='--')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Signal Power in counts')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.subplot(3, 1, 3)
    plt.axis('off')
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    result_text = (
    f"Observation Time: {timestamp}\n"
    f"Galactic Longitude (l): {L:.2f}°\n"
    f"Observed Frequency: {obs_freq_mhz:.6f} MHz\n"
    f"Vr: {results['Vr']} km/s\n"
    f"Ur: {results['Ur']} km/s\n"
    f"Ut: {results['Ut']} km/s\n"
    f"Total V: {results['V']} km/s\n"
    f"Galactic center Distance R: {results['R']} kpc"
    )

    plt.text(0.1, 0.5, result_text, fontfamily='monospace')
    plt.tight_layout()
    # Save plot with timestamp in filename, formatted as "YYYY-MM-DD at HH.MM.SS"
    plot_time = datetime.now().strftime("%Y-%m-%d at %H.%M.%S")
    filename = f'plot_trial_{num}_{plot_time}.png'
    plt.savefig(filename, dpi=100)
    plt.close()

    return filename, freq_hz, results, obs_freq_mhz

def save_data(avgps, avgps2, freq, results, trial_num, L, obs_freq, write_header=False, B=None, RA=None, DEC=None, filename=None):
    if filename is None:
        filename = csv_filename
    os.makedirs("data", exist_ok=True)
    mode = 'w' if write_header else 'a'
    with open(filename, mode, newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["Parameter", "Value", "Unit"])

        writer.writerow([f"Trial Number", trial_num, ""])
        writer.writerow(["Galactic Longitude (L)", L, "degrees"])
        if B is not None:
            writer.writerow(["Galactic Latitude (B)", B, "degrees"])
        if RA is not None:
            writer.writerow(["Right Ascension (RA)", RA, "degrees"])
        if DEC is not None:
            writer.writerow(["Declination (DEC)", DEC, "degrees"])
        writer.writerow(["observed Frequency", obs_freq, "MHz"])
        writer.writerow(["Observation Time", datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ""])
        writer.writerow(["Calculated Parameters"])
        writer.writerow(["Vr", results['Vr'], "km/s"])
        writer.writerow(["Ur", results['Ur'], "km/s"])
        writer.writerow(["Ut", results['Ut'], "km/s"])
        writer.writerow(["Total Velocity", results['V'], "km/s"])
        writer.writerow(["Distance from Sun", results['d'], "kpc"])
        writer.writerow(["Galactic center Distance", results['R'], "kpc"])
        writer.writerow(["Spectral Data"])
        writer.writerow(["Frequency (MHz)", "avgps", "avgps2"])

        for f_val, on, off in zip(freq/1e6, avgps.flatten(), avgps2.flatten()):
            writer.writerow([round(f_val, 4), round(on, 2), round(off, 2)])

        writer.writerow([])
        writer.writerow([])

def save_summary_csv(trial_num, obs_freq, results, L, B=None, RA=None, DEC=None, filename=None):
    if filename is None:
        filename = summary_filename
    os.makedirs("data", exist_ok=True)
    file_exists = os.path.isfile(filename)
    with open(filename, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow([
                "Trial Number", "Observed Frequency (MHz)", "Galactic Longitude (L)", "Galactic Latitude (B)",
                "Right Ascension (RA)", "Declination (DEC)",
                "Vr (km/s)", "Ur (km/s)", "Ut (km/s)", "Total V (km/s)", "Distance from Sun (kpc)", "Galactic Center Distance (kpc)"
            ])
        writer.writerow([
            trial_num, obs_freq, L,
            B if B is not None else "",
            RA if RA is not None else "",
            DEC if DEC is not None else "",
            results.get('Vr', ""), results.get('Ur', ""), results.get('Ut', ""), results.get('V', ""),
            results.get('d', ""), results.get('R', "")
        ])

st.title('Hydrogen Line Observation')
st.markdown("""
This application processes radio telescope data to observe the 21cm hydrogen line 
and calculate galactic rotation parameters.
""")
plot_spot = st.empty()
num_samples = st.sidebar.number_input('Samples per Capture', min_value=512000, max_value=20480000, value=20480000, step=512000)
sample_rate = 2.4e6
num_trials = st.sidebar.number_input('Number of Trials', 1, 20000, 1)
nfft = st.sidebar.selectbox('FFT Size', [64, 128, 256, 512], index=1)

# Replace this line:
# st.experimental_set_query_params()  # Ensures Streamlit session state is initialized
st.query_params.clear()  # Ensures Streamlit session state is initialized

# Add this near the top, before coordinate calculation
st_autorefresh = st.sidebar.checkbox("Auto-refresh coordinates", value=True)
if st_autorefresh:
    st.experimental_rerun = st.experimental_rerun if hasattr(st, "experimental_rerun") else lambda: None
    st_autorefresh_interval = st.sidebar.slider("Refresh interval (seconds)", 1, 60, 5)
    st_autorefresh_count = st.query_params.get("autorefresh_count", [0])
    st_autorefresh_count = int(st_autorefresh_count[0]) + 1
    import time
    if "last_autorefresh" not in st.session_state or \
       (datetime.now() - st.session_state.get("last_autorefresh", datetime.min)).total_seconds() > st_autorefresh_interval:
        st.session_state["last_autorefresh"] = datetime.now()
        st.query_params["autorefresh_count"] = st_autorefresh_count
        st.experimental_rerun()

def display_realtime_coordinates(alt_key="alt_slider", az_key="az_slider", show_sidebar=True):
    # Gauribidanur Radio Observatory coordinates
    latitude = 13.6029     # degrees North
    longitude = 77.4390    # degrees East
    elevation = 686.0      # meters

    if show_sidebar:
        st.sidebar.write(f"Latitude: {latitude}° N")
        st.sidebar.write(f"Longitude: {longitude}° E")
        st.sidebar.write(f"Elevation: {elevation} m")
        alt = st.sidebar.slider("Antenna Altitude (°)", 0.0, 90.0, 89.0, key=alt_key)
        az = st.sidebar.slider("Antenna Azimuth (°)", 0.0, 360.0, 180.0, key=az_key)
        # Sidebar live clock (auto-refresh with Streamlit rerun)
        clock_placeholder = st.sidebar.empty()
        now = datetime.now()  # Use local system time
        clock_placeholder.markdown(f"###  {now.strftime('%I:%M:%S %p')}")
        # Compute galactic longitude L from pointing
        obs_time = Time(datetime.utcnow()) # + timedelta(hours=5, minutes=30)
        obs_time_IST = obs_time + timedelta(hours=5, minutes=30)
        st.sidebar.write(obs_time_IST.strftime("%A, %B, %d, %Y, %I:%M%p IST"))
        location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=elevation*u.m)
        altaz = AltAz(obstime=obs_time, location=location)
        sky_coord = SkyCoord(alt=alt*u.deg, az=az*u.deg, frame=altaz)
        icrs_coord = sky_coord.transform_to('icrs')
        galactic_coord = sky_coord.transform_to('galactic')
        L = galactic_coord.l.deg
        B = galactic_coord.b.deg
        RA = icrs_coord.ra.deg
        DEC = icrs_coord.dec.deg
        RA_hms = icrs_coord.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
        st.sidebar.write(f"Calculated Right Ascension: **{RA:.2f}°** ({RA_hms} hms)")
        st.sidebar.write(f"Calculated Declination: **{DEC:.2f}°**")
        st.sidebar.write(f"Calculated Galactic Latitude: **{B:.2f}°**")
        st.sidebar.write(f"Calculated Galactic Longitude: **{L:.2f}°**")
    else:
        # Only compute, do not print to sidebar
        alt = st.session_state.get(alt_key, 89.0)
        az = st.session_state.get(az_key, 180.0)
        obs_time = Time(datetime.utcnow())
        location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=elevation*u.m)
        altaz = AltAz(obstime=obs_time, location=location)
        sky_coord = SkyCoord(alt=alt*u.deg, az=az*u.deg, frame=altaz)
        icrs_coord = sky_coord.transform_to('icrs')
        galactic_coord = sky_coord.transform_to('galactic')
        L = galactic_coord.l.deg
        B = galactic_coord.b.deg
        RA = icrs_coord.ra.deg
        DEC = icrs_coord.dec.deg
    return L, B, RA, DEC

# Initial sidebar display (static keys)
L, B, RA, DEC = display_realtime_coordinates()

if st.button('Start Observation Session'):
    status_text = st.empty()
    data_folder = "data"
    os.makedirs(data_folder, exist_ok=True)

    # Generate unique filenames with timestamp for each session, formatted as "YYYY-MM-DD at HH.MM.SS"
    session_time = datetime.now().strftime("%Y-%m-%d started at %H.%M.%S")
    csv_filename = os.path.join(data_folder, f"trials_results_{session_time}.csv")
    summary_filename = os.path.join(data_folder, f"trials_summary_{session_time}.csv")

    # Define save_data and save_summary_csv inside the button block so csv_filename and summary_filename are in scope
    def save_data(avgps, avgps2, freq, results, trial_num, L, obs_freq, write_header=False, B=None, RA=None, DEC=None, filename=None):
        if filename is None:
            filename = csv_filename
        os.makedirs("data", exist_ok=True)
        mode = 'w' if write_header else 'a'
        with open(filename, mode, newline='') as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow(["Parameter", "Value", "Unit"])

            writer.writerow([f"Trial Number", trial_num, ""])
            writer.writerow(["Galactic Longitude (L)", L, "degrees"])
            if B is not None:
                writer.writerow(["Galactic Latitude (B)", B, "degrees"])
            if RA is not None:
                writer.writerow(["Right Ascension (RA)", RA, "degrees"])
            if DEC is not None:
                writer.writerow(["Declination (DEC)", DEC, "degrees"])
            writer.writerow(["observed Frequency", obs_freq, "MHz"])
            writer.writerow(["Observation Time", datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ""])
            writer.writerow(["Calculated Parameters"])
            writer.writerow(["Vr", results['Vr'], "km/s"])
            writer.writerow(["Ur", results['Ur'], "km/s"])
            writer.writerow(["Ut", results['Ut'], "km/s"])
            writer.writerow(["Total Velocity", results['V'], "km/s"])
            writer.writerow(["Distance from Sun", results['d'], "kpc"])
            writer.writerow(["Galactic center Distance", results['R'], "kpc"])
            writer.writerow(["Spectral Data"])
            writer.writerow(["Frequency (MHz)", "avgps", "avgps2"])

            for f_val, on, off in zip(freq/1e6, avgps.flatten(), avgps2.flatten()):
                writer.writerow([round(f_val, 4), round(on, 2), round(off, 2)])

            writer.writerow([])
            writer.writerow([])

    def save_summary_csv(trial_num, obs_freq, results, L, B=None, RA=None, DEC=None, filename=None):
        if filename is None:
            filename = summary_filename
        os.makedirs("data", exist_ok=True)
        file_exists = os.path.isfile(filename)
        with open(filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow([
                    "Trial Number", "Observed Frequency (MHz)", "Galactic Longitude (L)", "Galactic Latitude (B)",
                    "Right Ascension (RA)", "Declination (DEC)",
                    "Vr (km/s)", "Ur (km/s)", "Ut (km/s)", "Total V (km/s)", "Distance from Sun (kpc)", "Galactic Center Distance (kpc)"
                ])
            writer.writerow([
                trial_num, obs_freq, L,
                B if B is not None else "",
                RA if RA is not None else "",
                DEC if DEC is not None else "",
                results.get('Vr', ""), results.get('Ur', ""), results.get('Ut', ""), results.get('V', ""),
                results.get('d', ""), results.get('R', "")
            ])

    for trial in range(1, num_trials + 1):
        try:
            # Use unique keys for each trial to avoid duplicate key error, and do not print to sidebar
            L, B, RA, DEC = display_realtime_coordinates(
                alt_key=f"alt_slider_{trial}", az_key=f"az_slider_{trial}", show_sidebar=False
            )
            status_text.text(f"Running Trial {trial}/{num_trials}...")
            aa = capture_rtl_data(1422000000, num_samples)
            bb = capture_rtl_data(1420405751, num_samples)
            avgps, avgps2 = process_data(aa, bb, nfft=nfft)
            plot_filename, freq_hz, results, obs_freq = create_plots(avgps, avgps2, trial, L, nfft)
            write_header = (trial == 1)
            # Pass avgps and avgps2 as numpy arrays (flatten if needed)
            save_data(avgps, avgps2, freq_hz, results, trial, L, obs_freq, write_header, B=B, RA=RA, DEC=DEC, filename=csv_filename)
            save_summary_csv(trial, obs_freq, results, L, B=B, RA=RA, DEC=DEC, filename=summary_filename)
            with plot_spot:
                st.image(plot_filename, caption=f'Trial {trial} Results')
            st.session_state.trial_num = trial
        except Exception as e:
            st.error(f"Error in trial {trial}: {str(e)}")
            break

    st.success("Observation session completed successfully!")

