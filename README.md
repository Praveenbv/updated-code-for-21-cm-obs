🌌 Hydrogen Line Observation Tool (21cm) with RTL-SDR
This Streamlit application allows real-time capture, processing, and analysis of 21 cm hydrogen line signals using an RTL-SDR dongle. It computes galactic rotation parameters based on spectral shifts, aiding studies of the Milky Way's dynamics.

📡 Features
🔭 Real-time data acquisition from RTL-SDR

📉 FFT-based power spectrum processing and comparison at 1422 MHz and 1420 MHz

🌌 Galactic coordinate and velocity calculation using astropy

📍 Interactive sky-pointing (azimuth & altitude) with live RA/DEC and galactic coordinates

📈 Visualization of hydrogen line data and spectral differences

📁 Automatic CSV logging of detailed trial data and summarized results

💾 Plot export with timestamped filenames

🧠 How It Works
Data Capture: Collects samples from RTL-SDR centered at two frequencies: 1422 MHz (baseline) and 1420.405751 MHz (expected hydrogen line).

Spectral Analysis: Performs FFT, averages spectra, and computes power differences.

Astro Calculations: Uses Doppler shift to calculate:

Radial velocity (Vr)

Total velocity (V)

Distance from Sun (d)

Distance to Galactic Center (R)

Visualization: Displays plots for both channels, difference spectrum, and computed astrophysical parameters.

Logging: Saves detailed and summary trial results in CSV format inside a data/ directory.

📦 Requirements
Install the following packages (preferably in a virtual environment):

bash
Copy
Edit
pip install streamlit numpy matplotlib pyrtlsdr astropy
You will also need:

A working RTL-SDR USB device

A suitable 21 cm (1420 MHz) antenna (e.g., horn, Yagi, or helical)

🚀 Running the App
bash
Copy
Edit
streamlit run astrop.py
🛠️ Sidebar Controls
Samples per Capture: Number of samples collected per frequency

FFT Size: Frequency resolution

Number of Trials: Total repetitions of observation cycle

Antenna Pointing: Altitude and Azimuth control

Auto-refresh Coordinates: Updates RA/DEC and Galactic coordinates continuously

📂 Outputs
All outputs are stored in the data/ folder:

trials_results_<timestamp>.csv – Raw and processed data for each trial

trials_summary_<timestamp>.csv – Summary of galactic parameters

plot_trial_<num>_<timestamp>.png – Spectral plots for each trial

🌍 Location Default
Observatory Location: Gauribidanur Radio Observatory (Latitude: 13.6029°N, Longitude: 77.4390°E, Elevation: 686m)

📚 References
Velocity and distance formulas are based on simplified galactic rotation models from low-cost 21cm observation literature

Uses the IAU 1958 recommended value of the rest frequency for hydrogen: 1420.405751 MH

🤝 Acknowledgments
Inspired by open-source efforts in low-cost radio astronomy

Built with ❤️ using Streamlit, Astropy, and RTL-SDR

