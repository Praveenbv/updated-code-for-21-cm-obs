Hydrogen Line Observation Streamlit App
This repository contains a Python Streamlit application designed to capture and analyze radio telescope data for observing the 21cm hydrogen line and performing galactic rotation calculations.

Features
Data Acquisition: Captures radio frequency data from RTL-SDR hardware.

Signal Processing: Averages power spectra for improved signal-to-noise ratio and detection.

Astrophysical Calculations: Computes galactic velocities and distances using measurements of the 21cm hydrogen line.

Data Visualization: Produces detailed plots of power spectra and a summary of calculated parameters.

Logging: Saves comprehensive data logs and summaries in CSV format for post-observation analysis.

Real-time Coordinates: Provides real-time astronomical coordinate calculations and visualizations based on the antenna's pointing direction.

Installation
Ensure you have Python installed on your system. You can then install the necessary dependencies using pip:

bash
pip install streamlit numpy matplotlib pyrtlsdr astropy
Usage
To launch the application, navigate to the repository's directory in your terminal and run the following command:

bash
streamlit run astrop.py
After running the command, a web interface will open in your browser. You can configure the observation parameters in the sidebar, including:

Number of samples per capture

Number of observation trials

FFT size

Antenna altitude and azimuth

Once configured, click the "Start Observation Session" button to begin data capture and analysis.

Detailed Functionality
Data Capture
The capture rtldata function interfaces with the connected RTL-SDR device to capture raw I/Q (in-phase and quadrature) samples based on the user-defined parameters.

Spectral Averaging
The returnaveragedspectras function takes the raw time-domain data, segments it, and computes an averaged FFT power spectrum. This process is crucial for reducing noise and making the faint hydrogen line signal detectable.

Galactic Calculations
The calculatevelocityanddistance function is the core of the scientific analysis. It uses the observed frequency of the hydrogen line and the galactic longitude of the observation to calculate key kinematic parameters, including radial velocities, distances from the Sun, and distances from the Galactic center.

Visualization
The createplots function generates a multi-panel plot for each observation trial. This plot displays the observed power spectra, the difference spectrum (to isolate the signal), and a text summary of the calculated galactic parameters.

Data Output
All data and plots are automatically saved to a data/ directory. For each observation trial, the application generates:

A CSV file containing the raw spectral data and all calculated parameters.

A PNG image of the plot for quick visual analysis.

Notes
The application is specifically configured to target the 1420 MHz neutral hydrogen line, which is used to study the structure and rotation of our Milky Way galaxy.

The script is pre-configured with the coordinates of the Gauribidanur Radio Observatory in Karnataka, India, for accurate astronomical calculations.
