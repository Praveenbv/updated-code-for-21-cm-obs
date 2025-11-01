# Hydrogen Line Observation Streamlit App

This repository contains a Streamlit-based Python application for capturing and analyzing radio telescope data to detect the 21 cm hydrogen line (1420 MHz) and compute galactic rotation parameters.

---

## Overview

This tool interfaces with an RTL-SDR receiver to capture raw I/Q radio data, perform FFT-based spectral averaging, and analyze the neutral hydrogen emission line to determine galactic velocities and distances.

Built using:
- Python
- Streamlit
- Astropy
- RTL-SDR

---

## Features

- Real-time Data Capture: Acquire radio signals from the connected RTL-SDR receiver.
- Spectral Averaging: Enhances signal-to-noise ratio for detecting the faint 21 cm line.
- Galactic Calculations: Computes radial velocity, galactocentric distance, and related kinematic parameters.
- Interactive UI: Configure sampling rate, FFT size, and antenna orientation via a simple web interface.
- Coordinate Tracking: Displays real-time Galactic and Equatorial coordinates based on antenna pointing.
- Data Logging: Saves all trials in organized CSV and PNG formats for later analysis.

---

## Installation

Ensure Python version 3.8 or higher is installed.

Install dependencies:

```bash
pip install streamlit numpy matplotlib pyrtlsdr astropy

Hydrogen Line Observation Streamlit App
Usage
-----
Clone this repository:
git clone https://github.com/yourusername/hydrogen-line-streamlit.git
cd hydrogen-line-streamlit
Run the app:
streamlit run astrop.py
Configure parameters in the Streamlit sidebar:
- Samples per Capture
- Number of Trials
- FFT Size
- Antenna Altitude and Azimuth
Click "Start Observation Session" to begin capturing and analyzing data.
Results are displayed in the main window and automatically saved in the /data/ folder.

Core Functionalities
--------------------
Function | Purpose
----------|----------
capture_rtl_data() | Acquires I/Q samples from the RTL-SDR.
return_averaged_spectras() | Computes averaged FFT spectra for noise reduction.
calculate_velocity_and_distance() | Derives galactic velocity, distance from Sun, and Galactic
center.
create_plots() | Generates and saves spectrum plots with overlaid results.
display_realtime_coordinates() | Computes real-time RA, DEC, Galactic L & B based on antenna
direction.

Output Data
------------
All outputs are saved in a timestamped folder under /data/:
- trials_results_YYYY-MM-DD.csv – Full spectral and computed data per trial
- trials_summary_YYYY-MM-DD.csv – Summary of galactic parameters
- plot_trial_X.png – Spectral plots with hydrogen line detections
Observatory Location
--------------------
Default site coordinates (Gauribidanur Radio Observatory, India):
- Latitude: 13.6029° N
- Longitude: 77.4390° E
- Elevation: 686 m
These coordinates are used for accurate astronomical calculations.
