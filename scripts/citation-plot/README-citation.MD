# Cloudy Citations Tracker

## Overview
This script tracks citations for various versions of the Cloudy code and Hazy documentation using the NASA ADS API. It retrieves citation counts per year and distinguishes between refereed and non-refereed sources.

## Author
Written by **Maryam Dehghanian** on **February 6, 2025**.

## How to Use
1. Ensure you have an **ADS API token**.
2. Replace `YOUR_ADS_API_TOKEN` in the script with your own ADS token.
3. Run the script using Python.

## ADS API Token
To access citation data, you must have an ADS API token:
- Visit [NASA ADS](https://ui.adsabs.harvard.edu/)
- Log in or create an account.
- Go to **Account Settings** > **API Tokens**.
- Generate a new token and replace `YOUR_ADS_API_TOKEN` in the script.

## Output
The script prints:
- **Total citations** for each version of Cloudy and Hazy.
- **Yearly refereed citations** across all versions.

## Dependencies
- `ads` (install via `pip install ads`)
- `matplotlib` (install via `pip install matplotlib`)


