# Data folder

This folder holds the raw data inputs. To keep the repository light, the two
large open-access OWID files are **not** bundled. Download them from the links
below and place them in this folder before running `22_empirical_application_2.R`.

## Files to download manually

| File | Size | Source page | Direct download |
|------|------|-------------|-----------------|
| `owid_co2_data.csv` | ~32 MB | https://github.com/owid/co2-data | https://nyc3.digitaloceanspaces.com/owid-public/data/co2/owid-co2-data.csv |
| `owid_energy_data.csv` | ~16 MB | https://github.com/owid/energy-data | https://nyc3.digitaloceanspaces.com/owid-public/data/energy/owid-energy-data.csv |

After downloading, rename them if necessary so the filenames match exactly:
`owid_co2_data.csv` and `owid_energy_data.csv`.

## Data fetched automatically at runtime

The World Bank World Development Indicators series used in both applications are
downloaded automatically by the `WDI` R package when the scripts run, so no
manual download is needed for those. The exact series codes and indicator links
are documented in the headers of `21_empirical_application_1.R` and
`22_empirical_application_2.R`.
