# Weather Data Fetcher

A simple R-based CLI tool to download Copernicus ERA5 weather data for specified months, years, and bounding box, and optionally create a hybrid weather dataset.

---

## Table of Contents

- [Requirements](#requirements)  
- [Installation](#installation)  
- [Input JSON Structure](#input-json-structure)  
- [Usage](#usage)  
- [Example](#example)  
- [Output](#output)  

---

## Requirements

- R (>= 4.0.0)  
- `ncdf4` R package  
- Internet connection (for Copernicus API)

---

## Installation

1. Clone this repository:  
   ```bash
   git clone https://github.com/yourusername/weather-fetcher.git
   cd weather-fetcher
   ```
2. Install required R packages (if not already installed):  
   ```r
   install.packages("ncdf4")
   ```

---

## Input JSON Structure

Your configuration JSON file should include the following fields:

```json
{
  "uid": "YOUR_COPERNICUS_UID",
  "key": "YOUR_COPERNICUS_KEY",
  "start_month": 1,
  "end_month": 12,
  "start_year": 2025,
  "end_year": 2025,
  "bounding_box": {
    "north": 40.613687151061505,
    "south": 35.99538225441254,
    "east": -89.09913230512305,
    "west": -95.76804054400543
  },
  "make_hybrid_weather": true
}
```

| Field                 | Type     | Description                                                                                           |
|-----------------------|----------|-------------------------------------------------------------------------------------------------------|
| `uid`                 | string   | Your Copernicus ERA5 API user ID.                                                                    |
| `key`                 | string   | Your Copernicus ERA5 API key.                                                                        |
| `start_month`         | integer  | Numeric month to start (1 = January, …, 12 = December).                                               |
| `end_month`           | integer  | Numeric month to end (inclusive).                                                                    |
| `start_year`          | integer  | Year (e.g., 2025) to start downloading data.                                                         |
| `end_year`            | integer  | Year to end downloading data.                                                                        |
| `bounding_box`        | object   | Geographic bounding box for data request.                                                             |
| └─ `north`            | numeric  | Maximum latitude (degrees north).                                                                    |
| └─ `south`            | numeric  | Minimum latitude (degrees north).                                                                    |
| └─ `east`             | numeric  | Maximum longitude (degrees east).                                                                    |
| └─ `west`             | numeric  | Minimum longitude (degrees east).                                                                    |
| `make_hybrid_weather` | boolean  | If `true`, create a combined/hybrid weather dataset from downloaded variables.                       |

---

## Usage

From the repository root, run:

```bash
Rscript main.R path/to/input.json
```

- **`main.R`**: Main script that parses the JSON configuration, fetches data via the Copernicus API, and (optionally) generates a hybrid weather file.
- **`path/to/input.json`**: Path to your JSON configuration file.

---

## Example

1. Create a file named `config.json`:

   ```json
   {
     "uid": "abc123def456",
     "key": "789ghi012jkl",
     "start_month": 1,
     "end_month": 12,
     "start_year": 2025,
     "end_year": 2025,
     "bounding_box": {
       "north": 40.613687151061505,
       "south": 35.99538225441254,
       "east": -89.09913230512305,
       "west": -95.76804054400543
     },
     "make_hybrid_weather": true
   }
   ```

2. Run the script:

   ```bash
   Rscript main.R config.json
   ```

---

## Output

- One NetCDF file per requested variable and year/month, saved in `./data/`.
- If `make_hybrid_weather = true`, a combined NetCDF named `hybrid_weather_{start_year}_{end_year}.nc` will be created in the project directory.

