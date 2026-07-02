# Toronto Suspected Non-Fatal Opioid Overdoses in Shelters: Geospatial Analysis and Spatial Regression

An address-level geospatial analysis of suspected non-fatal opioid overdoses in Toronto shelters and supportive housing sites, 2018–2025, with a spatial regression layer on the 158 official Toronto neighborhoods.

The analysis identifies which neighborhood-level features predict overdose incidence after accounting for spatial dependence, and where the model fails — a result that reinforces the descriptive finding that a small number of shelter addresses account for a disproportionate share of incidents.

---

**[▶ Explore the interactive map](https://ayokumo.github.io/toronto-overdose-spatial/overdose_map.html)** — click through overdose counts across Toronto shelter locations.

---

## Results at a glance

| Question | Finding |
|---|---|
| Is incidence concentrated? | Yes — a handful of addresses drive the majority of recorded incidents; several downtown neighborhoods have a concentration ratio of 1.0 (all incidents from a single address). |
| When did it intensify? | Most top-five addresses peak in 2020–2021, consistent with pandemic-era supply changes and service disruption. |
| Is the clustering "spatial"? | Largely not. Moran's I on OLS residuals = 0.057 (p = 0.086) — below significance. The clustering is a socioeconomic process that happens to be spatially organized. |
| What is the key predictor? | Low-income measure after-tax (LIM-AT) prevalence. Under HC3-robust inference, each +1 pp in LIM-AT ≈ +8% incidents (p = 0.023), controlling for population, area, downtown distance, and unemployment. |
| Which model was chosen? | Spatial lag (`ML_Lag`) — it eliminates residual autocorrelation (residual Moran's I = −0.024, p = 0.381); the spatial error model does not. |

---

## Key findings

**1. Concentration at the address level.** A handful of shelter and supportive housing addresses account for the majority of recorded incidents. The single highest-burden address contributes the largest share of overdoses in its neighborhood, and several downtown neighborhoods have concentration ratios of 1.0 — every recorded incident comes from a single address.

**2. The opioid crisis intensified sharply through 2020–2021.** Most of the top-five incident addresses show their peak in 2020 or 2021, consistent with the documented impact of pandemic-era supply changes (fentanyl analogs) and service disruption.

![Top overdose locations over time](figures/02_top_locations.png)

**3. After controlling for socioeconomic and exposure variables, most apparent geographic clustering is explained.** Moran's I on the OLS residuals is 0.057 (p = 0.086) — below conventional significance. The "spatial process" of overdose clustering is largely a *socioeconomic process that happens to be spatially organized*.

**4. Low-income concentration (LIM-AT prevalence) is the central socioeconomic predictor.** Under HC3 heteroskedasticity-consistent inference, each one-percentage-point increase in LIM-AT is associated with approximately 8% more incidents (p = 0.023), controlling for population, area, downtown distance, and unemployment.

![Coefficient comparison across models](figures/03_coefficient_comparison.png)

**5. The model's residuals confirm the address-level concentration finding.** The largest under-predictions concentrate in Moss Park, Harbourfront-CityPlace, and South Parkdale — exactly the neighborhoods where a single shelter address drives the count. No plausible set of neighborhood-level features can capture this, because the driver operates *below* the neighborhood scale.

![Observed incidents vs spatial error model residuals](figures/04_residual_choropleth.png)

---

## Methodology

The analysis proceeds in two parts:

### Descriptive geospatial analysis (Sections 1–8)

- Data cleaning of the City of Toronto Suspected Non-Fatal Opioid Overdoses in Shelters dataset (2018–2025), including treatment of suppressed counts (`< 5`) and address normalization
- Geocoding of unique shelter and supportive housing addresses using `geopy`, with results cached locally for reproducibility
- Year-over-year and address-level trend analysis
- Interactive [`folium` map](https://ayokumo.github.io/toronto-overdose-spatial/overdose_map.html) overlaying counts on Toronto's geography
- Targeted analysis of the Billy Bishop airport corridor, identifying a single address responsible for a disproportionate share of nearby incidents

### Spatial regression modeling (Section 8.5)

A diagnostic-then-correction sequence on the 158 official Toronto neighborhoods:

1. **Aggregation.** Address-level incident points spatially joined to neighborhood polygons via `geopandas`, projected to EPSG:32617 (UTM Zone 17N) for accurate distance and area calculations.
2. **Feature matrix.** Geometry-derived predictors (`log_area`, `dist_to_downtown_km`) combined with 2021 Census variables (`log_population`, `lim_at_prevalence`, `unemployment_rate`) from the City of Toronto Neighbourhood Profiles.
3. **Baseline OLS** with HC3 heteroskedasticity-consistent standard errors (`spreg.OLS` with `robust="white"`).
4. **Diagnostic battery.** Moran's I on residuals plus Lagrange Multiplier tests (lag, error, robust variants) using a queen-contiguity spatial weights matrix.
5. **Spatial extensions.** Both spatial lag (`ML_Lag`) and spatial error (`ML_Error`) models fit and compared on AIC, log-likelihood, and residual Moran's I.
6. **Model choice.** The spatial lag specification eliminates residual spatial autocorrelation (residual Moran's I = −0.024, p = 0.381) while the spatial error specification does not. The spatial lag is therefore preferred even though both models have similar AIC values.

A full discussion of findings, limitations, and natural next iterations (negative binomial GLM, H3 hex-level robustness check, panel spatial model with year effects) is in Section 8.5.10 of the notebook.

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── overdose_map.html                          # interactive folium map
├── notebooks/
│   └── overdose_analysis.ipynb                # main analysis notebook
├── data/
│   ├── README.md                              # data sources and download notes
│   ├── soois.csv                              # raw incident data (City of Toronto)
│   ├── geocoded_addresses.csv                 # cached geocoder output
│   └── toronto_neighbourhoods.geojson         # 158-neighborhood boundaries
│       # neighbourhood-profiles-2021-158-model.(xlsx|csv) is NOT committed —
│       # download it from the City of Toronto portal into data/ (see below).
└── figures/                                   # exported headline figures
```

---

## Reproducing the analysis

### Prerequisites

- Python 3.10+
- A working Python environment (`venv`, `conda`, or similar)

### Setup

```bash
git clone https://github.com/ayokumo/toronto-overdose-spatial.git
cd toronto-overdose-spatial
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### One manual data step

The 2021 Census profiles file is not redistributed in this repo. Download the
**Neighbourhood Profiles (2021, 158-model)** file from the
[City of Toronto Open Data portal](https://open.toronto.ca/dataset/neighbourhood-profiles/)
(the dataset is now marked *Retired*, so download it while it remains available)
and save it into `data/` as either:

- `data/neighbourhood-profiles-2021-158-model.xlsx`, or
- `data/neighbourhood-profiles-2021-158-model.csv`

The notebook's loader accepts either format.

### Run

Open `notebooks/overdose_analysis.ipynb` in Jupyter or VS Code and run cells
top-to-bottom. With the one data file above in place, the notebook runs end to
end. The geocoded address cache (`data/geocoded_addresses.csv`) is included, so
you do not need a Nominatim API key. The suppressed-count imputation is seeded
(`np.random.seed(42)`), so results are deterministic across runs.

Total runtime on a modern laptop is approximately 2–3 minutes, dominated by the
spatial weights matrix construction and the maximum-likelihood spatial
regression fits.

---

## Data sources

All inputs are publicly available and free to access:

| Dataset | Source | License |
|---|---|---|
| Suspected Non-Fatal Opioid Overdoses in Shelters | [City of Toronto Open Data](https://open.toronto.ca/dataset/suspected-non-fatal-opioid-overdoses-in-the-shelter-system/) | Open Government Licence – Toronto |
| Neighbourhoods (158-model) | [City of Toronto Open Data](https://open.toronto.ca/dataset/neighbourhoods/) | Open Government Licence – Toronto |
| Neighbourhood Profiles (2021 Census) | [City of Toronto Open Data](https://open.toronto.ca/dataset/neighbourhood-profiles/) | Open Government Licence – Toronto |
| Geocoding | OpenStreetMap via `geopy.Nominatim` | ODbL |

---

## Limitations and ethical considerations

This analysis works with publicly released, aggregated address-level counts. It does not use individual-level records and cannot identify individuals. Counts of fewer than 5 incidents at any address-year are suppressed by the City; this analysis imputes those values uniformly at random in the integer range **1–4** (a suppressed `< 5` count cannot be 5), using a fixed seed so results are reproducible. The imputation is documented in the notebook, and because the concentration finding is driven by high-count addresses, it is not sensitive to the imputed low counts.

Neighborhood-level associations do not translate to individual-level risk (ecological fallacy). The 158-neighborhood polygons are an administrative construct, and results may differ at finer (dissemination area, hex grid) or coarser scales — the Modifiable Areal Unit Problem is acknowledged in Section 8.5.10 as a limitation, and an H3 hex-resolution robustness check is identified as a natural next iteration.

The shelter-concentration finding should not be interpreted as a critique of shelter operators. Shelters serve populations at elevated risk, and a high overdose count at a shelter address is more plausibly read as evidence that the shelter is reaching the people who most need harm reduction services. Any policy reading of this analysis should center the people most affected.

---

## Author

**Ayokunmi Lawal** — Data Science (Hons.), Minor in Finance, York University

[GitHub](https://github.com/ayokumo) · Toronto, Ontario
