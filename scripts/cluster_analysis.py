"""Cluster analysis of Toronto shelter overdose incidents.

Two methods, run on the same data the notebook uses:

1. Local Moran's I (LISA) on neighbourhood log incidents, queen contiguity
   weights, 999 permutations, p < 0.05. Reports hotspot (HH) neighbourhoods
   and how many contiguous hotspot groups they form.
2. DBSCAN on the geocoded shelter addresses in metres (EPSG:32617).

The primary DBSCAN settings were fixed before any results were seen:
eps = 1000 m, min_samples = 3 addresses. A sensitivity grid is reported next
to it so the count can be read in context. Nothing is tuned toward a target.

Run from the repo root:
    python scripts/cluster_analysis.py

Writes results/*.csv, results/cluster_summary.md and figures/06_cluster_map.png.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from esda.moran import Moran_Local
from libpysal.weights import Queen
from sklearn.cluster import DBSCAN

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
OUT = ROOT / "results"
SEED = 42

# Fixed before looking at any result.
PRIMARY_EPS_M = 1000
PRIMARY_MIN_SAMPLES = 3
GRID_EPS = [500, 750, 1000, 1500]
GRID_MIN = [2, 3, 5]
LISA_P = 0.05


def load_address_totals() -> pd.DataFrame:
    """Reproduce the notebook's cleaning: impute '< 5' with 1 to 4 (seed 42),
    normalise addresses, total incidents per address, attach coordinates."""
    df = pd.read_csv(DATA / "soois.csv")
    df.columns = df.columns.str.lower().str.strip()
    col = "suspected_non_fatal_overdoses"

    np.random.seed(SEED)
    s = df[col].astype("object").map(lambda x: "" if pd.isna(x) else str(x))
    mask = s.str.contains(r"<\s*5", regex=True, na=False)
    s.loc[mask] = np.random.randint(1, 5, size=int(mask.sum())).astype(str)
    df[col] = s.str.extract(r"(\d+)")[0].fillna(0).astype(int)

    df["address"] = (
        df["address"]
        .str.upper()
        .str.strip()
        .str.replace("STREET", "ST", regex=False)
        .str.replace("AVENUE", "AVE", regex=False)
        .str.replace("ROAD", "RD", regex=False)
        .str.replace("BOULEVARD", "BLVD", regex=False)
    )
    totals = df.groupby("address")[col].sum().reset_index(name="incidents")

    geo = pd.read_csv(DATA / "geocoded_addresses.csv").drop_duplicates("address")
    merged = totals.merge(geo[["address", "lat", "lon"]], on="address", how="left")
    n_all, inc_all = len(merged), int(merged["incidents"].sum())
    merged = merged.dropna(subset=["lat", "lon"])
    print(f"Addresses with coordinates: {len(merged)} of {n_all}")
    print(f"Incidents covered: {int(merged['incidents'].sum())} of {inc_all}")
    merged.attrs["coverage"] = (len(merged), n_all, int(merged["incidents"].sum()), inc_all)
    return merged


def load_neighbourhoods() -> gpd.GeoDataFrame:
    g = gpd.read_file(DATA / "toronto_neighbourhoods.geojson").to_crs(epsg=32617)
    g = g[["AREA_NAME", "geometry"]].rename(columns={"AREA_NAME": "neighborhood"})
    f = pd.read_csv(DATA / "neighbourhood_features.csv")[
        ["neighborhood", "total_incidents", "log_incidents"]
    ]
    g = g.merge(f, on="neighborhood", how="left")
    assert g["log_incidents"].notna().all(), "name mismatch between geojson and features"
    return g.reset_index(drop=True)


def contiguous_groups(w, members: list[int]) -> list[list[int]]:
    """Connected components of `members` under the weights' adjacency."""
    members_set, seen, groups = set(members), set(), []
    for start in members:
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            comp.append(node)
            for nb in w.neighbors[node]:
                if nb in members_set and nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        groups.append(sorted(comp))
    return groups


def run_lisa(nb: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, dict]:
    w = Queen.from_dataframe(nb, ids=nb.index.tolist(), silence_warnings=True)
    w.transform = "r"
    lisa = Moran_Local(nb["log_incidents"].values, w, permutations=999, seed=SEED)
    nb = nb.copy()
    nb["lisa_q"] = lisa.q
    nb["lisa_p"] = lisa.p_sim
    label = {1: "HH", 2: "LH", 3: "LL", 4: "HL"}
    nb["lisa_class"] = [
        label[q] if p < LISA_P else "ns" for q, p in zip(nb["lisa_q"], nb["lisa_p"])
    ]
    hh = nb.index[nb["lisa_class"] == "HH"].tolist()
    groups = contiguous_groups(w, hh)
    hh_strict = nb.index[(nb["lisa_q"] == 1) & (nb["lisa_p"] < 0.01)].tolist()
    summary = {
        "counts": nb["lisa_class"].value_counts().to_dict(),
        "hh": hh,
        "groups": groups,
        "hh_p01": len(hh_strict),
        "n_zero": int((nb["total_incidents"] == 0).sum()),
    }
    return nb, summary


def run_dbscan(addr: pd.DataFrame, nb: gpd.GeoDataFrame):
    pts = gpd.GeoDataFrame(
        addr.copy(),
        geometry=gpd.points_from_xy(addr["lon"], addr["lat"]),
        crs="EPSG:4326",
    ).to_crs(epsg=32617)
    xy = np.column_stack([pts.geometry.x, pts.geometry.y])

    rows = []
    for eps in GRID_EPS:
        for ms in GRID_MIN:
            lab = DBSCAN(eps=eps, min_samples=ms).fit(xy).labels_
            k = len(set(lab)) - (1 if -1 in lab else 0)
            in_cl = pts.loc[lab != -1, "incidents"].sum()
            rows.append(
                {
                    "eps_m": eps,
                    "min_samples": ms,
                    "clusters": k,
                    "addresses_in_clusters": int((lab != -1).sum()),
                    "noise_addresses": int((lab == -1).sum()),
                    "incident_share_in_clusters": round(in_cl / pts["incidents"].sum(), 3),
                    "primary": eps == PRIMARY_EPS_M and ms == PRIMARY_MIN_SAMPLES,
                }
            )
    grid = pd.DataFrame(rows)

    pts["cluster"] = DBSCAN(eps=PRIMARY_EPS_M, min_samples=PRIMARY_MIN_SAMPLES).fit(xy).labels_
    joined = gpd.sjoin(pts, nb[["neighborhood", "geometry"]], how="left", predicate="within")
    pts["neighborhood"] = joined["neighborhood"].values
    return pts, grid


def describe_clusters(pts: gpd.GeoDataFrame) -> pd.DataFrame:
    rows = []
    total = pts["incidents"].sum()
    for cid in sorted(c for c in pts["cluster"].unique() if c != -1):
        sub = pts[pts["cluster"] == cid]
        top = sub.sort_values("incidents", ascending=False).iloc[0]
        by_nb = sub.groupby("neighborhood")["incidents"].sum().sort_values(ascending=False)
        rows.append(
            {
                "cluster": int(cid),
                "addresses": len(sub),
                "incidents": int(sub["incidents"].sum()),
                "share_of_all_incidents": round(sub["incidents"].sum() / total, 3),
                "top_address": top["address"],
                "top_address_incidents": int(top["incidents"]),
                "neighbourhoods": "; ".join(f"{n} ({int(v)})" for n, v in by_nb.items()),
            }
        )
    return pd.DataFrame(rows).sort_values("incidents", ascending=False).reset_index(drop=True)


def make_map(nb: gpd.GeoDataFrame, pts: gpd.GeoDataFrame, clusters: pd.DataFrame, run_date: str):
    fig, ax = plt.subplots(figsize=(11, 9))
    nb.plot(ax=ax, color="#f2f2f2", edgecolor="white", linewidth=0.4)
    hh = nb[nb["lisa_class"] == "HH"]
    if len(hh):
        hh.plot(ax=ax, color="#f4b6b6", edgecolor="#c0392b", linewidth=0.9)
    noise = pts[pts["cluster"] == -1]
    ax.scatter(noise.geometry.x, noise.geometry.y, s=np.sqrt(noise["incidents"]) * 3 + 4,
               color="#8c8c8c", alpha=0.55, linewidths=0, label="Not in a DBSCAN cluster")
    cmap = plt.get_cmap("tab10")
    rank = {int(c): i for i, c in enumerate(clusters["cluster"])}
    for cid, i in rank.items():
        sub = pts[pts["cluster"] == cid]
        ax.scatter(sub.geometry.x, sub.geometry.y, s=np.sqrt(sub["incidents"]) * 3 + 8,
                   color=cmap(i % 10), edgecolors="black", linewidths=0.5,
                   label=f"DBSCAN cluster {i + 1}: {len(sub)} addresses, {int(sub['incidents'].sum())} incidents")
    ax.plot([], [], marker="s", linestyle="", color="#f4b6b6", markeredgecolor="#c0392b",
            label="LISA hotspot neighbourhood (HH, p < 0.05)")
    ax.legend(loc="lower left", fontsize=8, frameon=True)
    ax.set_title(
        f"Overdose incident clusters in Toronto shelters, 2018 to 2025\n"
        f"DBSCAN (eps {PRIMARY_EPS_M} m, min samples {PRIMARY_MIN_SAMPLES}) and local Moran's I. Run {run_date}.",
        fontsize=11)
    ax.axis("off")
    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "06_cluster_map.png", dpi=150)
    plt.close(fig)


def main() -> None:
    run_date = date.today().isoformat()
    OUT.mkdir(exist_ok=True)
    addr = load_address_totals()
    nb = load_neighbourhoods()

    nb, lisa = run_lisa(nb)
    pts, grid = run_dbscan(addr, nb)
    clusters = describe_clusters(pts)
    clusters.insert(0, "rank", range(1, len(clusters) + 1))

    # Save tables
    nb.drop(columns="geometry").sort_values("lisa_p").to_csv(OUT / "lisa_neighbourhoods.csv", index=False)
    grid.to_csv(OUT / "dbscan_sensitivity.csv", index=False)
    clusters.to_csv(OUT / "dbscan_clusters_primary.csv", index=False)
    make_map(nb, pts, clusters.rename(columns={"rank": "r"}).assign(cluster=clusters["cluster"]), run_date)

    # Robustness only, not the headline: drop addresses with zero recorded incidents.
    pos = pts[pts["incidents"] > 0]
    pos_lab = DBSCAN(eps=PRIMARY_EPS_M, min_samples=PRIMARY_MIN_SAMPLES).fit(
        np.column_stack([pos.geometry.x, pos.geometry.y])
    ).labels_
    robust_k = len(set(pos_lab)) - (1 if -1 in pos_lab else 0)
    n_zero_addr = int((pts["incidents"] == 0).sum())

    n_addr, n_all, inc_cov, inc_all = addr.attrs["coverage"]
    hh_names = nb.loc[lisa["hh"], ["neighborhood", "total_incidents", "lisa_p"]].sort_values("total_incidents", ascending=False)
    group_names = [
        [nb.loc[i, "neighborhood"] for i in g] for g in lisa["groups"]
    ]

    lines = [
        f"Run date: {run_date}. Seed {SEED}.",
        f"Coverage: {n_addr} of {n_all} addresses have coordinates, covering {inc_cov} of {inc_all} incidents.",
        "",
        f"DBSCAN primary (eps {PRIMARY_EPS_M} m, min samples {PRIMARY_MIN_SAMPLES}, fixed in advance): {len(clusters)} clusters.",
        "",
        clusters.to_string(index=False),
        "",
        f"Robustness, addresses with at least one incident only ({len(pos)} of {len(pts)} addresses; {n_zero_addr} have zero): {robust_k} clusters at the primary settings.",
        "",
        "DBSCAN sensitivity grid:",
        grid.to_string(index=False),
        "",
        f"LISA classes at p < {LISA_P}: {lisa['counts']}",
        f"LISA HH hotspot neighbourhoods: {len(lisa['hh'])} in {len(lisa['groups'])} contiguous group(s). At p < 0.01: {lisa['hh_p01']} HH neighbourhoods.",
        f"Neighbourhoods with zero incidents: {lisa['n_zero']} of {len(nb)} (many ties at zero).",
        "",
        "HH neighbourhoods:",
        hh_names.to_string(index=False),
        "",
        "Contiguous HH groups:",
        *[f"  Group {i + 1}: {', '.join(g)}" for i, g in enumerate(group_names)],
    ]
    text = "\n".join(lines)
    (OUT / "cluster_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
