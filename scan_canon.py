# ============================================================
# NASA 5000-YEAR SOLAR ECLIPSE SCAN
# MECCA + JERUSALEM
# >= 70% SOLAR DISK AREA OBSCURATION
#
# Input:
#    C:\Users\...\Desktop\nasa.csv
#
# Output:
#    C:\Users\...\Desktop\mecca_jerusalem_70_EXACT.csv
#
# NASA Besselian elements:
# x, y, d, mu, l1, l2
#
# Obscuration:
#    overlap area of apparent Moon disk / apparent Sun disk
#
# IMPORTANT:
#    70% AREA obscuration != eclipse magnitude 0.70
# ============================================================

import os
import math
import time

import numpy as np
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = r"C:\Users\...\Desktop\nasa.csv"
OUTPUT_FILE = r"C:\Users\...\Desktop\mecca_jerusalem_70_EXACT.csv"

# Observer coordinates
MECCA_LAT = 21.4225
MECCA_LON = 39.8262

JERUSALEM_LAT = 31.7683
JERUSALEM_LON = 35.2137

THRESHOLD = 0.70

# Coarse scan around NASA Besselian t0
# NASA elements are fitted over approximately a 6-hour interval.
COARSE_POINTS = 241

# Fine scan around the coarse maximum
FINE_POINTS = 121

# Fine search half-width in hours
FINE_HALF_WIDTH = 0.08

# Search interval around t0
SEARCH_HALF_WIDTH = 3.0


# ============================================================
# EARTH PARAMETERS
# WGS84
# ============================================================

EARTH_FLATTENING = 1.0 / 298.257223563
EARTH_E2 = (
    EARTH_FLATTENING *
    (2.0 - EARTH_FLATTENING)
)


# ============================================================
# BASIC FUNCTIONS
# ============================================================

def poly(row, prefix, t):
    """
    Evaluate NASA Besselian polynomial:

        a = a0 + a1*t + a2*t^2 + a3*t^3

    where t is hours relative to t0.
    """

    a0 = float(row[f"{prefix}0"])
    a1 = float(row[f"{prefix}1"])
    a2 = float(row[f"{prefix}2"])

    a3_col = f"{prefix}3"

    if a3_col in row.index and pd.notna(row[a3_col]):
        a3 = float(row[a3_col])
    else:
        a3 = 0.0

    return (
        a0
        + a1 * t
        + a2 * t * t
        + a3 * t * t * t
    )


def circle_overlap_area(r1, r2, d):
    """
    Exact area of intersection of two circles.

    r1 = solar radius
    r2 = lunar radius
    d  = center separation

    Returns overlap area.
    """

    if r1 <= 0.0 or r2 <= 0.0:
        return 0.0

    # No overlap
    if d >= r1 + r2:
        return 0.0

    # One circle completely inside the other
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2

    # Partial overlap
    c1 = (
        d * d
        + r1 * r1
        - r2 * r2
    ) / (2.0 * d * r1)

    c2 = (
        d * d
        + r2 * r2
        - r1 * r1
    ) / (2.0 * d * r2)

    # Numerical protection
    c1 = max(-1.0, min(1.0, c1))
    c2 = max(-1.0, min(1.0, c2))

    alpha = math.acos(c1)
    beta = math.acos(c2)

    area = (
        r1 * r1 * alpha
        + r2 * r2 * beta
        - 0.5
        * math.sqrt(
            max(
                0.0,
                (-d + r1 + r2)
                * (d + r1 - r2)
                * (d - r1 + r2)
                * (d + r1 + r2)
            )
        )
    )

    return area


# ============================================================
# OBSERVER CONSTANTS
# ============================================================

def observer_constants(lat_deg):
    """
    Geocentric observer constants for WGS84.
    """

    phi = math.radians(lat_deg)

    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    denom = math.sqrt(
        1.0
        - EARTH_E2 * sin_phi * sin_phi
    )

    rho_sin_phi_prime = (
        (1.0 - EARTH_FLATTENING)
        * sin_phi
        / denom
    )

    rho_cos_phi_prime = (
        cos_phi
        / denom
    )

    return (
        rho_sin_phi_prime,
        rho_cos_phi_prime
    )


# ============================================================
# LOCAL ECLIPSE GEOMETRY
# ============================================================

def local_geometry(row, lat_deg, lon_deg, tau):
    """
    Calculate local eclipse geometry at a particular observer.

    tau:
        hours relative to NASA Besselian t0

    Returns:
        obscuration
        altitude
        separation
        solar_radius
        lunar_radius
        magnitude
    """

    # --------------------------------------------------------
    # Besselian elements
    # --------------------------------------------------------

    x = poly(row, "x", tau)
    y = poly(row, "y", tau)
    d_deg = poly(row, "d", tau)
    mu_deg = poly(row, "mu", tau)

    l1 = poly(row, "l1", tau)
    l2 = poly(row, "l2", tau)

    tan_f1 = float(row["tan_f1"])
    tan_f2 = float(row["tan_f2"])

    d = math.radians(d_deg)
    mu = math.radians(mu_deg)
    lon = math.radians(lon_deg)

    # --------------------------------------------------------
    # Observer constants
    # --------------------------------------------------------

    rho_sin, rho_cos = observer_constants(lat_deg)

    # Ephemeris hour angle
    theta = mu + lon

    # --------------------------------------------------------
    # Observer coordinates in fundamental-plane system
    #
    # xi, eta = observer position projected onto fundamental plane
    # zeta    = observer coordinate along shadow axis
    #
    # Standard Besselian observer geometry.
    # --------------------------------------------------------

    xi = (
        rho_cos
        * math.sin(theta)
    )

    eta = (
        rho_sin * math.cos(d)
        - rho_cos
        * math.sin(d)
        * math.cos(theta)
    )

    zeta = (
        rho_sin * math.sin(d)
        + rho_cos
        * math.cos(d)
        * math.cos(theta)
    )

    # --------------------------------------------------------
    # Shadow center relative to observer
    # --------------------------------------------------------

    u = x - xi
    v = y - eta

    separation = math.sqrt(
        u * u + v * v
    )

    # --------------------------------------------------------
    # Local penumbral / umbral radii
    #
    # L1_local = L1 - zeta*tan(f1)
    # L2_local = L2 - zeta*tan(f2)
    #
    # For apparent disks:
    #    L1 = R_sun + R_moon
    #    L2 = R_sun - R_moon
    #
    # Therefore:
    #    R_sun  = (L1 + L2)/2
    #    R_moon = (L1 - L2)/2
    # --------------------------------------------------------

    L1_local = (
        l1
        - zeta * tan_f1
    )

    L2_local = (
        l2
        - zeta * tan_f2
    )

    solar_radius = (
        L1_local + L2_local
    ) / 2.0

    lunar_radius = (
        L1_local - L2_local
    ) / 2.0

    # Invalid geometry
    if solar_radius <= 0.0 or lunar_radius <= 0.0:
        return {
            "obscuration": 0.0,
            "altitude": -90.0,
            "separation": separation,
            "solar_radius": solar_radius,
            "lunar_radius": lunar_radius,
            "magnitude": 0.0
        }

    # --------------------------------------------------------
    # Solar altitude
    # --------------------------------------------------------

    sin_alt = (
        rho_cos
        * math.cos(d)
        * math.cos(theta)
        + rho_sin
        * math.sin(d)
    )

    sin_alt = max(
        -1.0,
        min(1.0, sin_alt)
    )

    altitude = math.degrees(
        math.asin(sin_alt)
    )

    # --------------------------------------------------------
    # Eclipse magnitude
    # Reference calculation only (diameter-based metric)
    # --------------------------------------------------------

    magnitude = (
        L1_local - separation
    ) / (
        L1_local + L2_local
    )

    # --------------------------------------------------------
    # TRUE AREA OBSCURATION
    # Fraction of the SOLAR DISK covered by Moon.
    # --------------------------------------------------------

    overlap = circle_overlap_area(
        solar_radius,
        lunar_radius,
        separation
    )

    solar_area = (
        math.pi
        * solar_radius
        * solar_radius
    )

    if solar_area > 0.0:
        obscuration = (
            overlap / solar_area
        )
    else:
        obscuration = 0.0

    # Numerical protection
    obscuration = max(
        0.0,
        min(1.0, obscuration)
    )

    return {
        "obscuration": obscuration,
        "altitude": altitude,
        "separation": separation,
        "solar_radius": solar_radius,
        "lunar_radius": lunar_radius,
        "magnitude": magnitude
    }


# ============================================================
# MAXIMUM LOCAL OBSCURATION
# ============================================================

def maximum_obscuration(row, lat_deg, lon_deg):
    """
    Find the maximum observable solar-disk obscuration
    at a given location during one eclipse.

    Coarse scan followed by fine scan around peak.
    """

    # NASA Besselian t0
    t0 = float(row["t0"])

    # --------------------------------------------------------
    # Coarse scan
    # --------------------------------------------------------

    tau_values = np.linspace(
        -SEARCH_HALF_WIDTH,
        SEARCH_HALF_WIDTH,
        COARSE_POINTS
    )

    best = None

    for tau in tau_values:

        g = local_geometry(
            row,
            lat_deg,
            lon_deg,
            float(tau)
        )

        # Sun below horizon -> not observable
        if g["altitude"] <= 0.0:
            continue

        if best is None or g["obscuration"] > best[0]:
            best = (
                g["obscuration"],
                float(tau),
                g
            )

    if best is None:
        return {
            "obscuration": 0.0,
            "tau": np.nan,
            "clock": np.nan,
            "altitude": -90.0,
            "magnitude": 0.0
        }

    # --------------------------------------------------------
    # Fine scan
    # --------------------------------------------------------

    center = best[1]

    fine_values = np.linspace(
        center - FINE_HALF_WIDTH,
        center + FINE_HALF_WIDTH,
        FINE_POINTS
    )

    for tau in fine_values:

        g = local_geometry(
            row,
            lat_deg,
            lon_deg,
            float(tau)
        )

        if g["altitude"] <= 0.0:
            continue

        if g["obscuration"] > best[0]:
            best = (
                g["obscuration"],
                float(tau),
                g
            )

    obsc = best[0]
    tau = best[1]
    g = best[2]

    # Absolute Besselian decimal hour
    clock = t0 + tau

    return {
        "obscuration": obsc,
        "tau": tau,
        "clock": clock,
        "altitude": g["altitude"],
        "magnitude": g["magnitude"]
    }


# ============================================================
# LOAD DATA
# ============================================================

print()
print("=" * 70)
print("NASA 5000-YEAR SOLAR ECLIPSE SCAN")
print("MECCA + JERUSALEM")
print("=" * 70)
print()

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(
        f"CSV file not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(
    INPUT_FILE
)

print(
    f"Total NASA records: {len(df)}"
)

# ------------------------------------------------------------
# Required columns
# ------------------------------------------------------------

required = [
    "year",
    "month",
    "day",
    "eclipse_type",
    "cat_no",
    "t0",

    "x0", "x1", "x2", "x3",
    "y0", "y1", "y2", "y3",

    "d0", "d1", "d2",

    "mu0", "mu1", "mu2",

    "l10", "l11", "l12",
    "l20", "l21", "l22",

    "tan_f1",
    "tan_f2"
]

missing = [
    c for c in required
    if c not in df.columns
]

if missing:
    print()
    print("MISSING COLUMNS:")
    for c in missing:
        print("  ", c)

    print()
    print("Available columns:")
    print(df.columns.tolist())

    raise ValueError(
        "CSV layout does not match expected NASA Besselian elements format."
    )

print("Column validation: OK")


# ============================================================
# REMOVE INVALID ROWS
# ============================================================

numeric_columns = [
    "t0",
    "x0", "x1", "x2", "x3",
    "y0", "y1", "y2", "y3",
    "d0", "d1", "d2",
    "mu0", "mu1", "mu2",
    "l10", "l11", "l12",
    "l20", "l21", "l22",
    "tan_f1",
    "tan_f2"
]

for c in numeric_columns:
    df[c] = pd.to_numeric(
        df[c],
        errors="coerce"
    )

df = df.dropna(
    subset=numeric_columns
).reset_index(drop=True)

print(
    f"Valid NASA records: {len(df)}"
)

print()


# ============================================================
# SCAN
# ============================================================

results = []

start_time = time.time()

total = len(df)

for idx, row in df.iterrows():

    # --------------------------------------------------------
    # Mecca
    # --------------------------------------------------------

    mecca = maximum_obscuration(
        row,
        MECCA_LAT,
        MECCA_LON
    )

    # --------------------------------------------------------
    # Optimization:
    # Only calculate Jerusalem if Mecca meets the threshold
    # --------------------------------------------------------

    if mecca["obscuration"] < THRESHOLD:

        if (
            idx == 0
            or (idx + 1) % 500 == 0
        ):
            elapsed = time.time() - start_time
            speed = (idx + 1) / elapsed

            remaining = total - (idx + 1)

            eta = (
                remaining / speed
                if speed > 0
                else 0
            )

            print(
                f"Processed: {idx+1}/{total} "
                f"({100*(idx+1)/total:.1f}%) "
                f"| Speed: {speed:.1f}/s "
                f"| ETA: {eta:.1f} s",
                end="\r"
            )

        continue

    # --------------------------------------------------------
    # Jerusalem
    # --------------------------------------------------------

    jerusalem = maximum_obscuration(
        row,
        JERUSALEM_LAT,
        JERUSALEM_LON
    )

    # --------------------------------------------------------
    # Final criterion:
    # Mecca >= 70% area AND Jerusalem >= 70% area
    # --------------------------------------------------------

    if (
        mecca["obscuration"] >= THRESHOLD
        and
        jerusalem["obscuration"] >= THRESHOLD
    ):

        results.append({
            "year": int(row["year"]),
            "month": int(row["month"]),
            "day": int(row["day"]),
            "eclipse_type": row["eclipse_type"],
            "cat_no": row["cat_no"],

            "Mecca_obscuration_%":
                100.0 * mecca["obscuration"],

            "Jerusalem_obscuration_%":
                100.0 * jerusalem["obscuration"],

            "Mecca_magnitude":
                mecca["magnitude"],

            "Jerusalem_magnitude":
                jerusalem["magnitude"],

            "Mecca_altitude_deg":
                mecca["altitude"],

            "Jerusalem_altitude_deg":
                jerusalem["altitude"],

            "Mecca_tau_hours":
                mecca["tau"],

            "Jerusalem_tau_hours":
                jerusalem["tau"],

            "Mecca_t_decimal_hour":
                mecca["clock"],

            "Jerusalem_t_decimal_hour":
                jerusalem["clock"]
        })

        print()
        print(
            f"  >>> MATCH: "
            f"{int(row['year']):+05d}-"
            f"{int(row['month']):02d}-"
            f"{int(row['day']):02d} "
            f"{row['eclipse_type']} | "
            f"Mecca={100*mecca['obscuration']:.3f}% | "
            f"Jerusalem={100*jerusalem['obscuration']:.3f}%"
        )

    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------

    if (
        idx == 0
        or (idx + 1) % 500 == 0
    ):

        elapsed = time.time() - start_time
        speed = (idx + 1) / elapsed

        remaining = total - (idx + 1)

        eta = (
            remaining / speed
            if speed > 0
            else 0
        )

        print(
            f"Processed: {idx+1}/{total} "
            f"({100*(idx+1)/total:.1f}%) "
            f"| Speed: {speed:.1f}/s "
            f"| ETA: {eta:.1f} s",
            end="\r"
        )


# ============================================================
# SAVE
# ============================================================

print()
print()

result_df = pd.DataFrame(results)

if len(result_df) > 0:

    result_df = result_df.sort_values(
        by=[
            "year",
            "month",
            "day"
        ]
    )

    result_df.to_csv(
        OUTPUT_FILE,
        index=False,
        encoding="utf-8-sig"
    )

else:

    # Create an empty CSV with standard schema if no matches
    columns = [
        "year",
        "month",
        "day",
        "eclipse_type",
        "cat_no",
        "Mecca_obscuration_%",
        "Jerusalem_obscuration_%",
        "Mecca_magnitude",
        "Jerusalem_magnitude",
        "Mecca_altitude_deg",
        "Jerusalem_altitude_deg",
        "Mecca_tau_hours",
        "Jerusalem_tau_hours",
        "Mecca_t_decimal_hour",
        "Jerusalem_t_decimal_hour"
    ]

    result_df = pd.DataFrame(
        columns=columns
    )

    result_df.to_csv(
        OUTPUT_FILE,
        index=False,
        encoding="utf-8-sig"
    )


# ============================================================
# SUMMARY
# ============================================================

elapsed = time.time() - start_time

print("=" * 70)
print("SCAN COMPLETED")
print("=" * 70)

print(
    f"NASA records scanned : {total}"
)

print(
    f">= 70% dual matches  : {len(result_df)}"
)

print(
    f"Total elapsed time   : {elapsed:.2f} seconds"
)

print()
print(
    "Output file:"
)
print(
    OUTPUT_FILE
)

print()

if len(result_df) > 0:

    print("FIRST RESULTS:")
    print()

    display_columns = [
        "year",
        "month",
        "day",
        "eclipse_type",
        "Mecca_obscuration_%",
        "Jerusalem_obscuration_%",
        "Mecca_altitude_deg",
        "Jerusalem_altitude_deg"
    ]

    print(
        result_df[
            display_columns
        ].head(30).to_string(
            index=False
        )
    )

    print()
    print(
        f"Total of {len(result_df)} matching events found."
    )

else:

    print(
        "No eclipses reached >= 70% area obscuration at both locations."
    )

print()
print("=" * 70)
