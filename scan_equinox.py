# ============================================================
# NASA 5,000-YEAR VERNAL EQUINOX SOLAR ECLIPSE SCAN
# ONLY GLOBAL CANON (nasa.csv -> 11,898 RECORDS)
#
# Methodology:
#   Jean Meeus (Astronomical Algorithms, Ch. 27) Vernal Equinox JDE
#   versus NASA Five Millennium Canon of Solar Eclipses (-1999 to +3000)
#
# Input:
#   nasa.csv
#
# Output:
#   nasa_equinox_eclipses.csv
#
# Criteria:
#   Temporal proximity to Vernal Equinox zero-point: |Delta_t| <= 24 hours
# ============================================================

import os
import math
import numpy as np
import pandas as pd

# ============================================================
# SETTINGS & FILE PATHS
# ============================================================

INPUT_FILE = r"C:\Users\...\Desktop\nasa.csv"
OUTPUT_FILE = r"C:\Users\...\Desktop\nasa_equinox_eclipses.csv"

# Temporal proximity window around equinox (1.0 day = 24.0 hours)
EQUINOX_THRESHOLD_DAYS = 1.0


# ============================================================
# ASTRONOMICAL ALGORITHMS (JEAN MEEUS)
# ============================================================

def vernal_equinox_jde(year):
    """
    Compute the Julian Ephemeris Day (JDE) of the Vernal Equinox
    for a given astronomical year (-1999 to +3000).
    Reference: Jean Meeus, 'Astronomical Algorithms' (2nd ed., 1998), Chapter 27.
    """
    if -1000 <= year <= 1000:
        Y = year / 1000.0
        jde0 = (
            1721139.29189
            + 365242.13740 * Y
            + 0.06134 * (Y ** 2)
            + 0.00111 * (Y ** 3)
            - 0.00071 * (Y ** 4)
        )
    elif 1000 < year <= 3000:
        Y = (year - 2000) / 1000.0
        jde0 = (
            2451623.80984
            + 365242.37404 * Y
            + 0.05169 * (Y ** 2)
            - 0.00411 * (Y ** 3)
            - 0.00057 * (Y ** 4)
        )
    else:  # -2000 <= year < -1000
        Y = year / 1000.0
        jde0 = (
            1721139.29189
            + 365242.13740 * Y
            + 0.06134 * (Y ** 2)
            + 0.00111 * (Y ** 3)
            - 0.00071 * (Y ** 4)
        )

    # Periodic planetary perturbation terms (Meeus Table 27.C)
    T = (jde0 - 2451545.0) / 36525.0
    W = 35999.373 * T - 2.47
    delta_lambda = (
        1.0
        + 0.0334 * math.cos(math.radians(W))
        + 0.0007 * math.cos(math.radians(2 * W))
    )

    S = (
        485 * math.cos(math.radians(324.96 + 1934.136 * T))
        + 203 * math.cos(math.radians(276.11 + 72001.539 * T))
        + 199 * math.cos(math.radians(119.75 + 1318.515 * T))
        + 182 * math.cos(math.radians(263.36 + 34575.405 * T))
        + 156 * math.cos(math.radians(156.46 + 65928.934 * T))
        + 136 * math.cos(math.radians(41.85 + 20186.439 * T))
        + 77 * math.cos(math.radians(246.54 + 2924.669 * T))
        + 74 * math.cos(math.radians(320.61 + 71498.441 * T))
        + 70 * math.cos(math.radians(107.94 + 54994.431 * T))
        + 58 * math.cos(math.radians(278.57 + 6077.581 * T))
        + 52 * math.cos(math.radians(87.15 + 32243.061 * T))
        + 50 * math.cos(math.radians(26.80 + 73005.807 * T))
    )

    return jde0 + (0.00001 * S) / delta_lambda


def calendar_to_jd(year, month, day, decimal_hour=12.0):
    """
    Convert civil calendar date to Julian Day (JD).
    Adheres to the Julian calendar prior to 1582-10-15 and the
    proleptic Gregorian calendar thereafter (matching NASA Canon specifications).
    """
    if month <= 2:
        year -= 1
        month += 12

    A = int(year / 100)
    is_gregorian = (year > 1582) or (
        year == 1582 and (month > 10 or (month == 10 and day >= 15))
    )
    B = (2 - A + int(A / 4)) if is_gregorian else 0

    jd = (
        int(365.25 * (year + 4716))
        + int(30.6001 * (month + 1))
        + day
        + B
        - 1524.5
        + (decimal_hour / 24.0)
    )
    return jd


# ============================================================
# SCAN EXECUTION
# ============================================================

def scan_nasa_equinox():
    print()
    print("=" * 70)
    print("NASA 5,000-YEAR GLOBAL VERNAL EQUINOX SCAN (11,898 EVENTS)")
    print("=" * 70)
    print()

    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"File not found: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)
    total_records = len(df)
    print(f"Total NASA records loaded: {total_records}")

    equinox_matches = []

    for _, row in df.iterrows():
        y = int(row["year"])
        m = int(row["month"])
        d = int(row["day"])
        t0 = float(row["t0"]) if pd.notna(row["t0"]) else 12.0

        eclipse_jd = calendar_to_jd(y, m, d, t0)
        eq_jd = vernal_equinox_jde(y)

        diff_days = abs(eclipse_jd - eq_jd)

        if diff_days <= EQUINOX_THRESHOLD_DAYS:
            equinox_matches.append({
                "year": y,
                "month": m,
                "day": d,
                "eclipse_type": row["eclipse_type"],
                "cat_no": row.get("cat_no", np.nan),
                "saros": row.get("saros", np.nan),
                "magnitude": row.get("magnitude", np.nan),
                "t0_hour": t0,
                "diff_hours": round(diff_days * 24.0, 2),
            })

    results_df = pd.DataFrame(equinox_matches)
    matches_count = len(results_df)

    print()
    print("SCAN RESULTS:")
    print(f"Total NASA eclipses scanned : {total_records}")
    print(f"Equinox eclipses found      : {matches_count}")
    if matches_count > 0:
        freq = total_records / matches_count
        years_rec = 5000.0 / matches_count
        print(f"Empirical frequency         : 1 in every {freq:.2f} eclipses")
        print(f"Mean recurrence interval    : ~{years_rec:.2f} years")

    print()
    print("SAMPLE EQUINOX MATCHES (FIRST 15):")
    display_cols = ["year", "month", "day", "eclipse_type", "diff_hours", "magnitude"]
    print(results_df[display_cols].head(15).to_string(index=False))

    results_df.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")
    print()
    print(f"Full results saved to: {OUTPUT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    scan_nasa_equinox()
