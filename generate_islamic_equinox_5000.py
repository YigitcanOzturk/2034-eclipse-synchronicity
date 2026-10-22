# ============================================================
# 5,000-YEAR ISLAMIC CALENDAR (29 DHU AL-HIJJAH) vs. VERNAL EQUINOX SCAN
# -1999 to +3000 (Astronomical Years)
#
# Methodology:
#   1) Jean Meeus (Astronomical Algorithms, Ch. 27): Vernal Equinox JDE
#   2) Reingold & Dershowitz (Calendrical Calculations): Tabular Islamic Calendar
#
# Output:
#   islamic_equinox_5000.csv
# ============================================================

import os
import math
import numpy as np
import pandas as pd

OUTPUT_FILE = r"islamic_equinox_5000.csv"

# ============================================================
# 1. VERNAL EQUINOX (JEAN MEEUS)
# ============================================================

def vernal_equinox_jde(year):
    if -1000 <= year <= 1000:
        Y = year / 1000.0
        jde0 = (1721139.29189 + 365242.13740 * Y + 0.06134 * (Y ** 2) + 0.00111 * (Y ** 3) - 0.00071 * (Y ** 4))
    elif 1000 < year <= 3000:
        Y = (year - 2000) / 1000.0
        jde0 = (2451623.80984 + 365242.37404 * Y + 0.05169 * (Y ** 2) - 0.00411 * (Y ** 3) - 0.00057 * (Y ** 4))
    else:
        Y = year / 1000.0
        jde0 = (1721139.29189 + 365242.13740 * Y + 0.06134 * (Y ** 2) + 0.00111 * (Y ** 3) - 0.00071 * (Y ** 4))

    T = (jde0 - 2451545.0) / 36525.0
    W = 35999.373 * T - 2.47
    delta_lambda = 1.0 + 0.0334 * math.cos(math.radians(W)) + 0.0007 * math.cos(math.radians(2.0 * W))

    S = (
        485.0 * math.cos(math.radians(324.96 + 1934.136 * T))
        + 203.0 * math.cos(math.radians(276.11 + 72001.539 * T))
        + 199.0 * math.cos(math.radians(119.75 + 1318.515 * T))
        + 182.0 * math.cos(math.radians(263.36 + 34575.405 * T))
        + 156.0 * math.cos(math.radians(156.46 + 65928.934 * T))
        + 136.0 * math.cos(math.radians(41.85 + 20186.439 * T))
        + 77.0 * math.cos(math.radians(246.54 + 2924.669 * T))
        + 74.0 * math.cos(math.radians(320.61 + 71498.441 * T))
        + 70.0 * math.cos(math.radians(107.94 + 54994.431 * T))
        + 58.0 * math.cos(math.radians(278.57 + 6077.581 * T))
        + 52.0 * math.cos(math.radians(87.15 + 32243.061 * T))
        + 50.0 * math.cos(math.radians(26.80 + 73005.807 * T))
    )

    return jde0 + (0.00001 * S) / delta_lambda


# ============================================================
# 2. TABULAR ISLAMIC CALENDAR (REINGOLD & DERSHOWITZ)
# ============================================================

def islamic_leap_year(h_year):
    """30-year leap cycle (11 leap years: 355 days, 19 common: 354 days)."""
    return ((11 * h_year + 14) % 30) < 11

def islamic_new_year_rd(h_year):
    """Rata Die (RD) of 1 Muharram. Epoch: RD 227015 (Friday, 16 July 622 CE)."""
    return 227014 + (h_year - 1) * 354 + (3 + 11 * h_year) // 30

def islamic_dhu_al_hijjah_29_rd(h_year):
    """
    Rata Die of 29 Dhu al-Hijjah (Day 353 of a 354-day common year,
    or Day 353/354 before the new year).
    Equivalent to (1 Muharram of h_year + 1) - (1 or 2 days).
    """
    next_muharram = islamic_new_year_rd(h_year + 1)
    # If h_year is leap, Dhu al-Hijjah has 30 days; 29th is 2 days before next Muharram.
    # If common, 29th is 1 day before next Muharram.
    offset = 2 if islamic_leap_year(h_year) else 1
    return next_muharram - offset

def rd_to_jd(rd):
    return rd + 1721424.5


# ============================================================
# 3. SCAN & GENERATE 5,000-YEAR CSV
# ============================================================

def generate_islamic_scan():
    print("Generating 5,000-year Islamic calendar vs. Equinox dataset (-1999 to +3000)...")
    records = []

    for y in range(-1999, 3001):
        eq_jde = vernal_equinox_jde(y)
        eq_rd = eq_jde - 1721424.5

        # Estimate the Hijri year ending near the equinox of astronomical year y
        h_est = int((eq_rd - 227015) / 354.367) + 1

        best_h = None
        min_diff = 999999.0

        for h_cand in range(h_est - 2, h_est + 3):
            cand_rd = islamic_dhu_al_hijjah_29_rd(h_cand)
            diff = cand_rd - eq_rd
            if abs(diff) < abs(min_diff):
                min_diff = diff
                best_h = h_cand

        adar_rd = islamic_dhu_al_hijjah_29_rd(best_h)
        adar_jd = rd_to_jd(adar_rd)

        records.append({
            "astronomical_year": y,
            "hijri_year": best_h,
            "is_leap_year": islamic_leap_year(best_h),
            "equinox_jde": round(eq_jde, 4),
            "dhu_al_hijjah_29_jd": round(adar_jd, 4),
            "diff_days": round(min_diff, 3),
            "exact_same_day": abs(min_diff) <= 1.0
        })

    df = pd.DataFrame(records)
    df.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

    total = len(df)
    same_day_matches = df[df["exact_same_day"] == True]

    print("=" * 65)
    print(f"Total years generated       : {total}")
    print(f"Matches on Equinox (<= 1 day): {len(same_day_matches)}")
    print(f"Empirical recurrence        : 1 in every {total / len(same_day_matches):.2f} years")
    print(f"Output saved to             : {OUTPUT_FILE}")
    print("=" * 65)

if __name__ == "__main__":
    generate_islamic_scan()
