# ============================================================
# 5,000-YEAR HEBREW CALENDAR (29 ADAR) vs. VERNAL EQUINOX SCAN
# -1999 to +3000 (Astronomical Years)
#
# Methodology:
#   1) Jean Meeus (Astronomical Algorithms, Ch. 27): Vernal Equinox JDE
#   2) Reingold & Dershowitz (Calendrical Calculations): Hebrew Calendar
#
# Output:
#   hebrew_equinox_5000.csv
# ============================================================

import os
import math
import numpy as np
import pandas as pd

OUTPUT_FILE = r"hebrew_equinox_5000.csv"

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
# 2. HEBREW CALENDAR (REINGOLD & DERSHOWITZ)
# ============================================================

def hebrew_leap_year(h_year):
    return ((7 * h_year + 1) % 19) < 7

def hebrew_calendar_elapsed_days(h_year):
    months_elapsed = (235 * h_year - 234) // 19
    parts_elapsed = 12084 + 13753 * months_elapsed
    days = 29 * months_elapsed + parts_elapsed // 25920
    parts = parts_elapsed % 25920

    if (parts >= 19440) or \
       (((days % 7) == 2) and (parts >= 9924) and (not hebrew_leap_year(h_year))) or \
       (((days % 7) == 1) and (parts >= 16789) and hebrew_leap_year(h_year - 1)):
        alt_day = days + 1
    else:
        alt_day = days

    if (alt_day % 7) in [0, 3, 5]:
        return alt_day + 1
    return alt_day

def hebrew_new_year_rd(h_year):
    """Rata Die (RD) of 1 Tishrei."""
    return -1373427 + hebrew_calendar_elapsed_days(h_year)

def days_in_hebrew_year(h_year):
    return hebrew_new_year_rd(h_year + 1) - hebrew_new_year_rd(h_year)

def adar_end_rd(h_year):
    """
    Computes the RD of the terminal day of the Hebrew year:
    29 Adar in common years, or 29 Adar II in embolismic (leap) years.
    This corresponds to the day immediately preceding 1 Nisan.
    """
    tishrei_1 = hebrew_new_year_rd(h_year)
    year_len = days_in_hebrew_year(h_year)
    is_leap = hebrew_leap_year(h_year)

    # Days from 1 Tishrei to 1 Nisan:
    # Tishrei (30) + Cheshvan (29/30) + Kislev (29/30) + Tevet (29) + Shevat (30) + Adar(s) (29 or 30+29)
    cheshvan_len = 30 if (year_len % 10 == 5) else 29
    kislev_len = 29 if (year_len % 10 == 3) else 30
    tevet_len = 29
    shevat_len = 30
    adar_total_len = (30 + 29) if is_leap else 29

    # Day of 29 Adar (or 29 Adar II) is (1 Nisan - 1 day)
    days_to_adar_end = 30 + cheshvan_len + kislev_len + tevet_len + shevat_len + adar_total_len - 1
    return tishrei_1 + days_to_adar_end


def rd_to_jd(rd):
    """Converts Rata Die (RD) to Julian Day (JD at 12:00 UT)."""
    return rd + 1721424.5


# ============================================================
# 3. SCAN & GENERATE 5,000-YEAR CSV
# ============================================================

def generate_scan():
    print("Generating 5,000-year Hebrew calendar vs. Equinox dataset (-1999 to +3000)...")
    records = []

    # In AM reckoning: Year 1 AD is 3761 AM. Year 2034 AD is 5794 AM.
    # Astronomical year y maps to Hebrew year (y + 3760) for Nisan/Adar.
    for y in range(-1999, 3001):
        h_year = y + 3760

        eq_jde = vernal_equinox_jde(y)
        adar_rd = adar_end_rd(h_year)
        adar_jd = rd_to_jd(adar_rd)

        # Difference in days between 29 Adar (midday UT) and astronomical equinox
        diff_days = adar_jd - eq_jde
        is_leap = hebrew_leap_year(h_year)
        metonic_pos = h_year % 19

        records.append({
            "astronomical_year": y,
            "hebrew_year_AM": h_year,
            "is_leap_year": is_leap,
            "metonic_cycle_slot": metonic_pos,
            "equinox_jde": round(eq_jde, 4),
            "adar_end_jd": round(adar_jd, 4),
            "diff_days": round(diff_days, 3),
            "exact_same_day": abs(diff_days) <= 1.0
        })

    df = pd.DataFrame(records)
    df.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

    total = len(df)
    same_day_matches = df[df["exact_same_day"] == True]

    print("=" * 65)
    print(f"Total years generated       : {total}")
    print(f"Matches on Equinox (<= 1 day): {len(same_day_matches)}")
    print(f"Empirical frequency         : 1 in every {total / len(same_day_matches):.2f} years")
    print(f"Output saved to             : {OUTPUT_FILE}")
    print("=" * 65)


if __name__ == "__main__":
    generate_scan()
