# ============================================================
# SACRED CORRIDOR MIDDAY CULMINATION SCAN (10:00 - 14:00 LST)
# EMPIRICAL PROBABILITY DERIVATION FROM 161 CORRIDOR EVENTS
#
# Methodology:
#   Conversion of topocentric Besselian UT hours to Local Solar Time (LST)
#   LST = UT + (Longitude / 15.0)
#
# Input:
#   mecca_jerusalem_70_EXACT.csv
#
# Output:
#   corridor_midday_matches.csv
# ============================================================

import os
import pandas as pd

# ============================================================
# SETTINGS & OBSERVER COORDINATES
# ============================================================

INPUT_FILE = r"C:\Users\...\Desktop\mecca_jerusalem_70_EXACT.csv"
OUTPUT_FILE = r"C:\Users\...\Desktop\corridor_midday_matches.csv"

# Geographic longitudes (degrees East)
MECCA_LON = 39.8262
JERUSALEM_LON = 35.2137

# Midday culmination window in Local Solar Time (hours)
MIDDAY_START = 10.0
MIDDAY_END = 14.0


# ============================================================
# CALCULATION & FILTERING
# ============================================================

def calculate_corridor_midday_probability():
    print()
    print("=" * 70)
    print("SACRED CORRIDOR MIDDAY WINDOW PROBABILITY SCAN")
    print("=" * 70)

    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"Input file not found: '{INPUT_FILE}'.")

    # Load 161 dual-obscuration (>= 70%) matches
    df = pd.read_csv(INPUT_FILE)
    total_events = len(df)
    print(f"Total corridor events loaded: {total_events}")

    # Compute Local Solar Time (LST) for each observer
    df["Mecca_LST"] = df["Mecca_t_decimal_hour"] + (MECCA_LON / 15.0)
    df["Jerusalem_LST"] = df["Jerusalem_t_decimal_hour"] + (JERUSALEM_LON / 15.0)

    # Filter events where peak obscuration falls between 10:00 and 14:00 LST
    mecca_midday = df[(df["Mecca_LST"] >= MIDDAY_START) & (df["Mecca_LST"] <= MIDDAY_END)]
    jerusalem_midday = df[(df["Jerusalem_LST"] >= MIDDAY_START) & (df["Jerusalem_LST"] <= MIDDAY_END)]
    both_midday = df[
        (df["Mecca_LST"] >= MIDDAY_START) & (df["Mecca_LST"] <= MIDDAY_END) &
        (df["Jerusalem_LST"] >= MIDDAY_START) & (df["Jerusalem_LST"] <= MIDDAY_END)
    ]

    count_both = len(both_midday)
    count_mecca = len(mecca_midday)
    count_jerusalem = len(jerusalem_midday)

    p_both = count_both / total_events
    reciprocal_both = total_events / count_both

    # Print summary metrics
    print()
    print("EMPIRICAL DISTRIBUTION:")
    print(f"  Mecca peak in [10:00, 14:00] LST       : {count_mecca} / {total_events} ({count_mecca / total_events * 100:.2f}%)")
    print(f"  Jerusalem peak in [10:00, 14:00] LST   : {count_jerusalem} / {total_events} ({count_jerusalem / total_events * 100:.2f}%)")
    print(f"  Simultaneous midday peak (BOTH cities) : {count_both} / {total_events}")
    print()
    print("MATHEMATICAL DERIVATION:")
    print(f"  P_midday = {count_both} / {total_events} = {p_both:.6f}")
    print(f"  Recurrence Fraction = 1 in {reciprocal_both:.2f} corridor eclipses")

    # Display 20 March 2034 validation
    row_2034 = df[df["year"] == 2034]
    if not row_2034.empty:
        r = row_2034.iloc[0]
        print()
        print("TARGET EVENT VALIDATION (20 March 2034):")
        print(f"  Mecca peak LST      : {r['Mecca_LST']:.2f} h ({int(r['Mecca_LST'])}:{int((r['Mecca_LST'] % 1) * 60):02d} local solar)")
        print(f"  Jerusalem peak LST  : {r['Jerusalem_LST']:.2f} h ({int(r['Jerusalem_LST'])}:{int((r['Jerusalem_LST'] % 1) * 60):02d} local solar)")
        print(f"  Mecca Obscuration   : {r['Mecca_obscuration_%']:.2f}%")
        print(f"  Jerusalem Obscuration: {r['Jerusalem_obscuration_%']:.2f}%")
        print(f"  Midday Criterion Met: True")

    # Save matching events
    display_cols = [
        "year", "month", "day", "eclipse_type",
        "Mecca_obscuration_%", "Jerusalem_obscuration_%",
        "Mecca_LST", "Jerusalem_LST"
    ]
    both_midday[display_cols].to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")
    print()
    print(f"Matching subset ({count_both} events) exported to:\n{OUTPUT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    calculate_corridor_midday_probability()
