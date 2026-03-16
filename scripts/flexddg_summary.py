#!/usr/bin/env python3
"""
summarize_flexddg_local.py

Output:
    - flexddg_summary.csv
    - Console printout of ΔΔG table
"""

import sqlite3, glob, os
import pandas as pd

def extract_avg_scores(db_path):
    mut_tag = os.path.basename(db_path).replace("_ddG.db3", "")
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    try:
        total_id = cur.execute(
            "SELECT score_type_id FROM score_types WHERE score_type_name='total_score'"
        ).fetchone()[0]
    except TypeError:
        conn.close()
        return None

    query = """
        SELECT b.name AS state, s.score_value
        FROM structure_scores s
        JOIN batches b USING(batch_id)
        WHERE s.score_type_id = ?
    """
    rows = cur.execute(query, (total_id,)).fetchall()
    conn.close()

    score_lists = {}
    for r in rows:
        tag = r['state'].replace('_dbreport', '')
        score_lists.setdefault(tag, []).append(r['score_value'])

    required = ['bound_wt', 'unbound_wt', 'bound_mut', 'unbound_mut']
    if not all(k in score_lists for k in required):
        return None

    avg = {k: sum(v)/len(v) for k, v in score_lists.items()}
    ddg = (avg['bound_mut'] - avg['unbound_mut']) - (avg['bound_wt'] - avg['unbound_wt'])

    return {
        "mutation": mut_tag,
        "ddG_REU": round(ddg, 3),
        **{k: round(v, 3) for k, v in avg.items()}
    }

def main():
    db_files = glob.glob("*_ddG.db3")
    if not db_files:
        print("❌ No *_ddG.db3 files found.")
        return

    rows = []
    for db in db_files:
        result = extract_avg_scores(db)
        if result:
            rows.append(result)
        else:
            print(f"⚠️ Skipped: {db} (missing scores)")

    if not rows:
        print("❌ No valid ΔΔG values extracted.")
        return

    df = pd.DataFrame(rows).sort_values("ddG_REU")
    df.to_csv("flexddg_summary.csv", index=False)
    print("\n✅ Saved: flexddg_summary.csv\n")
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
