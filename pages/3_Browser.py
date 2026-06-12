import streamlit as st
from services.ensembl import get_sequence
from db.queries import run_query

st.title("Test Gene Browser")

genes = run_query("""
SELECT DISTINCT gene_symbol
FROM gene_annotation
WHERE gene_symbol IS NOT NULL
ORDER BY gene_symbol
""")

gene = st.sidebar.selectbox(
    "Gen",
    genes["gene_symbol"]
)


gene_info = run_query(f"""
SELECT *
FROM gene_annotation
WHERE gene_symbol = '{gene}'
LIMIT 1
""")
row = gene_info.iloc[0]
df_cpg = run_query(f"""
SELECT
    site_id,
    start_pos,
    distance_tss,
    cgi,
    cg_island
FROM cpg_annotation
WHERE gene = '{gene}'
ORDER BY start_pos
""")

df_cpg["label"] = (
    df_cpg["site_id"].astype(str)
    + " | "
    + df_cpg["start_pos"].astype(str)
)
if df_cpg.empty:
    st.warning("No hay CpGs asociados")
    st.stop()

selected_idx = st.sidebar.selectbox(
    "CpG",
    options=df_cpg.index,
    format_func=lambda i:
        f"{df_cpg.loc[i,'site_id']} | {int(df_cpg.loc[i,'start_pos']):,}"
)

selected_row = df_cpg.loc[selected_idx]

selected_pos = int(
    selected_row["start_pos"]
)
window = st.slider(
    "Ventana (bp)",
    200,
    5000,
    1000,
    100
)

seq_start = (
    selected_pos -
    window // 2
)

seq_end = (
    selected_pos +
    window // 2
)

st.info(
    f"{selected_row['site_id']} | "
    f"{selected_pos:,}"
)

chrom = str(row["chr"])
start = int(row["start_pos"])
end = int(row["end_pos"])

st.write(
    f"Cromosoma: {chrom} | "
    f"Start: {start:,} | "
    f"End: {end:,}"
)

sequence = get_sequence(
    chrom,
    seq_start,
    seq_end
)

import re

real_cpgs = []

for m in re.finditer(
    "CG",
    sequence
):

    real_cpgs.append(
        seq_start + m.start()
    )

st.code(sequence[:500])

st.metric(
    "CpGs reales",
    len(real_cpgs)
)

gcgc_sites = [
    seq_start + m.start()
    for m in re.finditer(
        "GCGC",
        sequence
    )
]

st.metric(
    "Motivos GCGC",
    len(gcgc_sites)
)
