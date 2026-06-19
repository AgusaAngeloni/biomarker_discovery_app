import streamlit as st
import plotly.express as px
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from db.queries import run_query

st.set_page_config(
    page_title="CpG Explorer",
    layout="wide"
)

st.title("CpG Explorer")

# ==========================================
# SIDEBAR
# ==========================================

st.sidebar.header("Filters")

tumor_map = {
    "CRC (COAD)": "COAD",
    "HCC (LIHC)": "LIHC",
    "LUAD": "LUAD",
    "LUSC": "LUSC"
}

tumor_label = st.sidebar.selectbox(
    "Tumor Type",
    list(tumor_map.keys())
)

tumor_type = tumor_map[tumor_label]

min_delta = st.sidebar.slider(
    "Minimum Delta",
    0.0,
    1.0,
    0.4
)

max_normal_median = st.sidebar.slider(
    "Max Median NT",
    0.0,
    1.0,
    0.08
)

max_pancan_t = st.sidebar.slider(
    "Max Median PanCancer T",
    0.0,
    1.0,
    0.08
)

max_leuco = st.sidebar.slider(
    "Max Median Leukocytes",
    0.0,
    1.0,
    0.08
)

min_hi = st.sidebar.slider(
    "Min HI",
    0.0,
    4.0,
    1.0
)
st.sidebar.markdown("---")

max_distance = st.sidebar.slider(
    "DMR max distance (bp)",
    100,
    5000,
    500,
    step=100
)

min_cpgs_region = st.sidebar.slider(
    "Minimum CpGs per DMR",
    2,
    20,
    3
)
# ==========================================
# QUERY
# ==========================================

query = f"""
SELECT
    ca.site_id,
    ca.gene,
    ca.start_pos,
    ca.chr,
    ca.distance_tss,
    ts.tumor_type,
    ts.delta_median,
    ts.tumor_median,
    ts.normal_median,
    ts.pantumor_median,
    ts.pannormal_median,
    ts.hi_index,
    cf.leukocyte_median

FROM tumor_summary ts

JOIN cpg_annotation ca
    ON ts.site_id = ca.site_id

JOIN cpg_features cf
    ON ts.site_id = cf.site_id

WHERE ts.tumor_type = '{tumor_type}'
AND ts.delta_median >= {min_delta}
AND ts.normal_median <= {max_normal_median}
AND ts.pantumor_median <= {max_pancan_t}
AND cf.leukocyte_median <= {max_leuco}
AND ts.hi_index >= {min_hi}

ORDER BY ts.delta_median DESC
LIMIT 500
"""

df = run_query(query)

# ==========================================
# RESULTS
# ==========================================

st.subheader("Biomarker Panel")

st.dataframe(
    df,
    use_container_width=True
)



st.subheader("Gene-level Biomarker Landscape")

gene_summary = (
    df
    .dropna(subset=["gene"])
    .groupby("gene")
    .agg(
        n_sites=("site_id", "count"),
        mean_delta=("delta_median", "mean"),
        mean_hi=("hi_index", "mean"),
        max_delta=("delta_median", "max"),
        max_hi=("hi_index", "max"),
        mean_leukocyte=("leukocyte_median", "mean")
    )
    .reset_index()
)

gene_summary["score"] = (
    gene_summary["n_sites"]
    * gene_summary["mean_delta"]
    * gene_summary["mean_hi"]
)

gene_summary = gene_summary.sort_values(
    "score",
    ascending=False
)

fig_bubble = px.scatter(
    gene_summary,

    x="mean_delta",
    y="mean_hi",

    size="n_sites",
    color="score",

    hover_name="gene",

    hover_data={
        "n_sites": True,
        "mean_delta": ":.3f",
        "mean_hi": ":.3f",
        "mean_leukocyte": ":.3f",
        "score": ":.3f"
    },

    labels={
        "mean_delta": "Mean Δβ",
        "mean_hi": "Mean HI Index",
        "n_sites": "Candidate CpGs"
    }
)

fig_bubble.update_layout(
    title=f"Gene-level methylation landscape ({tumor_type})",
    height=700
)

st.plotly_chart(
    fig_bubble,
    use_container_width=True
)




