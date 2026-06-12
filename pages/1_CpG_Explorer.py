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

tumor_type = st.sidebar.selectbox(
    "Tumor Type",
    ["COAD", "LUAD", "LUSC", "LIHC"]
)

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

min_HI = st.sidebar.slider(
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
    ts.pan_tumor_median,
    ts.pan_normal_median,
    ts.dispersion_index,
    cf.leukocyte_median

FROM tumor_summary ts

JOIN cpg_annotation ca
    ON ts.site_id = ca.site_id

JOIN cpg_features cf
    ON ts.site_id = cf.site_id

WHERE ts.tumor_type = '{tumor_type}'
AND ts.delta_median >= {min_delta}
AND ts.normal_median <= {max_normal_median}
AND ts.pan_tumor_median <= {max_pancan_t}
AND cf.leukocyte_median <= {max_leuco}
AND ts.dispersion_index >= {min_HI}

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

fig = px.scatter(
    df,
    x="delta_median",
    y="dispersion_index",
    hover_data=["site_id", "gene"]
)

st.plotly_chart(
    fig,
    use_container_width=True
)

st.subheader("CpG density heatmap")

# bins
delta_bins = np.arange(0, 1.05, 0.1)
hi_bins = np.arange(0, 4.5, 0.5)

# agrupar CpGs
heatmap_df = (
    df.assign(
        delta_bin=pd.cut(
            df["delta_median"],
            bins=delta_bins,
            include_lowest=True
        ),
        hi_bin=pd.cut(
            df["dispersion_index"],
            bins=hi_bins,
            include_lowest=True
        )
    )
    .groupby(
        ["hi_bin", "delta_bin"],
        observed=False
    )
    .size()
    .reset_index(name="n_sites")
)

# matriz
pivot = heatmap_df.pivot(
    index="hi_bin",
    columns="delta_bin",
    values="n_sites"
).fillna(0)

# etiquetas legibles
x_labels = [str(x) for x in pivot.columns]
y_labels = [str(y) for y in pivot.index]

# heatmap
fig_heat = go.Figure(
    data=go.Heatmap(
        z=pivot.values,
        x=x_labels,
        y=y_labels,
        hovertemplate=
            "Δβ: %{x}<br>" +
            "HI: %{y}<br>" +
            "CpGs: %{z}<extra></extra>"
    )
)

fig_heat.update_layout(
    title="Number of CpGs by Δβ and dispersion index",
    xaxis_title="Delta β",
    yaxis_title="Dispersion Index (HI)",
    height=600
)

st.plotly_chart(
    fig_heat,
    use_container_width=True
)

delta_thresholds = np.arange(0, 1.05, 0.05)
hi_thresholds = np.arange(0, 4.05, 0.1)

matrix = []

for hi in hi_thresholds:
    row = []

    for delta in delta_thresholds:

        n_sites = len(
            df[
                (df["delta_median"] >= delta)
                & (df["dispersion_index"] >= hi)
            ]
        )

        row.append(n_sites)

    matrix.append(row)

fig_heat = go.Figure(
    data=go.Heatmap(
        z=matrix,
        x=np.round(delta_thresholds, 2),
        y=np.round(hi_thresholds, 2),
        hovertemplate=
            "Δβ ≥ %{x}<br>" +
            "HI ≥ %{y}<br>" +
            "CpGs: %{z}<extra></extra>"
    )
)

fig_heat.update_layout(
    title="CpGs retained after threshold filtering",
    xaxis_title="Minimum Δβ",
    yaxis_title="Minimum HI"
)

st.plotly_chart(
    fig_heat,
    use_container_width=True
)

st.subheader("Genes with most candidate CpGs")

gene_counts = (
    df.groupby("gene")
      .size()
      .reset_index(name="n_sites")
      .sort_values(
          "n_sites",
          ascending=False
      )
      .head(30)
)

fig_genes = px.bar(
    gene_counts,
    x="gene",
    y="n_sites",
    labels={
        "gene": "Gene",
        "n_sites": "Candidate CpGs"
    }
)

st.plotly_chart(
    fig_genes,
    use_container_width=True
)

st.dataframe(
    gene_counts,
    use_container_width=True
)

gene_score = (
    df.groupby("gene")
      .agg(
          n_sites=("site_id","count"),
          mean_delta=("delta_median","mean"),
          mean_HI=("dispersion_index","mean")
      )
      .reset_index()
)

gene_score["score"] = (
    gene_score["n_sites"]
    * gene_score["mean_delta"]
    * gene_score["mean_HI"]
)

gene_score = gene_score.sort_values(
    "score",
    ascending=False
)

st.subheader("DMR Analysis")

if len(df) > 0:
    dmr_df = (
        df
        .sort_values(
            ["chr", "start_pos"]
        )
        .copy()
    )

    region_ids = []

    region = 0
    prev_chr = None
    prev_pos = None

    for row in dmr_df.itertuples():

        new_region = False

        if prev_chr is None:
            new_region = True

        elif row.chr != prev_chr:
            new_region = True

        elif (row.start_pos - prev_pos) > max_distance:
            new_region = True

        if new_region:
            region += 1

        region_ids.append(region)

        prev_chr = row.chr
        prev_pos = row.start_pos

    dmr_df["region_id"] = region_ids

    dmrs = (
        dmr_df
        .groupby("region_id")
        .agg(
            chr=("chr", "first"),
            start=("start_pos", "min"),
            end=("start_pos", "max"),
            n_cpgs=("site_id", "count"),
            mean_delta=("delta_median", "mean"),
            max_delta=("delta_median", "max"),
            mean_hi=("dispersion_index", "mean"),
            mean_leukocyte=("leukocyte_median", "mean"),
            genes=("gene",
                lambda x: ";".join(
                    sorted(set(x.dropna()))
                )
            )
        )   
       .reset_index()
    )

    dmrs["region_size"] = (
        dmrs["end"] - dmrs["start"]
    )

    dmrs = dmrs[
        dmrs["n_cpgs"] >= min_cpgs_region
    ]

    dmrs["score"] = (
        dmrs["mean_delta"]
        * dmrs["mean_hi"]
        * np.log2(dmrs["n_cpgs"] + 1)
    )

    dmrs = dmrs.sort_values(
        "score",
        ascending=False
    )

    st.metric(
        "Candidate DMRs",
        len(dmrs)
    )

    st.dataframe(
        dmrs,
        use_container_width=True
    )

    fig_dmrs = px.bar(
        dmrs.head(20),
        x="genes",
        y="score",
        hover_data=[
            "n_cpgs",
            "mean_delta",
            "mean_hi",
            "region_size"
        ]
    )

    fig_dmrs.update_layout(
        title="Top Candidate DMRs"
    )

    st.plotly_chart(
        fig_dmrs,
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
        mean_HI=("dispersion_index", "mean"),
        max_delta=("delta_median", "max"),
        max_HI=("dispersion_index", "max"),
        mean_leukocyte=("leukocyte_median", "mean")
    )
    .reset_index()
)

gene_summary["score"] = (
    gene_summary["n_sites"]
    * gene_summary["mean_delta"]
    * gene_summary["mean_HI"]
)

gene_summary = gene_summary.sort_values(
    "score",
    ascending=False
)
fig_bubble = px.scatter(
    gene_summary,

    x="mean_delta",
    y="mean_HI",

    size="n_sites",
    color="score",

    hover_name="gene",

    hover_data={
        "n_sites": True,
        "mean_delta": ":.3f",
        "mean_HI": ":.3f",
        "mean_leukocyte": ":.3f",
        "score": ":.3f"
    },

    labels={
        "mean_delta": "Mean Δβ",
        "mean_HI": "Mean Dispersion Index",
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

st.write(f"Rows returned: {len(df)}")
