import streamlit as st
import plotly.graph_objects as go

from db.queries import run_query

st.set_page_config(
    page_title="Gene Explorer",
    layout="wide"
)

st.title("Gene Explorer")

# ==========================================
# SIDEBAR
# ==========================================

st.sidebar.header("Filters")

tumor_type = st.sidebar.selectbox(
    "Tumor Type",
    ["COAD", "LUAD", "LUSC", "LIHC"]
)

gene = st.sidebar.text_input(
    "Gene Symbol",
    "SEPTIN9"
)

# ==========================================
# QUERY
# ==========================================

query_curve = f"""
SELECT
    ca.site_id,
    ca.gene,
    ca.start_pos,

    ts.delta_median,
    ts.dispersion_index,
    ts.tumor_median,
    ts.normal_median,

    ts.pan_tumor_median,
    ts.pan_normal_median,

    cf.leukocyte_median,

    ec.spearman_r

FROM tumor_summary ts

JOIN cpg_annotation ca
    ON ts.site_id = ca.site_id

JOIN cpg_features cf
    ON ts.site_id = cf.site_id

LEFT JOIN expression_correlation ec
    ON ts.site_id = ec.site_id
    AND ec.tumor_type = ts.tumor_type


WHERE ca.gene ILIKE '%{gene}%'
AND ts.tumor_type = '{tumor_type}'

ORDER BY ca.start_pos
LIMIT 300
"""

graph = run_query(query_curve)

#------------------- Grafico Curvas------------------------
site = graph["site_id"]
positions = graph["start_pos"]
median_tumor = graph["tumor_median"]
median_normal = graph["normal_median"]
pan_tumor = graph["pan_tumor_median"]
pan_normal = graph["pan_normal_median"]
delta_median = graph["delta_median"]
HI = graph["dispersion_index"]
leukocytes = graph["leukocyte_median"]
expression = graph["spearman_r"]

hover_text = (
    graph["site_id"] + " | " +
    graph["start_pos"].astype(str)
)

fig2 = go.Figure()
# =========================================================
# LEYENDAS
# =========================================================
curves = [
    ("Median Type NT", median_normal, "rgba(255,255,0,1)", "median_normal", "line"),
    ("Median PanCan T", pan_tumor, "rgba(0,255,0,1)", "pan_tumor", "line"),
    ("Median PanCan NT", pan_normal, "rgba(0,128,0,1)", "pan_normal", "line"),
    ("Median Leukocytes", leukocytes, "rgba(130,168,240,1)", "leukocytes", "line"),
  #  ("Expression", expression, "rgba(127,127,127,1)", "expression", "line"),
    ("Delta", delta_median, "rgba(236,0,139,1)", "delta", "marker"),
    ("HI", HI, "rgba(255,0,0,1)", "HI", "marker")
]

# =========================================================
# CURVAS
# =========================================================

fig2.add_trace(
    go.Scatter(
        x=positions,
        y=median_tumor,
        mode="lines",
        name="Median Type T",
        customdata=graph["site_id"],
        line=dict(
            width=2,
            color="rgba(255, 165, 0,1)",
            shape="spline",
            smoothing=1.3
        ),
        fill="tozeroy",
        fillcolor="rgba(255, 165, 0,1)"
    )
)

for name, values, color, group, trace_type in curves:
    if trace_type == "line":
        fig2.add_trace(
            go.Scatter(
                x=positions,
                y=values,
                mode="lines",
                name=name,
                line=dict(
                    color=color,
                    width=2,
                    shape="spline",
                    smoothing=1.3
                )
            )
        )

    else:
        fig2.add_trace(
            go.Scatter(
                x=positions,
                y=values,
                mode="markers",
                yaxis="y2",
                name=name,
                marker=dict(
                    color=color,
                    size=8
                )
            )
        )

# =========================================================
# LAYOUT
# =========================================================
show_ih = st.checkbox("Mostrar IH", value=True)

if show_ih:
    fig2.add_trace(
        go.Scatter(
            x=positions,
            y=HI,
            mode="markers",
            marker=dict(
                color="rgba(255,0,0,1)",
                size=8
            ),
            name="HI",
            yaxis="y2"
        )
    )
    y2_max = HI.max()
else:
    y2_max = 1

fig2.update_traces(
    hovertext=hover_text,
    hovertemplate=
        "%{hovertext}<br>" +
        "%{fullData.name}: %{y:.3f}" +
        "<extra></extra>"
)

fig2.update_layout(
    template="plotly_white",
    height=500,
    hovermode="x unified",
    margin=dict(
        l=40,
        r=40,
        t=50,
        b=50
    ),

    legend=dict(
        orientation="v",
        y=1,
        x=1.1,
        bgcolor="rgba(255,255,255,0)",
        borderwidth=0,
        groupclick="togglegroup",
        tracegroupgap=1,
        font=dict(
            family="Arial Narrow",
            size=18,
            color="black"
        )
    ),

    title=dict(
        text=f"Aberrant tumor hypermethylation area for <b>{gene}</b>  - {tumor_type}",
        x=0,
        xanchor="left",

        font=dict(
            family="Arial Narrow",
            size=20,
            color="black"
        )
    ),

    xaxis=dict(
        title=dict(
            text="Genomic Position",
            font=dict(
                family="Arial Narrow",
                color="black",
                size=18
            )
        ),
        tickfont=dict(
            family="Arial Narrow",
            color="black",
            size=18
        ),
        tickformat=",d",  
        rangeslider=dict(visible=True),
        showgrid=True,
        unifiedhovertitle=dict(
            text="<b>%{customdata[0]}</b> | %{x:,d}"
        )
    ),

    yaxis=dict(
        range=[-0.2, 1],
        title=dict(
            text="Methylation",
            font=dict(
                family="Arial Narrow",
                color="black",
                size=18
            )
        ),
        tickfont=dict(
            family="Arial Narrow",
            color="black",
            size=18
        ),
        showgrid=True
    ),
    yaxis2=dict(
        title=dict(
            text="Delta/HI",
            font=dict(
                family="Arial Narrow",
                color="black",
                size=18
            )
        ),
        overlaying="y",
        side="right",
        range=[-0.2, y2_max],
        tickfont=dict(
            family="Arial Narrow",
            color="black",
            size=18
        ),
        )
)

config = {
    "toImageButtonOptions": {
        "format": "svg",
        "filename": f"{gene}_{tumor_type}",
        "height": 600,
        "width": 1400,
        "scale": 1
    }
}

st.plotly_chart(
    fig2,
    use_container_width=True,
    config=config
)

st.dataframe(
    graph,
    use_container_width=True
)

#st.write(f"Rows returned: {len(df)}")
