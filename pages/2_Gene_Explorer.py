import html
import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from db.queries import run_query
from services.ensembl import get_sequence


# ============================================================
# Page config
# ============================================================

st.set_page_config(
    page_title="Gene Explorer",
    layout="wide",
)

st.title("Gene Explorer")


# ============================================================
# Constants
# ============================================================

BASES_PER_ROW = 40

DEFAULT_SITE_WINDOW = 500
MIN_SITE_WINDOW = 200
MAX_SITE_WINDOW = 1500


# ============================================================
# Sequence browser helper functions
# ============================================================

def normalize_gene_symbol(gene: str) -> str:
    """
    Normalize gene symbols used in cpg_gene_map and gene_annotation.
    """
    return str(gene).strip().upper()


def normalize_cpg_c_position(raw_pos: int, sequence: str, seq_start: int) -> int | None:
    """
    Convert a raw CpG coordinate into the actual C coordinate of the CpG
    dinucleotide inside the displayed sequence.
    """
    sequence = sequence.upper()
    idx = raw_pos - seq_start

    def is_cg_at(i: int) -> bool:
        return (
            0 <= i < len(sequence) - 1
            and sequence[i] == "C"
            and sequence[i + 1] == "G"
        )

    if is_cg_at(idx):
        return raw_pos

    for shift in [-1, 1, -2, 2]:
        if is_cg_at(idx + shift):
            return raw_pos + shift

    return None


def build_gcgc_roles(sequence: str, seq_start: int) -> dict[int, list[str]]:
    """
    Detect GCGC motifs and assign CSS classes by genomic position.
    """
    roles: dict[int, list[str]] = {}
    sequence = sequence.upper()

    for match in re.finditer("GCGC", sequence):
        motif_start = seq_start + match.start()

        motif_roles = [
            "gcgc-left",
            "gcgc-middle",
            "gcgc-middle",
            "gcgc-right",
        ]

        for offset, role in enumerate(motif_roles):
            genomic_pos = motif_start + offset
            roles.setdefault(genomic_pos, []).append("gcgc-site")
            roles[genomic_pos].append(role)

    return roles


def render_sequence_matrix(
    sequence: str,
    seq_start: int,
    manifest_cpg_positions: set[int],
    selected_c_pos: int | None,
    gcgc_roles: dict[int, list[str]],
    bases_per_row: int = BASES_PER_ROW,
) -> str:
    """
    Render the DNA sequence as a nucleotide matrix.
    """
    sequence = sequence.upper()

    html_parts = [
        '<div class="sequence-card">',
        '<div class="sequence-title">Sequence — Nucleotide View</div>',
    ]

    for row_start in range(0, len(sequence), bases_per_row):
        row_seq = sequence[row_start:row_start + bases_per_row]
        row_genomic_start = seq_start + row_start

        html_parts.append('<div class="sequence-row">')
        html_parts.append(f'<span class="coord">{row_genomic_start:,}</span>')

        for i, base in enumerate(row_seq):
            genomic_pos = row_genomic_start + i
            safe_base = html.escape(base)

            classes = ["base"]

            if base == "A":
                classes.append("base-a")
            elif base == "C":
                classes.append("base-c")
            elif base == "G":
                classes.append("base-g")
            elif base == "T":
                classes.append("base-t")
            else:
                classes.append("base-n")

            if genomic_pos in gcgc_roles:
                classes.extend(gcgc_roles[genomic_pos])

            if genomic_pos in manifest_cpg_positions:
                classes.append("manifest-cpg")

            if selected_c_pos is not None and genomic_pos == selected_c_pos:
                classes.append("selected-cpg")

            class_str = " ".join(classes)

            html_parts.append(
                f'<span class="{class_str}" '
                f'title="pos {genomic_pos:,} | base {safe_base}">'
                f'{safe_base}</span>'
            )

        html_parts.append("</div>")

    html_parts.append(
        """
        <div class="legend">
            <span class="legend-item"><span class="dot dot-green"></span>CpG del gráfico/manifest</span>
            <span class="legend-item"><span class="dot dot-orange"></span>CpG seleccionado</span>
            <span class="legend-item"><span class="box-black"></span>Motivo GCGC</span>
        </div>
        </div>
        """
    )

    return "".join(html_parts)


def extract_selected_table_row(table_event) -> int | None:
    """
    Extract selected row index from st.dataframe(..., on_select='rerun').
    """
    if table_event is None:
        return None

    selection = getattr(table_event, "selection", None)

    if selection is None and isinstance(table_event, dict):
        selection = table_event.get("selection")

    if selection is None:
        return None

    rows = getattr(selection, "rows", None)

    if rows is None and isinstance(selection, dict):
        rows = selection.get("rows", [])

    if not rows:
        return None

    return int(rows[0])


# ============================================================
# Data loading helpers
# ============================================================

@st.cache_data(ttl=600, show_spinner=False)
def load_gene_options(tumor_type: str, min_sites: int) -> pd.DataFrame:
    """
    Load only genes with CpG sites available for the selected tumor type.

    Important:
    This uses cpg_gene_map instead of cpg_annotation.gene because the manifest
    can contain compound annotations such as 'AC111182.2;SEPTIN9'.
    """
    query = """
    SELECT
        cgm.gene_symbol,
        COUNT(DISTINCT cgm.site_id) AS n_sites
    FROM cpg_gene_map cgm
    JOIN tumor_summary ts
        ON cgm.site_id = ts.site_id
    WHERE
        ts.tumor_type = :tumor_type
        AND cgm.gene_symbol IS NOT NULL
        AND cgm.gene_symbol <> ''
    GROUP BY cgm.gene_symbol
    HAVING COUNT(DISTINCT cgm.site_id) >= :min_sites
    ORDER BY cgm.gene_symbol
    """

    return run_query(
        query,
        params={
            "tumor_type": tumor_type,
            "min_sites": int(min_sites),
        },
    )


@st.cache_data(ttl=600, show_spinner=False)
def load_gene_graph(gene: str, tumor_type: str) -> pd.DataFrame:
    """
    Load CpG-level methylation summaries for one gene and tumor type.

    The gene match is exact against cpg_gene_map.gene_symbol, avoiding false
    matches caused by compound manifest gene fields.
    """
    query = """
    SELECT DISTINCT
        ca.site_id,
        cgm.gene_symbol AS gene,
        ca.chr,
        ca.start_pos,
        ca.end_pos,
        ts.delta_median,
        ts.hi_index,
        ts.tumor_median,
        ts.normal_median,
        ts.pantumor_median AS pan_tumor_median,
        ts.pannormal_median AS pan_normal_median,
        cf.leukocyte_median,
        ec.spearman_r
    FROM tumor_summary ts
    JOIN cpg_annotation ca
        ON ts.site_id = ca.site_id
    JOIN cpg_gene_map cgm
        ON ts.site_id = cgm.site_id
    LEFT JOIN cpg_features cf
        ON ts.site_id = cf.site_id
    LEFT JOIN expression_correlation ec
        ON ts.site_id = ec.site_id
       AND ec.tumor_type = ts.tumor_type
    WHERE
        cgm.gene_symbol = :gene
        AND ts.tumor_type = :tumor_type
    ORDER BY ca.start_pos
    LIMIT 1000
    """

    df = run_query(
        query,
        params={
            "gene": normalize_gene_symbol(gene),
            "tumor_type": tumor_type,
        },
    )

    if df.empty:
        return df

    numeric_cols = [
        "start_pos",
        "end_pos",
        "delta_median",
        "hi_index",
        "tumor_median",
        "normal_median",
        "pan_tumor_median",
        "pan_normal_median",
        "leukocyte_median",
        "spearman_r",
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["start_pos"]).copy()
    df["start_pos"] = df["start_pos"].astype(int)

    return df


@st.cache_data(ttl=600, show_spinner=False)
def load_gene_annotation(gene: str) -> pd.DataFrame:
    """
    Load clean gene-level coordinates, if available.
    """
    query = """
    SELECT
        gene_symbol,
        ensembl_id,
        chr,
        start_pos,
        end_pos,
        strand,
        tss
    FROM gene_annotation
    WHERE gene_symbol = :gene
    LIMIT 1
    """

    return run_query(
        query,
        params={"gene": normalize_gene_symbol(gene)},
    )


# ============================================================
# CSS
# ============================================================

st.markdown(
    """
    <style>
    .sequence-card {
        background-color: white;
        border: 1px solid #d9e2ec;
        border-radius: 8px;
        padding: 16px;
        overflow-x: auto;
        margin-top: 12px;
    }

    .sequence-title {
        font-weight: 700;
        font-size: 16px;
        margin-bottom: 14px;
        color: #102a43;
    }

    .sequence-row {
        white-space: nowrap;
        line-height: 1.9;
    }

    .coord {
        display: inline-block;
        width: 105px;
        color: #9aa5b1;
        font-family: monospace;
        font-size: 13px;
        text-align: right;
        margin-right: 12px;
    }

    .base {
        display: inline-block;
        width: 17px;
        height: 19px;
        line-height: 19px;
        text-align: center;
        font-family: monospace;
        font-size: 14px;
        font-weight: 700;
        margin: 0 1px;
        border-radius: 3px;
        box-sizing: border-box;
        border: 1px solid transparent;
        background-color: transparent;
    }

    .base-a { color: #16a34a; }
    .base-c { color: #ef4444; }
    .base-g { color: #2563eb; }
    .base-t { color: #9333ea; }
    .base-n { color: #6b7280; }

    .gcgc-site {
        color: #ffffff !important;
        font-weight: 900;
        border-top: 2px solid #000000;
        border-bottom: 2px solid #000000;
        background-color: #000000;
    }

    .gcgc-left {
        border-left: 2px solid #000000;
        border-top-left-radius: 4px;
        border-bottom-left-radius: 4px;
    }

    .gcgc-right {
        border-right: 2px solid #000000;
        border-top-right-radius: 4px;
        border-bottom-right-radius: 4px;
    }

    .manifest-cpg {
        color: #000000 !important;
        border: 2px solid #22c55e !important;
        background-color: #22ff5e80;
    }

    .selected-cpg {
        color: #000000 !important;
        border: 2px solid #ffa500 !important;
        background-color: #ffa500 !important;
        font-weight: 900;
    }

    .legend {
        margin-top: 14px;
        font-size: 13px;
        color: #52616b;
    }

    .legend-item {
        display: inline-block;
        margin-right: 22px;
    }

    .dot {
        display: inline-block;
        width: 12px;
        height: 12px;
        border-radius: 50%;
        margin-right: 5px;
        vertical-align: middle;
    }

    .dot-green { background-color: #22c55e; }
    .dot-orange { background-color: #ffa500; }

    .box-black {
        display: inline-block;
        width: 14px;
        height: 12px;
        border: 2px solid black;
        margin-right: 5px;
        vertical-align: middle;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# Sidebar
# ============================================================

st.sidebar.header("Filters")

tumor_map = {
    "CRC (COAD)": "COAD",
    "HCC (LIHC)": "LIHC",
    "LUAD": "LUAD",
    "LUSC": "LUSC"
}

tumor_label = st.sidebar.selectbox(
    "Tumor Type",
    list(tumor_map.keys()),
)

tumor_type = tumor_map[tumor_label]

min_gene_sites_browser = st.sidebar.slider(
    "Minimum CpGs per gene in browser",
    min_value=1,
    max_value=50,
    value=1,
    step=1,
)

gene_options = load_gene_options(
    tumor_type=tumor_type,
    min_sites=min_gene_sites_browser,
)

if gene_options.empty:
    st.warning(
        "No genes with available CpGs were found for that tumor type"
        "and that minimum number of sites."
    )
    st.stop()

gene_options["gene_symbol"] = gene_options["gene_symbol"].astype(str).str.upper()
gene_options = gene_options.drop_duplicates(subset=["gene_symbol"]).sort_values("gene_symbol")

n_sites_by_gene = dict(
    zip(
        gene_options["gene_symbol"],
        gene_options["n_sites"].astype(int),
    )
)

gene_symbols = gene_options["gene_symbol"].tolist()

default_gene = "SEPTIN9"
default_index = gene_symbols.index(default_gene) if default_gene in gene_symbols else 0

gene = st.sidebar.selectbox(
    "Gene Symbol",
    gene_symbols,
    index=default_index,
    format_func=lambda g: f"{g} ({n_sites_by_gene.get(g, 0)} CpGs)",
)

gene = normalize_gene_symbol(gene)


# ============================================================
# Load graph data
# ============================================================

graph = load_gene_graph(
    gene=gene,
    tumor_type=tumor_type,
)

if graph.empty:
    st.warning(
        f"No CpGs were found for {gene} in {tumor_type}."
        "Check that cpg_gene_map is loaded and that gene_symbol is normalized."
    )
    st.stop()

for optional_col in [
    "pan_tumor_median",
    "pan_normal_median",
    "leukocyte_median",
    "spearman_r",
    "hi_index",
]:
    if optional_col not in graph.columns:
        graph[optional_col] = np.nan


# ============================================================
# Gene annotation metadata
# ============================================================

gene_annotation = load_gene_annotation(gene)

if gene_annotation.empty:
    st.info(
        "{gene} has CpGs in cpg_gene_map, but I didn't find coordinates"
        "in gene_annotation. The browser will use the annotated chromosome for each CpG."
    )
else:
    gene_row = gene_annotation.iloc[0]

    with st.expander("Gene annotation"):
        st.write(
            {
                "gene_symbol": gene_row.get("gene_symbol"),
                "ensembl_id": gene_row.get("ensembl_id"),
                "chr": gene_row.get("chr"),
                "start_pos": gene_row.get("start_pos"),
                "end_pos": gene_row.get("end_pos"),
                "strand": gene_row.get("strand"),
                "tss": gene_row.get("tss"),
            }
        )


# ============================================================
# Plot
# ============================================================

site = graph["site_id"]
positions = graph["start_pos"]

median_tumor = graph["tumor_median"]
median_normal = graph["normal_median"]
pan_tumor = graph["pan_tumor_median"]
pan_normal = graph["pan_normal_median"]
delta_median = graph["delta_median"]
hi_index = graph["hi_index"]
leukocytes = graph["leukocyte_median"]

hover_text = (
    graph["site_id"].astype(str)
    + " | "
    + graph["start_pos"].astype(str)
)

show_hi = st.checkbox("Show HI. Checking this box will modify the right axis.", value=False)

fig2 = go.Figure()

fig2.add_trace(
    go.Scatter(
        x=positions,
        y=median_tumor,
        mode="lines",
        name="Median Type T",
        customdata=site,
        line=dict(
            width=2,
            color="rgba(255,165,0,1)",
            shape="spline",
            smoothing=1.3,
        ),
        fill="tozeroy",
        fillcolor="rgba(255,165,0,1)",
    )
)

curves = [
    ("Median Type NT", median_normal, "rgba(255,255,0,1)", "line", "y"),
    ("Median PanCan T", pan_tumor, "rgba(0,255,0,1)", "line", "y"),
    ("Median PanCan NT", pan_normal, "rgba(0,128,0,1)", "line", "y"),
    ("Median Leukocytes", leukocytes, "rgba(130,168,240,1)", "line", "y"),
  # ("Expression", expression, "rgba(127,127,127,1)", "expression", "line"),
    ("Delta", delta_median, "rgba(236,0,139,1)", "marker", "y2"),
]

for name, values, color, trace_type, yaxis in curves:
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
                    smoothing=1.3,
                ),
            )
        )
    else:
        fig2.add_trace(
            go.Scatter(
                x=positions,
                y=values,
                mode="markers",
                yaxis=yaxis,
                name=name,
                marker=dict(
                    color=color,
                    size=8,
                ),
            )
        )

if show_hi:
    fig2.add_trace(
        go.Scatter(
            x=positions,
            y=hi_index,
            mode="markers",
            marker=dict(
                color="rgba(255,0,0,1)",
                size=8,
            ),
            name="HI",
            yaxis="y2",
        )
    )
    y2_max = hi_index.max()
else:
    y2_max = 1

finite_y2 = pd.concat(
    [
        pd.Series(delta_median).replace([np.inf, -np.inf], np.nan),
        pd.Series(hi_index).replace([np.inf, -np.inf], np.nan),
    ],
    ignore_index=True,
).dropna()

fig2.update_traces(
    hovertext=hover_text,
    hovertemplate="%{hovertext}<br>%{fullData.name}: %{y:.3f}<extra></extra>",
)

fig2.update_layout(
    template="plotly_white",
    height=500,
    hovermode="x unified",
    margin=dict(l=40, r=40, t=50, b=50),
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
            color="black",
        ),
    ),
    title=dict(
        text=f"Aberrant tumor hypermethylation area for <b>{gene}</b> - {tumor_type}",
        x=0,
        xanchor="left",
        font=dict(
            family="Arial Narrow",
            size=20,
            color="black",
        ),
    ),
    xaxis=dict(
        title=dict(
            text="Genomic Position",
            font=dict(
                family="Arial Narrow",
                color="black",
                size=18,
            ),
        ),
        tickfont=dict(
            family="Arial Narrow",
            color="black",
            size=18,
        ),
        tickformat=",d",
        rangeslider=dict(visible=True),
        showgrid=True,
    ),
    yaxis=dict(
        range=[-0.2, 1],
        title=dict(
            text="Methylation",
            font=dict(
                family="Arial Narrow",
                color="black",
                size=18,
            ),
        ),
        tickfont=dict(
            family="Arial Narrow",
            color="black",
            size=18,
        ),
        showgrid=True,
    ),
    yaxis2=dict(
        title=dict(
            text="Delta/HI",
            font=dict(
                family="Arial Narrow",
                color="black",
                size=18,
            ),
        ),
        overlaying="y",
        side="right",
        range=[-0.2, y2_max],
        tickfont=dict(
            family="Arial Narrow",
            color="black",
            size=18,
        ),
    ),
)

config = {
    "toImageButtonOptions": {
        "format": "svg",
        "filename": f"{gene}_{tumor_type}",
        "height": 600,
        "width": 1400,
        "scale": 1,
    }
}

st.plotly_chart(
    fig2,
    use_container_width=True,
    config=config,
)


# ============================================================
# CpG table controlling the sequence browser
# ============================================================

st.subheader("CpG sites in the chart")
st.caption(
    "Select a row in the table to open a sequence window centered on that site."
    "The browser window will remain empty until you select a site."
)

site_window = st.slider(
    "Window around the selected site (bp)",
    min_value=MIN_SITE_WINDOW,
    max_value=MAX_SITE_WINDOW,
    value=DEFAULT_SITE_WINDOW,
    step=500,
)

graph_table = graph.reset_index(drop=True).copy()

table_event = st.dataframe(
    graph_table,
    use_container_width=True,
    on_select="rerun",
    selection_mode="single-row",
    hide_index=True,
)

selected_table_idx = extract_selected_table_row(table_event)


# ============================================================
# Sequence browser
# ============================================================

st.subheader("Browser")

if selected_table_idx is None:
    st.info("No site selected. Choose a row from the table to see the sequence.")
else:
    selected_row = graph_table.iloc[selected_table_idx]

    selected_pos = int(selected_row["start_pos"])
    selected_site_id = str(selected_row["site_id"])
    selected_gene = normalize_gene_symbol(selected_row["gene"])

    raw_start = max(1, selected_pos - site_window // 2)
    raw_end = selected_pos + site_window // 2

    # For a site-centered sequence browser, the most robust chromosome source
    # is the CpG annotation row itself. This avoids failures when gene_annotation
    # is incomplete or when gene symbols were compound in the original manifest.
    chrom = str(selected_row["chr"]).replace("chr", "").replace("CHR", "").strip()

    if chrom in {"", "nan", "None"}:
        st.error("No gene for CpG was found.")
        st.stop()

    seq_start = max(1, int(raw_start))
    seq_end = int(raw_end)

    try:
        sequence = get_sequence(chrom, seq_start, seq_end).upper()
    except Exception as exc:
        st.error(f"The sequence could not be obtained from Ensembl: {exc}")
        sequence = ""

    if sequence:
        window_cpgs = graph_table[
            (graph_table["start_pos"] >= seq_start)
            & (graph_table["start_pos"] <= seq_end)
        ].copy()

        if not window_cpgs.empty:
            window_cpgs["highlight_c_pos"] = window_cpgs["start_pos"].apply(
                lambda p: normalize_cpg_c_position(
                    raw_pos=int(p),
                    sequence=sequence,
                    seq_start=seq_start,
                )
            )

            window_cpgs_display = window_cpgs.dropna(subset=["highlight_c_pos"]).copy()
            window_cpgs_display["highlight_c_pos"] = (
                window_cpgs_display["highlight_c_pos"].astype(int)
            )

            manifest_cpg_positions = set(
                window_cpgs_display["highlight_c_pos"].astype(int).tolist()
            )
        else:
            window_cpgs_display = window_cpgs
            manifest_cpg_positions = set()

        selected_c_pos = normalize_cpg_c_position(
            raw_pos=selected_pos,
            sequence=sequence,
            seq_start=seq_start,
        )

        gcgc_roles = build_gcgc_roles(sequence, seq_start)
        col_a, col_c, col_d = st.columns(3)

        with col_a:
            st.metric("Site", selected_site_id)

        with col_c:
            st.metric("Longitud", f"{len(sequence):,} bp")

        with col_d:
            st.metric("Motivos GCGC", len(re.findall("GCGC", sequence)))

        col_b1, col_b2, _ = st.columns([1, 2, 1])

        with col_b1:
            st.metric("Chromosome", f"chr{chrom}")

        with col_b2:
            st.metric("Region", f"{seq_start}-{seq_end}")

        if selected_c_pos is not None:
            st.info(
                f"Selected site: {selected_site_id} | "
                f"gen: {selected_gene} | "
                f"raw position: {selected_pos:,} | "
                f"corrected C: {selected_c_pos:,} | "
                f"ventana: {site_window:,} bp"
            )
        else:
            st.warning(
                f"Selected site: {selected_site_id} | gene: {selected_gene} | "
                f"Crude position: {selected_pos:,}."
                f"Could not be normalized to a CpG within the sequence."
            )

        sequence_html = render_sequence_matrix(
            sequence=sequence,
            seq_start=seq_start,
            manifest_cpg_positions=manifest_cpg_positions,
            selected_c_pos=selected_c_pos,
            gcgc_roles=gcgc_roles,
            bases_per_row=BASES_PER_ROW,
        )

        st.markdown(sequence_html, unsafe_allow_html=True)

