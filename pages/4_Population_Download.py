#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import copy
import re
from typing import Iterable

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from db.queries import run_query


# ============================================================
# Page configuration
# ============================================================

st.set_page_config(
    page_title="Population Bubble Comparison — Download",
    page_icon="",
    layout="wide",
)

st.title("Population Bubble Comparison — Download")
st.caption(
    "Visual comparison of region-level methylation evidence in the full "
    "cohort and after excluding samples annotated as Asian."
)


# ============================================================
# Constants
# ============================================================

TUMOR_LABELS = {
    "COAD": "CRC (COAD)",
    "LIHC": "HCC (LIHC)",
    "LUAD": "LUAD",
    "LUSC": "LUSC",
}

MODE_LABELS = {
    "full": "Full cohort",
    "asian_excluded": "Asian-excluded cohort",
    "asian_only": "Asian-only cohort",
}

MODE_COLORS = {
    "full": "#2F80ED",            # stronger blue
    "asian_excluded": "#5F6368",  # darker gray
    "asian_only": "#173F6B",      # dark blue
}

# Same tumor / non-tumor palette used in Gene Explorer.
SAMPLE_TUMOR_COLOR = "rgba(18,144,152,1)"
SAMPLE_NON_TUMOR_COLOR = "rgba(240,145,62,1)"
SAMPLE_TUMOR_LABEL = "Tumor (T)"
SAMPLE_NON_TUMOR_LABEL = "Non-Tumor (NT)"

REGION_LINK_COLOR = "rgba(75,85,99,1)"
ASIAN_REGION_LINK_COLOR = "rgba(23,63,107,0.28)"

PLOT_HEIGHT = 700
PLOT_WIDTH = 1400


# ============================================================
# Generic helpers
# ============================================================

@st.cache_data(ttl=600, show_spinner=False)
def get_table_columns(table_name: str) -> set[str]:
    query = """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = 'public'
          AND table_name = :table_name
        ORDER BY ordinal_position
    """
    try:
        df = run_query(query, params={"table_name": table_name})
    except Exception:
        return set()

    if df.empty or "column_name" not in df.columns:
        return set()
    return set(df["column_name"].astype(str).tolist())


def sql_col(
    alias: str,
    columns: set[str],
    candidates: Iterable[str],
    fallback: str = "NULL",
) -> str:
    """Return the first available compatible SQL column."""
    for column in candidates:
        if column in columns:
            return f'{alias}."{column}"'
    return fallback


def clean_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def split_gene_text(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []

    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "na"}:
        return []

    genes = [
        gene.strip().upper()
        for gene in text.replace(",", ";").split(";")
    ]
    return list(
        dict.fromkeys(
            gene
            for gene in genes
            if gene and gene.lower() not in {"nan", "none", "null", "na"}
        )
    )


def join_gene_values(values: pd.Series) -> str:
    genes: list[str] = []
    for value in values:
        genes.extend(split_gene_text(value))
    return ";".join(sorted(set(genes)))


def choose_main_gene(value: object) -> str:
    genes = split_gene_text(value)
    if not genes:
        return ""

    technical_prefixes = ("AC", "AL", "AP", "LOC", "LINC", "MIR")

    def priority(gene: str) -> tuple[int, int, str]:
        technical = gene.startswith(technical_prefixes) or "." in gene
        return (1 if technical else 0, len(gene), gene)

    return sorted(genes, key=priority)[0]


def physical_region_id(chrom: object, start: object, end: object) -> str:
    chrom_clean = str(chrom).replace("chr", "").replace("CHR", "").strip()
    return f"chr{chrom_clean}:{int(start)}-{int(end)}"


def get_lung_cross_tumor_type(tumor_type: str) -> str | None:
    if tumor_type == "LUAD":
        return "LUSC"
    if tumor_type == "LUSC":
        return "LUAD"
    return None


def plot_svg_config(
    filename: str,
    height: int = PLOT_HEIGHT,
    width: int = PLOT_WIDTH,
) -> dict:
    """Configure the Plotly camera button to export the figure as SVG."""
    return {
        "toImageButtonOptions": {
            "format": "svg",
            "filename": filename,
            "height": height,
            "width": width,
            "scale": 1,
        }
    }


def apply_plot_style(
    fig: go.Figure,
    *,
    height: int = PLOT_HEIGHT,
    width: int = PLOT_WIDTH,
) -> go.Figure:
    font_family = "Arial Narrow"

    fig.update_layout(
        template="plotly_white",
        autosize=False,
        width=width,
        height=height,
        paper_bgcolor="white",
        plot_bgcolor="white",
        hovermode="closest",
        margin=dict(l=60, r=60, t=80, b=70),
        font=dict(family=font_family, color="black"),
        title=dict(
            x=0,
            xanchor="left",
            font=dict(family=font_family, size=24, color="black"),
        ),
        legend=dict(
            orientation="v",
            y=1,
            x=1.02,
            bgcolor="rgba(255,255,255,0)",
            borderwidth=0,
            font=dict(family=font_family, size=18, color="black"),
        ),
    )

    fig.update_xaxes(
        title_font=dict(family=font_family, size=24, color="black"),
        tickfont=dict(family=font_family, size=18, color="black"),
        showgrid=True,
        zeroline=True,
    )
    fig.update_yaxes(
        title_font=dict(family=font_family, size=24, color="black"),
        tickfont=dict(family=font_family, size=18, color="black"),
        showgrid=True,
        zeroline=True,
    )
    return fig


def make_large_font_export_figure(
    fig: go.Figure,
    title_size: int = 30,
    axis_title_size: int = 28,
    axis_tick_size: int = 26,
    legend_size: int = 24,
    annotation_size: int = 20,
    trace_text_size: int = 22,
) -> go.Figure:
    """Create a copy of a Plotly figure with larger fonts for SVG export."""
    fig_export = copy.deepcopy(fig)

    fig_export.update_layout(
        title=dict(font=dict(size=title_size)),
        legend=dict(font=dict(size=legend_size)),
    )

    for axis_name in ("xaxis", "yaxis", "xaxis2", "yaxis2"):
        axis = getattr(fig_export.layout, axis_name, None)
        if axis is not None:
            axis.update(
                title=dict(font=dict(size=axis_title_size)),
                tickfont=dict(size=axis_tick_size),
            )

    if getattr(fig_export.layout, "annotations", None):
        annotations = []
        for annotation in fig_export.layout.annotations:
            annotation_json = annotation.to_plotly_json()
            font = annotation_json.get("font", {}) or {}
            font["size"] = annotation_size
            annotation_json["font"] = font
            annotations.append(annotation_json)
        fig_export.update_layout(annotations=annotations)

    # Enlarge labels drawn directly by bar traces.
    for trace in fig_export.data:
        if isinstance(trace, go.Bar):
            trace.update(textfont=dict(size=trace_text_size))

    return fig_export


def show_large_font_svg_export(
    fig: go.Figure,
    filename: str,
    expander_label: str,
    *,
    height: int,
    width: int,
) -> None:
    """
    Display a second browser-rendered Plotly figure for SVG export.

    This follows the same approach as Gene_Explorer_Download and
    Region_Explorer_Download: the user opens the preview and downloads the SVG
    with Plotly's camera icon. No TSV or Kaleido export is used.
    """
    fig_export = make_large_font_export_figure(fig)
    export_key = re.sub(r"[^A-Za-z0-9_-]+", "_", filename)

    with st.expander(expander_label, expanded=False):
        st.caption(
            "Open this preview and use the camera icon in the Plotly toolbar. "
            "The SVG keeps the same colors, markers, lines and dimensions; "
            "only the typography is enlarged for figure export."
        )
        st.plotly_chart(
            fig_export,
            use_container_width=True,
            key=f"svg_export_{export_key}",
            config=plot_svg_config(
                filename=filename,
                height=height,
                width=width,
            ),
        )


def common_marker_sizeref(
    region_frames: dict[str, pd.DataFrame],
    max_px: float = 44.0,
) -> float:
    """Use one area scale across all cohorts so bubble sizes are comparable."""
    maxima: list[float] = []
    for frame in region_frames.values():
        if frame.empty or "n_qualifying_sites" not in frame.columns:
            continue
        values = pd.to_numeric(
            frame["n_qualifying_sites"], errors="coerce"
        ).dropna()
        if not values.empty:
            maxima.append(float(values.max()))

    global_max = max(maxima, default=1.0)
    return max(2.0 * global_max / (max_px**2), 1e-6)


def add_physical_region_id(cpgs: pd.DataFrame) -> pd.DataFrame:
    """Add the physical-region key used to pair matching cohort bubbles."""
    if cpgs.empty:
        return cpgs.copy()
    out = cpgs.copy()
    out["physical_region_id"] = [
        physical_region_id(row.chr, row.browser_start, row.browser_end)
        for row in out.itertuples(index=False)
    ]
    return out



# ============================================================
# Availability and methylation sample counts
# ============================================================

@st.cache_data(ttl=600, show_spinner=False)
def load_available_population_data() -> pd.DataFrame:
    query = """
        SELECT DISTINCT
            tumor_type,
            population_mode,
            pan_reference_mode
        FROM tumor_summary_population
        ORDER BY tumor_type, population_mode, pan_reference_mode
    """
    return run_query(query)


@st.cache_data(ttl=600, show_spinner=False)
def load_methylation_sample_counts(
    tumor_type: str,
    pan_reference_mode: str,
) -> pd.DataFrame:
    """Return sample numbers actually represented in the methylation matrix."""
    columns = get_table_columns("tumor_summary_population")
    required = {
        "tumor_type",
        "population_mode",
        "pan_reference_mode",
        "tumor_total_samples",
        "normal_total_samples",
    }
    if not required.issubset(columns):
        return pd.DataFrame()

    query = """
        SELECT
            population_mode,
            MAX(tumor_total_samples)::INTEGER AS tumor_n,
            MAX(normal_total_samples)::INTEGER AS normal_n
        FROM tumor_summary_population
        WHERE tumor_type = :tumor_type
          AND pan_reference_mode = :pan_reference_mode
          AND population_mode IN (
              'full',
              'asian_excluded',
              'asian_only'
          )
        GROUP BY population_mode
        ORDER BY population_mode
    """
    counts = run_query(
        query,
        params={
            "tumor_type": tumor_type,
            "pan_reference_mode": pan_reference_mode,
        },
    )
    if counts.empty:
        return counts

    for column in ["tumor_n", "normal_n"]:
        counts[column] = pd.to_numeric(
            counts[column], errors="coerce"
        ).fillna(0).astype(int)

    by_mode = counts.set_index("population_mode")
    rows: list[dict[str, object]] = []

    for mode in ["full", "asian_excluded"]:
        if mode in by_mode.index:
            rows.append(
                {
                    "mode": mode,
                    "cohort": MODE_LABELS[mode],
                    "tumor_n": int(by_mode.loc[mode, "tumor_n"]),
                    "normal_n": int(by_mode.loc[mode, "normal_n"]),
                    "count_source": "Directly stored by pipeline 04b",
                }
            )

    if "asian_only" in by_mode.index:
        asian_tumor_n = int(by_mode.loc["asian_only", "tumor_n"])
        asian_normal_n = int(by_mode.loc["asian_only", "normal_n"])
        asian_source = "Directly stored by pipeline 04b"
    elif {"full", "asian_excluded"}.issubset(by_mode.index):
        asian_tumor_n = max(
            0,
            int(by_mode.loc["full", "tumor_n"])
            - int(by_mode.loc["asian_excluded", "tumor_n"]),
        )
        asian_normal_n = max(
            0,
            int(by_mode.loc["full", "normal_n"])
            - int(by_mode.loc["asian_excluded", "normal_n"]),
        )
        asian_source = "Calculated as Full minus Asian-excluded"
    else:
        asian_tumor_n = 0
        asian_normal_n = 0
        asian_source = "Unavailable"

    rows.append(
        {
            "mode": "asian_only",
            "cohort": MODE_LABELS["asian_only"],
            "tumor_n": asian_tumor_n,
            "normal_n": asian_normal_n,
            "count_source": asian_source,
        }
    )

    return pd.DataFrame(rows)


# ============================================================
# CpG query with the filters used by Population Sensitivity
# ============================================================

@st.cache_data(ttl=600, show_spinner=True)
def load_population_mode_cpgs(
    tumor_type: str,
    population_mode: str,
    pan_reference_mode: str,
    min_delta: float,
    max_normal_median: float,
    max_pan_normal_median: float,
    max_pan_tumor_median: float,
    max_pb: float,
    min_hi: float,
    cross_tumor_type: str | None,
    max_cross_tumor_median: float | None,
) -> pd.DataFrame:
    """Load qualifying CpGs for one population mode."""
    population_columns = get_table_columns("tumor_summary_population")
    region_columns = get_table_columns("biomarker_region")
    cpg_feature_columns = get_table_columns("cpg_features")
    expression_columns = get_table_columns("expression_correlation")
    tumor_summary_columns = get_table_columns("tumor_summary")

    if not population_columns:
        return pd.DataFrame()

    genes_all_expr = sql_col(
        "r",
        region_columns,
        ["gene_symbols_all"],
        'r."gene_symbol"',
    )
    pan_tumor_expr = sql_col(
        "s",
        population_columns,
        ["pan_tumor_median", "pantumor_median"],
        "NULL",
    )
    pan_normal_expr = sql_col(
        "s",
        population_columns,
        ["pan_normal_median", "pannormal_median"],
        "NULL",
    )
    pb_expr = sql_col(
        "cf",
        cpg_feature_columns,
        ["pb_median", "leukocyte_median"],
        "NULL",
    )

    # Expression evidence is optional. When the table/columns do not exist,
    # return NULL and keep the page usable with the expression filter disabled.
    expression_join = ""
    spearman_expr = "NULL"
    if {"site_id", "tumor_type", "spearman_r"}.issubset(expression_columns):
        expression_join = """
            LEFT JOIN (
                SELECT
                    site_id,
                    tumor_type,
                    MIN(spearman_r) AS spearman_r
                FROM expression_correlation
                GROUP BY site_id, tumor_type
            ) eb
              ON eb.site_id = s.site_id
             AND eb.tumor_type = s.tumor_type
        """
        spearman_expr = "eb.spearman_r"

    cross_join = ""
    cross_select = "NULL AS cross_tumor_median"
    cross_condition = "TRUE"
    use_cross_filter = (
        cross_tumor_type in {"LUAD", "LUSC"}
        and max_cross_tumor_median is not None
        and {"site_id", "tumor_type", "tumor_median"}.issubset(
            tumor_summary_columns
        )
    )
    if use_cross_filter:
        cross_join = """
            LEFT JOIN tumor_summary ts_cross
              ON ts_cross.site_id = s.site_id
             AND ts_cross.tumor_type = :cross_tumor_type
        """
        cross_select = "ts_cross.tumor_median AS cross_tumor_median"
        cross_condition = (
            "COALESCE(ts_cross.tumor_median, 1) "
            "<= :max_cross_tumor_median"
        )

    query = f"""
        SELECT
            brc.region_id,
            r.gene_symbol AS gene,
            {genes_all_expr} AS genes_all,
            r.chr,
            r.browser_start,
            r.browser_end,
            r.browser_length,
            brc.site_id,
            brc.start_pos,
            s.delta_median,
            s.hi_index,
            s.tumor_median,
            s.normal_median,
            {pan_tumor_expr} AS pan_tumor_median,
            {pan_normal_expr} AS pan_normal_median,
            {pb_expr} AS pb_median,
            {spearman_expr} AS spearman_r,
            {cross_select}
        FROM tumor_summary_population s
        JOIN biomarker_region_cpg brc
          ON brc.site_id = s.site_id
        JOIN biomarker_region r
          ON r.region_id = brc.region_id
        LEFT JOIN cpg_features cf
          ON cf.site_id = s.site_id
        {expression_join}
        {cross_join}
        WHERE s.tumor_type = :tumor_type
          AND s.population_mode = :population_mode
          AND s.pan_reference_mode = :pan_reference_mode
          AND COALESCE(s.delta_median, -1) >= :min_delta
          AND COALESCE(s.normal_median, 1) <= :max_normal_median
          AND COALESCE({pan_normal_expr}, 1) <= :max_pan_normal_median
          AND COALESCE({pan_tumor_expr}, 1) <= :max_pan_tumor_median
          AND COALESCE({pb_expr}, 1) <= :max_pb
          AND COALESCE(s.hi_index, 0) >= :min_hi
          AND {cross_condition}
        ORDER BY r.chr, r.browser_start, brc.start_pos
    """

    params: dict[str, object] = {
        "tumor_type": tumor_type,
        "population_mode": population_mode,
        "pan_reference_mode": pan_reference_mode,
        "min_delta": float(min_delta),
        "max_normal_median": float(max_normal_median),
        "max_pan_normal_median": float(max_pan_normal_median),
        "max_pan_tumor_median": float(max_pan_tumor_median),
        "max_pb": float(max_pb),
        "min_hi": float(min_hi),
    }
    if use_cross_filter:
        params["cross_tumor_type"] = str(cross_tumor_type)
        params["max_cross_tumor_median"] = float(max_cross_tumor_median)

    df = run_query(query, params=params)
    if df.empty:
        return df

    return clean_numeric(
        df,
        [
            "browser_start",
            "browser_end",
            "browser_length",
            "start_pos",
            "delta_median",
            "hi_index",
            "tumor_median",
            "normal_median",
            "pan_tumor_median",
            "pan_normal_median",
            "pb_median",
            "spearman_r",
            "cross_tumor_median",
        ],
    )


def aggregate_regions(cpgs: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Aggregate qualifying CpGs into the app's physical regions."""
    if cpgs.empty:
        return pd.DataFrame()

    out = cpgs.copy()
    out["physical_region_id"] = [
        physical_region_id(row.chr, row.browser_start, row.browser_end)
        for row in out.itertuples(index=False)
    ]

    # Avoid counting the same CpG twice when a physical interval has more than
    # one gene-oriented region_id.
    out = out.drop_duplicates(subset=["physical_region_id", "site_id"])

    regions = (
        out.groupby(
            [
                "physical_region_id",
                "chr",
                "browser_start",
                "browser_end",
                "browser_length",
            ],
            dropna=False,
        )
        .agg(
            genes_all=("genes_all", join_gene_values),
            genes_from_rows=("gene", join_gene_values),
            n_qualifying_sites=("site_id", "nunique"),
            mean_delta=("delta_median", "mean"),
            mean_hi=("hi_index", "mean"),
            mean_tumor_median=("tumor_median", "mean"),
            mean_normal_median=("normal_median", "mean"),
            mean_pan_tumor_median=("pan_tumor_median", "mean"),
            mean_pan_normal_median=("pan_normal_median", "mean"),
            mean_pb_median=("pb_median", "mean"),
            mean_cross_tumor_median=("cross_tumor_median", "mean"),
            mean_spearman_r=("spearman_r", "mean"),
        )
        .reset_index()
    )

    regions["genes_all"] = np.where(
        regions["genes_all"].astype(str).str.strip() != "",
        regions["genes_all"],
        regions["genes_from_rows"],
    )
    regions["gene_main"] = regions["genes_all"].apply(choose_main_gene)
    regions["population_mode"] = mode
    regions["cohort_label"] = MODE_LABELS[mode]

    return regions.replace([np.inf, -np.inf], np.nan)


# ============================================================
# Plot builders
# ============================================================

def make_sample_count_plot(counts: pd.DataFrame) -> go.Figure | None:
    if counts.empty:
        return None

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=counts["cohort"],
            y=counts["tumor_n"],
            name=SAMPLE_TUMOR_LABEL,
            marker=dict(color=SAMPLE_TUMOR_COLOR),
            text=counts["tumor_n"],
            textposition="outside",
            hovertemplate=(
                "%{x}<br>" + SAMPLE_TUMOR_LABEL + ": %{y}<extra></extra>"
            ),
        )
    )
    fig.add_trace(
        go.Bar(
            x=counts["cohort"],
            y=counts["normal_n"],
            name=SAMPLE_NON_TUMOR_LABEL,
            marker=dict(color=SAMPLE_NON_TUMOR_COLOR),
            text=counts["normal_n"],
            textposition="outside",
            hovertemplate=(
                "%{x}<br>" + SAMPLE_NON_TUMOR_LABEL + ": %{y}<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Methylation samples in each population definition",
        xaxis_title="Population definition",
        yaxis_title="Methylation samples (N)",
        barmode="group",
        uniformtext_minsize=12,
        uniformtext_mode="hide",
    )
    return apply_plot_style(fig, height=520)


def add_region_connection_trace(
    fig: go.Figure,
    source: pd.DataFrame,
    target: pd.DataFrame,
    *,
    target_mode: str,
    line_color: str,
    line_dash: str = "solid",
) -> None:
    """Connect the same physical region across two displayed cohorts."""
    if source.empty or target.empty:
        return

    source_points = source[
        ["physical_region_id", "mean_hi", "mean_delta"]
    ].rename(
        columns={
            "mean_hi": "source_hi",
            "mean_delta": "source_delta",
        }
    )
    target_points = target[
        ["physical_region_id", "mean_hi", "mean_delta"]
    ].rename(
        columns={
            "mean_hi": "target_hi",
            "mean_delta": "target_delta",
        }
    )
    shared = source_points.merge(
        target_points, on="physical_region_id", how="inner"
    ).dropna(
        subset=["source_hi", "source_delta", "target_hi", "target_delta"]
    )
    if shared.empty:
        return

    line_x: list[float | None] = []
    line_y: list[float | None] = []
    hover: list[str | None] = []
    for row in shared.itertuples(index=False):
        line_x.extend([row.source_hi, row.target_hi, None])
        line_y.extend([row.source_delta, row.target_delta, None])
        text = (
            f"{row.physical_region_id}<br>"
            f"Full: HI={row.source_hi:.3f}, Δβ={row.source_delta:.3f}<br>"
            f"{MODE_LABELS[target_mode]}: HI={row.target_hi:.3f}, "
            f"Δβ={row.target_delta:.3f}"
        )
        hover.extend([text, text, None])

    fig.add_trace(
        go.Scatter(
            x=line_x,
            y=line_y,
            mode="lines",
            name=f"Same region",
            line=dict(color=line_color, width=1.2, dash=line_dash),
            hovertext=hover,
            hovertemplate="%{hovertext}<extra></extra>",
            connectgaps=False,
            showlegend=True,
        )
    )


def make_bubble_plot(
    region_frames: dict[str, pd.DataFrame],
    *,
    connect_matching_regions: bool,
) -> go.Figure:
    fig = go.Figure()
    sizeref = common_marker_sizeref(region_frames)

    if connect_matching_regions:
        full = region_frames.get("full", pd.DataFrame())
        excluded = region_frames.get("asian_excluded", pd.DataFrame())
        asian = region_frames.get("asian_only", pd.DataFrame())
        add_region_connection_trace(
            fig,
            full,
            excluded,
            target_mode="asian_excluded",
            line_color=REGION_LINK_COLOR,
        )
        if not asian.empty:
            add_region_connection_trace(
                fig,
                full,
                asian,
                target_mode="asian_only",
                line_color=ASIAN_REGION_LINK_COLOR,
                line_dash="dot",
            )

    for mode in ["full", "asian_excluded", "asian_only"]:
        regions = region_frames.get(mode, pd.DataFrame())
        if regions.empty:
            continue

        sizes = (
            pd.to_numeric(
                regions["n_qualifying_sites"], errors="coerce"
            )
            .fillna(1)
            .clip(lower=1)
        )

        customdata = np.column_stack(
            [
                regions["gene_main"].fillna(""),
                regions["physical_region_id"].astype(str),
                regions["n_qualifying_sites"],
                regions["mean_tumor_median"],
                regions["mean_normal_median"],
                regions["mean_pan_tumor_median"],
                regions["mean_pan_normal_median"],
                regions["mean_pb_median"],
                regions["mean_cross_tumor_median"],
                regions["mean_spearman_r"],
            ]
        )

        fig.add_trace(
            go.Scatter(
                x=regions["mean_hi"],
                y=regions["mean_delta"],
                mode="markers",
                name=MODE_LABELS[mode],
                marker=dict(
                    symbol="circle",
                    size=sizes,
                    sizemode="area",
                    sizeref=sizeref,
                    sizemin=6,
                    color=MODE_COLORS[mode],
                    opacity=0.78 if mode != "full" else 0.86,
                    line=dict(width=1, color="white"),
                ),
                customdata=customdata,
                hovertemplate=(
                    "%{customdata[0]} | %{customdata[1]}<br>"
                    f"Cohort: {MODE_LABELS[mode]}<br>"
                    "Mean HI: %{x:.3f}<br>"
                    "Mean Δβ: %{y:.3f}<br>"
                    "Qualifying CpGs: %{customdata[2]:.0f}<br>"
                    "Mean tumor β: %{customdata[3]:.3f}<br>"
                    "Mean non-tumor β: %{customdata[4]:.3f}<br>"
                    "Mean PanCan tumor β: %{customdata[5]:.3f}<br>"
                    "Mean PanCan non-tumor β: %{customdata[6]:.3f}<br>"
                    "Mean PB β: %{customdata[7]:.3f}<br>"
                    "Mean cross-tumor β: %{customdata[8]:.3f}<br>"
                    "Mean Spearman r: %{customdata[9]:.3f}"
                    "<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title="Population comparison of region-level methylation candidates",
        xaxis_title="Mean HI in region",
        yaxis_title="Mean Δβ in region",
    )
    return apply_plot_style(fig)


# ============================================================
# Validate data and build controls
# ============================================================

population_columns = get_table_columns("tumor_summary_population")
if not population_columns:
    st.error(
        "The table `tumor_summary_population` is not available. This page uses "
        "the same population summary table as Population Sensitivity Explorer."
    )
    st.stop()

try:
    available = load_available_population_data()
except Exception as exc:
    st.error(f"Could not read population summary data: {exc}")
    st.stop()

if available.empty:
    st.warning("`tumor_summary_population` exists but contains no rows.")
    st.stop()

required_modes = {"full", "asian_excluded"}
mode_sets = (
    available.groupby("tumor_type")["population_mode"]
    .apply(lambda values: set(values.dropna().astype(str)))
)
available_tumors = sorted(
    tumor for tumor, modes in mode_sets.items() if required_modes.issubset(modes)
)

if not available_tumors:
    st.error(
        "No tumor contains both `full` and `asian_excluded` population modes."
    )
    st.stop()

preferred_order = [
    tumor
    for tumor in ["COAD", "LIHC", "LUAD", "LUSC"]
    if tumor in available_tumors
]
preferred_order.extend(
    tumor for tumor in available_tumors if tumor not in preferred_order
)

st.sidebar.header("Population comparison")

tumor_type = st.sidebar.selectbox(
    "Tumor Type",
    preferred_order,
    format_func=lambda value: TUMOR_LABELS.get(value, value),
)

tumor_available = available[
    available["tumor_type"].astype(str) == tumor_type
].copy()

full_pan_modes = set(
    tumor_available.loc[
        tumor_available["population_mode"].astype(str) == "full",
        "pan_reference_mode",
    ].dropna().astype(str)
)
excluded_pan_modes = set(
    tumor_available.loc[
        tumor_available["population_mode"].astype(str) == "asian_excluded",
        "pan_reference_mode",
    ].dropna().astype(str)
)

# Keep the same fixed full pan-cancer reference used in the original
# population-sensitivity analysis, without exposing it as a sidebar control.
pan_reference_mode = "full"
if (
    pan_reference_mode not in full_pan_modes
    or pan_reference_mode not in excluded_pan_modes
):
    st.error(
        "The selected tumor does not contain both `full` and "
        "`asian_excluded` summaries under `pan_reference_mode = full`."
    )
    st.stop()

asian_only_available = bool(
    (
        tumor_available["population_mode"].astype(str).eq("asian_only")
        & tumor_available["pan_reference_mode"].astype(str).eq(
            pan_reference_mode
        )
    ).any()
)

st.sidebar.header("CpG filters")
min_delta = st.sidebar.slider("Δβ", 0.0, 1.0, 0.45, 0.01)
max_normal_median = st.sidebar.slider(
    "Median NT β", 0.0, 1.0, 0.06, 0.01
)
max_pan_normal_median = st.sidebar.slider(
    "Median PanCan NT β", 0.0, 1.0, 0.06, 0.01
)
max_pan_tumor_median = st.sidebar.slider(
    "Median PanCan T β", 0.0, 1.0, 0.06, 0.01
)
max_pb = st.sidebar.slider("Median PB β", 0.0, 1.0, 0.04, 0.01)
min_hi = st.sidebar.slider("HI", 0.0, 5.0, 1.5, 0.05)

cross_tumor_type = get_lung_cross_tumor_type(tumor_type)
if cross_tumor_type is not None:
    max_cross_tumor_median = st.sidebar.slider(
        f"Median {cross_tumor_type} T β",
        0.0,
        1.0,
        1.00,
        0.01,
        help=(
            "Uses the original full-cohort tumor_summary as a fixed "
            "cross-subtype specificity reference. Set to 1.00 to disable."
        ),
    )
else:
    max_cross_tumor_median = None

min_region_sites = st.sidebar.slider(
    "Minimum qualifying CpGs per region", 1, 20, 1, 1
)

# Fixed display behavior for the download page.
# Matching regions are always connected and all filtered regions are shown.
connect_matching_regions = True



# ============================================================
# Methylation sample counts
# ============================================================

methylation_sample_counts = load_methylation_sample_counts(
    tumor_type=tumor_type,
    pan_reference_mode=pan_reference_mode,
)

if methylation_sample_counts.empty:
    st.warning(
        "Methylation sample counts could not be read from "
        "`tumor_summary_population`. Confirm that it contains "
        "`tumor_total_samples` and `normal_total_samples`."
    )
else:
    sample_fig = make_sample_count_plot(methylation_sample_counts)
    if sample_fig is not None:
        sample_filename = (
            f"{tumor_type}_population_methylation_sample_counts"
        )
        st.plotly_chart(
            sample_fig,
            use_container_width=True,
            config=plot_svg_config(
                filename=sample_filename,
                height=520,
                width=PLOT_WIDTH,
            ),
        )
        show_large_font_svg_export(
            sample_fig,
            filename=f"{sample_filename}_large_font",
            expander_label="SVG export — methylation sample-count chart",
            height=520,
            width=PLOT_WIDTH,
        )

    asian_source = methylation_sample_counts.loc[
        methylation_sample_counts["mode"].eq("asian_only"),
        "count_source",
    ]
    if not asian_source.empty and asian_source.iloc[0].startswith("Calculated"):
        st.caption(
            "Asian-only sample counts were calculated as Full minus "
            "Asian-excluded. The plotted counts refer only to samples "
            "represented in the methylation matrix."
        )


# ============================================================
# Bubble data
# ============================================================

modes_to_load = ["full", "asian_excluded"]
if asian_only_available:
    modes_to_load.append("asian_only")

region_frames_all: dict[str, pd.DataFrame] = {}
region_frames: dict[str, pd.DataFrame] = {}

try:
    for mode in modes_to_load:
        cpgs = load_population_mode_cpgs(
            tumor_type=tumor_type,
            population_mode=mode,
            pan_reference_mode=pan_reference_mode,
            min_delta=min_delta,
            max_normal_median=max_normal_median,
            max_pan_normal_median=max_pan_normal_median,
            max_pan_tumor_median=max_pan_tumor_median,
            max_pb=max_pb,
            min_hi=min_hi,
            cross_tumor_type=cross_tumor_type,
            max_cross_tumor_median=max_cross_tumor_median,
        )
        cpgs = add_physical_region_id(cpgs)
        regions = aggregate_regions(cpgs, mode)

        if not regions.empty:
            regions = regions[
                regions["n_qualifying_sites"] >= min_region_sites
            ].copy()

            regions = regions.sort_values(
                ["n_qualifying_sites", "mean_delta", "mean_hi"],
                ascending=[False, False, False],
            ).reset_index(drop=True)

            valid_region_ids = set(regions["physical_region_id"].astype(str))
            cpgs = cpgs[
                cpgs["physical_region_id"].astype(str).isin(valid_region_ids)
            ].copy()

        region_frames_all[mode] = regions
        region_frames[mode] = regions.copy()
except Exception as exc:
    st.error(f"Could not build the population bubble plot: {exc}")
    st.stop()

if all(frame.empty for frame in region_frames_all.values()):
    st.warning("No regions passed the selected filters.")
    st.stop()

# Compact visual summary without calculating scores or rankings.
summary_columns = st.columns(len(modes_to_load))
for column, mode in zip(summary_columns, modes_to_load):
    with column:
        full_frame = region_frames_all.get(mode, pd.DataFrame())
        st.metric(
            MODE_LABELS[mode],
            f"{len(full_frame):,} regions",
        )

bubble_fig = make_bubble_plot(
    region_frames,
    connect_matching_regions=connect_matching_regions,
)
bubble_filename = (
    f"{tumor_type}_population_bubble_comparison"
    f"_delta-{min_delta:.2f}_hi-{min_hi:.2f}"
)
st.plotly_chart(
    bubble_fig,
    use_container_width=True,
    config=plot_svg_config(
        filename=bubble_filename,
        height=PLOT_HEIGHT,
        width=PLOT_WIDTH,
    ),
)
show_large_font_svg_export(
    bubble_fig,
    filename=f"{bubble_filename}_large_font",
    expander_label="SVG export — population bubble comparison",
    height=PLOT_HEIGHT,
    width=PLOT_WIDTH,
)




