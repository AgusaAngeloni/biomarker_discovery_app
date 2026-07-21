#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import copy
import re
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from db.queries import run_query

# ============================================================
# Page config
# ============================================================

st.set_page_config(
    page_title="Region Candidate Scatter",
    page_icon="🧬",
    layout="wide",
)
st.title("Region Candidate Explorer")

# ============================================================
# Helpers
# ============================================================

TUMOR_MAP = {
    "COAD (CRC)": "COAD",
    "LIHC (HCC)": "LIHC",
    "LUAD": "LUAD",
    "LUSC": "LUSC",
}


def get_lung_cross_tumor_type(tumor_type: str) -> str | None:
    """Return the opposite NSCLC subtype used for cross-tumor filtering."""
    if tumor_type == "LUAD":
        return "LUSC"
    if tumor_type == "LUSC":
        return "LUAD"
    return None


@st.cache_data(ttl=600, show_spinner=False)
def get_table_columns(table_name: str) -> set[str]:
    query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema = 'public'
      AND table_name = :table_name
    """
    df = run_query(query, params={"table_name": table_name})
    if df.empty:
        return set()
    return set(df["column_name"].astype(str).tolist())


def sql_col(alias: str, columns: set[str], candidates: list[str], fallback: str = "NULL") -> str:
    for col in candidates:
        if col in columns:
            return f'{alias}."{col}"'
    return fallback


def clean_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def make_coordinate_region_id(chrom: object, start: object, end: object) -> str:
    """
    Coordinate-derived ID used only as an auxiliary label.

    It does not replace the original biomarker_region.region_id.
    """
    chrom_clean = str(chrom).replace("chr", "").replace("CHR", "").strip()
    return f"chr{chrom_clean}:{int(start)}-{int(end)}"


def choose_primary_region_id(values: pd.Series) -> str:
    """
    Preserve the original gene region ID from biomarker_region.

    If the same gene coordinates contain multiple original IDs because of
    gene annotations, keep all of them in region_ids and use the first one
    alphabetically as the displayed gene_region_id.
    """
    ids = sorted({str(v) for v in values if str(v).strip()})
    return ids[0] if ids else ""


def split_gene_text(value: object) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "na"}:
        return []
    genes = [g.strip().upper() for g in text.replace(",", ";").split(";")]
    genes = [g for g in genes if g and g.lower() not in {"nan", "none", "null", "na"}]
    return list(dict.fromkeys(genes))


def join_gene_values(values: pd.Series) -> str:
    genes: list[str] = []
    for value in values:
        genes.extend(split_gene_text(value))
    return ";".join(sorted(set(genes)))


def choose_main_gene(genes: list[str]) -> str:
    """Prefer simple protein-coding-like symbols over AC/LOC-style annotations."""
    clean = [str(g).strip().upper() for g in genes if str(g).strip()]
    clean = list(dict.fromkeys(clean))
    if not clean:
        return ""

    def priority(gene: str) -> tuple[int, int, str]:
        technical_prefixes = ("AC", "AL", "AP", "LOC", "LINC", "MIR")
        is_technical = gene.startswith(technical_prefixes) or "." in gene
        return (1 if is_technical else 0, len(gene), gene)

    return sorted(clean, key=priority)[0]


def alphabetic_category_orders(df: pd.DataFrame) -> dict[str, list[str]]:
    genes = sorted(df["gene_main"].dropna().astype(str).unique().tolist()) if "gene_main" in df.columns else []
    return {"gene_main": genes}


def plot_svg_config(filename: str, height: int = 700, width: int = 1400) -> dict:
    """
    Plotly toolbar SVG export configuration.

    These dimensions match the fixed Plotly layout applied by
    apply_gene_explorer_plot_style(), so downloaded SVGs keep the same
    typography, color and text size used on screen.
    """
    return {
        "toImageButtonOptions": {
            "format": "svg",
            "filename": filename,
            "height": height,
            "width": width,
            "scale": 1,
        }
    }


def apply_gene_explorer_plot_style(
    fig: go.Figure,
    height: int = 700,
    width: int = 1400,
    title_size: int = 20,
    axis_title_size: int = 18,
    axis_tick_size: int = 18,
    legend_size: int = 18,
    colorbar_title_size: int = 18,
    colorbar_tick_size: int = 18,
) -> go.Figure:
    gene_font = "Arial Narrow"

    fig.update_layout(
        template="plotly_white",
        autosize=False,
        width=width,
        height=height,
        paper_bgcolor="white",
        plot_bgcolor="white",
        hovermode="closest",
        margin=dict(l=40, r=40, t=50, b=50),
        font=dict(
            family=gene_font,
            color="black",
        ),
        title=dict(
            x=0,
            xanchor="left",
            font=dict(
                family=gene_font,
                size=title_size,
                color="black",
            ),
        ),
        legend=dict(
            orientation="v",
            y=1,
            x=1.06,
            bgcolor="rgba(255,255,255,0)",
            borderwidth=0,
            groupclick="togglegroup",
            tracegroupgap=1,
            font=dict(
                family=gene_font,
                size=legend_size,
                color="black",
            ),
        ),
        coloraxis_colorbar=dict(
            title=dict(
                font=dict(
                    family=gene_font,
                    size=colorbar_title_size,
                    color="black",
                ),
            ),
            tickfont=dict(
                family=gene_font,
                size=colorbar_tick_size,
                color="black",
            ),
        ),
    )

    fig.update_xaxes(
        title_font=dict(
            family=gene_font,
            size=axis_title_size,
            color="black",
        ),
        tickfont=dict(
            family=gene_font,
            size=axis_tick_size,
            color="black",
        ),
        showgrid=True,
    )

    fig.update_yaxes(
        title_font=dict(
            family=gene_font,
            size=axis_title_size,
            color="black",
        ),
        tickfont=dict(
            family=gene_font,
            size=axis_tick_size,
            color="black",
        ),
        showgrid=True,
    )

    return fig


def make_large_font_export_figure(
    fig: go.Figure,
    title_size: int = 30,
    axis_title_size: int = 28,
    axis_tick_size: int = 26,
    legend_size: int = 22,
    colorbar_title_size: int = 24,
    colorbar_tick_size: int = 22,
) -> go.Figure:
    """
    Copy the on-screen Plotly figure and change only font sizes.

    The copied figure keeps the original traces, marker colors, marker sizes,
    color scales, dimensions, margins and legend position.
    """
    fig_export = copy.deepcopy(fig)

    fig_export.update_layout(
        title=dict(font=dict(size=title_size)),
        legend=dict(font=dict(size=legend_size)),
    )

    fig_export.update_xaxes(
        title_font=dict(size=axis_title_size),
        tickfont=dict(size=axis_tick_size),
    )

    fig_export.update_yaxes(
        title_font=dict(size=axis_title_size),
        tickfont=dict(size=axis_tick_size),
    )

    if "coloraxis" in fig_export.layout:
        fig_export.update_layout(
            coloraxis_colorbar=dict(
                title=dict(font=dict(size=colorbar_title_size)),
                tickfont=dict(size=colorbar_tick_size),
            )
        )

    return fig_export


def show_large_font_svg_export(
    fig: go.Figure,
    filename: str,
    expander_label: str,
) -> None:
    """
    Show a second browser-rendered Plotly figure for SVG export.

    The download is performed from Plotly's camera button, not through Kaleido.
    This preserves the same browser-rendered colors as the interactive chart.
    """
    fig_export = make_large_font_export_figure(fig)
    width = int(fig.layout.width or 1400)
    height = int(fig.layout.height or 700)
    export_key = re.sub(r"[^A-Za-z0-9_-]+", "_", filename)

    with st.expander(expander_label, expanded=False):
        st.caption(
            "Open this preview and use the camera icon in the Plotly toolbar. "
            "The SVG keeps the same dimensions, colors, markers and layout; "
            "only the font is larger."
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



def extract_selected_table_row(table_event) -> int | None:
    """Extract selected row index from st.dataframe(..., on_select='rerun')."""
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


def normalize_chrom_for_sql(chrom: object) -> str:
    """Return chromosome without chr/CHR prefix for SQL filtering."""
    return str(chrom).replace("chr", "").replace("CHR", "").strip()


def parse_semicolon_values(value: object) -> list[str]:
    """Parse semicolon-separated strings such as cpg_sites or region_ids."""
    if value is None or pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


@st.cache_data(ttl=600, show_spinner=True)
def load_complete_gene_profile(
    tumor_type: str,
    gene: str,
) -> pd.DataFrame:
    """
    Load the complete Gene Explorer-like methylation profile for the selected gene.

    This intentionally does NOT restrict CpGs to the selected region. The selected
    filtered region is added later as a light-blue region track on top of the
    full-gene profile.
    """
    ts_cols = get_table_columns("tumor_summary")

    hi_expr = sql_col("ts", ts_cols, ["hi_index", "dispersion_index", "HI_index"], "NULL")
    pan_tumor_expr = sql_col(
        "ts",
        ts_cols,
        ["pan_tumor_median", "pantumor_median", "panTumor_median"],
        "NULL",
    )
    pan_normal_expr = sql_col(
        "ts",
        ts_cols,
        ["pan_normal_median", "pannormal_median", "panNormal_median"],
        "NULL",
    )

    query = f"""
    WITH expr_best AS (
        SELECT
            site_id,
            tumor_type,
            MIN(spearman_r) AS spearman_r
        FROM expression_correlation
        GROUP BY site_id, tumor_type
    )
    SELECT DISTINCT
        ca.site_id,
        cgm.gene_symbol AS gene,
        ca.chr,
        ca.start_pos,
        ca.end_pos,
        ts.tumor_type,
        ts.delta_median,
        ts.tumor_median,
        ts.normal_median,
        {pan_tumor_expr} AS pan_tumor_median,
        {pan_normal_expr} AS pan_normal_median,
        {hi_expr} AS hi_index,
        cf.pb_median,
        eb.spearman_r
    FROM tumor_summary ts
    JOIN cpg_annotation ca
        ON ts.site_id = ca.site_id
    JOIN cpg_gene_map cgm
        ON ts.site_id = cgm.site_id
    LEFT JOIN cpg_features cf
        ON ts.site_id = cf.site_id
    LEFT JOIN expr_best eb
        ON eb.site_id = ts.site_id
       AND eb.tumor_type = ts.tumor_type
    WHERE
        ts.tumor_type = :tumor_type
        AND cgm.gene_symbol = :gene
    ORDER BY ca.start_pos
    LIMIT 5000
    """

    df = run_query(
        query,
        params={
            "tumor_type": tumor_type,
            "gene": str(gene).strip().upper(),
        },
    )

    if df.empty:
        return df

    numeric_cols = [
        "start_pos",
        "end_pos",
        "delta_median",
        "tumor_median",
        "normal_median",
        "pan_tumor_median",
        "pan_normal_median",
        "hi_index",
        "pb_median",
        "spearman_r",
    ]
    df = clean_numeric(df, numeric_cols)
    df = df.dropna(subset=["site_id", "start_pos"]).copy()
    df["start_pos"] = df["start_pos"].astype(int)

    return df.sort_values("start_pos").reset_index(drop=True)

@st.cache_data(ttl=600, show_spinner=True)
def load_region_candidate_cpgs(
    tumor_type: str,
    min_delta: float,
    max_normal_median: float,
    max_pan_normal_median: float,
    max_pan_tumor_median: float,
    max_pb: float,
    min_hi: float,
    cross_tumor_type: str | None = None,
    max_cross_tumor_median: float | None = None,
) -> pd.DataFrame:
    ts_cols = get_table_columns("tumor_summary")
    r_cols = get_table_columns("biomarker_region")
    seq_cols = get_table_columns("biomarker_region_sequence_score")

    hi_expr = sql_col("ts", ts_cols, ["hi_index", "dispersion_index", "HI_index"], "NULL")
    pan_tumor_expr = sql_col(
        "ts",
        ts_cols,
        ["pan_tumor_median", "pantumor_median", "panTumor_median"],
        "NULL",
    )
    pan_normal_expr = sql_col(
        "ts",
        ts_cols,
        ["pan_normal_median", "pannormal_median", "panNormal_median"],
        "NULL",
    )
    genes_all_expr = sql_col("r", r_cols, ["gene_symbols_all"], 'r."gene_symbol"')
    seq_score_expr = sql_col("seq", seq_cols, ["sequence_score"], "0")
    gcgc_density_expr = sql_col("seq", seq_cols, ["gcgc_density_per_100bp"], "NULL")
    gc_fraction_expr = sql_col("seq", seq_cols, ["gc_fraction"], "NULL")

    cross_join = ""
    cross_select = "NULL AS cross_tumor_type, NULL AS cross_tumor_median"
    cross_filter = ""
    if cross_tumor_type in {"LUAD", "LUSC"} and max_cross_tumor_median is not None:
        cross_join = """
    LEFT JOIN tumor_summary ts_cross
        ON ts_cross.site_id = ts.site_id
       AND ts_cross.tumor_type = :cross_tumor_type
        """
        cross_select = "ts_cross.tumor_type AS cross_tumor_type, ts_cross.tumor_median AS cross_tumor_median"
        cross_filter = "AND COALESCE(ts_cross.tumor_median, 1) <= :max_cross_tumor_median"

    query = f"""
    WITH expr_best AS (
        SELECT
            site_id,
            tumor_type,
            MIN(spearman_r) AS spearman_r
        FROM expression_correlation
        GROUP BY site_id, tumor_type
    )
    SELECT
        brc.region_id,
        r.gene_symbol AS gene,
        {genes_all_expr} AS genes_all,
        r.chr,
        r.core_start,
        r.core_end,
        r.browser_start,
        r.browser_end,
        r.browser_length,
        r.core_length,
        r.n_manifest_cpgs,
        COALESCE({seq_score_expr}, 0) AS sequence_score,
        {gcgc_density_expr} AS gcgc_density_per_100bp,
        {gc_fraction_expr} AS gc_fraction,
        brc.site_id,
        brc.start_pos,
        ts.tumor_type,
        ts.delta_median,
        ts.normal_median,
        {cross_select},
        {pan_tumor_expr} AS pan_tumor_median,
        {pan_normal_expr} AS pan_normal_median,
        {hi_expr} AS hi_index,
        cf.pb_median,
        eb.spearman_r
    FROM tumor_summary ts
    JOIN biomarker_region_cpg brc
        ON brc.site_id = ts.site_id
    JOIN biomarker_region r
        ON r.region_id = brc.region_id
    LEFT JOIN biomarker_region_sequence_score seq
        ON seq.region_id = r.region_id
    LEFT JOIN cpg_features cf
        ON cf.site_id = ts.site_id
    {cross_join}
    LEFT JOIN expr_best eb
        ON eb.site_id = ts.site_id
       AND eb.tumor_type = ts.tumor_type
    WHERE
        ts.tumor_type = :tumor_type
        AND ts.delta_median >= :min_delta
        AND ts.normal_median <= :max_normal_median
        AND COALESCE({pan_normal_expr}, 1) <= :max_pan_normal_median
        AND COALESCE({pan_tumor_expr}, 1) <= :max_pan_tumor_median
        AND COALESCE(cf.pb_median, 1) <= :max_pb
        AND COALESCE({hi_expr}, 0) >= :min_hi
        {cross_filter}
    ORDER BY ts.delta_median DESC
    """

    params = {
        "tumor_type": tumor_type,
        "min_delta": float(min_delta),
        "max_normal_median": float(max_normal_median),
        "max_pan_normal_median": float(max_pan_normal_median),
        "max_pan_tumor_median": float(max_pan_tumor_median),
        "max_pb": float(max_pb),
        "min_hi": float(min_hi),
    }
    if cross_tumor_type in {"LUAD", "LUSC"} and max_cross_tumor_median is not None:
        params["cross_tumor_type"] = cross_tumor_type
        params["max_cross_tumor_median"] = float(max_cross_tumor_median)

    return run_query(query, params=params)


@st.cache_data(ttl=600, show_spinner=True)
def load_unfiltered_region_universe_by_tumor(tumor_type: str) -> pd.DataFrame:
    """
    Count the available region-model universe before applying candidate filters.

    This is the shared CpG-to-region screening universe used as starting point,
    restricted to the selected tumor type only because the app works tumor by tumor.
    """
    query = """
    SELECT
        ts.tumor_type,
        COUNT(DISTINCT brc.site_id) AS total_cpg_sites_before_filters,
        COUNT(DISTINCT NULLIF(r.gene_symbol, '')) AS total_genes_before_filters,
        COUNT(DISTINCT brc.region_id) AS total_regions_before_filters,
        COUNT(*) AS total_region_cpg_links
    FROM tumor_summary ts
    JOIN biomarker_region_cpg brc
        ON brc.site_id = ts.site_id
    JOIN biomarker_region r
        ON r.region_id = brc.region_id
    WHERE ts.tumor_type = :tumor_type
    GROUP BY ts.tumor_type
    """

    df = run_query(query, params={"tumor_type": tumor_type})

    if df.empty:
        return df

    numeric_cols = [
        "total_cpg_sites_before_filters",
        "total_genes_before_filters",
        "total_regions_before_filters",
        "total_region_cpg_links",
    ]
    df = clean_numeric(df, numeric_cols)

    return df.reset_index(drop=True)



def aggregate_gene_regions(cpgs: pd.DataFrame, apply_expression_filter: bool) -> pd.DataFrame:
    if cpgs.empty:
        return pd.DataFrame()

    out = cpgs.copy()
    out = clean_numeric(
        out,
        [
            "core_start",
            "core_end",
            "browser_start",
            "browser_end",
            "browser_length",
            "core_length",
            "n_manifest_cpgs",
            "sequence_score",
            "gcgc_density_per_100bp",
            "gc_fraction",
            "start_pos",
            "delta_median",
            "normal_median",
            "cross_tumor_median",
            "pan_tumor_median",
            "pan_normal_median",
            "hi_index",
            "pb_median",
            "spearman_r",
        ],
    )

    gene_cols = ["chr", "browser_start", "browser_end"]
    out = out.drop_duplicates(subset=gene_cols + ["site_id"]).copy()

    group_cols = gene_cols + ["browser_length"]

    region = (
        out.groupby(group_cols, dropna=False)
        .agg(
            gene_region_id=("region_id", choose_primary_region_id),
            region_ids=("region_id", lambda x: ";".join(sorted(set(map(str, x))))),
            genes_all=("genes_all", join_gene_values),
            genes_from_rows=("gene", join_gene_values),
            n_qualifying_sites=("site_id", "nunique"),
            n_manifest_cpgs=("n_manifest_cpgs", "max"),
            core_start=("core_start", "max"),
            core_end=("core_end", "max"),
            sequence_size_core=("core_length", "max"),
            gcgc_density_per_100bp=("gcgc_density_per_100bp", "max"),
            gc_fraction=("gc_fraction", "max"),
            mean_delta=("delta_median", "mean"),
            mean_hi=("hi_index", "mean"),
            mean_normal_median=("normal_median", "mean"),
            mean_cross_tumor_median=("cross_tumor_median", "mean"),
            mean_pan_tumor_median=("pan_tumor_median", "mean"),
            mean_pan_normal_median=("pan_normal_median", "mean"),
            mean_pb_median=("pb_median", "mean"),
            mean_spearman_r=("spearman_r", "mean"),
            sequence_score=("sequence_score", "max"),
            cpg_sites=("site_id", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )

    # Keep gene_region_id as the original biomarker_region.region_id.
    # coordinate_region_id is only an auxiliary coordinate label.
    region["coordinate_region_id"] = [
        make_coordinate_region_id(row.chr, int(row.browser_start), int(row.browser_end))
        for row in region.itertuples(index=False)
    ]

    region["genes_all"] = np.where(
        region["genes_all"].astype(str).str.strip() != "",
        region["genes_all"],
        region["genes_from_rows"],
    )
    region["gene_main"] = region["genes_all"].apply(lambda x: choose_main_gene(split_gene_text(x)))

    region["fraction_qualifying_sites"] = (
        region["n_qualifying_sites"]
        / pd.to_numeric(region["n_manifest_cpgs"], errors="coerce").replace(0, np.nan)
    ).fillna(0).clip(lower=0, upper=1)

    region["expression_signal"] = (
        -pd.to_numeric(region["mean_spearman_r"], errors="coerce")
    ).clip(lower=0).fillna(0)
    region["expression_size"] = (region["expression_signal"] ** 2 + 0.01).astype(float)

    region["sequence_site_score"] = (
        pd.to_numeric(region["sequence_score"], errors="coerce").fillna(0).clip(lower=0)
        * pd.to_numeric(region["n_qualifying_sites"], errors="coerce").fillna(0).clip(lower=0)
    )

    if apply_expression_filter:
        region["expression_score_component"] = 100 * region["expression_signal"] * region["n_qualifying_sites"]
    else:
        region["expression_score_component"] = 0.0

    region["final_region_score"] = region["sequence_site_score"] + region["expression_score_component"]
    region["final_region_score_size"] = region["final_region_score"].clip(lower=0.01)

    region["region_label"] = (
        region["gene_main"].astype(str)
        + " | chr"
        + region["chr"].astype(str)
        + ":"
        + region["browser_start"].astype(int).astype(str)
        + "-"
        + region["browser_end"].astype(int).astype(str)
    )

    return region.replace([np.inf, -np.inf], np.nan)


def make_region_scatter(
    df: pd.DataFrame,
    size_col: str,
    title: str,
    size_label: str,
    size_max: int = 45,
    color_col: str = "gene_main",
    color_label: str = "Gene",
):
    """
    Region-level bubble plot with focused hover.

    Hover shows:
    - Gene | region
    - X axis: mean_hi
    - Y axis: mean_delta
    - Bubble diameter: size_col
    """
    category_orders = alphabetic_category_orders(df)

    labels = {
        "mean_hi": "Mean HI in region",
        "mean_delta": "Mean Δβ in region",
        "gene_main": "Gene",
        size_col: size_label,
        color_col: color_label,
    }

    plot_df = df.copy()

    plot_df["_hover_size_value"] = pd.to_numeric(
        plot_df[size_col],
        errors="coerce"
    )

    # For the methylation-expression bubble plot, marker size uses the
    # biologically interpretable value:
    #     expression_signal = max(0, -mean_spearman_r)
    plot_df["_hover_mean_spearman_r"] = (
        pd.to_numeric(plot_df.get("mean_spearman_r", np.nan), errors="coerce")
        if "mean_spearman_r" in plot_df.columns
        else np.nan
    )
    plot_df["_hover_expression_signal"] = (
        pd.to_numeric(plot_df.get("expression_signal", np.nan), errors="coerce")
        if "expression_signal" in plot_df.columns
        else np.nan
    )
    plot_df["_hover_expression_size"] = (
        pd.to_numeric(plot_df.get("expression_size", np.nan), errors="coerce")
        if "expression_size" in plot_df.columns
        else np.nan
    )

    # -------------------------------
    # Hover gene
    # -------------------------------
    if "gene_main" in plot_df.columns:
        plot_df["_hover_gene"] = plot_df["gene_main"].astype(str)
    elif "genes_all" in plot_df.columns:
        plot_df["_hover_gene"] = plot_df["genes_all"].astype(str)
    else:
        plot_df["_hover_gene"] = ""

    # -------------------------------
    # Hover region
    # -------------------------------
    if "physical_region_id" in plot_df.columns:
        plot_df["_hover_region"] = plot_df["physical_region_id"].astype(str)
    elif "gene_region_id" in plot_df.columns:
        plot_df["_hover_region"] = plot_df["gene_region_id"].astype(str)
    elif "coordinate_region_id" in plot_df.columns:
        plot_df["_hover_region"] = plot_df["coordinate_region_id"].astype(str)
    elif "region_ids" in plot_df.columns:
        plot_df["_hover_region"] = plot_df["region_ids"].astype(str)
    else:
        plot_df["_hover_region"] = ""

    fig = px.scatter(
        plot_df,
        x="mean_hi",
        y="mean_delta",
        size=size_col,
        color=color_col,
        custom_data=[
            "_hover_gene",
            "_hover_region",
            "_hover_size_value",
            "_hover_mean_spearman_r",
            "_hover_expression_signal",
            "_hover_expression_size",
        ],
        labels=labels,
        category_orders=category_orders,
        title=title,
        size_max=size_max,
    )

    if size_col == "expression_signal":
        fig.update_traces(
            hovertemplate=(
                "<b>%{customdata[0]} | %{customdata[1]}</b><br>"
                "Mean HI in region: %{x:.3f}<br>"
                "Mean Δβ in region: %{y:.3f}<br>"
                "Mean Spearman r: %{customdata[3]:.3f}<br>"
                "Methylation-expression signal: %{customdata[4]:.3f}<br>"
                "<extra></extra>"
            )
        )
    else:
        fig.update_traces(
            hovertemplate=(
                "<b>%{customdata[0]} | %{customdata[1]}</b><br>"
                "Mean HI in region: %{x:.3f}<br>"
                "Mean Δβ in region: %{y:.3f}<br>"
                f"{size_label}: " + "%{customdata[2]:.3f}"
                "<extra></extra>"
            )
        )

    fig.update_layout(
        legend_traceorder="normal",
    )

    fig = apply_gene_explorer_plot_style(fig, height=700, width=1400)

    return fig


# ============================================================
# Sidebar
# ============================================================

st.sidebar.header("CpG filters")

tumor_label = st.sidebar.selectbox("Tumor Type", list(TUMOR_MAP.keys()))
tumor_type = TUMOR_MAP[tumor_label]

min_delta = st.sidebar.slider("Min Δβ", 0.0, 1.0, 0.55, 0.01)
max_normal_median = st.sidebar.slider("Max Median NT β", 0.0, 1.0, 0.06, 0.01)

cross_tumor_type = get_lung_cross_tumor_type(tumor_type)
if cross_tumor_type is not None:
    max_cross_tumor_median = st.sidebar.slider(
        f"Max Median {cross_tumor_type} T β",
        0.0,
        1.0,
        0.06,
        0.01,
        help=(
            f"Cross-subtype filter. When analyzing {tumor_type}, this keeps only CpGs "
            f"with low tumor methylation in {cross_tumor_type}. Set to 1.00 to disable."
        ),
    )
else:
    max_cross_tumor_median = None

max_pan_normal_median = st.sidebar.slider("Max Median PanCan NT β", 0.0, 1.0, 0.06, 0.01)
max_pan_tumor_median = st.sidebar.slider("Max Median PanCan T β", 0.0, 1.0, 0.06, 0.01)
max_pb = st.sidebar.slider("Max Median PB β", 0.0, 1.0, 0.04, 0.01)
min_hi = st.sidebar.slider("Min HI", 0.0, 5.0, 2.4, 0.05)

st.sidebar.header("Methylation-Expression Association Filter")
apply_expression_filter = st.sidebar.checkbox("Apply Methylation-Expression Association Filter", value=False)
max_mean_spearman_r = st.sidebar.slider(
    "Max Mean Methylation-Expression Association",
    -1.0,
    1.0,
    -0.16,
    0.01,
    disabled=not apply_expression_filter,
)

st.sidebar.header("Export")
enable_large_font_svg_export = st.sidebar.checkbox(
    "Show SVG export previews with larger font",
    value=False,
    help=(
        "Keeps the main plots unchanged and adds browser-rendered export previews. "
        "Download each SVG from the Plotly camera icon."
    ),
)

# ============================================================
# General region universe overview
# ============================================================

unfiltered_universe = load_unfiltered_region_universe_by_tumor(tumor_type=tumor_type)

if not unfiltered_universe.empty:
    selected_universe_row = unfiltered_universe.iloc[0]

    st.subheader("Starting region universe before filters")
    st.caption(
        "Counts before applying Δβ, normal methylation, PanCan, PB and HI filters. "
        "This is the shared CpG-to-region screening universe used as the starting point."
    )

    u1, u2, u3, u4 = st.columns(4)
    with u1:
        st.metric(
            "CpGs mapped",
            f"{int(selected_universe_row['total_cpg_sites_before_filters']):,}",
        )
    with u2:
        st.metric(
            "Mapped Genes",
            f"{int(selected_universe_row['total_genes_before_filters']):,}",
        )
    with u3:
        st.metric(
            "Mapped Regions",
            f"{int(selected_universe_row['total_regions_before_filters']):,}",
        )
    with u4:
        st.metric(
            "Mapped Region-CpG links",
            f"{int(selected_universe_row['total_region_cpg_links']):,}",
        )
else:
    selected_universe_row = None


# ============================================================
# Load and aggregate
# ============================================================

cpgs = load_region_candidate_cpgs(
    tumor_type=tumor_type,
    min_delta=min_delta,
    max_normal_median=max_normal_median,
    max_pan_normal_median=max_pan_normal_median,
    max_pan_tumor_median=max_pan_tumor_median,
    max_pb=max_pb,
    min_hi=min_hi,
    cross_tumor_type=cross_tumor_type,
    max_cross_tumor_median=max_cross_tumor_median,
)

if cpgs.empty:
    st.warning("No CpGs passed the current filters.")
    st.stop()

regions = aggregate_gene_regions(cpgs, apply_expression_filter=apply_expression_filter)

if apply_expression_filter:
    regions = regions[
        regions["mean_spearman_r"].notna()
        & (regions["mean_spearman_r"] <= max_mean_spearman_r)
    ].copy()

regions = regions.sort_values(
    ["final_region_score", "n_qualifying_sites", "mean_delta", "mean_hi"],
    ascending=False,
).reset_index(drop=True)

if regions.empty:
    st.warning("No regions passed the current filters.")
    st.stop()


# ============================================================
# Metrics
# ============================================================

qualifying_cpg_count = int(cpgs["site_id"].nunique())
candidate_gene_count = int(regions["gene_main"].nunique()) if "gene_main" in regions.columns else 0
candidate_region_count = int(regions["gene_region_id"].nunique())

m1, m2, m3, m4, m5 = st.columns(5)
with m1:
    st.metric("Qualifying CpGs", f"{qualifying_cpg_count:,}")
with m2:
    st.metric("Candidate genes", f"{candidate_gene_count:,}")
with m3:
    st.metric("Gene regions", f"{candidate_region_count:,}")
with m4:
    st.metric("Max CpGs in a region", f"{int(regions['n_qualifying_sites'].max()):,}")
with m5:
    st.metric("Best final score", f"{regions['final_region_score'].max():.2f}")

# ============================================================
# Plots
# ============================================================

st.subheader("Region Candidates")

plot_df = regions.head(1000).copy()
plot_df = plot_df.sort_values("gene_main", ascending=True).reset_index(drop=True)

fig_sites = make_region_scatter(
    plot_df,
    size_col="n_qualifying_sites",
    title=f"{tumor_type} - DMR Plot",
    size_label="Qualifying CpGs in region",
    size_max=45,
)
st.plotly_chart(
    fig_sites,
    use_container_width=True,
    config=plot_svg_config(
        filename=f"{tumor_type} - DMR Plot",
        height=700,
        width=1400,
    ),
)
if enable_large_font_svg_export:
    show_large_font_svg_export(
        fig_sites,
        filename=f"{tumor_type}_DMR_plot_larger_font",
        expander_label="DMR plot — SVG export with larger font",
    )

fig_expression = make_region_scatter(
    plot_df,
    size_col="expression_signal",
    title=f"{tumor_type} - Methylation-expression association Plot",
    size_label="Methylation-expression signal",
    size_max=45,
)
st.plotly_chart(
    fig_expression,
    use_container_width=True,
    config=plot_svg_config(
        filename=f"{tumor_type} - Methylation-expression association Plot",
        height=700,
        width=1400,
    ),
)
if enable_large_font_svg_export:
    show_large_font_svg_export(
        fig_expression,
        filename=f"{tumor_type}_methylation_expression_plot_larger_font",
        expander_label="Methylation-expression plot — SVG export with larger font",
    )

score_plot_df = plot_df[pd.to_numeric(plot_df["final_region_score"], errors="coerce").fillna(0) > 0].copy()
if score_plot_df.empty:
    st.warning("No positive final_region_score values available for score-sized plot.")
else:
    fig_score = make_region_scatter(
        score_plot_df,
        size_col="final_region_score",
        title=f"{tumor_type} - Region-level candidates sized by region score",
        size_label="Final score",
        size_max=55,
        color_col="final_region_score",
        color_label="Final score",
    )
    st.plotly_chart(
        fig_score,
        use_container_width=True,
        config=plot_svg_config(
            filename=f"{tumor_type} - Region-level candidates sized by region score",
            height=700,
            width=1400,
        ),
    )
    if enable_large_font_svg_export:
        show_large_font_svg_export(
            fig_score,
            filename=f"{tumor_type}_region_score_plot_larger_font",
            expander_label="Region score plot — SVG export with larger font",
        )

sequence_size_plot_df = plot_df[
    pd.to_numeric(plot_df["sequence_size_core"], errors="coerce").fillna(0) > 0
].copy()


# ============================================================
# Region table
# ============================================================

st.subheader("Region Candidate Table")

region_cols = [
    "gene_main",
    "chr",
    "final_region_score",
    "sequence_score",
    "gcgc_density_per_100bp",
    "gc_fraction",
    "sequence_size_core",
    "core_start",
    "core_end",
    "n_qualifying_sites",
    "n_manifest_cpgs",
    "fraction_qualifying_sites",
    "mean_delta",
    "mean_hi",
    "mean_normal_median",
    "mean_cross_tumor_median",
    "mean_pan_tumor_median",
    "mean_pan_normal_median",
    "mean_pb_median",
    "sequence_site_score",
    "cpg_sites",
    "region_ids",    
    "mean_spearman_r",
    "expression_signal",
    "gene_region_id",
]
visible_cols = [c for c in region_cols if c in regions.columns]

region_table = regions[visible_cols].head(1000).reset_index(drop=True).copy()

st.dataframe(
    region_table,
    use_container_width=True,
    hide_index=True,
)

st.download_button(
    "Download region candidates TSV",
    data=regions.to_csv(index=False, sep="\t").encode("utf-8"),
    file_name=f"region_candidate_scatter_{tumor_type.lower()}.tsv",
    mime="text/tab-separated-values",
)

st.info(
    "Score: sequence_site_score = sequence_score × n_qualifying_sites. "
    "If the expression filter is enabled, final_region_score also includes "
    "100 × max(0, -mean methylation-expression association) × n_qualifying_sites. "
)


