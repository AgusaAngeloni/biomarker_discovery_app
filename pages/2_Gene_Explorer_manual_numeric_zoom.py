import html
import io
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
    page_icon="🧬",
    layout="wide",
)

st.title("Gene Explorer — Candidate Regions")

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
    Detect true GCGC motifs and assign CSS classes by genomic position.

    Important:
    Do not iterate over the string "GCGC", because that iterates over
    the characters G, C, G and C and ends up marking almost every G/C
    base as a motif, which makes the browser look black.
    """
    roles: dict[int, list[str]] = {}
    sequence = sequence.upper()

    for match in re.finditer("GCGC", sequence):
        motif_start = int(seq_start) + match.start()
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
            <span class="legend-item"><span class="dot dot-green"></span>CpG manifest</span>
            <span class="legend-item"><span class="box-black"></span>GCGC motif</span>
        </div>
        """
    )

    return "".join(html_parts)



def render_sequence_svg(
    sequence: str,
    seq_start: int,
    manifest_cpg_positions: set[int],
    gcgc_roles: dict[int, list[str]],
    title: str,
    subtitle: str = "",
    bases_per_row: int = BASES_PER_ROW,
) -> str:
    """
    Render the displayed sequence browser as SVG for download.
    """
    sequence = sequence.upper()

    cell_w = 18
    cell_h = 22
    row_h = 28
    coord_w = 115
    margin_x = 24
    margin_y = 70
    title_h = 42
    legend_h = 42

    n_rows = int(np.ceil(len(sequence) / bases_per_row))
    width = margin_x * 2 + coord_w + bases_per_row * cell_w + 20
    height = margin_y + title_h + n_rows * row_h + legend_h + 24

    base_colors = {
        "A": "#16a34a",
        "C": "#ef4444",
        "G": "#2563eb",
        "T": "#9333ea",
    }
    
    svg_font = "'Arial Narrow', Arial, sans-serif"
    svg_mono_font = "'Arial Narrow', Arial, sans-serif" 

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin_x}" y="30" font-family="{svg_font}" font-size="18" font-weight="700"     fill="#102a43">{html.escape(title)}</text>',
        f'<text x="{margin_x}" y="52" font-family="{svg_font}" font-size="13" fill="#52616b">{html.escape(subtitle)}</text>',
    ]

    y0 = margin_y + title_h

    for row_start in range(0, len(sequence), bases_per_row):
        row_seq = sequence[row_start:row_start + bases_per_row]
        row_idx = row_start // bases_per_row
        y = y0 + row_idx * row_h
        row_genomic_start = seq_start + row_start

        parts.append(
            f'<text x="{margin_x + coord_w - 8}" y="{y + 16}" '
            f'text-anchor="end" font-family="{svg_font}" font-size="12" '
            f'fill="#9aa5b1">{row_genomic_start:,}</text>'
        )

        for i, base in enumerate(row_seq):
            genomic_pos = row_genomic_start + i
            x = margin_x + coord_w + i * cell_w
            safe_base = html.escape(base)

            fill = "transparent"
            stroke = "transparent"
            stroke_w = 1
            text_fill = base_colors.get(base, "#6b7280")
            font_weight = 700

            if genomic_pos in manifest_cpg_positions:
                fill = "#bbf7d0"
                stroke = "#22c55e"
                stroke_w = 2
                text_fill = "#000000"

            if genomic_pos in gcgc_roles:
                fill = "#000000"
                stroke = "#000000"
                stroke_w = 2
                text_fill = "#ffffff"
                font_weight = 900

            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h}" '
                f'rx="3" ry="3" fill="{fill}" stroke="{stroke}" '
                f'stroke-width="{stroke_w}"/>'
            )

            parts.append(
                f'<text x="{x + (cell_w - 2) / 2}" y="{y + 16}" '
                f'text-anchor="middle" font-family="{svg_mono_font}" '
                f'font-size="14" font-weight="{font_weight}" '
                f'fill="{text_fill}">{safe_base}</text>'
            )

    legend_y = y0 + n_rows * row_h + 26

    parts.extend([
    f'<circle cx="{margin_x + 8}" cy="{legend_y - 4}" r="6" fill="#22c55e"/>',
    f'<text x="{margin_x + 22}" y="{legend_y}" font-family="{svg_font}" font-size="13" fill="#52616b">CpG manifest</text>',

    f'<rect x="{margin_x + 210}" y="{legend_y - 13}" width="16" height="12" fill="#000000"/>',
    f'<text x="{margin_x + 232}" y="{legend_y}" font-family="{svg_font}" font-size="13" fill="#52616b">GCGC motif</text>',

    '</svg>',
    ])

    return "".join(parts)


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


@st.cache_data(ttl=600, show_spinner=False)
def get_table_columns(table_name: str) -> set[str]:
    """
    Return PostgreSQL column names for a public table.

    Used to keep the page compatible with slightly different versions of
    tumor_summary, cpg_features, biomarker_cpg_score and sequence-score tables.
    """
    query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema = 'public'
      AND table_name = :table_name
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
    candidates: list[str],
    fallback: str = "NULL",
) -> str:
    """
    Return the first available SQL column expression from a list of candidates.
    """
    for col in candidates:
        if col in columns:
            return f'{alias}."{col}"'
    return fallback


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
        ts.pan_tumor_median AS pan_tumor_median,
        ts.pan_normal_median AS pan_normal_median,
        cf.leukocyte_median,
        ec.spearman_r
    FROM tumor_summary ts
    JOIN cpg_annotation ca
        ON ts.site_id = ca.site_id
    JOIN cpg_gene_map cgm
        ON ts.site_id = cgm.site_id
    JOIN cpg_features cf
        ON ts.site_id = cf.site_id
    JOIN expression_correlation ec
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
            "gene_pattern": f"%{normalize_gene_symbol(gene)}%",
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



@st.cache_data(ttl=600, show_spinner=False)
def load_gene_regions_for_profile(gene: str) -> pd.DataFrame:
    """
    Load physical CpG regions early so they can be drawn as a light-blue
    region track inside the methylation gene profile plot.
    """
    query = """
    SELECT
        r.region_id,
        r.gene_symbol,
        r.chr,
        r.core_start,
        r.core_end,
        r.core_length,
        r.browser_start,
        r.browser_end,
        r.browser_length,
        r.flank_bp,
        r.n_manifest_cpgs,
        r.n_manifest_c,
        r.cpg_density_per_100bp,
        r.region_rank_by_density,
        COALESCE(seq.sequence_score, 0) AS sequence_score,
        seq.sequence_available,
        seq.sequence_length,
        seq.n_cg_sequence,
        seq.n_gcgc,
        seq.gc_fraction
    FROM biomarker_region r
    LEFT JOIN biomarker_region_sequence_score seq
        ON r.region_id = seq.region_id
    WHERE r.gene_symbol = :gene
    ORDER BY r.chr, r.core_start, r.core_end
    """

    df = run_query(query, params={"gene": normalize_gene_symbol(gene)})

    if df.empty:
        return df

    numeric_cols = [
        "core_start", "core_end", "core_length", "browser_start", "browser_end",
        "browser_length", "flank_bp", "n_manifest_cpgs", "n_manifest_c",
        "cpg_density_per_100bp", "region_rank_by_density", "sequence_score",
        "sequence_length", "n_cg_sequence", "n_gcgc", "gc_fraction",
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["gene_symbol"] = df["gene_symbol"].map(normalize_gene_symbol)
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", case=False, regex=False).str.strip()
    df = df.dropna(subset=["region_id", "gene_symbol", "chr", "core_start", "core_end"]).copy()

    return df.reset_index(drop=True)




@st.cache_data(ttl=600, show_spinner=False)
def load_filtered_candidate_region_cpgs_for_gene(
    gene: str,
    tumor_type: str,
    min_delta: float,
    max_normal_median: float,
    max_pan_normal_median: float,
    max_pan_tumor_median: float,
    max_leukocyte: float,
    min_hi: float,
) -> pd.DataFrame:
    """
    Load CpGs inside physical regions that pass the selected candidate filters.

    This uses the existing PostgreSQL region model directly and does not need
    biomarker_candidate_region_curve or any precomputed parquet/table.
    """
    ts_cols = get_table_columns("tumor_summary")
    r_cols = get_table_columns("biomarker_region")
    seq_cols = get_table_columns("biomarker_region_sequence_score")
    cf_cols = get_table_columns("cpg_features")
    bcs_cols = get_table_columns("biomarker_cpg_score")
    ec_cols = get_table_columns("expression_correlation")

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

    seq_join = ""
    seq_select = """
        NULL AS sequence_score,
        NULL AS sequence_available,
        NULL AS sequence_length,
        NULL AS n_c_sequence,
        NULL AS n_g_sequence,
        NULL AS n_cg_sequence,
        NULL AS n_gcgc,
        NULL AS gc_fraction
    """
    if seq_cols:
        seq_join = """
        LEFT JOIN biomarker_region_sequence_score seq
            ON seq.region_id = r.region_id
        """
        seq_select = f"""
            {sql_col("seq", seq_cols, ["sequence_score"], "NULL")} AS sequence_score,
            {sql_col("seq", seq_cols, ["sequence_available"], "NULL")} AS sequence_available,
            {sql_col("seq", seq_cols, ["sequence_length"], "NULL")} AS sequence_length,
            {sql_col("seq", seq_cols, ["n_c_sequence", "n_c", "n_cytosine"], "NULL")} AS n_c_sequence,
            {sql_col("seq", seq_cols, ["n_g_sequence", "n_g", "n_guanine"], "NULL")} AS n_g_sequence,
            {sql_col("seq", seq_cols, ["n_cg_sequence", "n_cg", "cg_count"], "NULL")} AS n_cg_sequence,
            {sql_col("seq", seq_cols, ["n_gcgc", "gcgc_count"], "NULL")} AS n_gcgc,
            {sql_col("seq", seq_cols, ["gc_fraction"], "NULL")} AS gc_fraction
        """

    cf_join = ""
    leukocyte_select = "NULL AS leukocyte_median"
    if cf_cols:
        cf_join = """
        LEFT JOIN cpg_features cf
            ON cf.site_id = brc.site_id
        """
        leukocyte_select = f'{sql_col("cf", cf_cols, ["leukocyte_median"], "NULL")} AS leukocyte_median'

    bcs_join = ""
    bcs_select = """
        NULL AS biological_score,
        NULL AS delta_score,
        NULL AS normal_low_score,
        NULL AS leukocyte_low_score,
        NULL AS hi_score,
        NULL AS expression_score,
        NULL AS passes_loose_seed,
        NULL AS passes_default_filter,
        NULL AS passes_strict_filter
    """
    if bcs_cols:
        if "gene_symbol" in bcs_cols:
            bcs_join = """
            LEFT JOIN biomarker_cpg_score bcs
                ON bcs.site_id = brc.site_id
               AND bcs.tumor_type = ts.tumor_type
               AND (
                    bcs.gene_symbol = brc.gene_symbol
                    OR bcs.gene_symbol IS NULL
                    OR bcs.gene_symbol = ''
               )
            """
        else:
            bcs_join = """
            LEFT JOIN biomarker_cpg_score bcs
                ON bcs.site_id = brc.site_id
               AND bcs.tumor_type = ts.tumor_type
            """
        bcs_select = f"""
            {sql_col("bcs", bcs_cols, ["biological_score"], "NULL")} AS biological_score,
            {sql_col("bcs", bcs_cols, ["delta_score"], "NULL")} AS delta_score,
            {sql_col("bcs", bcs_cols, ["normal_low_score"], "NULL")} AS normal_low_score,
            {sql_col("bcs", bcs_cols, ["leukocyte_low_score"], "NULL")} AS leukocyte_low_score,
            {sql_col("bcs", bcs_cols, ["hi_score"], "NULL")} AS hi_score,
            {sql_col("bcs", bcs_cols, ["expression_score"], "NULL")} AS expression_score,
            {sql_col("bcs", bcs_cols, ["passes_loose_seed"], "NULL")} AS passes_loose_seed,
            {sql_col("bcs", bcs_cols, ["passes_default_filter"], "NULL")} AS passes_default_filter,
            {sql_col("bcs", bcs_cols, ["passes_strict_filter"], "NULL")} AS passes_strict_filter
        """

    expr_cte = ""
    expr_join = ""
    spearman_select = "NULL AS spearman_r"
    if ec_cols:
        expr_cte = """
        expr_best AS (
            SELECT
                site_id,
                tumor_type,
                MIN(spearman_r) AS spearman_r
            FROM expression_correlation
            GROUP BY site_id, tumor_type
        ),
        """
        expr_join = """
        LEFT JOIN expr_best eb
            ON eb.site_id = brc.site_id
           AND eb.tumor_type = ts.tumor_type
        """
        spearman_select = "eb.spearman_r AS spearman_r"

    query = f"""
    WITH
        {expr_cte}
        candidate_cpgs AS (
            SELECT DISTINCT
                r.region_id,
                r.gene_symbol AS region_gene_symbol,
                CAST({genes_all_expr} AS TEXT) AS region_genes_all,
                brc.gene_symbol AS cpg_gene_symbol,
                r.chr,
                r.core_start,
                r.core_end,
                r.core_length,
                r.browser_start,
                r.browser_end,
                r.browser_length,
                r.flank_bp,
                r.n_manifest_cpgs,
                r.n_manifest_c,
                r.cpg_density_per_100bp,
                r.region_rank_by_density,
                {seq_select},
                brc.site_id,
                brc.start_pos,
                brc.cpg_order,
                ts.tumor_type,
                ts.delta_median,
                ts.tumor_median,
                ts.normal_median,
                {pan_tumor_expr} AS pan_tumor_median,
                {pan_normal_expr} AS pan_normal_median,
                {hi_expr} AS hi_index,
                {leukocyte_select},
                {spearman_select},
                {bcs_select}
            FROM biomarker_region r
            JOIN biomarker_region_cpg brc
                ON brc.region_id = r.region_id
            JOIN tumor_summary ts
                ON ts.site_id = brc.site_id
               AND ts.tumor_type = :tumor_type
            {seq_join}
            {cf_join}
            {expr_join}
            {bcs_join}
            WHERE
                (
                    r.gene_symbol = :gene
                    OR brc.gene_symbol = :gene
                    OR CAST({genes_all_expr} AS TEXT) ILIKE :gene_pattern
                )
                AND ts.delta_median >= :min_delta
                AND ts.normal_median <= :max_normal_median
                AND COALESCE({pan_normal_expr}, 1) <= :max_pan_normal_median
                AND COALESCE({pan_tumor_expr}, 1) <= :max_pan_tumor_median
                AND COALESCE({leukocyte_select.split(" AS ")[0]}, 1) <= :max_leukocyte
                AND COALESCE({hi_expr}, 0) >= :min_hi
        )
    SELECT *
    FROM candidate_cpgs
    ORDER BY chr, core_start, core_end, cpg_order, start_pos
    """

    df = run_query(
        query,
        params={
            "gene": normalize_gene_symbol(gene),
            "gene_pattern": f"%{normalize_gene_symbol(gene)}%",
            "tumor_type": tumor_type,
            "min_delta": float(min_delta),
            "max_normal_median": float(max_normal_median),
            "max_pan_normal_median": float(max_pan_normal_median),
            "max_pan_tumor_median": float(max_pan_tumor_median),
            "max_leukocyte": float(max_leukocyte),
            "min_hi": float(min_hi),
        },
    )

    if df.empty:
        return df

    numeric_cols = [
        "core_start", "core_end", "core_length", "browser_start", "browser_end",
        "browser_length", "flank_bp", "n_manifest_cpgs", "n_manifest_c",
        "cpg_density_per_100bp", "region_rank_by_density", "sequence_score",
        "sequence_length", "n_c_sequence", "n_g_sequence", "n_cg_sequence", "n_gcgc",
        "gc_fraction", "start_pos", "cpg_order", "delta_median", "tumor_median",
        "normal_median", "pan_tumor_median", "pan_normal_median", "hi_index",
        "leukocyte_median", "spearman_r", "biological_score", "delta_score",
        "normal_low_score", "leukocyte_low_score", "hi_score", "expression_score",
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["chr"] = (
        df["chr"]
        .astype(str)
        .str.replace("chr", "", case=False, regex=False)
        .str.strip()
    )

    df = df.dropna(subset=["region_id", "site_id", "start_pos", "core_start", "core_end"]).copy()
    df["start_pos"] = df["start_pos"].astype(int)

    return df.reset_index(drop=True)


def aggregate_filtered_candidate_regions(
    candidate_cpgs: pd.DataFrame,
    apply_expression_filter: bool,
    max_mean_spearman_r: float,
) -> pd.DataFrame:
    """
    Aggregate candidate CpGs to one row per physical region.

    The region_id is preserved from biomarker_region. No hash/REGION ID is generated.
    """
    if candidate_cpgs.empty:
        return pd.DataFrame()

    df = candidate_cpgs.copy()

    group_cols = [
        "region_id",
        "region_gene_symbol",
        "region_genes_all",
        "chr",
        "core_start",
        "core_end",
        "core_length",
        "browser_start",
        "browser_end",
        "browser_length",
    ]

    region_summary = (
        df.groupby(group_cols, dropna=False)
        .agg(
            n_qualifying_cpgs=("site_id", "nunique"),
            qualifying_cpg_sites=("site_id", lambda x: ";".join(sorted(set(map(str, x))))),
            qualifying_cpg_positions=("start_pos", lambda x: ";".join(map(str, sorted(set(map(int, x)))))),
            n_manifest_cpgs=("n_manifest_cpgs", "max"),
            n_manifest_c=("n_manifest_c", "max"),
            cpg_density_per_100bp=("cpg_density_per_100bp", "max"),
            region_rank_by_density=("region_rank_by_density", "max"),
            sequence_score=("sequence_score", "max"),
            sequence_available=("sequence_available", "first"),
            sequence_length=("sequence_length", "max"),
            n_c_sequence=("n_c_sequence", "max"),
            n_g_sequence=("n_g_sequence", "max"),
            n_cg_sequence=("n_cg_sequence", "max"),
            n_gcgc=("n_gcgc", "max"),
            gc_fraction=("gc_fraction", "max"),
            mean_delta=("delta_median", "mean"),
            max_delta=("delta_median", "max"),
            mean_hi=("hi_index", "mean"),
            max_hi=("hi_index", "max"),
            mean_tumor_median=("tumor_median", "mean"),
            max_tumor_median=("tumor_median", "max"),
            mean_normal_median=("normal_median", "mean"),
            min_normal_median=("normal_median", "min"),
            mean_pan_tumor_median=("pan_tumor_median", "mean"),
            mean_pan_normal_median=("pan_normal_median", "mean"),
            mean_leukocyte_median=("leukocyte_median", "mean"),
            min_leukocyte_median=("leukocyte_median", "min"),
            mean_spearman_r=("spearman_r", "mean"),
            min_spearman_r=("spearman_r", "min"),
            max_biological_score=("biological_score", "max"),
            mean_biological_score=("biological_score", "mean"),
            max_delta_score=("delta_score", "max"),
            max_normal_low_score=("normal_low_score", "max"),
            max_leukocyte_low_score=("leukocyte_low_score", "max"),
            max_hi_score=("hi_score", "max"),
            max_expression_score=("expression_score", "max"),
        )
        .reset_index()
    )

    if apply_expression_filter:
        region_summary = region_summary[
            region_summary["mean_spearman_r"].notna()
            & (region_summary["mean_spearman_r"] <= float(max_mean_spearman_r))
        ].copy()

        if not region_summary.empty:
            valid_region_ids = set(region_summary["region_id"].astype(str))
            candidate_cpgs.drop(
                candidate_cpgs.index[
                    ~candidate_cpgs["region_id"].astype(str).isin(valid_region_ids)
                ],
                inplace=True,
            )

    region_summary["fraction_qualifying_cpgs"] = (
        pd.to_numeric(region_summary["n_qualifying_cpgs"], errors="coerce")
        / pd.to_numeric(region_summary["n_manifest_cpgs"], errors="coerce").replace(0, np.nan)
    ).fillna(0).clip(0, 1)

    region_summary["region_label"] = (
        region_summary["region_id"].astype(str)
        + " | chr"
        + region_summary["chr"].astype(str)
        + ":"
        + region_summary["core_start"].astype(int).astype(str)
        + "-"
        + region_summary["core_end"].astype(int).astype(str)
    )

    return region_summary.replace([np.inf, -np.inf], np.nan).sort_values(
        ["max_tumor_median", "n_qualifying_cpgs", "mean_delta"],
        ascending=False,
    ).reset_index(drop=True)


def build_filtered_region_track_trace(
    filtered_regions: pd.DataFrame,
    x_min: int,
    x_max: int,
    y_level: float = -0.16,
) -> tuple[list, list, list]:
    """
    Build x/y/hover vectors for Plotly line segments marking filtered regions.
    """
    region_x: list = []
    region_y: list = []
    region_hover: list = []

    if filtered_regions.empty:
        return region_x, region_y, region_hover

    region_track = filtered_regions.copy()
    region_track = region_track[
        (pd.to_numeric(region_track["core_end"], errors="coerce") >= x_min)
        & (pd.to_numeric(region_track["core_start"], errors="coerce") <= x_max)
    ].copy()

    for _, region_row in region_track.iterrows():
        r_start = int(region_row["core_start"])
        r_end = int(region_row["core_end"])
        region_id = str(region_row["region_id"])

        region_x.extend([r_start, r_end, None])
        region_y.extend([y_level, y_level, None])
        region_hover.extend([region_id, region_id, None])

    return region_x, region_y, region_hover




def make_editable_profile_svg(
    graph: pd.DataFrame,
    filtered_regions: pd.DataFrame,
    tumor_type: str,
    gene: str,
    x_range: tuple[int, int],
    show_hi: bool,
    y2_max: float,
) -> bytes:
    """
    Build an Illustrator-friendly SVG version of the methylation profile.

    This export uses the same visual configuration as the Plotly curve plot:
    - smoothed methylation curves
    - teal tumor curve with light fill
    - orange NT curve
    - green PanCancer T / PanCancer curves
    - purple leukocyte curve
    - black delta markers
    - red HI markers when enabled
    - crimson candidate-region borders and bottom track
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except Exception as exc:
        raise RuntimeError(
            "matplotlib is required to generate the editable SVG export. "
            "Add matplotlib to the environment if it is not installed."
        ) from exc

    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["font.family"] = "Arial"

    plot_graph = graph.copy().sort_values("start_pos")
    x = pd.to_numeric(plot_graph["start_pos"], errors="coerce")

    x_start, x_end = int(x_range[0]), int(x_range[1])
    if x_end <= x_start:
        raise ValueError("Invalid x-axis range for SVG export.")

    def numeric_series(col: str) -> pd.Series:
        if col not in plot_graph.columns:
            return pd.Series(np.nan, index=plot_graph.index)
        return pd.to_numeric(plot_graph[col], errors="coerce")

    def smooth_curve_xy(
        x_values: pd.Series,
        y_values: pd.Series,
        points_per_segment: int = 24,
        clip_y: tuple[float, float] | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Return smoothed x/y arrays for SVG export.

        Preference:
        1. scipy PCHIP if available, because it preserves local shape well.
        2. Cubic Hermite fallback using finite-difference slopes.

        The fallback avoids adding scipy as a hard dependency.
        """
        data = pd.DataFrame(
            {
                "x": pd.to_numeric(x_values, errors="coerce"),
                "y": pd.to_numeric(y_values, errors="coerce"),
            }
        ).replace([np.inf, -np.inf], np.nan).dropna()

        data = data.sort_values("x").drop_duplicates(subset=["x"])
        data = data[(data["x"] >= x_start) & (data["x"] <= x_end)]

        if data.shape[0] < 2:
            return np.array([], dtype=float), np.array([], dtype=float)

        x_arr = data["x"].to_numpy(dtype=float)
        y_arr = data["y"].to_numpy(dtype=float)

        if data.shape[0] < 4:
            if clip_y is not None:
                y_arr = np.clip(y_arr, clip_y[0], clip_y[1])
            return x_arr, y_arr

        try:
            from scipy.interpolate import PchipInterpolator

            n_points = max(250, int((len(x_arr) - 1) * points_per_segment))
            xs = np.linspace(float(x_arr.min()), float(x_arr.max()), n_points)
            ys = PchipInterpolator(x_arr, y_arr)(xs)

            if clip_y is not None:
                ys = np.clip(ys, clip_y[0], clip_y[1])

            return xs, ys
        except Exception:
            pass

        n = len(x_arr)
        slopes = np.zeros(n, dtype=float)

        dx = np.diff(x_arr)
        dy = np.diff(y_arr)
        safe_dx = np.where(dx == 0, np.nan, dx)

        segment_slopes = dy / safe_dx
        segment_slopes = np.nan_to_num(segment_slopes, nan=0.0, posinf=0.0, neginf=0.0)

        slopes[0] = segment_slopes[0]
        slopes[-1] = segment_slopes[-1]

        for i in range(1, n - 1):
            denom = x_arr[i + 1] - x_arr[i - 1]
            slopes[i] = (y_arr[i + 1] - y_arr[i - 1]) / denom if denom != 0 else 0.0

        xs_parts: list[np.ndarray] = []
        ys_parts: list[np.ndarray] = []

        for i in range(n - 1):
            x0 = x_arr[i]
            x1 = x_arr[i + 1]
            h = x1 - x0

            if h <= 0:
                continue

            t = np.linspace(0, 1, points_per_segment, endpoint=False)
            h00 = 2 * t**3 - 3 * t**2 + 1
            h10 = t**3 - 2 * t**2 + t
            h01 = -2 * t**3 + 3 * t**2
            h11 = t**3 - t**2

            xs = x0 + t * h
            ys = (
                h00 * y_arr[i]
                + h10 * h * slopes[i]
                + h01 * y_arr[i + 1]
                + h11 * h * slopes[i + 1]
            )

            xs_parts.append(xs)
            ys_parts.append(ys)

        xs_parts.append(np.array([x_arr[-1]], dtype=float))
        ys_parts.append(np.array([y_arr[-1]], dtype=float))

        xs = np.concatenate(xs_parts)
        ys = np.concatenate(ys_parts)

        if clip_y is not None:
            ys = np.clip(ys, clip_y[0], clip_y[1])

        return xs, ys

    median_tumor_svg = numeric_series("tumor_median")
    median_normal_svg = numeric_series("normal_median")
    pan_tumor_svg = numeric_series("pan_tumor_median")
    pan_normal_svg = numeric_series("pan_normal_median")
    leukocyte_svg = numeric_series("leukocyte_median")
    delta_svg = numeric_series("delta_median")
    hi_svg = numeric_series("hi_index")

    fig, ax = plt.subplots(figsize=(14, 6))

    x_tumor, y_tumor = smooth_curve_xy(x, median_tumor_svg, clip_y=(-0.2, 1.0))
    x_normal, y_normal = smooth_curve_xy(x, median_normal_svg, clip_y=(-0.2, 1.0))
    x_pan_tumor, y_pan_tumor = smooth_curve_xy(x, pan_tumor_svg, clip_y=(-0.2, 1.0))
    x_pan_normal, y_pan_normal = smooth_curve_xy(x, pan_normal_svg, clip_y=(-0.2, 1.0))
    x_leukocyte, y_leukocyte = smooth_curve_xy(x, leukocyte_svg, clip_y=(-0.2, 1.0))

    if len(x_tumor) > 0:
        ax.plot(
            x_tumor,
            y_tumor,
            color="#129098",
            linewidth=3,
            label="Median Type T",
            solid_joinstyle="round",
            solid_capstyle="round",
        )
        ax.fill_between(
            x_tumor,
            np.zeros(len(x_tumor), dtype=float),
            np.nan_to_num(y_tumor, nan=0.0),
            color="#c3e2ea",
            alpha=1.0,
            linewidth=0,
        )

    if len(x_normal) > 0:
        ax.plot(
            x_normal,
            y_normal,
            color="#f0913e",
            linewidth=3,
            label="Median Type NT",
            solid_joinstyle="round",
            solid_capstyle="round",
        )

    if len(x_pan_tumor) > 0:
        ax.plot(
            x_pan_tumor,
            y_pan_tumor,
            color="#57ac3a",
            linewidth=3,
            label="Median PanCancer T",
            solid_joinstyle="round",
            solid_capstyle="round",
        )

    if len(x_pan_normal) > 0:
        ax.plot(
            x_pan_normal,
            y_pan_normal,
            color="#22672e",
            linewidth=3,
            label="Median PanCancer",
            solid_joinstyle="round",
            solid_capstyle="round",
        )

    if len(x_leukocyte) > 0:
        ax.plot(
            x_leukocyte,
            y_leukocyte,
            color="#a04589",
            linewidth=3,
            label="Median Leukocytes",
            solid_joinstyle="round",
            solid_capstyle="round",
        )

    first_region_label = True
    if filtered_regions is not None and not filtered_regions.empty:
        region_track = filtered_regions.copy()
        region_track["core_start"] = pd.to_numeric(region_track["core_start"], errors="coerce")
        region_track["core_end"] = pd.to_numeric(region_track["core_end"], errors="coerce")
        region_track = region_track.dropna(subset=["core_start", "core_end"])
        region_track = region_track[
            (region_track["core_end"] >= x_start)
            & (region_track["core_start"] <= x_end)
        ]

        for _, region_row in region_track.iterrows():
            r_start = int(region_row["core_start"])
            r_end = int(region_row["core_end"])
            label = "Candidate region" if first_region_label else None
            first_region_label = False

            ax.axvline(r_start, color="#dc143c", linestyle=":", linewidth=2.5)
            ax.axvline(r_end, color="#dc143c", linestyle=":", linewidth=2.5)
            ax.plot(
                [r_start, r_end],
                [-0.16, -0.16],
                color="#dc143c",
                linewidth=10,
                solid_capstyle="butt",
                label=label,
            )

    ax2 = ax.twinx()

    marker_data = pd.DataFrame(
        {
            "x": x,
            "delta": delta_svg,
            "hi": hi_svg,
        }
    ).replace([np.inf, -np.inf], np.nan)
    marker_data = marker_data[
        (marker_data["x"] >= x_start) & (marker_data["x"] <= x_end)
    ]

    ax2.scatter(
        marker_data["x"],
        marker_data["delta"],
        color="#000000",
        s=26,
        label="Delta",
        zorder=5,
    )

    if show_hi:
        ax2.scatter(
            marker_data["x"],
            marker_data["hi"],
            color="#912523",
            s=26,
            label="HI",
            zorder=5,
        )

    ax.set_xlim(x_start, x_end)
    ax.set_ylim(-0.2, 1)
    ax2.set_ylim(-0.2, max(1.0, float(y2_max) if pd.notna(y2_max) else 1.0))

    ax.set_title(
        f"{tumor_type} - {gene} Methylation gene profile",
        loc="left",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_xlabel("Genomic Position", fontsize=13)
    ax.set_ylabel("Methylation", fontsize=13)
    ax2.set_ylabel("Delta/HI", fontsize=13)

    ax.grid(True, linewidth=0.5, alpha=0.35)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))

    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(
        handles1 + handles2,
        labels1 + labels2,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        frameon=False,
        fontsize=10,
    )

    fig.tight_layout()

    buffer = io.BytesIO()
    fig.savefig(buffer, format="svg", bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return buffer.getvalue()



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

    .base-a { color: #000000; }
    .base-c { color: #ef4444; }
    .base-g { color: #2563eb; }
    .base-t { color: #000000; }
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

st.sidebar.header("Candidate region filters")

show_filtered_candidate_regions = st.sidebar.checkbox(
    "Show filtered candidate regions",
    value=True,
    help=(
        "Show only physical regions containing CpGs that pass the selected candidate filters. "
        "This does not require a precomputed curve table."
    ),
)

candidate_min_delta = st.sidebar.slider(
    "Candidate minimum Delta",
    0.0,
    1.0,
    0.50,
    0.01,
)

candidate_max_normal_median = st.sidebar.slider(
    "Candidate max Median NT",
    0.0,
    1.0,
    0.06,
    0.01,
)

candidate_max_pan_normal_median = st.sidebar.slider(
    "Candidate max PanCancer",
    0.0,
    1.0,
    0.08,
    0.01,
)

candidate_max_pan_tumor_median = st.sidebar.slider(
    "Candidate max PanCancer T",
    0.0,
    1.0,
    0.08,
    0.01,
)

candidate_max_leukocyte = st.sidebar.slider(
    "Candidate max Leukocytes",
    0.0,
    1.0,
    0.05,
    0.01,
)

candidate_min_hi = st.sidebar.slider(
    "Candidate min HI",
    0.0,
    5.0,
    1.60,
    0.05,
)

candidate_apply_expression_filter = st.sidebar.checkbox(
    "Apply candidate expression filter",
    value=False,
)

candidate_max_mean_spearman_r = st.sidebar.slider(
    "Candidate max mean Spearman r",
    -1.0,
    1.0,
    -0.06,
    0.01,
    disabled=not candidate_apply_expression_filter,
)


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


filtered_candidate_cpgs = pd.DataFrame()
filtered_candidate_regions = pd.DataFrame()
filtered_candidate_error = None

if show_filtered_candidate_regions:
    try:
        filtered_candidate_cpgs = load_filtered_candidate_region_cpgs_for_gene(
            gene=gene,
            tumor_type=tumor_type,
            min_delta=candidate_min_delta,
            max_normal_median=candidate_max_normal_median,
            max_pan_normal_median=candidate_max_pan_normal_median,
            max_pan_tumor_median=candidate_max_pan_tumor_median,
            max_leukocyte=candidate_max_leukocyte,
            min_hi=candidate_min_hi,
        )

        filtered_candidate_regions = aggregate_filtered_candidate_regions(
            candidate_cpgs=filtered_candidate_cpgs,
            apply_expression_filter=candidate_apply_expression_filter,
            max_mean_spearman_r=candidate_max_mean_spearman_r,
        )
    except Exception as exc:
        filtered_candidate_cpgs = pd.DataFrame()
        filtered_candidate_regions = pd.DataFrame()
        filtered_candidate_error = exc


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
expression = graph["spearman_r"]

hover_text = (
    graph["site_id"].astype(str)
    + " | "
    + graph["start_pos"].astype(str)
)

# ============================================================
# Manual Plotly zoom/export
# ============================================================

position_values = pd.to_numeric(positions, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()

if position_values.empty:
    st.warning("No valid genomic positions are available for plotting.")
    st.stop()

plot_x_min = int(position_values.min())
plot_x_max = int(position_values.max())

st.sidebar.header("Plot zoom / SVG export")

apply_manual_plot_zoom = st.sidebar.checkbox(
    "Apply manual genomic zoom",
    value=False,
    help=(
        "Enter genomic coordinates manually. The Plotly plot and the SVG "
        "download from the Plotly toolbar will use this same visible range."
    ),
)

zoom_start = st.sidebar.number_input(
    "Zoom start coordinate",
    min_value=plot_x_min,
    max_value=plot_x_max,
    value=plot_x_min,
    step=1,
    disabled=not apply_manual_plot_zoom,
)

zoom_end = st.sidebar.number_input(
    "Zoom end coordinate",
    min_value=plot_x_min,
    max_value=plot_x_max,
    value=plot_x_max,
    step=1,
    disabled=not apply_manual_plot_zoom,
)

if apply_manual_plot_zoom:
    zoom_start = int(zoom_start)
    zoom_end = int(zoom_end)

    if zoom_end <= zoom_start:
        st.sidebar.error("Zoom end coordinate must be greater than zoom start coordinate.")
        st.stop()

    plotly_export_suffix = f"{zoom_start}_{zoom_end}"
    plotly_xaxis_range = [zoom_start, zoom_end]
else:
    zoom_start = plot_x_min
    zoom_end = plot_x_max
    plotly_export_suffix = "full_gene"
    plotly_xaxis_range = None

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
            width=3,
            color="rgba(18,144,152,1)",
            shape="spline",
            smoothing=1.3,
        ),
        fill="tozeroy",
        fillcolor="rgba(195,226,234, 1)",
    )
)

curves = [
    ("Median Type NT", median_normal, "rgba(240,145,62,1)", "line", "y"),
    ("Median PanCan T", pan_tumor, "rgba(87,172,58,1)", "line", "y"),
    ("Median PanCan NT", pan_normal, "rgba(34,103,46,1)", "line", "y"),
    ("Median Leukocytes", leukocytes, "rgba(160,69,137,1)", "line", "y"),
    #("Expression", expression, "rgba(127,127,127,1)", "expression", "line"),
    ("Delta", delta_median, "rgba(0,0,0,1)", "marker", "y2"),
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
                    width=3,
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
                color="rgba(145,37,35,1)",
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

for _, row in filtered_candidate_regions.iterrows():
    region_id = str(row["region_id"])
    region_start = int(row["core_start"])
    region_end = int(row["core_end"])

    # Mark only the candidate-region borders.
    # No background band is used, so the methylation curves remain fully visible.
    for x_value, boundary_name in [
        (region_start, "Candidate region start"),
        (region_end, "Candidate region end"),
    ]:
        fig2.add_vline(
            x=x_value,
            line_width=2.5,
            line_dash="dot",
            line_color="rgba(220,20,60,0.95)",
            annotation_text="",
        )

    fig2.add_trace(
        go.Scatter(
            x=[region_start, region_end],
            y=[-0.16, -0.16],
            mode="lines",
            name="Candidate region",
            line=dict(
                color="rgba(220,20,60,1)",
                width=10,
            ),
            hovertext=[
                f"{region_id} | {region_start:,}-{region_end:,}",
                f"{region_id} | {region_start:,}-{region_end:,}",
            ],
            hovertemplate="%{hovertext}<extra></extra>",
            showlegend=False,
        )
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
        text=f"{tumor_type} - <b>{gene}</b> Methylation gene profile",
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
        range=plotly_xaxis_range,
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
        "filename": f"{gene}_{tumor_type}_{plotly_export_suffix}",
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

if apply_manual_plot_zoom:
    st.info(
        f"Manual Plotly zoom applied: chr region {zoom_start}-{zoom_end}. "
        "Use the camera/download button in the Plotly toolbar to export this exact visible range as SVG."
    )

# ============================================================
# Editable SVG export
# ============================================================

x_min_available = int(pd.to_numeric(positions, errors="coerce").min())
x_max_available = int(pd.to_numeric(positions, errors="coerce").max())

with st.expander("Editable SVG export for Illustrator"):
    st.caption(
        "Use this export when the default Plotly SVG is difficult to edit in Illustrator. "
        "The exported SVG uses smoothed curves and the same numeric zoom coordinates defined in the sidebar."
    )

    st.caption(
        f"Editable SVG range: {int(zoom_start)}-{int(zoom_end)}. "
        "This uses the same numeric coordinates defined in the sidebar."
    )

    editable_svg_x_range = (int(zoom_start), int(zoom_end))

    try:
        editable_curve_svg = make_editable_profile_svg(
            graph=graph,
            filtered_regions=filtered_candidate_regions,
            tumor_type=tumor_type,
            gene=gene,
            x_range=(int(editable_svg_x_range[0]), int(editable_svg_x_range[1])),
            show_hi=show_hi,
            y2_max=y2_max,
        )

        st.download_button(
            label="Download editable smooth curve SVG for Illustrator",
            data=editable_curve_svg,
            file_name=f"{gene}_{tumor_type}_methylation_profile_editable_smooth.svg",
            mime="image/svg+xml",
        )
    except Exception as exc:
        st.warning(f"Editable SVG export could not be generated: {exc}")


# ----------------------------
# Expression Graph
# ----------------------------
fig3 = go.Figure()

fig3.add_trace(
    go.Scatter(
        x=positions,
        y=median_tumor,
        mode="lines",
        name="Median Type T",
        customdata=site,
        line=dict(
            width=3,
            color="rgba(18,144,152,1)",
            shape="spline",
            smoothing=1.3,
        ),
        fill="tozeroy",
        fillcolor="rgba(195,226,234, 1)",
    )
)

curves = [
    ("Median Type NT", median_normal, "rgba(240,145,62,1)"),
    ("Median PanCan T", pan_tumor, "rgba(87,172,58,1)"),
    ("Median PanCan NT", pan_normal, "rgba(34,103,46,1)"),
    ("Median Leukocytes", leukocytes, "rgba(160,69,137,1)"),
    ("Expression", expression, "rgba(0,0,0,1)"),
]

for name, values, color in curves:
    fig3.add_trace(
        go.Scatter(
            x=positions,
            y=values,
            mode="lines",
            name=name,
            line=dict(
                color=color,
                width=3,
                shape="spline",
                smoothing=1.3,
            ),
        )
    )


fig3.update_traces(
    hovertext=hover_text,
    hovertemplate="%{hovertext}<br>%{fullData.name}: %{y:.3f}<extra></extra>",
)

fig3.update_layout(
    template="plotly_white",
    height=560,
    hovermode="x unified",
    margin=dict(l=40, r=40, t=50, b=50),
    legend=dict(
        orientation="v",
        y=1,
        x=1.02,
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
        text=f"{tumor_type} - <b>{gene}</b> Methylation gene profile expression",
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
        range=plotly_xaxis_range,
        rangeslider=dict(visible=True),
        showgrid=True,
    ),
    yaxis=dict(
        # Upper methylation space fixed at 1.
        # Lower space free enough for Spearman r.
        range=[-0.6, 1],
        title=dict(
            text="Methylation / Expression Spearman r",
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
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor="rgba(0,0,0,0.45)",
    ),
)

config = {
    "toImageButtonOptions": {
        "format": "svg",
        "filename": f"{gene}_{tumor_type}_methylation_expression_{plotly_export_suffix}",
        "height": 650,
        "width": 1400,
        "scale": 1,
    }
}

st.plotly_chart(
    fig3,
    use_container_width=True,
    config=config,
)

# ============================================================
# Filtered candidate region tables
# ============================================================

selected_candidate_region_for_browser = None

st.subheader("Filtered candidate regions")
st.caption(
    "These regions are calculated live from the selected filters. "
    "No biomarker_candidate_region_curve table or parquet is required."
)

if not show_filtered_candidate_regions:
    st.info("Enable 'Show filtered candidate regions' in the sidebar to calculate this table.")
elif filtered_candidate_error is not None:
    st.warning("Filtered candidate regions could not be calculated.")
elif filtered_candidate_regions.empty:
    st.info("No filtered candidate regions passed the current filters.")
else:
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Candidate regions", f"{filtered_candidate_regions['region_id'].nunique():,}")
    with c2:
        st.metric("Candidate CpGs", f"{filtered_candidate_cpgs['site_id'].nunique():,}")
    with c3:
        st.metric("Max tumor median", f"{filtered_candidate_regions['max_tumor_median'].max():.3f}")
    with c4:
        st.metric("Max biological score", f"{filtered_candidate_regions['max_biological_score'].max():.2f}")

    region_table_cols = [
        "region_id",
        "region_gene_symbol",
        "chr",
        "core_start",
        "core_end",
        "core_length",
        "browser_start",
        "browser_end",
        "browser_length",
        "sequence_length",
        "n_qualifying_cpgs",
        "n_manifest_cpgs",
        "fraction_qualifying_cpgs",
        "n_manifest_c",
        "cpg_density_per_100bp",
        "n_c_sequence",
        "n_g_sequence",
        "n_cg_sequence",
        "n_gcgc",
        "gc_fraction",
        "sequence_score",
        "mean_delta",
        "max_delta",
        "mean_hi",
        "max_hi",
        "mean_tumor_median",
        "mean_normal_median",
        "mean_pan_tumor_median",
        "mean_pan_normal_median",
        "mean_leukocyte_median",
        "mean_spearman_r",
        "min_spearman_r",
        "qualifying_cpg_sites",
        "qualifying_cpg_positions",
    ]
    region_table_cols = [col for col in region_table_cols if col in filtered_candidate_regions.columns]

    st.markdown("**Region-level table**")
    st.caption(
        "Select one candidate region here to open it in the sequence browser below."
    )

    candidate_region_table = (
        filtered_candidate_regions[region_table_cols]
        .reset_index(drop=True)
        .copy()
    )

    candidate_region_event = st.dataframe(
        candidate_region_table,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
    )

    selected_candidate_region_idx = extract_selected_table_row(candidate_region_event)

    if selected_candidate_region_idx is not None:
        selected_candidate_region_for_browser = candidate_region_table.iloc[
            selected_candidate_region_idx
        ].copy()

        st.session_state["selected_region_id"] = str(
            selected_candidate_region_for_browser.get("region_id", "")
        )
        st.session_state["selected_physical_region_id"] = str(
            selected_candidate_region_for_browser.get("region_id", "")
        )
        st.session_state["selected_gene_symbol"] = normalize_gene_symbol(
            selected_candidate_region_for_browser.get("region_gene_symbol", gene)
        )
        st.session_state["selected_tumor_type"] = tumor_type
        st.session_state["selected_region_chr"] = str(
            selected_candidate_region_for_browser.get("chr", "")
        ).replace("chr", "").replace("CHR", "").strip()
        st.session_state["selected_region_start"] = int(
            selected_candidate_region_for_browser.get("browser_start", 0)
        )
        st.session_state["selected_region_end"] = int(
            selected_candidate_region_for_browser.get("browser_end", 0)
        )
        st.session_state["selected_browser_start"] = int(
            selected_candidate_region_for_browser.get("browser_start", 0)
        )
        st.session_state["selected_browser_end"] = int(
            selected_candidate_region_for_browser.get("browser_end", 0)
        )
        st.session_state["selected_region_source"] = "Region-level table"

        st.success(
            "Selected candidate region will be opened in the sequence browser below: "
            f"{selected_candidate_region_for_browser.get('region_id', '')}"
        )

    st.download_button(
        label="Download filtered candidate regions CSV",
        data=filtered_candidate_regions[region_table_cols].to_csv(index=False).encode("utf-8"),
        file_name=f"{gene}_{tumor_type}_filtered_candidate_regions.csv",
        mime="text/csv",
    )

    cpg_table_cols = [
        "region_id",
        "region_gene_symbol",
        "region_genes_all",
        "cpg_gene_symbol",
        "chr",
        "core_start",
        "core_end",
        "browser_start",
        "browser_end",
        "site_id",
        "start_pos",
        "cpg_order",
        "tumor_type",
        "delta_median",
        "tumor_median",
        "normal_median",
        "pan_tumor_median",
        "pan_normal_median",
        "hi_index",
        "leukocyte_median",
        "spearman_r",
        "biological_score",
        "delta_score",
        "normal_low_score",
        "leukocyte_low_score",
        "hi_score",
        "expression_score",
        "passes_loose_seed",
        "passes_default_filter",
        "passes_strict_filter",
        "n_manifest_cpgs",
        "n_manifest_c",
        "cpg_density_per_100bp",
        "sequence_score",
        "sequence_available",
        "sequence_length",
        "n_c_sequence",
        "n_g_sequence",
        "n_cg_sequence",
        "n_gcgc",
        "gc_fraction",
    ]
    cpg_table_cols = [col for col in cpg_table_cols if col in filtered_candidate_cpgs.columns]

    st.markdown("**CpG-level table inside filtered candidate regions**")
    st.dataframe(
        filtered_candidate_cpgs[cpg_table_cols],
        use_container_width=True,
        hide_index=True,
    )

    st.download_button(
        label="Download filtered candidate CpGs CSV",
        data=filtered_candidate_cpgs[cpg_table_cols].to_csv(index=False).encode("utf-8"),
        file_name=f"{gene}_{tumor_type}_filtered_candidate_cpgs.csv",
        mime="text/csv",
    )

# ============================================================
# Region browser helpers
# ============================================================

@st.cache_data(ttl=600, show_spinner=False)
def load_gene_regions(gene: str) -> pd.DataFrame:
    """
    Load physical regions for the selected gene from the normalized region model.
    This reads biomarker_region and optional sequence scores from PostgreSQL.
    """
    query = """
    SELECT
        r.region_id,
        r.gene_symbol,
        r.chr,
        r.core_start,
        r.core_end,
        r.core_length,
        r.browser_start,
        r.browser_end,
        r.browser_length,
        r.flank_bp,
        r.n_manifest_cpgs,
        r.n_manifest_c,
        r.cpg_density_per_100bp,
        r.region_rank_by_density,
        COALESCE(seq.sequence_score, 0) AS sequence_score,
        seq.sequence_available,
        seq.sequence_length,
        seq.n_cg_sequence,
        seq.n_gcgc,
        seq.gc_fraction
    FROM biomarker_region r
    LEFT JOIN biomarker_region_sequence_score seq
        ON r.region_id = seq.region_id
    WHERE r.gene_symbol = :gene
    ORDER BY r.chr, r.core_start, r.core_end
    """
    df = run_query(query, params={"gene": normalize_gene_symbol(gene)})

    if df.empty:
        return df

    numeric_cols = [
        "core_start", "core_end", "core_length", "browser_start", "browser_end",
        "browser_length", "flank_bp", "n_manifest_cpgs", "n_manifest_c",
        "cpg_density_per_100bp", "region_rank_by_density", "sequence_score",
        "sequence_length", "n_cg_sequence", "n_gcgc", "gc_fraction",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["gene_symbol"] = df["gene_symbol"].map(normalize_gene_symbol)
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", case=False, regex=False).str.strip()
    df = df.dropna(subset=["region_id", "gene_symbol", "chr", "core_start", "core_end"]).copy()
    return df.reset_index(drop=True)


@st.cache_data(ttl=600, show_spinner=False)
def load_region_cpgs(region_id: str, tumor_type: str) -> pd.DataFrame:
    """
    Load CpGs belonging to one region with tumor-specific summaries and CpG scores.
    """
    query = """
    SELECT DISTINCT
        rc.region_id,
        rc.site_id,
        rc.gene_symbol,
        rc.chr,
        rc.start_pos,
        rc.cpg_order,
        ts.tumor_median,
        ts.normal_median,
        ts.delta_median,
        ts.hi_index,
        ts.pan_tumor_median,
        ts.pan_normal_median,
        cf.leukocyte_median,
        ec.spearman_r,
        bcs.biological_score,
        bcs.delta_score,
        bcs.normal_low_score,
        bcs.leukocyte_low_score,
        bcs.hi_score,
        bcs.expression_score,
        bcs.passes_loose_seed,
        bcs.passes_default_filter,
        bcs.passes_strict_filter
    FROM biomarker_region_cpg rc
    LEFT JOIN tumor_summary ts
        ON rc.site_id = ts.site_id
        AND ts.tumor_type = :tumor_type
    LEFT JOIN cpg_features cf
        ON rc.site_id = cf.site_id
    LEFT JOIN expression_correlation ec
        ON rc.site_id = ec.site_id
        AND ec.tumor_type = :tumor_type
        AND (
            ec.gene_symbol = rc.gene_symbol
            OR ec.gene_symbol IS NULL
            OR ec.gene_symbol = ''
        )
    LEFT JOIN biomarker_cpg_score bcs
        ON rc.site_id = bcs.site_id
        AND bcs.tumor_type = :tumor_type
        AND bcs.gene_symbol = rc.gene_symbol
    WHERE rc.region_id = :region_id
    ORDER BY rc.cpg_order, rc.start_pos
    """
    df = run_query(
        query,
        params={
            "region_id": str(region_id),
            "tumor_type": tumor_type,
        },
    )

    if df.empty:
        return df

    numeric_cols = [
        "start_pos", "cpg_order", "tumor_median", "normal_median", "delta_median",
        "hi_index", "pan_tumor_median", "pan_normal_median", "leukocyte_median",
        "spearman_r", "biological_score", "delta_score", "normal_low_score",
        "leukocyte_low_score", "hi_score", "expression_score",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["site_id", "start_pos"]).copy()
    df["start_pos"] = df["start_pos"].astype(int)
    return df.reset_index(drop=True)


@st.cache_data(ttl=3600, show_spinner=False)
def cached_get_sequence(chrom: str, start: int, end: int) -> str:
    """Cached wrapper around the app Ensembl sequence service."""
    chrom = str(chrom).replace("chr", "").replace("CHR", "").strip()
    return str(get_sequence(chrom, int(start), int(end))).upper()


def count_overlapping_motif(sequence: str, motif: str) -> int:
    """Count overlapping motif occurrences."""
    sequence = str(sequence).upper()
    motif = str(motif).upper()
    return len(re.findall(f"(?={re.escape(motif)})", sequence))


def render_region_sequence_matrix(
    sequence: str,
    seq_start: int,
    manifest_cpg_positions: set[int],
    core_start: int,
    core_end: int,
    gcgc_roles: dict[int, list[str]],
    bases_per_row: int = BASES_PER_ROW,
) -> str:
    """Render the region browser using the selectable sequence layout from Region Browser."""
    sequence = sequence.upper()

    html_parts = [
        '<div class="region-sequence-card">',
        '<div class="sequence-title">Sequence — Region View</div>',
    ]

    for row_start in range(0, len(sequence), bases_per_row):
        row_seq = sequence[row_start:row_start + bases_per_row]
        row_genomic_start = int(seq_start) + row_start

        html_parts.append('<div class="region-sequence-row">')
        html_parts.append(f'<span class="region-coord" aria-hidden="true">{row_genomic_start:,}</span>')
        html_parts.append('<span class="region-sequence-bases">')

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

            if int(core_start) <= genomic_pos <= int(core_end):
                classes.append("region-core")

            if genomic_pos in gcgc_roles:
                classes.extend(gcgc_roles[genomic_pos])

            if genomic_pos in manifest_cpg_positions:
                classes.append("manifest-cpg")

            class_str = " ".join(classes)
            html_parts.append(
                f'<span class="{class_str}" '
                f'title="pos {genomic_pos:,} | base {safe_base}">'
                f'{safe_base}</span>'
            )

        html_parts.append("</span>")
        html_parts.append("</div>")

    html_parts.append(
        """
        <div class="legend">
            <span class="legend-item"><span class="dot dot-green"></span>CpG manifest</span>
            <span class="legend-item"><span class="box-black"></span>GCGC motif</span>
            <span class="legend-item"><span class="line-orange"></span>Core region</span>
        </div>
        """
    )
    html_parts.append("</div>")
    return "".join(html_parts)


def render_region_sequence_svg(
    sequence: str,
    seq_start: int,
    manifest_cpg_positions: set[int],
    core_start: int,
    core_end: int,
    gcgc_roles: dict[int, list[str]],
    title: str,
    subtitle: str = "",
    bases_per_row: int = BASES_PER_ROW,
) -> str:
    """Render the region browser as SVG for download."""
    sequence = sequence.upper()

    cell_w = 18
    cell_h = 22
    row_h = 28
    coord_w = 115
    margin_x = 24
    margin_y = 70
    title_h = 42
    legend_h = 42

    n_rows = int(np.ceil(len(sequence) / bases_per_row))
    width = margin_x * 2 + coord_w + bases_per_row * cell_w + 20
    height = margin_y + title_h + n_rows * row_h + legend_h + 24

    base_colors = {
        "A": "#000000",
        "C": "#ef4444",
        "G": "#2563eb",
        "T": "#000000",
    }

    svg_font = "'Arial Narrow', Arial, sans-serif"
    svg_mono_font = "'Arial Narrow', Arial, sans-serif"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin_x}" y="30" font-family="{svg_font}" font-size="18" font-weight="700" fill="#102a43">{html.escape(title)}</text>',
        f'<text x="{margin_x}" y="52" font-family="{svg_font}" font-size="13" fill="#52616b">{html.escape(subtitle)}</text>',
    ]

    y0 = margin_y + title_h

    for row_start in range(0, len(sequence), bases_per_row):
        row_seq = sequence[row_start:row_start + bases_per_row]
        row_idx = row_start // bases_per_row
        y = y0 + row_idx * row_h
        row_genomic_start = int(seq_start) + row_start

        parts.append(
            f'<text x="{margin_x + coord_w - 8}" y="{y + 16}" '
            f'text-anchor="end" font-family="{svg_font}" font-size="12" '
            f'fill="#9aa5b1">{row_genomic_start:,}</text>'
        )

        for i, base in enumerate(row_seq):
            genomic_pos = row_genomic_start + i
            x = margin_x + coord_w + i * cell_w
            safe_base = html.escape(base)

            fill = "transparent"
            stroke = "transparent"
            stroke_w = 1
            text_fill = base_colors.get(base, "#6b7280")
            font_weight = 700

            if int(core_start) <= genomic_pos <= int(core_end):
                fill = "#fff7ed"
                stroke = "#fb923c"
                stroke_w = 1

            if genomic_pos in manifest_cpg_positions:
                fill = "#bbf7d0"
                stroke = "#22c55e"
                stroke_w = 2
                text_fill = "#000000"
                font_weight = 900

            if genomic_pos in gcgc_roles:
                fill = "#000000"
                stroke = "#000000"
                stroke_w = 2
                text_fill = "#ffffff"
                font_weight = 900

            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h}" '
                f'rx="3" ry="3" fill="{fill}" stroke="{stroke}" '
                f'stroke-width="{stroke_w}"/>'
            )
            parts.append(
                f'<text x="{x + (cell_w - 2) / 2}" y="{y + 16}" '
                f'text-anchor="middle" font-family="{svg_mono_font}" '
                f'font-size="14" font-weight="{font_weight}" '
                f'fill="{text_fill}">{safe_base}</text>'
            )

    legend_y = y0 + n_rows * row_h + 26
    parts.extend([
        f'<circle cx="{margin_x + 8}" cy="{legend_y - 4}" r="6" fill="#22c55e"/>',
        f'<text x="{margin_x + 22}" y="{legend_y}" font-family="{svg_font}" font-size="13" fill="#52616b">CpG manifest</text>',
        f'<rect x="{margin_x + 210}" y="{legend_y - 13}" width="16" height="12" fill="#000000"/>',
        f'<text x="{margin_x + 232}" y="{legend_y}" font-family="{svg_font}" font-size="13" fill="#52616b">GCGC motif</text>',
        f'<rect x="{margin_x + 370}" y="{legend_y - 13}" width="16" height="12" fill="#fff7ed" stroke="#fb923c"/>',
        f'<text x="{margin_x + 392}" y="{legend_y}" font-family="{svg_font}" font-size="13" fill="#52616b">Core region</text>',
        '</svg>',
    ])
    return "".join(parts)


st.markdown(
    """
    <style>
    .region-sequence-card {
        background-color: white;
        border: 1px solid #d9e2ec;
        border-radius: 8px;
        padding: 16px;
        overflow-x: auto;
        margin-top: 12px;
    }

    .region-sequence-row {
        display: flex;
        align-items: center;
        white-space: nowrap;
        line-height: 1.9;
    }

    .region-coord {
        flex: 0 0 105px;
        color: #9aa5b1;
        font-family: monospace;
        font-size: 13px;
        text-align: right;
        margin-right: 12px;
        user-select: none;
        -webkit-user-select: none;
        -moz-user-select: none;
        -ms-user-select: none;
        pointer-events: none;
    }

    .region-sequence-bases {
        display: inline-block;
        user-select: text;
        -webkit-user-select: text;
    }

    .region-core {
        background-color: #fff7ed;
        border-bottom: 2px solid #fb923c;
    }

    .line-orange {
        display: inline-block;
        width: 18px;
        height: 10px;
        border-bottom: 3px solid #fb923c;
        background: #fff7ed;
        margin-right: 5px;
        vertical-align: middle;
    }

    .box-black {
        background: black;
    }

    /* GCGC motifs must override core and manifest-CpG styling. */
    .sequence-card .gcgc-site,
    .region-sequence-card .gcgc-site {
        color: #ffffff !important;
        background-color: #000000 !important;
        border-color: #000000 !important;
        border-top: 2px solid #000000 !important;
        border-bottom: 2px solid #000000 !important;
        font-weight: 900 !important;
    }

    .sequence-card .gcgc-left,
    .region-sequence-card .gcgc-left {
        border-left: 2px solid #000000 !important;
        border-top-left-radius: 4px !important;
        border-bottom-left-radius: 4px !important;
    }

    .sequence-card .gcgc-right,
    .region-sequence-card .gcgc-right {
        border-right: 2px solid #000000 !important;
        border-top-right-radius: 4px !important;
        border-bottom-right-radius: 4px !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# Browser controls
# ============================================================

st.subheader("Browser")
st.caption(
    "Select one row in the Region-level table above to open that candidate region here. "
    "The Physical CpG regions table was removed to keep the browser focused on filtered candidate regions."
)

if show_filtered_candidate_regions and not filtered_candidate_regions.empty:
    st.info(
        "Filtered candidate regions are shown as a blue track inside the "
        "Methylation gene profile plot above."
    )

sequence_view = st.radio(
    "Sequence window",
    ["browser", "core", "custom flank"],
    horizontal=True,
    index=0,
    help=(
        "browser uses browser_start/browser_end; core uses core_start/core_end; "
        "custom flank adds flanks around the core."
    ),
)

custom_flank_bp = 0
if sequence_view == "custom flank":
    custom_flank_bp = st.slider(
        "Custom flank around core region (bp)",
        min_value=0,
        max_value=2000,
        value=500,
        step=100,
    )


# ============================================================
# Region sequence browser
# ============================================================

if selected_candidate_region_for_browser is not None:
    selected_region = selected_candidate_region_for_browser.copy()

    if "gene_symbol" not in selected_region.index:
        selected_region["gene_symbol"] = selected_region.get("region_gene_symbol", gene)

    browser_selection_source = "Region-level table"
else:
    st.info("No candidate region selected. Choose one row from the Region-level table above to see the sequence.")
    st.stop()

region_id = str(selected_region["region_id"])
region_gene = normalize_gene_symbol(selected_region.get("gene_symbol", selected_region.get("region_gene_symbol", gene)))
chrom = str(selected_region["chr"]).replace("chr", "").replace("CHR", "").strip()
core_start = int(selected_region["core_start"])
core_end = int(selected_region["core_end"])

if sequence_view == "browser":
    seq_start = int(selected_region["browser_start"])
    seq_end = int(selected_region["browser_end"])
elif sequence_view == "core":
    seq_start = core_start
    seq_end = core_end
else:
    seq_start = max(1, core_start - int(custom_flank_bp))
    seq_end = core_end + int(custom_flank_bp)

seq_start = max(1, int(seq_start))
seq_end = int(seq_end)

if seq_end <= seq_start:
    st.error("Invalid sequence coordinates for the selected region.")
    st.stop()

region_cpgs = load_region_cpgs(region_id=region_id, tumor_type=tumor_type)

if region_cpgs.empty:
    st.warning("No CpGs were found for this region in biomarker_region_cpg.")
    manifest_positions_raw: list[int] = []
else:
    manifest_positions_raw = region_cpgs["start_pos"].dropna().astype(int).tolist()

try:
    sequence = cached_get_sequence(chrom=chrom, start=seq_start, end=seq_end).upper()
except Exception as exc:
    st.error(f"The sequence could not be obtained from Ensembl: {exc}")
    st.stop()

if not sequence:
    st.warning("No sequence was returned for this region.")
    st.stop()

highlight_positions: list[int] = []
not_found_positions: list[int] = []

for pos in manifest_positions_raw:
    if seq_start <= int(pos) <= seq_end:
        c_pos = normalize_cpg_c_position(raw_pos=int(pos), sequence=sequence, seq_start=seq_start)
        if c_pos is None:
            not_found_positions.append(int(pos))
        else:
            highlight_positions.append(int(c_pos))

manifest_cpg_positions = set(highlight_positions)
gcgc_roles = build_gcgc_roles(sequence, seq_start)

n_c = sequence.count("C")
n_g = sequence.count("G")
n_cg = count_overlapping_motif(sequence, "CG")
n_gcgc = count_overlapping_motif(sequence, "GCGC")
gc_fraction = (n_c + n_g) / len(sequence) if sequence else np.nan

col_a, col_b, col_c, col_d = st.columns(4)
with col_a:
    st.metric("Region ID", region_id)
with col_b:
    st.metric("Gene", region_gene)
with col_c:
    st.metric("Length", f"{len(sequence):,} bp")
with col_d:
    st.metric("GCGC motifs", f"{n_gcgc:,}")

col_e, col_f, col_g, col_h = st.columns(4)
with col_e:
    st.metric("Chromosome", f"chr{chrom}")
with col_f:
    st.metric("Core region", f"{core_start}-{core_end}")
with col_g:
    st.metric("Displayed region", f"{seq_start}-{seq_end}")
with col_h:
    st.metric("Manifest CpGs", f"{len(manifest_cpg_positions)}/{len(manifest_positions_raw)}")

metadata = (
    f"Region: {region_id} | gene: {region_gene} | chr{chrom}:{seq_start}-{seq_end} | "
    f"core: {core_start}-{core_end} | manifest CpGs: {len(manifest_cpg_positions)}/{len(manifest_positions_raw)} | "
    f"CG: {n_cg:,} | GCGC: {n_gcgc:,} | GC fraction: {gc_fraction:.3f}"
)
st.info(f"Browser source: {browser_selection_source} | {metadata}")

if not_found_positions:
    st.warning(
        "Some manifest positions could not be normalized to a CpG C inside "
        f"the displayed sequence: {not_found_positions[:10]}"
        + ("..." if len(not_found_positions) > 10 else "")
    )

# CpG table intentionally removed. CpGs are still used internally to mark
# manifest CpG positions in the sequence browser.

sequence_html = render_region_sequence_matrix(
    sequence=sequence,
    seq_start=seq_start,
    manifest_cpg_positions=manifest_cpg_positions,
    core_start=core_start,
    core_end=core_end,
    gcgc_roles=gcgc_roles,
    bases_per_row=BASES_PER_ROW,
)

st.markdown(sequence_html, unsafe_allow_html=True)
st.caption("Tip: select only the sequence letters. The left genomic coordinates are separated and should not be copied.")

sequence_svg = render_region_sequence_svg(
    sequence=sequence,
    seq_start=seq_start,
    manifest_cpg_positions=manifest_cpg_positions,
    core_start=core_start,
    core_end=core_end,
    gcgc_roles=gcgc_roles,
    title=f"{region_id} | {region_gene} | chr{chrom}:{seq_start}-{seq_end}",
    subtitle=metadata,
    bases_per_row=BASES_PER_ROW,
)

st.download_button(
    label="Download sequence browser SVG",
    data=sequence_svg.encode("utf-8"),
    file_name=f"{region_gene}_{region_id}_chr{chrom}_{seq_start}_{seq_end}.svg",
    mime="image/svg+xml",
)
