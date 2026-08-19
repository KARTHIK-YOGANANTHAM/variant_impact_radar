"""
Variant Impact Radar — Workbench Edition
----------------------------------------
A multi-page Streamlit application for rapid computational interpretation of
dbSNP / COSMIC / HGMD-style variant identifiers using the Ensembl VEP REST API.

Run:
    pip install -r requirements.txt
    streamlit run app.py

Core data source:
    Ensembl REST VEP API
    https://rest.ensembl.org/documentation/info/vep_id_get

Important:
    The "Impact Radar Score" is a transparent heuristic for demonstration and
    research communication. It is NOT a clinical pathogenicity classifier and
    must not be used for diagnosis or patient management.
"""

from __future__ import annotations

import base64
import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import plotly.graph_objects as go
import requests
import streamlit as st


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

APP_TITLE = "Variant Impact Radar"
APP_VERSION = "2.0.0"

ENSEMBL_BASE = "https://rest.ensembl.org"
VEP_ENDPOINT = f"{ENSEMBL_BASE}/vep/human/id"
LOOKUP_SYMBOL_ENDPOINT = f"{ENSEMBL_BASE}/lookup/symbol/homo_sapiens"
ENSEMBL_DOCS_URL = "https://rest.ensembl.org/documentation/info/vep_id_get"

REQUEST_TIMEOUT = 30

EXAMPLE_VARIANTS = [
    "rs7903146",
    "rs1801133",
    "rs1801131",
    "rs1800795",
]

# ---------------------------------------------------------------------------
# Developer profile — rendered on the "About the developer" page.
#
# EDIT ME: every value below is a placeholder. Swap in your real name, bio,
# links, stats, and timeline. The portrait is loaded from the first existing
# file in DEV_PHOTO_CANDIDATES — drop your photo at assets/developer.png
# (next to assets/logo.svg) and it appears automatically; until then a
# neutral "blank profile" silhouette is rendered inline.
# ---------------------------------------------------------------------------

DEV_PHOTO_CANDIDATES = (
    "developer.png",
    "developer.jpg",
    "developer.jpeg",
    "developer.webp",
)

DEVELOPER_PROFILE = {
    "name": "Karthik Yoganantham",
    "role": "Bioinformatics Analyst · Computational Biologist",
    "location": "Kanchipuram, India",
    "availability": "Open to collaborations",
    "tagline": (
        "Turning dense genomic annotation into interfaces people actually "
        "want to read."
    ),
    "bio": (
        "I am a bioinformatics engineer and data analyst working at "
        "the intersection of genomics and software: REST-API-driven analysis "
        "tools, reproducible pipelines, and dashboards that keep every number "
        "traceable back to its source. Variant Impact Radar grew out of that "
        "habit — a small, transparent workbench that shows exactly how a raw "
        "Ensembl VEP response becomes an evidence summary. Away from the "
        "keyboard I am usually reading population-genetics papers, sketching "
        "interface ideas, or fine-tuning espresso ratios."
    ),
    "email": "karthikyoganantham@gmail.com",
    "linkedin": "www.linkedin.com/in/karthik-yoganantham",
    "github": "https://github.com/KARTHIK-YOGANANTHAM",
    "portfolio": "https://karthikyoganantham.framer.website/",
    "skills": [
        "Python", "Streamlit", "Plotly", "pandas", "REST APIs",
        "Ensembl VEP", "R / Bioconductor", "SQL", "Docker", "Git & CI",
        "Linux", "NGS pipelines", "Data visualization", "HTML / CSS", "ML", "Scikit-Learn", "RNN"
    ],
    # (label, value, sub-caption) — rendered as metric cards.
    "stats": [
        ("Experience", "1+ years", "Bioinformatics & ML Analyst"),
        ("Projects shipped", "12", "Tools, pipelines & dashboards"),
        ("Open-source repos", "8", "Python · R · JavaScript"),
        ("Publications", "1", "Co-authored, peer-reviewed"),
    ],
    # (year, title, description) — rendered as a vertical timeline.
    "timeline": [
        (
            "2026", "Variant Impact Radar v2",
            "Rebuilt the workbench as a multi-page Streamlit application "
            "with a transparent scoring heuristic and shareable reports.",
        ),
        (
            "2025", "Gene Express analysis of Diabetic Complications",
            "Gene Express analysis and visualization of Diabetic Retinopathy and its related complications using R, "
            "A comprehensive Study",
        ),
        (
            "2025", "MSc Bioinformatics",
            "Thesis on regulatory-variant prioritization using public "
            "annotation resources.",
        ),
        (
            "2024", "Multiple Disease Prediction Project",
            "Machine Learning based Project on Multiple Disease Prediction "
            "using Patient MetaData",
        ),
    ],
}

# Transparent heuristic weights.
CONSEQUENCE_WEIGHTS = {
    "transcript_ablation": 10,
    "splice_acceptor_variant": 9,
    "splice_donor_variant": 9,
    "stop_gained": 9,
    "frameshift_variant": 9,
    "start_lost": 8,
    "stop_lost": 8,
    "missense_variant": 6,
    "inframe_insertion": 5,
    "inframe_deletion": 5,
    "protein_altering_variant": 5,
    "regulatory_region_variant": 3,
    "TF_binding_site_variant": 3,
    "promoter_variant": 3,
    "5_prime_UTR_variant": 2,
    "3_prime_UTR_variant": 2,
    "non_coding_transcript_exon_variant": 2,
    "intron_variant": 1,
    "synonymous_variant": 1,
    "intergenic_variant": 0,
}

MAX_SCORE = 10.0

# ---------------------------------------------------------------------------
# Design tokens — "gel under UV": warm near-black, amber luminescence,
# chromatogram base colours, serif display type.
# ---------------------------------------------------------------------------

INK = "#14110C"
SURFACE = "#1B1712"
SURFACE_2 = "#211C15"
LINE = "#2C261C"
PORCELAIN = "#EDE7DD"
MUTED = "#A2988A"
AMBER = "#E8A33D"
AMBER_DIM = "rgba(232,163,61,0.14)"

TIER_COLORS = {
    "HIGH": "#E2593B",
    "MODERATE": "#E8A33D",
    "LOW–MODERATE": "#B4AC5C",
    "LOW": "#7E8B79",
}

BASE_COLORS = {  # Sanger chromatogram convention, tuned for dark ground
    "A": "#7FB069",
    "C": "#6FA8DC",
    "G": "#E8A33D",
    "T": "#E2593B",
}

SERIES_COLORS = ["#E8A33D", "#7FB069", "#6FA8DC", "#E2593B", "#B08BC9", "#8B8577"]

POPULATION_LABELS = {
    "af": "1000G · global",
    "afr": "1000G · African",
    "amr": "1000G · American",
    "eas": "1000G · East Asian",
    "eur": "1000G · European",
    "sas": "1000G · South Asian",
    "gnomade": "gnomAD exomes · all",
    "gnomadg": "gnomAD genomes · all",
    "aa": "ESP · African American",
    "ea": "ESP · European American",
}


# ---------------------------------------------------------------------------
# Page config + CSS
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=":material/biotech:",
    layout="wide",
    initial_sidebar_state="expanded",
)


def inject_css() -> None:
    st.markdown(
        """
<style>
@import url('https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght@0,9..144,400;0,9..144,500;0,9..144,600;1,9..144,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500;600&display=swap');

:root {
    --ink: #14110C;
    --surface: #1B1712;
    --surface-2: #211C15;
    --line: #2C261C;
    --porcelain: #EDE7DD;
    --muted: #A2988A;
    --amber: #E8A33D;
    --amber-dim: rgba(232,163,61,0.14);
    --serif: "Fraunces", Georgia, serif;
    --sans: "IBM Plex Sans", -apple-system, sans-serif;
    --mono: "IBM Plex Mono", ui-monospace, monospace;
}

.stApp { font-family: var(--sans); }
div.block-container { max-width: 1180px; padding-top: 1.4rem; padding-bottom: 3.5rem; }

/* Quiet the Streamlit chrome */
#MainMenu, footer { visibility: hidden; }
div[data-testid="stToolbar"], div[data-testid="stDecoration"] { display: none; }
header[data-testid="stHeader"] { background: transparent; }

/* Sidebar console */
section[data-testid="stSidebar"] {
    background: #171310;
    border-right: 1px solid var(--line);
}
section[data-testid="stSidebar"] hr { border-color: var(--line); }

/* Inputs */
div[data-testid="stTextInput"] input {
    font-family: var(--mono);
    letter-spacing: .02em;
    background: var(--surface);
    border: 1px solid var(--line);
    border-radius: 8px;
}
div[data-testid="stTextInput"] input:focus { border-color: var(--amber); }

/* Buttons */
.stButton button, .stDownloadButton button, .stLinkButton a,
div[data-testid="stFormSubmitButton"] button {
    border-radius: 8px;
    border: 1px solid var(--line);
    font-weight: 600;
}
div[data-testid="stFormSubmitButton"] button[kind="primary"],
.stButton button[kind="primary"] {
    background: var(--amber);
    color: #1A150D;
    border-color: var(--amber);
}
div[data-testid="stFormSubmitButton"] button[kind="primary"]:hover {
    background: #F0B254;
    border-color: #F0B254;
    color: #1A150D;
}

/* Tabs read like instrument channel labels */
button[data-baseweb="tab"] p {
    font-family: var(--mono);
    font-size: .78rem;
    text-transform: uppercase;
    letter-spacing: .10em;
}

/* Expanders */
div[data-testid="stExpander"] details {
    border: 1px solid var(--line);
    border-radius: 10px;
    background: var(--surface);
}

/* ---------------- Custom components ---------------- */

.vr-eyebrow {
    font-family: var(--mono);
    font-size: .72rem;
    letter-spacing: .16em;
    text-transform: uppercase;
    color: var(--amber);
}

.vr-section { margin: 2.1rem 0 0.9rem; }
.vr-section .vr-eyebrow { color: var(--muted); }
.vr-section-title {
    font-family: var(--serif);
    font-size: 1.5rem;
    font-weight: 500;
    color: var(--porcelain);
    margin: .1rem 0 .15rem;
    letter-spacing: .01em;
}
.vr-section-desc { color: var(--muted); font-size: .9rem; margin: 0; }
.vr-section-rule {
    height: 1px;
    background: linear-gradient(90deg, var(--amber) 0, var(--amber) 42px, var(--line) 42px);
    margin-top: .65rem;
}

/* Hero: the locus strip. Ruler ticks along the bottom edge echo a genome
   browser coordinate track. */
.vr-hero {
    position: relative;
    background: linear-gradient(180deg, var(--surface-2) 0%, var(--surface) 100%);
    border: 1px solid var(--line);
    border-radius: 14px;
    padding: 1.6rem 1.8rem 2.0rem;
    overflow: hidden;
}
.vr-hero::after {
    content: "";
    position: absolute;
    left: 0; right: 0; bottom: 0;
    height: 14px;
    background:
        repeating-linear-gradient(90deg,
            rgba(232,163,61,.55) 0 1px, transparent 1px 12px),
        repeating-linear-gradient(90deg,
            rgba(232,163,61,.28) 0 1px, transparent 1px 60px);
    background-size: auto 7px, auto 14px;
    background-repeat: repeat-x;
    background-position: bottom left;
    opacity: .8;
}
.vr-hero-top {
    display: flex; justify-content: space-between; align-items: center;
    flex-wrap: wrap; gap: .6rem;
}
.vr-hero-id {
    font-family: var(--serif);
    font-size: clamp(2.1rem, 4.4vw, 3.3rem);
    font-weight: 500;
    color: var(--porcelain);
    line-height: 1.08;
    margin: .5rem 0 .8rem;
}
.vr-hero-locus {
    display: flex; align-items: center; flex-wrap: wrap; gap: .9rem;
    font-family: var(--mono);
    font-size: 1.02rem;
    color: var(--porcelain);
}
.vr-hero-meta {
    margin-top: .9rem;
    font-family: var(--mono);
    font-size: .74rem;
    letter-spacing: .06em;
    color: var(--muted);
}
.vr-gene-chip {
    font-family: var(--mono);
    font-size: .82rem;
    padding: .18rem .6rem;
    border-radius: 999px;
    background: var(--amber-dim);
    border: 1px solid rgba(232,163,61,.35);
    color: var(--amber);
}
.vr-base { font-weight: 600; }
.vr-allele-sep { color: var(--muted); padding: 0 .1rem; }

.vr-tier {
    display: inline-flex; align-items: center; gap: .45rem;
    font-family: var(--mono);
    font-size: .74rem;
    letter-spacing: .12em;
    text-transform: uppercase;
    padding: .3rem .7rem;
    border-radius: 999px;
    border: 1px solid;
}
.vr-tier-dot {
    width: 7px; height: 7px; border-radius: 50%;
    box-shadow: 0 0 8px 1px currentColor;
}
@media (prefers-reduced-motion: no-preference) {
    .vr-tier-dot { animation: vr-glow 2.6s ease-in-out infinite; }
    @keyframes vr-glow { 50% { box-shadow: 0 0 3px 0 currentColor; } }
}

/* Metric cards */
.vr-card {
    background: var(--surface);
    border: 1px solid var(--line);
    border-radius: 12px;
    padding: 1rem 1.1rem;
    height: 100%;
}
.vr-card-label {
    font-family: var(--mono);
    color: var(--muted);
    font-size: .70rem;
    text-transform: uppercase;
    letter-spacing: .12em;
    margin-bottom: .45rem;
}
.vr-card-value {
    font-size: 1.32rem;
    font-weight: 600;
    color: var(--porcelain);
    line-height: 1.25;
    overflow-wrap: anywhere;
}
.vr-card-value.mono { font-family: var(--mono); font-size: 1.08rem; }
.vr-card-sub { color: var(--muted); font-size: .78rem; margin-top: .3rem; }

/* Biological flow pipeline — a true sequence, so the numbering is earned */
.vr-flow { display: flex; align-items: stretch; flex-wrap: wrap; gap: .5rem; }
.vr-flow-step {
    flex: 1 1 150px;
    background: var(--surface);
    border: 1px solid var(--line);
    border-radius: 12px;
    padding: .85rem .95rem;
    position: relative;
}
.vr-flow-num {
    font-family: var(--mono);
    font-size: .68rem;
    letter-spacing: .14em;
    color: var(--amber);
    margin-bottom: .35rem;
}
.vr-flow-label {
    font-family: var(--mono);
    font-size: .68rem;
    text-transform: uppercase;
    letter-spacing: .12em;
    color: var(--muted);
    margin-bottom: .3rem;
}
.vr-flow-value {
    font-size: .92rem; font-weight: 600; color: var(--porcelain);
    overflow-wrap: anywhere;
}
.vr-flow-arrow {
    align-self: center;
    color: var(--muted);
    font-family: var(--mono);
    padding: 0 .1rem;
}

/* Narrative interpretation — the one serif prose moment */
.vr-story {
    font-family: var(--serif);
    font-size: 1.04rem;
    line-height: 1.8;
    color: var(--porcelain);
    background: var(--surface);
    border: 1px solid var(--line);
    border-left: 3px solid var(--amber);
    border-radius: 0 12px 12px 0;
    padding: 1.3rem 1.5rem;
}

/* Small chips (clinical significance etc.) */
.vr-chip {
    display: inline-block;
    font-family: var(--mono);
    font-size: .76rem;
    padding: .22rem .6rem;
    border-radius: 999px;
    border: 1px solid var(--line);
    background: var(--surface-2);
    color: var(--porcelain);
    margin: .12rem .18rem .12rem 0;
}
.vr-chip a { color: inherit; text-decoration: none; }

/* Empty state */
.vr-empty {
    border: 1px dashed var(--line);
    border-radius: 14px;
    padding: 2.6rem 2rem;
    text-align: center;
    background: rgba(255,255,255,0.012);
}
.vr-empty-title {
    font-family: var(--serif);
    font-size: 1.35rem;
    color: var(--porcelain);
    margin: .55rem 0 .4rem;
}
.vr-empty-body { color: var(--muted); font-size: .92rem; max-width: 480px; margin: 0 auto; }

/* Sidebar console labels */
.vr-side-label {
    font-family: var(--mono);
    font-size: .68rem;
    text-transform: uppercase;
    letter-spacing: .14em;
    color: var(--muted);
    margin: .9rem 0 .3rem;
}
.vr-side-note {
    font-size: .78rem;
    color: var(--muted);
    border-top: 1px solid var(--line);
    padding-top: .8rem;
    margin-top: 1rem;
    line-height: 1.55;
}

.vr-footer {
    margin-top: 3rem;
    padding-top: .9rem;
    border-top: 1px solid var(--line);
    font-family: var(--mono);
    font-size: .72rem;
    letter-spacing: .05em;
    color: var(--muted);
}

/* ------------- About-the-developer page ------------- */

.vr-dev-hero {
    position: relative;
    display: flex;
    align-items: center;
    flex-wrap: wrap;
    gap: 1.7rem;
    background: linear-gradient(180deg, var(--surface-2) 0%, var(--surface) 100%);
    border: 1px solid var(--line);
    border-radius: 14px;
    padding: 1.8rem 1.9rem 2.2rem;
    overflow: hidden;
}
.vr-dev-hero::after {
    content: "";
    position: absolute;
    left: 0; right: 0; bottom: 0;
    height: 14px;
    background:
        repeating-linear-gradient(90deg,
            rgba(232,163,61,.55) 0 1px, transparent 1px 12px),
        repeating-linear-gradient(90deg,
            rgba(232,163,61,.28) 0 1px, transparent 1px 60px);
    background-size: auto 7px, auto 14px;
    background-repeat: repeat-x;
    background-position: bottom left;
    opacity: .8;
}
.vr-dev-photo {
    width: 152px;
    height: 152px;
    border-radius: 50%;
    object-fit: cover;
    border: 2px solid var(--amber);
    box-shadow: 0 0 0 6px var(--amber-dim), 0 0 26px rgba(232,163,61,.22);
    background: var(--surface-2);
    flex: 0 0 auto;
}
.vr-dev-hero-body { flex: 1 1 320px; min-width: 260px; }
.vr-dev-name {
    font-family: var(--serif);
    font-size: clamp(1.9rem, 3.6vw, 2.7rem);
    font-weight: 500;
    color: var(--porcelain);
    line-height: 1.12;
    margin: .25rem 0 .35rem;
}
.vr-dev-loc {
    font-family: var(--mono);
    font-size: .78rem;
    letter-spacing: .06em;
    color: var(--muted);
    margin-bottom: .7rem;
}
.vr-dev-tagline {
    font-family: var(--serif);
    font-style: italic;
    font-size: 1.04rem;
    color: var(--porcelain);
    opacity: .92;
    max-width: 560px;
    line-height: 1.55;
}
.vr-dev-timeline {
    border-left: 1px solid var(--line);
    margin: .5rem 0 .2rem .4rem;
    padding-left: 1.25rem;
}
.vr-dev-tl-item { position: relative; padding-bottom: 1.2rem; }
.vr-dev-tl-item:last-child { padding-bottom: .2rem; }
.vr-dev-tl-item::before {
    content: "";
    position: absolute;
    left: calc(-1.25rem - 4.5px);
    top: .32rem;
    width: 8px; height: 8px;
    border-radius: 50%;
    background: var(--amber);
    box-shadow: 0 0 8px 1px rgba(232,163,61,.55);
}
.vr-dev-tl-year {
    font-family: var(--mono);
    font-size: .70rem;
    letter-spacing: .12em;
    text-transform: uppercase;
    color: var(--amber);
    margin-bottom: .15rem;
}
.vr-dev-tl-title { font-weight: 600; color: var(--porcelain); font-size: .95rem; }
.vr-dev-tl-desc {
    color: var(--muted);
    font-size: .84rem;
    margin-top: .18rem;
    line-height: 1.55;
    max-width: 640px;
}
</style>
        """,
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalize_variant_id(value: str) -> str:
    value = (value or "").strip()
    value = re.sub(r"\s+", "", value)
    return value


def valid_variant_id(value: str) -> bool:
    # Common dbSNP/COSMIC/HGMD-style identifiers. We intentionally allow
    # broader alphanumeric IDs rather than rejecting valid future IDs.
    return bool(re.fullmatch(r"[A-Za-z0-9_.:\-]+", value))


def safe_get(mapping: Any, *keys: str, default: Any = None) -> Any:
    if not isinstance(mapping, dict):
        return default
    for key in keys:
        if key in mapping:
            return mapping[key]
    return default


def flatten(value: Any) -> List[Any]:
    """Flatten nested lists/tuples/sets while keeping scalars."""
    if isinstance(value, (list, tuple, set)):
        out: List[Any] = []
        for item in value:
            out.extend(flatten(item))
        return out
    return [value]


def unique_strings(values: Iterable[Any]) -> List[str]:
    seen = set()
    result = []
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if not s or s.lower() in {"none", "nan"}:
            continue
        if s not in seen:
            seen.add(s)
            result.append(s)
    return result


def format_number(value: Any, decimals: int = 3) -> str:
    try:
        n = float(value)
        if pd.isna(n):
            return "—"
        return f"{n:.{decimals}f}"
    except Exception:
        return "—"


def evidence_level(score: float) -> Tuple[str, str]:
    if score >= 8:
        return "HIGH", "Potentially high-impact computational evidence"
    if score >= 5:
        return "MODERATE", "Moderate computational impact"
    if score >= 2.5:
        return "LOW–MODERATE", "Limited computational impact"
    return "LOW", "Low computational impact"


def score_badge(score: float) -> str:
    label, _ = evidence_level(score)
    return label


def esc(value: Any) -> str:
    return html.escape(str(value))


def humanize(term: str) -> str:
    return str(term).replace("_", " ")


# ---------------------------------------------------------------------------
# Ensembl API
# ---------------------------------------------------------------------------

@st.cache_data(ttl=900, show_spinner=False)
def fetch_vep(variant_id: str) -> Dict[str, Any]:
    """
    Fetch one variant annotation from Ensembl VEP.

    Ensembl VEP supports identifier-based queries and optional annotations
    including CADD and AlphaMissense.
    """
    params = {
        "content-type": "application/json",
        "CADD": "1",
        "AlphaMissense": "1",
    }
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": f"VariantImpactRadar/{APP_VERSION} (Streamlit demo)",
    }

    response = requests.get(
        f"{VEP_ENDPOINT}/{variant_id}",
        params=params,
        headers=headers,
        timeout=REQUEST_TIMEOUT,
    )

    if response.status_code == 429:
        raise RuntimeError(
            "Ensembl rate limit reached (HTTP 429). Please retry shortly."
        )

    if response.status_code == 404:
        raise RuntimeError(
            f"Variant '{variant_id}' was not found by the Ensembl VEP endpoint."
        )

    if not response.ok:
        message = response.text[:500].strip()
        raise RuntimeError(
            f"Ensembl VEP returned HTTP {response.status_code}. {message}"
        )

    try:
        payload = response.json()
    except ValueError as exc:
        raise RuntimeError("Ensembl returned a non-JSON response.") from exc

    if isinstance(payload, list):
        if not payload:
            raise RuntimeError("Ensembl returned no annotation for this variant.")
        # Usually one object for one identifier.
        return payload[0]

    if isinstance(payload, dict):
        return payload

    raise RuntimeError("Unexpected response format from Ensembl VEP.")


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_gene_summary(gene_symbol: str) -> Dict[str, Any]:
    """Optional gene-level context using Ensembl symbol lookup."""
    if not gene_symbol:
        return {}

    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": f"VariantImpactRadar/{APP_VERSION} (Streamlit demo)",
    }
    params = {
        "content-type": "application/json",
    }

    try:
        response = requests.get(
            f"{LOOKUP_SYMBOL_ENDPOINT}/{gene_symbol}",
            headers=headers,
            params=params,
            timeout=REQUEST_TIMEOUT,
        )
        if not response.ok:
            return {}
        return response.json()
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Parsing VEP output
# ---------------------------------------------------------------------------

def get_transcript_consequences(record: Dict[str, Any]) -> List[Dict[str, Any]]:
    tc = record.get("transcript_consequences")
    if isinstance(tc, list):
        return [x for x in tc if isinstance(x, dict)]
    return []


def get_regulatory_consequences(record: Dict[str, Any]) -> List[Dict[str, Any]]:
    rc = record.get("regulatory_feature_consequences")
    if isinstance(rc, list):
        return [x for x in rc if isinstance(x, dict)]
    return []


def get_all_consequence_terms(record: Dict[str, Any]) -> List[str]:
    values: List[Any] = []

    for key in ("most_severe_consequence",):
        if record.get(key):
            values.append(record[key])

    for item in get_transcript_consequences(record):
        values.append(item.get("consequence_terms", []))

    for item in get_regulatory_consequences(record):
        values.append(item.get("consequence_terms", []))

    return unique_strings(flatten(values))


def choose_primary_transcript(
    record: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    transcripts = get_transcript_consequences(record)
    if not transcripts:
        return None

    # Priority: canonical > MANE > most detailed consequence > first item.
    def rank(item: Dict[str, Any]) -> Tuple[int, int, int]:
        is_canonical = int(bool(item.get("canonical")))
        mane = int(bool(item.get("mane_select") or item.get("mane_plus_clinical")))
        consequence_count = len(item.get("consequence_terms") or [])
        return (is_canonical, mane, consequence_count)

    return sorted(transcripts, key=rank, reverse=True)[0]


def choose_gene_symbol(record: Dict[str, Any], tx: Optional[Dict[str, Any]]) -> str:
    candidates = []
    if tx:
        candidates += [
            tx.get("gene_symbol"),
            tx.get("gene_id"),
        ]
    candidates += [
        record.get("gene_symbol"),
        record.get("gene_symbol_source"),
        record.get("gene_id"),
    ]
    values = unique_strings(candidates)
    return values[0] if values else "Unknown"


def extract_scores(tx: Optional[Dict[str, Any]]) -> Dict[str, Optional[float]]:
    """
    CADD and AlphaMissense fields can appear in plugin/annotation structures
    depending on VEP configuration. This parser intentionally checks several
    likely representations.
    """
    cadd: Optional[float] = None
    alphamissense: Optional[float] = None

    if not tx:
        return {"cadd": cadd, "alphamissense": alphamissense}

    # CADD may be represented as a scalar in some VEP deployments.
    direct_cadd = tx.get("CADD_PHRED", tx.get("CADD"))
    if isinstance(direct_cadd, (int, float)):
        cadd = float(direct_cadd)

    # AlphaMissense may be scalar or nested.
    direct_am = tx.get("am_pathogenicity")
    if direct_am is None:
        direct_am = tx.get("AlphaMissense")
    if isinstance(direct_am, (int, float)):
        alphamissense = float(direct_am)

    # Try nested "colocated_variants" / plugin-style structures.
    for container in (
        tx.get("colocated_variants"),
        tx.get("transcript_consequences"),
        tx.get("plugins"),
    ):
        if isinstance(container, list):
            for item in container:
                if not isinstance(item, dict):
                    continue
                for key in ("CADD_PHRED", "CADD"):
                    v = item.get(key)
                    if cadd is None and isinstance(v, (int, float)):
                        cadd = float(v)
                for key in ("am_pathogenicity", "AlphaMissense"):
                    v = item.get(key)
                    if alphamissense is None and isinstance(v, (int, float)):
                        alphamissense = float(v)

    return {"cadd": cadd, "alphamissense": alphamissense}


def derive_consequence_summary(record: Dict[str, Any], tx: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    terms = []
    if tx:
        terms = unique_strings(tx.get("consequence_terms", []))
    if not terms:
        terms = get_all_consequence_terms(record)

    most_severe = record.get("most_severe_consequence")
    if not most_severe and terms:
        most_severe = max(
            terms,
            key=lambda t: CONSEQUENCE_WEIGHTS.get(t, 0),
        )

    return {
        "terms": terms,
        "primary": most_severe or "Unknown",
        "protein_change": (
            tx.get("amino_acids")
            or tx.get("protein_end")
            or tx.get("protein_id")
            or ""
        ) if tx else "",
        "codons": tx.get("codons") if tx else None,
    }


def extract_population_frequencies(record: Dict[str, Any]) -> pd.DataFrame:
    """
    Population allele frequencies from colocated dbSNP records (when the VEP
    response includes them). Returned as a tidy dataframe.
    """
    rows: List[Dict[str, Any]] = []
    for cv in record.get("colocated_variants") or []:
        if not isinstance(cv, dict):
            continue
        freqs = cv.get("frequencies")
        if not isinstance(freqs, dict):
            continue
        for allele, pops in freqs.items():
            if not isinstance(pops, dict):
                continue
            for pop_code, value in pops.items():
                try:
                    freq = float(value)
                except (TypeError, ValueError):
                    continue
                label = POPULATION_LABELS.get(pop_code)
                if label is None:
                    code = str(pop_code)
                    if code.startswith("gnomade_"):
                        label = f"gnomAD exomes · {code.split('_', 1)[1].upper()}"
                    elif code.startswith("gnomadg_"):
                        label = f"gnomAD genomes · {code.split('_', 1)[1].upper()}"
                    else:
                        label = code.upper()
                rows.append(
                    {
                        "Population": label,
                        "Allele": str(allele),
                        "Frequency": freq,
                    }
                )
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates().sort_values(
            ["Allele", "Frequency"], ascending=[True, False]
        )
    return df


def extract_clinical_context(record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Clinical-significance terms, PubMed references, and phenotype names found
    in the colocated-variant annotations of the VEP response.
    """
    clin_sig: List[Any] = []
    pubmed: List[Any] = []
    minor_allele = None
    minor_allele_freq = None

    for cv in record.get("colocated_variants") or []:
        if not isinstance(cv, dict):
            continue
        clin_sig.extend(flatten(cv.get("clin_sig") or []))
        pubmed.extend(flatten(cv.get("pubmed") or []))
        if minor_allele is None and cv.get("minor_allele"):
            minor_allele = cv.get("minor_allele")
        if minor_allele_freq is None and cv.get("minor_allele_freq") is not None:
            minor_allele_freq = cv.get("minor_allele_freq")

    phenotypes: List[Any] = []
    for item in record.get("phenotypes") or []:
        if isinstance(item, dict):
            phenotypes.append(item.get("phenotype") or item.get("trait"))
        else:
            phenotypes.append(item)

    return {
        "clin_sig": unique_strings(clin_sig),
        "pubmed": unique_strings(pubmed),
        "phenotypes": unique_strings(phenotypes),
        "minor_allele": minor_allele,
        "minor_allele_freq": minor_allele_freq,
    }


# ---------------------------------------------------------------------------
# Transparent scoring
# ---------------------------------------------------------------------------

def calculate_impact_score(
    record: Dict[str, Any],
    tx: Optional[Dict[str, Any]],
    gene_summary: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Demonstration score:
      - consequence class contributes most
      - CADD / AlphaMissense add modest evidence
      - canonical/MANE transcript adds reproducibility context
      - existing clinical/phenotype annotations add context
    This is deliberately NOT a clinical pathogenicity score.
    """
    consequence = derive_consequence_summary(record, tx)
    primary = consequence["primary"]

    score = float(CONSEQUENCE_WEIGHTS.get(primary, 0))
    reasons = [
        f"Primary consequence: {primary} (+{CONSEQUENCE_WEIGHTS.get(primary, 0):g})"
    ]

    scores = extract_scores(tx)
    cadd = scores["cadd"]
    am = scores["alphamissense"]

    if cadd is not None:
        if cadd >= 30:
            score += 2.0
            reasons.append("CADD-Phred ≥30 (+2)")
        elif cadd >= 20:
            score += 1.0
            reasons.append("CADD-Phred ≥20 (+1)")
        elif cadd >= 10:
            score += 0.5
            reasons.append("CADD-Phred ≥10 (+0.5)")

    if am is not None:
        if am >= 0.9:
            score += 1.5
            reasons.append("AlphaMissense ≥0.90 (+1.5)")
        elif am >= 0.7:
            score += 1.0
            reasons.append("AlphaMissense ≥0.70 (+1)")
        elif am >= 0.5:
            score += 0.5
            reasons.append("AlphaMissense ≥0.50 (+0.5)")

    if tx:
        if tx.get("canonical"):
            score += 0.5
            reasons.append("Canonical transcript (+0.5)")
        if tx.get("mane_select") or tx.get("mane_plus_clinical"):
            score += 0.5
            reasons.append("MANE transcript annotation (+0.5)")

    # Annotation/context bonuses should be small, because the score is not a
    # clinical classifier.
    if record.get("colocated_variants"):
        score += 0.5
        reasons.append("Known colocated-variant context (+0.5)")

    if record.get("phenotypes"):
        score += 0.5
        reasons.append("Phenotype annotation present (+0.5)")

    if gene_summary:
        score += 0.25
        reasons.append("Gene-level Ensembl context available (+0.25)")

    score = max(0.0, min(score, MAX_SCORE))

    return {
        "score": round(score, 2),
        "label": score_badge(score),
        "reasons": reasons,
        "cadd": cadd,
        "alphamissense": am,
    }


# ---------------------------------------------------------------------------
# Interpretation
# ---------------------------------------------------------------------------

def biological_interpretation(
    record: Dict[str, Any],
    tx: Optional[Dict[str, Any]],
    gene_symbol: str,
    consequence: Dict[str, Any],
    score: Dict[str, Any],
    gene_summary: Dict[str, Any],
) -> str:
    primary = consequence["primary"].replace("_", " ")
    chrom = record.get("seq_region_name", "unknown")
    start = record.get("start", "unknown")
    allele = record.get("allele_string", "unknown")
    assembly = record.get("assembly_name", "GRCh38")

    tx_id = tx.get("transcript_id", "selected transcript") if tx else "available transcript"
    protein = ""
    if tx:
        aa = tx.get("amino_acids")
        codons = tx.get("codons")
        if aa:
            protein = f" The transcript-level annotation reports the amino-acid change '{aa}'."
        elif codons:
            protein = f" The reported codon context is '{codons}'."

    severity_phrase = {
        "transcript ablation": "a potentially severe loss-of-transcript event",
        "splice acceptor variant": "a variant at a splice-acceptor boundary",
        "splice donor variant": "a variant at a splice-donor boundary",
        "stop gained": "a premature stop-gain consequence",
        "frameshift variant": "a reading-frame-disrupting consequence",
        "start lost": "a start-codon-disrupting consequence",
        "stop lost": "a stop-codon-disrupting consequence",
        "missense variant": "a protein-altering missense consequence",
        "inframe insertion": "an in-frame insertion",
        "inframe deletion": "an in-frame deletion",
        "protein altering variant": "a protein-altering consequence",
        "regulatory region variant": "a regulatory-region consequence",
        "tf binding site variant": "a transcription-factor-binding-site consequence",
        "promoter variant": "a promoter-associated consequence",
        "5 prime utr variant": "a 5' UTR consequence",
        "3 prime utr variant": "a 3' UTR consequence",
        "intron variant": "an intronic consequence",
        "synonymous variant": "a synonymous consequence",
        "intergenic variant": "an intergenic consequence",
    }.get(primary, f"a {primary} consequence")

    cadd_text = (
        f"CADD-Phred is approximately {score['cadd']:.2f}."
        if score.get("cadd") is not None
        else "No CADD value was available in the returned annotation."
    )
    am_text = (
        f"AlphaMissense is approximately {score['alphamissense']:.3f}."
        if score.get("alphamissense") is not None
        else "No AlphaMissense value was available for the selected annotation."
    )

    gene_description = safe_get(
        gene_summary,
        "description",
        default=None,
    )
    gene_context = (
        f" Ensembl describes {gene_symbol} as: {gene_description}."
        if gene_description
        else ""
    )

    return (
        f"This analysis places {record.get('id', 'the queried variant')} at "
        f"{chrom}:{start} on {assembly}, with allele representation "
        f"{allele}. The selected transcript is {tx_id}, where the primary "
        f"annotated consequence is {severity_phrase}.{protein} "
        f"{cadd_text} {am_text}{gene_context} "
        f"The current heuristic impact score is {score['score']:.2f}/10 "
        f"({score['label']}). This score is an evidence-aggregation demo and "
        f"should not be interpreted as a clinical pathogenicity verdict."
    )


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def make_transcript_dataframe(record: Dict[str, Any]) -> pd.DataFrame:
    rows = []
    for tx in get_transcript_consequences(record):
        terms = ", ".join(unique_strings(tx.get("consequence_terms", [])))
        rows.append(
            {
                "Transcript": tx.get("transcript_id", "—"),
                "Gene": tx.get("gene_symbol", tx.get("gene_id", "—")),
                "Consequence": terms or "—",
                "Canonical": "Yes" if tx.get("canonical") else "No",
                "MANE": (
                    "Yes"
                    if tx.get("mane_select") or tx.get("mane_plus_clinical")
                    else "No"
                ),
                "Amino acids": tx.get("amino_acids", "—"),
                "Codons": tx.get("codons", "—"),
            }
        )
    return pd.DataFrame(rows)


def make_evidence_dataframe(
    consequence: Dict[str, Any],
    score: Dict[str, Any],
) -> pd.DataFrame:
    rows = [
        {
            "Evidence": "Primary consequence",
            "Value": consequence["primary"],
            "Contribution": CONSEQUENCE_WEIGHTS.get(consequence["primary"], 0),
        },
        {
            "Evidence": "CADD-Phred",
            "Value": (
                score["cadd"] if score["cadd"] is not None else "Not available"
            ),
            "Contribution": "see thresholds",
        },
        {
            "Evidence": "AlphaMissense",
            "Value": (
                score["alphamissense"]
                if score["alphamissense"] is not None
                else "Not available"
            ),
            "Contribution": "see thresholds",
        },
    ]
    if score["reasons"]:
        # Add all non-primary reasons without trying to duplicate thresholds.
        for reason in score["reasons"][1:]:
            rows.append(
                {
                    "Evidence": "Additional context",
                    "Value": reason,
                    "Contribution": "",
                }
            )
    df = pd.DataFrame(rows)
    # Mixed numeric/text cells are display-only; cast for clean Arrow tables.
    df["Value"] = df["Value"].astype(str)
    df["Contribution"] = df["Contribution"].astype(str)
    return df


# ---------------------------------------------------------------------------
# Charts — themed to the gel-glow palette
# ---------------------------------------------------------------------------

def _theme_figure(fig: go.Figure, height: int) -> go.Figure:
    fig.update_layout(
        height=height,
        margin=dict(l=28, r=28, t=30, b=24),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="IBM Plex Mono, monospace", size=12, color=MUTED),
    )
    return fig


def radar_axis_values(
    score: float, cadd: Optional[float], am: Optional[float]
) -> List[float]:
    """Normalize external evidence to comparable 0–10 display values only."""
    cadd_norm = min((cadd / 30.0) * 10.0, 10.0) if cadd is not None else 0.0
    am_norm = min(max(am, 0.0) * 10.0, 10.0) if am is not None else 0.0
    consequence_norm = min(score, 10.0)
    return [score, cadd_norm, am_norm, consequence_norm]


RADAR_CATEGORIES = ["Overall", "CADD", "AlphaMissense", "Consequence"]


def make_radar_chart(score: float, cadd: Optional[float], am: Optional[float]) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatterpolar(
            r=radar_axis_values(score, cadd, am),
            theta=RADAR_CATEGORIES,
            fill="toself",
            name="Evidence profile",
            line=dict(color=AMBER, width=2),
            fillcolor="rgba(232,163,61,0.22)",
        )
    )
    fig.update_layout(
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            radialaxis=dict(
                visible=True,
                range=[0, 10],
                gridcolor="rgba(237,231,221,0.10)",
                linecolor="rgba(237,231,221,0.14)",
                tickfont=dict(size=10),
            ),
            angularaxis=dict(
                gridcolor="rgba(237,231,221,0.10)",
                linecolor="rgba(237,231,221,0.14)",
            ),
        ),
        showlegend=False,
    )
    return _theme_figure(fig, 340)


def make_score_gauge(score: float) -> go.Figure:
    fig = go.Figure(
        go.Indicator(
            mode="gauge+number",
            value=score,
            number={
                "suffix": " / 10",
                "font": {
                    "family": "Fraunces, Georgia, serif",
                    "size": 46,
                    "color": PORCELAIN,
                },
            },
            gauge={
                "axis": {
                    "range": [0, 10],
                    "tickwidth": 1,
                    "tickcolor": LINE,
                    "tickfont": {"size": 10, "color": MUTED},
                },
                "bar": {"color": AMBER, "thickness": 0.30},
                "bgcolor": "rgba(0,0,0,0)",
                "borderwidth": 0,
                "steps": [
                    {"range": [0, 2.5], "color": "rgba(232,163,61,0.05)"},
                    {"range": [2.5, 5], "color": "rgba(232,163,61,0.10)"},
                    {"range": [5, 7.5], "color": "rgba(232,163,61,0.17)"},
                    {"range": [7.5, 10], "color": "rgba(232,163,61,0.26)"},
                ],
            },
        )
    )
    return _theme_figure(fig, 300)


def make_frequency_chart(freq_df: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    alleles = sorted(freq_df["Allele"].unique())
    for i, allele in enumerate(alleles):
        sub = freq_df[freq_df["Allele"] == allele].sort_values("Frequency")
        fig.add_trace(
            go.Bar(
                y=sub["Population"],
                x=sub["Frequency"],
                orientation="h",
                name=str(allele),
                marker_color=BASE_COLORS.get(
                    str(allele), SERIES_COLORS[i % len(SERIES_COLORS)]
                ),
                marker_line_width=0,
            )
        )
    n_pops = freq_df["Population"].nunique()
    fig.update_layout(
        barmode="group",
        xaxis=dict(
            tickformat=".1%",
            gridcolor="rgba(237,231,221,0.08)",
            zerolinecolor="rgba(237,231,221,0.14)",
        ),
        yaxis=dict(gridcolor="rgba(0,0,0,0)"),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            title_text="Allele",
        ),
    )
    return _theme_figure(fig, max(280, 26 * n_pops + 120))


def make_compare_bar(rows: List[Dict[str, Any]]) -> go.Figure:
    fig = go.Figure(
        go.Bar(
            y=[r["Variant"] for r in rows],
            x=[r["Impact score"] for r in rows],
            orientation="h",
            marker_color=[
                TIER_COLORS.get(r["Evidence tier"], AMBER) for r in rows
            ],
            marker_line_width=0,
            text=[f'{r["Impact score"]:.2f}' for r in rows],
            textposition="outside",
            textfont=dict(family="IBM Plex Mono, monospace"),
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, 11],
            gridcolor="rgba(237,231,221,0.08)",
            title="Impact score / 10",
        ),
        yaxis=dict(autorange="reversed"),
        showlegend=False,
    )
    return _theme_figure(fig, max(240, 64 * len(rows) + 110))


def make_compare_radar(items: List[Dict[str, Any]]) -> go.Figure:
    fig = go.Figure()
    for i, item in enumerate(items):
        color = SERIES_COLORS[i % len(SERIES_COLORS)]
        fig.add_trace(
            go.Scatterpolar(
                r=radar_axis_values(
                    item["score"]["score"],
                    item["score"].get("cadd"),
                    item["score"].get("alphamissense"),
                ),
                theta=RADAR_CATEGORIES,
                fill="toself",
                name=item["variant_id"],
                line=dict(color=color, width=2),
                opacity=0.85,
            )
        )
    fig.update_layout(
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            radialaxis=dict(
                visible=True,
                range=[0, 10],
                gridcolor="rgba(237,231,221,0.10)",
                linecolor="rgba(237,231,221,0.14)",
                tickfont=dict(size=10),
            ),
            angularaxis=dict(
                gridcolor="rgba(237,231,221,0.10)",
                linecolor="rgba(237,231,221,0.14)",
            ),
        ),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=-0.18),
    )
    return _theme_figure(fig, 420)


# ---------------------------------------------------------------------------
# HTML report (downloadable)
# ---------------------------------------------------------------------------

def build_html_report(
    variant_id: str,
    record: Dict[str, Any],
    gene_symbol: str,
    consequence: Dict[str, Any],
    score: Dict[str, Any],
    interpretation: str,
    transcript_df: pd.DataFrame,
) -> str:
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    transcript_html = (
        transcript_df.to_html(
            index=False,
            escape=True,
            classes="tbl",
            border=0,
        )
        if not transcript_df.empty
        else "<p>No transcript consequences returned.</p>"
    )

    raw_json = html.escape(json.dumps(record, indent=2, ensure_ascii=False))

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Variant Impact Radar — {html.escape(variant_id)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
body {{
    font-family: "IBM Plex Sans", Inter, Arial, sans-serif;
    margin: 0;
    background: #FAF7F0;
    color: #1C1712;
}}
.container {{
    max-width: 1100px;
    margin: 30px auto;
    padding: 0 20px 40px;
}}
.hero {{
    background: white;
    border-radius: 16px;
    padding: 28px;
    border: 1px solid #E8E2D4;
    border-top: 3px solid #B8791C;
}}
h1 {{ margin: 0 0 6px; font-family: Georgia, "Times New Roman", serif; font-weight: 600; }}
h2 {{ margin-top: 30px; font-family: Georgia, "Times New Roman", serif; font-weight: 600; }}
.muted {{ color: #7A7264; }}
.mono {{ font-family: "IBM Plex Mono", ui-monospace, monospace; }}
.grid {{
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 14px;
    margin-top: 20px;
}}
.card {{
    background: white;
    border: 1px solid #E8E2D4;
    border-radius: 12px;
    padding: 16px;
}}
.label {{
    font-family: "IBM Plex Mono", ui-monospace, monospace;
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: .10em;
    color: #7A7264;
}}
.value {{ font-size: 24px; font-weight: 700; margin-top: 6px; }}
.story {{
    background: white;
    border: 1px solid #E8E2D4;
    border-left: 3px solid #B8791C;
    border-radius: 8px;
    padding: 18px;
    line-height: 1.7;
}}
.tbl {{
    width: 100%;
    border-collapse: collapse;
    background: white;
}}
.tbl th, .tbl td {{
    border-bottom: 1px solid #EFEAE0;
    padding: 9px;
    text-align: left;
    font-size: 13px;
}}
pre {{
    background: #1B1712;
    color: #EDE7DD;
    padding: 16px;
    border-radius: 12px;
    overflow: auto;
    font-size: 12px;
}}
.note {{
    margin-top: 25px;
    padding: 14px;
    background: #FBF3E1;
    border: 1px solid #EBD8AC;
    border-radius: 10px;
}}
@media (max-width: 800px) {{
    .grid {{ grid-template-columns: repeat(2, 1fr); }}
}}
</style>
</head>
<body>
<div class="container">
<div class="hero">
    <h1>Variant Impact Radar</h1>
    <div class="muted">Computational variant evidence report · {html.escape(generated)}</div>
    <h2 class="mono">{html.escape(variant_id)}</h2>
</div>

<div class="grid">
    <div class="card"><div class="label">Gene</div><div class="value">{html.escape(gene_symbol)}</div></div>
    <div class="card"><div class="label">Consequence</div><div class="value">{html.escape(consequence["primary"])}</div></div>
    <div class="card"><div class="label">Impact score</div><div class="value">{score["score"]:.2f}/10</div></div>
    <div class="card"><div class="label">Category</div><div class="value">{html.escape(score["label"])}</div></div>
</div>

<h2>Biological interpretation</h2>
<div class="story">{html.escape(interpretation)}</div>

<h2>Evidence summary</h2>
<div class="card">
    <p><strong>Genome position:</strong> <span class="mono">{html.escape(str(record.get("seq_region_name", "—")))}:{html.escape(str(record.get("start", "—")))}</span></p>
    <p><strong>Assembly:</strong> {html.escape(str(record.get("assembly_name", "GRCh38")))} </p>
    <p><strong>Allele string:</strong> <span class="mono">{html.escape(str(record.get("allele_string", "—")))}</span></p>
    <p><strong>CADD-Phred:</strong> {html.escape(format_number(score.get("cadd"), 3))}</p>
    <p><strong>AlphaMissense:</strong> {html.escape(format_number(score.get("alphamissense"), 3))}</p>
</div>

<h2>Transcript consequences</h2>
{transcript_html}

<h2>Raw Ensembl VEP response</h2>
<pre>{raw_json}</pre>

<div class="note">
<strong>Scientific-use note:</strong>
The impact score shown here is a transparent computational heuristic for
education, rapid exploration, and portfolio demonstration. It is not a
validated clinical pathogenicity classifier and should not be used for
diagnosis, treatment, or patient management.
</div>
</div>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Reusable UI fragments
# ---------------------------------------------------------------------------

def html_block(markup: str) -> str:
    """Collapse indentation and blank lines so Markdown never interprets
    indented HTML as a code block."""
    return "".join(line.strip() for line in markup.splitlines() if line.strip())


def section_header(eyebrow: str, title: str, desc: Optional[str] = None) -> None:
    desc_html = f'<p class="vr-section-desc">{esc(desc)}</p>' if desc else ""
    st.markdown(
        html_block(
            f"""
            <div class="vr-section">
                <div class="vr-eyebrow">{esc(eyebrow)}</div>
                <div class="vr-section-title">{esc(title)}</div>
                {desc_html}
                <div class="vr-section-rule"></div>
            </div>
            """
        ),
        unsafe_allow_html=True,
    )


def metric_card(
    label: str, value: str, sub: Optional[str] = None,
    mono: bool = False, raw: bool = False,
) -> str:
    value_html = value if raw else esc(value)
    sub_html = f'<div class="vr-card-sub">{esc(sub)}</div>' if sub else ""
    mono_class = " mono" if mono else ""
    return (
        f'<div class="vr-card">'
        f'<div class="vr-card-label">{esc(label)}</div>'
        f'<div class="vr-card-value{mono_class}">{value_html}</div>'
        f"{sub_html}</div>"
    )


def metric_row(cards: List[str]) -> None:
    cols = st.columns(len(cards))
    for col, card in zip(cols, cards):
        with col:
            st.markdown(card, unsafe_allow_html=True)


def tier_chip(label: str, description: str) -> str:
    color = TIER_COLORS.get(label, AMBER)
    return (
        f'<span class="vr-tier" title="{esc(description)}" '
        f'style="color:{color}; border-color:{color};">'
        f'<span class="vr-tier-dot" style="background:{color};"></span>'
        f"{esc(label)} evidence</span>"
    )


def colored_alleles(allele_string: str) -> str:
    parts = []
    for ch in str(allele_string):
        if ch.upper() in BASE_COLORS:
            parts.append(
                f'<span class="vr-base" style="color:{BASE_COLORS[ch.upper()]}">{esc(ch)}</span>'
            )
        elif ch == "/":
            parts.append('<span class="vr-allele-sep">/</span>')
        else:
            parts.append(esc(ch))
    return "".join(parts)


def render_hero(result: Dict[str, Any]) -> None:
    record = result["record"]
    score = result["score"]
    chrom = record.get("seq_region_name", "—")
    start = record.get("start", "—")
    try:
        position = f"{int(start):,}"
    except (TypeError, ValueError):
        position = str(start)
    assembly = record.get("assembly_name", "GRCh38")
    allele = record.get("allele_string", "—")
    _, tier_desc = evidence_level(score["score"])

    analyzed = result.get("analyzed_at_utc", "")
    try:
        analyzed_display = datetime.fromisoformat(analyzed).strftime(
            "%Y-%m-%d %H:%M UTC"
        )
    except (TypeError, ValueError):
        analyzed_display = str(analyzed)

    gene_chip = (
        f'<span class="vr-gene-chip">{esc(result["gene_symbol"])}</span>'
        if result["gene_symbol"] != "Unknown"
        else ""
    )

    st.markdown(
        html_block(f"""
        <div class="vr-hero">
            <div class="vr-hero-top">
                <span class="vr-eyebrow">Ensembl VEP · {esc(assembly)}</span>
                {tier_chip(score["label"], tier_desc)}
            </div>
            <div class="vr-hero-id">{esc(result["variant_id"])}</div>
            <div class="vr-hero-locus">
                <span>chr{esc(chrom)} : {esc(position)}</span>
                <span>{colored_alleles(allele)}</span>
                {gene_chip}
            </div>
            <div class="vr-hero-meta">
                Annotation retrieved · {esc(analyzed_display)} ·
                response cached for 15 minutes
            </div>
        </div>
        """),
        unsafe_allow_html=True,
    )


def render_flow(result: Dict[str, Any]) -> None:
    tx = result["tx"]
    steps = [
        ("01", "Variant", result["variant_id"]),
        ("02", "Gene", result["gene_symbol"]),
        ("03", "Transcript", tx.get("transcript_id", "—") if tx else "—"),
        ("04", "Consequence", humanize(result["consequence"]["primary"])),
        ("05", "Evidence", result["score"]["label"]),
    ]
    cells = []
    for i, (num, label, value) in enumerate(steps):
        cells.append(
            f'<div class="vr-flow-step">'
            f'<div class="vr-flow-num">{num}</div>'
            f'<div class="vr-flow-label">{esc(label)}</div>'
            f'<div class="vr-flow-value">{esc(value)}</div>'
            f"</div>"
        )
        if i < len(steps) - 1:
            cells.append('<div class="vr-flow-arrow">→</div>')
    st.markdown(
        f'<div class="vr-flow">{"".join(cells)}</div>', unsafe_allow_html=True
    )


def render_empty_state() -> None:
    st.markdown(
        html_block("""
        <div class="vr-empty">
            <div class="vr-eyebrow">No analysis loaded</div>
            <div class="vr-empty-title">Query a variant to begin</div>
            <div class="vr-empty-body">
                Enter a dbSNP identifier such as <strong>rs7903146</strong> in the
                console on the left, or pick one of the example variants.
                The report is generated from the current Ensembl VEP response.
            </div>
        </div>
        """),
        unsafe_allow_html=True,
    )


def page_footer() -> None:
    st.markdown(
        f'<div class="vr-footer">{esc(APP_TITLE)} v{esc(APP_VERSION)} · '
        f"Ensembl VEP-powered research prototype · "
        f"Not for clinical decision-making.</div>",
        unsafe_allow_html=True,
    )


def current_result() -> Optional[Dict[str, Any]]:
    return st.session_state.get("current")


def require_result() -> Optional[Dict[str, Any]]:
    result = current_result()
    if result is None:
        render_empty_state()
        page_footer()
    return result


# ---------------------------------------------------------------------------
# Analysis pipeline + session state
# ---------------------------------------------------------------------------

def init_state() -> None:
    st.session_state.setdefault("variant_query", "rs7903146")
    st.session_state.setdefault("current", None)
    st.session_state.setdefault("history", {})
    st.session_state.setdefault("last_error", None)
    st.session_state.setdefault("boot_done", False)


def perform_analysis(raw_value: str) -> bool:
    variant_id = normalize_variant_id(raw_value)

    if not variant_id:
        st.session_state["last_error"] = (
            "Enter a variant identifier such as rs7903146."
        )
        return False

    if not valid_variant_id(variant_id):
        st.session_state["last_error"] = (
            "Variant ID contains unsupported characters."
        )
        return False

    with st.spinner(f"Querying Ensembl VEP for {variant_id}…"):
        try:
            record = fetch_vep(variant_id)
        except Exception as exc:
            st.session_state["last_error"] = str(exc)
            return False

    tx = choose_primary_transcript(record)
    gene_symbol = choose_gene_symbol(record, tx)
    gene_summary = (
        fetch_gene_summary(gene_symbol) if gene_symbol != "Unknown" else {}
    )
    consequence = derive_consequence_summary(record, tx)
    score = calculate_impact_score(record, tx, gene_summary)
    interpretation = biological_interpretation(
        record, tx, gene_symbol, consequence, score, gene_summary,
    )
    transcript_df = make_transcript_dataframe(record)
    evidence_df = make_evidence_dataframe(consequence, score)
    freq_df = extract_population_frequencies(record)
    clinical = extract_clinical_context(record)

    result = {
        "variant_id": variant_id,
        "record": record,
        "tx": tx,
        "gene_symbol": gene_symbol,
        "gene_summary": gene_summary,
        "consequence": consequence,
        "score": score,
        "interpretation": interpretation,
        "transcript_df": transcript_df,
        "evidence_df": evidence_df,
        "freq_df": freq_df,
        "clinical": clinical,
        "analyzed_at_utc": datetime.now(timezone.utc).isoformat(),
    }

    st.session_state["current"] = result
    history = st.session_state["history"]
    history.pop(variant_id, None)  # move to most-recent position
    history[variant_id] = result
    st.session_state["last_error"] = None
    st.session_state["flash"] = f"Annotation retrieved for {variant_id}"
    try:
        st.query_params["v"] = variant_id
    except Exception:
        pass
    return True


def _on_example_pick() -> None:
    value = st.session_state.get("example_pick")
    st.session_state["example_pick"] = None
    if value:
        st.session_state["variant_query"] = value
        st.session_state["pending_query"] = value
        st.session_state["user_query"] = True


def _on_history_pick() -> None:
    value = st.session_state.get("history_pick")
    st.session_state["history_pick"] = None
    history = st.session_state.get("history", {})
    if value and value in history:
        st.session_state["current"] = history[value]
        st.session_state["variant_query"] = value
        st.session_state["flash"] = f"Loaded {value} from session history"


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------

PAGES: Dict[str, Any] = {}


def overview_page() -> None:
    result = require_result()
    if result is None:
        return

    render_hero(result)

    consequence = result["consequence"]
    score = result["score"]

    st.markdown('<div style="height:.9rem"></div>', unsafe_allow_html=True)
    metric_row(
        [
            metric_card("Gene", result["gene_symbol"]),
            metric_card(
                "Primary consequence",
                humanize(consequence["primary"]),
                sub=f'{len(consequence["terms"])} annotated term(s)',
            ),
            metric_card("Impact score", f'{score["score"]:.2f} / 10'),
            metric_card(
                "Evidence category",
                score["label"],
                sub=evidence_level(score["score"])[1],
            ),
        ]
    )

    section_header(
        "Pipeline", "Biological flow",
        "How the identifier resolves into an evidence category.",
    )
    render_flow(result)

    section_header(
        "Signals", "Impact profile & evidence radar",
        "The heuristic score alongside its normalized evidence axes.",
    )
    left, right = st.columns(2)
    with left:
        st.plotly_chart(
            make_score_gauge(score["score"]),
            width="stretch",
            config={"displayModeBar": False},
        )
        with st.popover("Why this score?", width="stretch"):
            for reason in score["reasons"]:
                st.markdown(f"- {reason}")
            st.page_link(
                PAGES["methodology"],
                label="Full scoring methodology",
                icon=":material/menu_book:",
            )
    with right:
        st.plotly_chart(
            make_radar_chart(
                score["score"], score.get("cadd"), score.get("alphamissense"),
            ),
            width="stretch",
            config={"displayModeBar": False},
        )
        st.page_link(
            PAGES["evidence"],
            label="Open the full evidence breakdown",
            icon=":material/rule:",
        )

    section_header(
        "Narrative", "Biological interpretation",
        "A plain-language reading of the returned annotation.",
    )
    st.markdown(
        f'<div class="vr-story">{esc(result["interpretation"])}</div>',
        unsafe_allow_html=True,
    )

    page_footer()


def evidence_page() -> None:
    result = require_result()
    if result is None:
        return

    record = result["record"]
    score = result["score"]
    clinical = result["clinical"]

    section_header(
        "Locus", "Variant metadata",
        "Genomic coordinates as reported by the Ensembl VEP response.",
    )
    metric_row(
        [
            metric_card(
                "Chromosome", str(record.get("seq_region_name", "—")), mono=True
            ),
            metric_card("Position", str(record.get("start", "—")), mono=True),
            metric_card(
                "Assembly", str(record.get("assembly_name", "GRCh38")), mono=True
            ),
            metric_card(
                "Alleles",
                colored_alleles(record.get("allele_string", "—")),
                mono=True,
                raw=True,
            ),
        ]
    )

    section_header(
        "Channels", "Evidence detail",
        "Score construction, population frequencies, and clinical annotations.",
    )
    tab_score, tab_freq, tab_clin = st.tabs(
        ["Score breakdown", "Population frequencies", "Clinical annotations"]
    )

    with tab_score:
        e_left, e_right = st.columns([1, 1])
        with e_left:
            st.dataframe(
                result["evidence_df"],
                width="stretch",
                hide_index=True,
            )
        with e_right:
            st.markdown("**Scoring logic**")
            for reason in score["reasons"]:
                st.markdown(f"- {reason}")
        st.caption(
            "The score combines consequence class and selected annotation "
            "evidence using fixed, visible rules. It is intentionally not a "
            "clinical classifier."
        )

    with tab_freq:
        freq_df = result["freq_df"]
        if freq_df.empty:
            st.info(
                "No population allele frequencies were included in this VEP "
                "response."
            )
        else:
            st.plotly_chart(
                make_frequency_chart(freq_df),
                width="stretch",
                config={"displayModeBar": False},
            )
            st.dataframe(
                freq_df,
                width="stretch",
                hide_index=True,
                column_config={
                    "Frequency": st.column_config.NumberColumn(format="%.4f"),
                },
            )

    with tab_clin:
        has_any = (
            clinical["clin_sig"]
            or clinical["phenotypes"]
            or clinical["pubmed"]
            or clinical["minor_allele"]
        )
        if not has_any:
            st.info(
                "No clinical-significance or phenotype annotations were "
                "included in this VEP response."
            )
        else:
            if clinical["minor_allele"]:
                maf = clinical.get("minor_allele_freq")
                m1, m2 = st.columns(2)
                with m1:
                    st.markdown(
                        metric_card(
                            "Minor allele",
                            str(clinical["minor_allele"]),
                            mono=True,
                        ),
                        unsafe_allow_html=True,
                    )
                with m2:
                    st.markdown(
                        metric_card(
                            "Minor allele frequency",
                            format_number(maf, 4) if maf is not None else "—",
                            mono=True,
                        ),
                        unsafe_allow_html=True,
                    )
                st.markdown("")
            if clinical["clin_sig"]:
                st.markdown("**Reported clinical significance**")
                chips = "".join(
                    f'<span class="vr-chip">{esc(humanize(s))}</span>'
                    for s in clinical["clin_sig"]
                )
                st.markdown(chips, unsafe_allow_html=True)
            if clinical["phenotypes"]:
                st.markdown("**Associated phenotypes**")
                for name in clinical["phenotypes"][:12]:
                    st.markdown(f"- {name}")
                if len(clinical["phenotypes"]) > 12:
                    st.caption(
                        f'+ {len(clinical["phenotypes"]) - 12} more in the raw '
                        "response (Reports & export)."
                    )
            if clinical["pubmed"]:
                st.markdown(
                    f'**Literature** · {len(clinical["pubmed"])} linked '
                    "PubMed record(s)"
                )
                chips = "".join(
                    f'<span class="vr-chip"><a href="https://pubmed.ncbi.nlm.nih.gov/{esc(pmid)}/" '
                    f'target="_blank" rel="noopener">PMID {esc(pmid)}</a></span>'
                    for pmid in clinical["pubmed"][:8]
                )
                st.markdown(chips, unsafe_allow_html=True)

    page_footer()


def transcripts_page() -> None:
    result = require_result()
    if result is None:
        return

    tx = result["tx"]
    gene_summary = result["gene_summary"]
    transcript_df = result["transcript_df"]

    section_header(
        "Selection", "Selected transcript",
        "Chosen by canonical status, then MANE annotation, then detail.",
    )
    if tx:
        metric_row(
            [
                metric_card(
                    "Transcript", str(tx.get("transcript_id", "—")), mono=True
                ),
                metric_card(
                    "Canonical", "Yes" if tx.get("canonical") else "No"
                ),
                metric_card(
                    "MANE",
                    "Yes"
                    if tx.get("mane_select") or tx.get("mane_plus_clinical")
                    else "No",
                ),
                metric_card(
                    "Protein change",
                    str(tx.get("amino_acids", "—")),
                    mono=True,
                ),
            ]
        )
    else:
        st.info("No transcript consequence objects were returned.")

    section_header(
        "Consequences", "Transcript consequences",
        "Every transcript-level annotation in the VEP response.",
    )
    if transcript_df.empty:
        st.info("No transcript consequence objects were returned.")
    else:
        view = st.segmented_control(
            "Filter",
            ["All", "Canonical only", "MANE only"],
            default="All",
            label_visibility="collapsed",
        )
        filtered = transcript_df
        if view == "Canonical only":
            filtered = transcript_df[transcript_df["Canonical"] == "Yes"]
        elif view == "MANE only":
            filtered = transcript_df[transcript_df["MANE"] == "Yes"]

        st.dataframe(
            filtered,
            width="stretch",
            hide_index=True,
            column_config={
                "Canonical": st.column_config.TextColumn(width="small"),
                "MANE": st.column_config.TextColumn(width="small"),
            },
        )
        st.caption(
            f"{len(filtered)} of {len(transcript_df)} transcript "
            "consequence(s) shown."
        )

    section_header(
        "Context", "Gene context",
        "Gene-level summary from the Ensembl symbol lookup.",
    )
    if gene_summary:
        metric_row(
            [
                metric_card(
                    "Gene ID", str(gene_summary.get("id", "—")), mono=True
                ),
                metric_card("Biotype", str(gene_summary.get("biotype", "—"))),
                metric_card("Species", str(gene_summary.get("species", "—"))),
                metric_card(
                    "Assembly", str(gene_summary.get("assembly_name", "—")),
                    mono=True,
                ),
            ]
        )
        if gene_summary.get("description"):
            st.markdown("")
            st.markdown(f'**Description** — {gene_summary["description"]}')
    else:
        st.info("Gene-level lookup was unavailable for this symbol.")

    section_header(
        "Resources", "External evidence links",
        "Primary databases for deeper manual review.",
    )
    ensembl_variant_url = (
        f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={result['variant_id']}"
    )
    dbsnp_search_url = (
        f"https://www.ncbi.nlm.nih.gov/snp/?term={result['variant_id']}"
    )
    ensembl_gene_url = (
        f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={result['gene_symbol']}"
        if result["gene_symbol"] != "Unknown"
        else None
    )

    l1, l2, l3 = st.columns(3)
    with l1:
        st.link_button(
            "Ensembl variant page",
            ensembl_variant_url,
            width="stretch",
        )
    with l2:
        st.link_button(
            "NCBI dbSNP search",
            dbsnp_search_url,
            width="stretch",
        )
    with l3:
        if ensembl_gene_url:
            st.link_button(
                "Ensembl gene page",
                ensembl_gene_url,
                width="stretch",
            )

    page_footer()


def compare_page() -> None:
    history: Dict[str, Dict[str, Any]] = st.session_state.get("history", {})

    section_header(
        "Session", "Compare variants",
        "Side-by-side evidence for variants analyzed in this session.",
    )

    if not history:
        render_empty_state()
        page_footer()
        return

    if len(history) == 1:
        st.info(
            "One variant is in the session so far. Analyze another from the "
            "console to unlock side-by-side comparison."
        )

    ids = list(history.keys())
    default = ids[-4:]
    selected = st.multiselect(
        "Variants to compare", ids, default=default,
    )

    if not selected:
        st.info("Select at least one analyzed variant.")
        page_footer()
        return

    rows = []
    for vid in selected:
        r = history[vid]
        rows.append(
            {
                "Variant": vid,
                "Gene": r["gene_symbol"],
                "Primary consequence": humanize(r["consequence"]["primary"]),
                "Impact score": r["score"]["score"],
                "Evidence tier": r["score"]["label"],
                "CADD-Phred": r["score"].get("cadd"),
                "AlphaMissense": r["score"].get("alphamissense"),
            }
        )

    st.dataframe(
        pd.DataFrame(rows),
        width="stretch",
        hide_index=True,
        column_config={
            "Impact score": st.column_config.ProgressColumn(
                min_value=0, max_value=10, format="%.2f",
            ),
            "CADD-Phred": st.column_config.NumberColumn(format="%.2f"),
            "AlphaMissense": st.column_config.NumberColumn(format="%.3f"),
        },
    )

    c_left, c_right = st.columns(2)
    with c_left:
        st.plotly_chart(
            make_compare_bar(rows),
            width="stretch",
            config={"displayModeBar": False},
        )
    with c_right:
        st.plotly_chart(
            make_compare_radar([history[vid] for vid in selected]),
            width="stretch",
            config={"displayModeBar": False},
        )

    st.caption(
        "Scores are the same transparent heuristic shown on the Overview "
        "page; comparison does not imply clinical ranking."
    )
    page_footer()


def reports_page() -> None:
    result = require_result()
    if result is None:
        return

    record = result["record"]
    score = result["score"]
    consequence = result["consequence"]
    tx = result["tx"]

    section_header(
        "Deliverables", "Export report",
        "Portable artifacts built from the current analysis.",
    )

    export_payload = {
        "app": APP_TITLE,
        "version": APP_VERSION,
        "analyzed_at_utc": result["analyzed_at_utc"],
        "variant_id": result["variant_id"],
        "gene_symbol": result["gene_symbol"],
        "primary_consequence": consequence["primary"],
        "consequence_terms": consequence["terms"],
        "impact_score": score["score"],
        "impact_category": score["label"],
        "cadd_phred": score.get("cadd"),
        "alphamissense": score.get("alphamissense"),
        "interpretation": result["interpretation"],
        "variant_metadata": {
            "assembly": record.get("assembly_name"),
            "chromosome": record.get("seq_region_name"),
            "start": record.get("start"),
            "end": record.get("end"),
            "allele_string": record.get("allele_string"),
        },
        "selected_transcript": tx,
        "gene_summary": result["gene_summary"],
        "scoring_reasons": score["reasons"],
        "raw_ensembl_vep": record,
    }

    json_bytes = json.dumps(
        export_payload, indent=2, ensure_ascii=False, default=str,
    ).encode("utf-8")

    csv_bytes = result["transcript_df"].to_csv(index=False).encode("utf-8")

    html_report = build_html_report(
        variant_id=result["variant_id"],
        record=record,
        gene_symbol=result["gene_symbol"],
        consequence=consequence,
        score=score,
        interpretation=result["interpretation"],
        transcript_df=result["transcript_df"],
    )

    x1, x2, x3 = st.columns(3)
    with x1:
        st.markdown(
            metric_card(
                "Machine-readable", "JSON",
                sub="Full payload incl. raw VEP response",
                mono=True,
            ),
            unsafe_allow_html=True,
        )
        st.download_button(
            "Download JSON",
            data=json_bytes,
            file_name=f"{result['variant_id']}_variant_report.json",
            mime="application/json",
            width="stretch",
            icon=":material/download:",
        )
    with x2:
        st.markdown(
            metric_card(
                "Tabular", "CSV",
                sub="Transcript consequence table",
                mono=True,
            ),
            unsafe_allow_html=True,
        )
        st.download_button(
            "Download transcript CSV",
            data=csv_bytes,
            file_name=f"{result['variant_id']}_transcripts.csv",
            mime="text/csv",
            width="stretch",
            icon=":material/download:",
        )
    with x3:
        st.markdown(
            metric_card(
                "Shareable", "HTML",
                sub="Self-contained printable report",
                mono=True,
            ),
            unsafe_allow_html=True,
        )
        st.download_button(
            "Download HTML report",
            data=html_report.encode("utf-8"),
            file_name=f"{result['variant_id']}_variant_report.html",
            mime="text/html",
            width="stretch",
            icon=":material/download:",
        )

    section_header(
        "Source", "Raw Ensembl VEP response",
        "The unmodified annotation this analysis was built from.",
    )
    with st.expander("Raw Ensembl VEP response", expanded=False):
        st.json(record)

    page_footer()


def methodology_page() -> None:
    section_header(
        "Reference", "What this app does",
        "Scope, data sources, and the exact scoring rules.",
    )
    st.markdown(
        "Variant Impact Radar queries **Ensembl VEP** for a single variant "
        "identifier, summarizes the returned consequence evidence, adds a "
        "transparent heuristic score, and builds a shareable report. It is a "
        "research and education prototype: every number shown can be traced "
        "back to a visible rule on this page or to the raw VEP response."
    )
    st.link_button(
        "Open Ensembl REST documentation",
        ENSEMBL_DOCS_URL,
        icon=":material/open_in_new:",
    )

    section_header(
        "Rules", "Consequence weights",
        "The base contribution of each consequence class (0–10).",
    )
    weights_df = pd.DataFrame(
        [
            {"Consequence": humanize(k), "Weight": v}
            for k, v in sorted(
                CONSEQUENCE_WEIGHTS.items(), key=lambda kv: -kv[1]
            )
        ]
    )
    st.dataframe(
        weights_df,
        width="stretch",
        hide_index=True,
        column_config={
            "Weight": st.column_config.ProgressColumn(
                min_value=0, max_value=10, format="%d",
            ),
        },
    )

    section_header(
        "Rules", "Annotation evidence & context bonuses",
        "Fixed additive adjustments applied on top of the consequence weight.",
    )
    b_left, b_right = st.columns(2)
    with b_left:
        st.markdown("**Annotation thresholds**")
        st.markdown(
            "- CADD-Phred ≥ 30 → +2.0\n"
            "- CADD-Phred ≥ 20 → +1.0\n"
            "- CADD-Phred ≥ 10 → +0.5\n"
            "- AlphaMissense ≥ 0.90 → +1.5\n"
            "- AlphaMissense ≥ 0.70 → +1.0\n"
            "- AlphaMissense ≥ 0.50 → +0.5"
        )
    with b_right:
        st.markdown("**Context bonuses**")
        st.markdown(
            "- Canonical transcript → +0.5\n"
            "- MANE transcript annotation → +0.5\n"
            "- Known colocated-variant context → +0.5\n"
            "- Phenotype annotation present → +0.5\n"
            "- Gene-level Ensembl context available → +0.25"
        )
    st.caption(f"The total is clamped to the 0–{MAX_SCORE:g} range.")

    section_header(
        "Legend", "Evidence tiers",
        "How the clamped score maps to a category.",
    )
    tiers = [
        ("HIGH", "score ≥ 8", "Potentially high-impact computational evidence"),
        ("MODERATE", "score ≥ 5", "Moderate computational impact"),
        ("LOW–MODERATE", "score ≥ 2.5", "Limited computational impact"),
        ("LOW", "score < 2.5", "Low computational impact"),
    ]
    for label, rule, desc in tiers:
        c1, c2 = st.columns([1, 3])
        with c1:
            st.markdown(tier_chip(label, desc), unsafe_allow_html=True)
        with c2:
            st.markdown(f"`{rule}` — {desc}")

    section_header("Scope", "Scientific-use note")
    st.warning(
        "Demo / research use only. The score is not a clinical pathogenicity "
        "assessment and must not be used for diagnosis, treatment, or "
        "patient management.",
        icon=":material/science:",
    )
    st.caption(
        f"{APP_TITLE} v{APP_VERSION} · Ensembl VEP responses cached for "
        "15 minutes · gene lookups cached for 60 minutes."
    )
    page_footer()


# ---------------------------------------------------------------------------
# About-the-developer page
# ---------------------------------------------------------------------------

def _developer_photo_src() -> str:
    """
    Return a data-URI for the developer portrait.

    The first existing file in DEV_PHOTO_CANDIDATES is embedded, so adding a
    real photo is a pure asset drop (e.g. assets/developer.png) with no code
    change. Until then a neutral "blank profile" silhouette — matching the
    classic placeholder avatar — is rendered inline so the page never shows
    a broken image.
    """
    base = Path(__file__).parent
    for rel in DEV_PHOTO_CANDIDATES:
        photo = base / rel
        if photo.exists():
            mime = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".webp": "image/webp",
            }.get(photo.suffix.lower(), "image/png")
            encoded = base64.b64encode(photo.read_bytes()).decode("ascii")
            return f"data:{mime};base64,{encoded}"

    # Placeholder silhouette (grey bust), tuned to the dark theme.
    placeholder_svg = (
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 200 200'>"
        "<rect width='200' height='200' fill='#211C15'/>"
        "<circle cx='100' cy='76' r='36' fill='#8B8577'/>"
        "<path d='M100 120 c-43 0 -63 27 -67 56 a100 100 0 0 0 134 0 "
        "c-4 -29 -24 -56 -67 -56 z' fill='#8B8577'/>"
        "</svg>"
    )
    encoded = base64.b64encode(placeholder_svg.encode("utf-8")).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def about_developer_page() -> None:
    dev = DEVELOPER_PROFILE

    # Hero — portrait, name, role, and tagline in the locus-strip style.
    st.markdown(
        html_block(f"""
        <div class="vr-dev-hero">
            <img class="vr-dev-photo" src="{_developer_photo_src()}"
                 alt="Portrait of {esc(dev['name'])}">
            <div class="vr-dev-hero-body">
                <div class="vr-eyebrow">{esc(dev['role'])}</div>
                <div class="vr-dev-name">{esc(dev['name'])}</div>
                <div class="vr-dev-loc">
                    {esc(dev['location'])} · {esc(dev['availability'])}
                </div>
                <div class="vr-dev-tagline">“{esc(dev['tagline'])}”</div>
            </div>
        </div>
        """),
        unsafe_allow_html=True,
    )

    st.markdown('<div style="height:.9rem"></div>', unsafe_allow_html=True)
    l1, l2, l3, l4 = st.columns(4)
    with l1:
        st.link_button(
            "LinkedIn", dev["linkedin"],
            width="stretch", icon=":material/open_in_new:",
        )
    with l2:
        st.link_button(
            "GitHub", dev["github"],
            width="stretch", icon=":material/code:",
        )
    with l3:
        st.link_button(
            "Portfolio", dev["portfolio"],
            width="stretch", icon=":material/language:",
        )
    with l4:
        st.link_button(
            "Email", f"mailto:{dev['email']}",
            width="stretch", icon=":material/mail:",
        )

    section_header(
        "Profile", "Biography",
        "Who built Variant Impact Radar, and why.",
    )
    st.markdown(
        f'<div class="vr-story">{esc(dev["bio"])}</div>',
        unsafe_allow_html=True,
    )

    section_header(
        "Snapshot", "At a glance",
        "A quick read of experience and output.",
    )
    metric_row(
        [metric_card(label, value, sub=sub) for label, value, sub in dev["stats"]]
    )

    section_header(
        "Toolbox", "Skills & technologies",
        "The day-to-day stack behind projects like this one.",
    )
    skill_chips = "".join(
        f'<span class="vr-chip">{esc(skill)}</span>' for skill in dev["skills"]
    )
    st.markdown(skill_chips, unsafe_allow_html=True)

    section_header(
        "Journey", "Experience highlights",
        "Selected milestones — newest first.",
    )
    timeline_items = "".join(
        f'<div class="vr-dev-tl-item">'
        f'<div class="vr-dev-tl-year">{esc(year)}</div>'
        f'<div class="vr-dev-tl-title">{esc(title)}</div>'
        f'<div class="vr-dev-tl-desc">{esc(desc)}</div>'
        f"</div>"
        for year, title, desc in dev["timeline"]
    )
    st.markdown(
        f'<div class="vr-dev-timeline">{timeline_items}</div>',
        unsafe_allow_html=True,
    )

    section_header(
        "Build notes", "Behind this project",
        "How this workbench came together.",
    )
    st.markdown(
        "Variant Impact Radar doubles as a portfolio piece: a worked example "
        "of a clean, API-driven analytical application — input validation, "
        "cached HTTP calls, defensive parsing of a deeply nested VEP "
        "response, and a fixed, fully visible scoring scheme. Everything on "
        "screen can be audited against the raw Ensembl payload on the "
        "*Reports & export* page."
    )
    st.page_link(
        PAGES["methodology"],
        label="Read the full scoring methodology",
        icon=":material/menu_book:",
    )

    section_header(
        "Connect", "Find me online",
        "Collaborations, feedback, or a friendly hello.",
    )
    contact_links = [
        ("LinkedIn", dev["linkedin"]),
        ("GitHub", dev["github"]),
        ("Portfolio", dev["portfolio"]),
        (dev["email"], f"mailto:{dev['email']}"),
    ]
    contact_chips = "".join(
        f'<span class="vr-chip"><a href="{esc(url)}" target="_blank" '
        f'rel="noopener">{esc(label)}</a></span>'
        for label, url in contact_links
    )
    st.markdown(contact_chips, unsafe_allow_html=True)

    page_footer()


# ---------------------------------------------------------------------------
# Sidebar console
# ---------------------------------------------------------------------------

def render_sidebar() -> bool:
    with st.sidebar:
        st.markdown(
            '<div class="vr-side-label">Query console</div>',
            unsafe_allow_html=True,
        )
        with st.form("query_form", border=False):
            st.text_input(
                "dbSNP / variant identifier",
                key="variant_query",
                placeholder="e.g. rs7903146",
                help=(
                    "Start with a dbSNP rsID. Ensembl VEP also supports other "
                    "identifier classes (COSMIC, HGMD-style)."
                ),
            )
            submitted = st.form_submit_button(
                "Analyze variant",
                type="primary",
                width="stretch",
            )

        if st.session_state.get("last_error"):
            st.error(st.session_state["last_error"], icon=":material/error:")

        st.markdown(
            '<div class="vr-side-label">Example variants</div>',
            unsafe_allow_html=True,
        )
        st.pills(
            "Example variants",
            EXAMPLE_VARIANTS,
            selection_mode="single",
            key="example_pick",
            on_change=_on_example_pick,
            label_visibility="collapsed",
        )

        history = st.session_state.get("history", {})
        if history:
            st.markdown(
                '<div class="vr-side-label">Session history</div>',
                unsafe_allow_html=True,
            )
            recent = list(history.keys())[::-1][:6]
            st.pills(
                "Session history",
                recent,
                selection_mode="single",
                key="history_pick",
                on_change=_on_history_pick,
                label_visibility="collapsed",
            )

        st.markdown(
            '<div class="vr-side-note">'
            "Queries Ensembl VEP, summarizes the returned consequence "
            "evidence, adds a transparent heuristic score, and builds a "
            "shareable report."
            "<br><br>"
            "<strong>Research demo only.</strong> The score is not a clinical "
            "pathogenicity assessment."
            "</div>",
            unsafe_allow_html=True,
        )
        st.link_button(
            "Ensembl REST documentation",
            ENSEMBL_DOCS_URL,
            width="stretch",
            icon=":material/open_in_new:",
        )
        st.caption(f"{APP_TITLE} v{APP_VERSION}")

    return submitted


# ---------------------------------------------------------------------------
# App entry
# ---------------------------------------------------------------------------

def main() -> None:
    init_state()
    inject_css()

    logo_path = Path(__file__).parent / "logo.svg"
    if logo_path.exists():
        st.logo(str(logo_path), size="large")

    PAGES["overview"] = st.Page(
        overview_page, title="Overview", url_path="overview",
        icon=":material/radar:", default=True,
    )
    PAGES["evidence"] = st.Page(
        evidence_page, title="Evidence & scoring", url_path="evidence",
        icon=":material/rule:",
    )
    PAGES["transcripts"] = st.Page(
        transcripts_page, title="Transcripts & gene", url_path="transcripts",
        icon=":material/account_tree:",
    )
    PAGES["compare"] = st.Page(
        compare_page, title="Compare variants", url_path="compare",
        icon=":material/compare_arrows:",
    )
    PAGES["reports"] = st.Page(
        reports_page, title="Reports & export", url_path="reports",
        icon=":material/download:",
    )
    PAGES["methodology"] = st.Page(
        methodology_page, title="Methodology", url_path="methodology",
        icon=":material/menu_book:",
    )
    PAGES["about"] = st.Page(
        about_developer_page, title="About the developer", url_path="about",
        icon=":material/person:",
    )

    pg = st.navigation(
        {
            "Workbench": [
                PAGES["overview"], PAGES["evidence"], PAGES["transcripts"],
            ],
            "Tools": [PAGES["compare"], PAGES["reports"]],
            "Reference": [PAGES["methodology"], PAGES["about"]],
        }
    )

    # First load: honour a shared ?v= link, else auto-run the default example
    # so the workbench opens populated (same behaviour as v1). This runs
    # before the sidebar so the input widget can still be synced safely.
    if st.session_state["current"] is None and not st.session_state["boot_done"]:
        st.session_state["boot_done"] = True
        try:
            shared = st.query_params.get("v")
        except Exception:
            shared = None
        boot_query = (
            shared or st.session_state.get("variant_query") or "rs7903146"
        )
        st.session_state["variant_query"] = boot_query
        st.session_state["pending_query"] = boot_query

    submitted = render_sidebar()

    if submitted:
        st.session_state["pending_query"] = st.session_state.get(
            "variant_query", ""
        )

    # Example-pill picks (set via callback) also count as user-initiated.
    user_initiated = submitted or st.session_state.pop("user_query", False)
    pending = st.session_state.pop("pending_query", None)
    if pending:
        ok = perform_analysis(pending)
        # Jump to the Overview only for analyses the user just launched;
        # the silent boot analysis must not hijack deep links to other pages.
        if ok and user_initiated and pg is not PAGES["overview"]:
            st.switch_page(PAGES["overview"])

    flash = st.session_state.pop("flash", None)
    if flash:
        st.toast(flash, icon=":material/check_circle:")

    pg.run()


if __name__ == "__main__":
    main()
