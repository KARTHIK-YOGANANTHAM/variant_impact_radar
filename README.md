# variant_impact_radar
An API based Variant Analysis and Interpretation Tool

# Variant Impact Radar — Workbench Edition (v2.0.0)

Multi-page Streamlit workbench that turns a dbSNP / COSMIC / HGMD-style
variant identifier into a structured computational evidence report using the
Ensembl VEP REST API.

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

Keep the folder structure intact — `.streamlit/config.toml` carries the theme
and `assets/logo.svg` is the sidebar wordmark:

```
variant-impact-radar/
├── app.py
├── requirements.txt
├── assets/
│   └── logo.svg
└── .streamlit/
    └── config.toml
```

## Pages

- **Overview** — locus hero, summary cards, biological flow, impact gauge,
  evidence radar, narrative interpretation.
- **Evidence & scoring** — variant metadata, score breakdown, population
  allele frequencies, clinical annotations (when present in the response).
- **Transcripts & gene** — selected transcript, filterable transcript
  consequence table, gene context, external database links.
- **Compare variants** — side-by-side scores for variants analyzed this
  session.
- **Reports & export** — JSON / CSV / HTML downloads and the raw VEP response.
- **Methodology** — the full, transparent scoring rules and disclaimer.

Analyses are shareable via URL (`?v=rs7903146`).

> Research / education demo only. The impact score is a transparent heuristic,
> not a clinical pathogenicity classifier.
