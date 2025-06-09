import json
import os
import pysam
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

# File paths
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68Samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/88_samples/1000g-ONT-88Samples.vcf.gz"
VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/503_samples/1000g-ONT-STRchive-503Samples.vcf.gz"

JSON_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json"

# OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Tandem_Repeats_Plots/68Samples_tr_plots"
# OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Tandem_Repeats_Plots/88Samples_tr_plots"
OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Tandem_Repeats_Plots/503Samples_tr_plots"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load JSON metadata for loci
with open(JSON_PATH, "r") as file:
    loci_data = json.load(file)

# Function to annotate the plot with range lines and labels
def add_range_marker_or_line(fig, x0, x1, label_id, label_positions,
                              chart_width_px, x_span,
                              base_y_line=1.02, base_y_label=1.1,
                              bump_y=0.08, bump_px_threshold=60):
    
    # Adds either a horizontal line or vertical line (if a single point) to the plot
    # with a label above it. Bumps labels if overlapping.
    
    x_center = (x0 + x1) / 2
    pixels_per_unit = chart_width_px / x_span

    # Bump logic based on visual overlap in pixel space
    level = 0
    for other_x in label_positions.values():
        if abs(other_x - x_center) * pixels_per_unit < bump_px_threshold:
            level += 1

    label_positions[label_id] = x_center
    y_line = base_y_line + bump_y * level
    y_label = base_y_label + bump_y * level

    label = f"{label_id.capitalize()}" if x0 == x1 else f"{label_id.capitalize()}"

    # Add label above the line
    fig.add_annotation(
        x=x_center, y=y_label,
        text=label, showarrow=False,
        xref="x", yref="paper",
        font=dict(color="black", size=12)
    )

    # Add line (dashed if a single point)
    if x0 == x1:
        fig.add_shape(
            type="line", x0=x0, x1=x0, y0=0, y1=1,
            xref="x", yref="paper",
            line=dict(color="black", width=1, dash="dot")
        )
    else:
        fig.add_shape(
            type="line", x0=x0, x1=x1, y0=y_line, y1=y_line,
            xref="x", yref="paper",
            line=dict(color="black", width=1)
        )

# Read the VCF
vcf_in = pysam.VariantFile(VCF_PATH)

for record in vcf_in.fetch():
    chrom = record.chrom
    pos = record.pos

    matched_locus = next(
        (entry for entry in loci_data
         if chrom == entry["chrom"] and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
        None
    )
    if not matched_locus:
        print(f"Skipped: {chrom}:{pos} — no matching locus in JSON")
        continue

    # Calculate motif length as # of bases
    motif = matched_locus.get("reference_motif_reference_orientation")
    # If motif is a list (e.g., ['GCC']), join to string and strip quotes/brackets
    if isinstance(motif, list):
        motif_str = "".join(motif)
    else:
        motif_str = str(motif)
    motif_length = len(motif_str)
    motif = motif_str

    # print(f"Processing {chrom}:{pos} - {motif} ({motif_length})")

    allele_lengths = []
    for sample in record.samples.values():
        al_lengths = sample.get("AL")
        if al_lengths:
            allele_lengths.extend([(length) for length in al_lengths if length is not None])
        
    repeat_counts = []
    for sample in record.samples.values():
        al_lengths = sample.get("AL")
        if al_lengths:
            repeat_counts.extend([(length) // motif_length for length in al_lengths if length is not None])

    if not repeat_counts:
        continue

    # print(f"{gene} ({chrom}:{pos}): {al_lengths}")

    gene = matched_locus.get("gene", "Unknown")
    disease = matched_locus.get("disease", "Unknown")
    benign_min = matched_locus.get("benign_min")
    benign_max = matched_locus.get("benign_max")
    pathogenic_min = matched_locus.get("pathogenic_min")
    pathogenic_max = matched_locus.get("pathogenic_max")

    # --- SAVE DEBUG INFO TO CSV ---
    debug_info = {
        "Chromosome": chrom,
        "Position": pos,
        "Gene": gene,
        "Disease": disease,
        "Motif": motif,
        "Allele lengths": ";".join(map(str, allele_lengths)),
        "Motif length": motif_length,
        "Repeat counts": ";".join(map(str, repeat_counts)),
        "Benign min": benign_min,
        "Benign max": benign_max,
        "Pathogenic min": pathogenic_min,
        "Pathogenic max": pathogenic_max,
        "Output file": os.path.join(OUTPUT_DIR, f"{gene}_{chrom}_{pos}_allele_dist.html"),
        "Number of alleles": len(repeat_counts),
        "Min repeat count": min(repeat_counts),
        "Max repeat count": max(repeat_counts),
        "Repeat counts (unique)": len(set(repeat_counts))
    }

    # # Uncomment the following lines to save debug info to a CSV file:
    # OUTPUT_DIR2 = "/Users/annelisethorn/Documents/Anschutz/Code/Plotting/CSVs"
    # debug_csv_path = os.path.join(OUTPUT_DIR2, "debug_info.csv")
    # write_header = not os.path.exists(debug_csv_path)
    # pd.DataFrame([debug_info]).to_csv(debug_csv_path, mode='a', header=write_header, index=False)
    # # --- END DEBUG INFO ---

    df = pd.DataFrame({"Repeat Count": repeat_counts})
    min_bin, max_bin = min(repeat_counts), max(repeat_counts)

    fig = px.histogram(
        df, x="Repeat Count",
        nbins=(max_bin - min_bin + 1),
        title=f"Allele Size Distribution<br><sup>{chrom.title()} - {pos} | {gene} - {disease.title()}</sup>",
    )

    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=160),
        title_font_size=18,
        title_font=dict(weight="bold"),
        title_font_color="black",
        xaxis_title="Repeat Count",
        yaxis_title="Allele Count",
        bargap=0.7,
        plot_bgcolor='white'
    )

    # Add zero lines and styling
    fig.update_xaxes(
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        tickcolor='black',
        showline=True,
        linecolor='black',
        linewidth=1,
        zeroline=False
    )

    fig.update_yaxes(
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        tickcolor='black',
        showline=True,
        linecolor='black',
        linewidth=1,
        zeroline=False
    )

    fig.update_traces(
        hovertemplate="Repeat Count=%{x}<br>Allele Count=%{y}<extra></extra>"
    )

    all_x = repeat_counts + [
        benign_min or 0, benign_max or 0,
        pathogenic_min or 0, pathogenic_max or 0
    ]
    fig.update_layout(xaxis_range=[0, max(all_x) + 10])

    x_range = fig.layout.xaxis.range
    x_span = x_range[1] - x_range[0]

    label_positions = {}
    if benign_min is not None and benign_max is not None:
        add_range_marker_or_line(
            fig, benign_min, benign_max,
            label_id="normal",
            label_positions=label_positions,
            chart_width_px=fig.layout.width,
            x_span=x_span
        )
    if pathogenic_min is not None and pathogenic_max is not None:
        add_range_marker_or_line(
            fig, pathogenic_min, pathogenic_max,
            label_id="pathogenic",
            label_positions=label_positions,
            chart_width_px=fig.layout.width,
            x_span=x_span
        )

    # Save plots as HTML
    plot_filename = f"{gene}_{chrom}_{pos}_allele_dist.html"
    fig.write_html(os.path.join(OUTPUT_DIR, plot_filename))
print(f"Saved plots in {OUTPUT_DIR}")
