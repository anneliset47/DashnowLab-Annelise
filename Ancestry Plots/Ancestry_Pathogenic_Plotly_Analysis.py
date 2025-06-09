import pandas as pd
import re
import plotly.express as px
import os

# Load raw allele data
df = pd.read_csv('/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/locus_sample_population_map_expanded_AllAleles.csv')

# Ensure 'is_pathogenic' is boolean
if 'is_pathogenic' in df.columns:
    df['is_pathogenic'] = df['is_pathogenic'].astype(str).str.lower().map({'true': True, 'false': False})

# Group by disease, gene, and population_description
agg_df = df.groupby(['disease', 'gene', 'population_description'])

# Count of pathogenic
pathogenic_count = agg_df['is_pathogenic'].sum().reset_index(name='pathogenic_count')

# Count of total
total_count = agg_df.size().reset_index(name='total_count')

# Merge counts
df_agg = pd.merge(pathogenic_count, total_count, on=['disease', 'gene', 'population_description'])

# Define a mapping to rename population descriptions for consistency
population_rename_map = {
    'Unknown': 'Unknown',
    'Finnish in Finland': 'Finnish',
    'Han Chinese South, China': 'East Asian',
    'Puerto Rican in Puerto Rico': 'Admixed American',
    'Colombian in Medellin, Colombia': 'Admixed American',
    'African Caribbean in Barbados': 'African/African American',
    'Peruvian in Lima, Peru': 'Admixed American',
    'Kinh in Ho Chi Minh City, Vietnam': 'East Asian',
    'Gambian in Western Division Ð Mandinka': 'African/African American',
    'Punjabi in Lahore, Pakistan': 'South Asian',
    'Esan in Nigeria': 'African/African American',
    'Mende in Sierra Leone': 'African/African American',
    'Sri Lankan Tamil in the UK': 'South Asian',
    'Bengali in Bangladesh': 'South Asian',
    'Han Chinese in Beijing, China': 'East Asian',
    'Japanese in Tokyo, Japan': 'East Asian',
    'Luhya in Webuye, Kenya': 'African/African American',
    'Toscani in Italia': 'European (non Finnish)'
}

# Clean and map population names
df_agg['population_description'] = df_agg['population_description'].map(
    population_rename_map
).fillna(df_agg['population_description'])

# Group again in case multiple entries mapped to same population
df_agg = df_agg.groupby(['disease', 'gene', 'population_description'])[['pathogenic_count', 'total_count']].sum().reset_index()

# Add percentage
df_agg['percentage'] = df_agg['pathogenic_count'] / df_agg['total_count'] * 100

# Save processed data (optional)
# df_agg.to_csv('/Users/annelisethorn/Documents/Anschutz/Code/Ancestry/CSVs/PathogenicTotalCounts_UpdatedNames_PercentColumn.csv', index=False)

# File paths
OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Ancestry_Plots/Updated_C2_2"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Plotting function
def create_horizontal_bar_plot(df, gene, disease):
    filtered_df = df[(df['gene'] == gene) & (df['disease'] == disease)].copy()
    if filtered_df.empty:
        print(f"No data available for gene: {gene} and disease: {disease}")
        return None

    # Round display percentage
    filtered_df['percentage_display'] = filtered_df['percentage'].round(1).apply(lambda x: int(x) if x == 0.0 else x)

    fig = px.bar(
        filtered_df,
        x='percentage',
        y='population_description',
        orientation='h',
        color_discrete_sequence=["blue"],
        title=f"Pathogenic Genotype Distribution<br><sup>{gene} - {disease}</sup>",
        labels={'percentage': 'Pathogenic Genotype (%)', 'population_description': 'Population'},
        text='percentage_display'
    )

    fig.update_traces(
        texttemplate='%{text}',
        textposition='outside',
        cliponaxis=False,
        customdata=filtered_df[['pathogenic_count', 'total_count']].values,
        hovertemplate=(
            "Population: %{y}<br>"
            "Pathogenic Genotype (%): %{x}<br>"
            "# Pathogenic: %{customdata[0]}<br>"
            "# Total: %{customdata[1]}<extra></extra>"
        )
    )

    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=160),
        title_font_size=18,
        title_font=dict(weight="bold"),
        title_font_color="black",
        xaxis_title="Pathogenic Genotypes (%)",
        yaxis_title="",
        bargap=0.5,
        plot_bgcolor='white',
        xaxis=dict(range=[0, 100], showline=True, linecolor='black'),
        yaxis=dict(ticks='outside', showline=True, linecolor='black'),
        showlegend=False
    )

    return fig

# Loop and generate plots
for gene in df_agg['gene'].unique():
    diseases = df_agg[df_agg['gene'] == gene]['disease'].unique()
    for disease in diseases:
        filtered_df = df_agg[(df_agg['gene'] == gene) & (df_agg['disease'] == disease)].copy()
        fig = create_horizontal_bar_plot(filtered_df, gene, disease)
        if fig is None:
            continue
        safe_gene = re.sub(r'[\\/]', '_', gene)
        safe_disease = re.sub(r'[\\/]', '_', disease)
        plot_filename = f"{safe_gene}_{safe_disease}_ancestry_plot.html"
        fig.write_html(os.path.join(OUTPUT_DIR, plot_filename))

print("Saved plots")

