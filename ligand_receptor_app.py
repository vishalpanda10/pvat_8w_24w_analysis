

import streamlit as st
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from plotnine import ggplot, theme, element_text

import liana as li
import tempfile

# Define the source and target labels
source_labels = ['Adipocytes', 'Endothelial_Cells', 'Fibroblasts', 'Immune_Cells',
                 'Mesothelial_Cells', 'Neuronal_Cells', 'Pericytes', 'SMCs']
target_labels = ['Adipocytes', 'Endothelial_Cells', 'Fibroblasts', 'Immune_Cells',
                 'Mesothelial_Cells', 'Neuronal_Cells', 'Pericytes', 'SMCs']

# Define the file paths
file_paths = {
    "pvat_8weeks_control_male": "pvat_8weeks_control_male.h5ad",
    "pvat_8weeks_control_female": "pvat_8weeks_control_female.h5ad",
    "pvat_8weeks_hf_male": "pvat_8weeks_hf_male.h5ad",
    "pvat_8weeks_hf_female": "pvat_8weeks_hf_female.h5ad",
    "pvat_24weeks_control_male": "pvat_24weeks_control_male.h5ad",
    "pvat_24weeks_control_female": "pvat_24weeks_control_female.h5ad",
    "pvat_24weeks_hf_male": "pvat_24weeks_hf_male.h5ad",
    "pvat_24weeks_hf_female": "pvat_24weeks_hf_female.h5ad"
}

# Load the datasets
datasets = {
    "8w_ctrl_male": sc.read_h5ad(file_paths["pvat_8weeks_control_male"]),
    "8w_ctrl_female": sc.read_h5ad(file_paths["pvat_8weeks_control_female"]),
    "8w_hf_male": sc.read_h5ad(file_paths["pvat_8weeks_hf_male"]),
    "8w_hf_female": sc.read_h5ad(file_paths["pvat_8weeks_hf_female"]),
    "24w_ctrl_male": sc.read_h5ad(file_paths["pvat_24weeks_control_male"]),
    "24w_ctrl_female": sc.read_h5ad(file_paths["pvat_24weeks_control_female"]),
    "24w_hf_male": sc.read_h5ad(file_paths["pvat_24weeks_hf_male"]),
    "24w_hf_female": sc.read_h5ad(file_paths["pvat_24weeks_hf_female"])
}

# Streamlit app
st.title("Ligand-Receptor Interactions")

# Layout for select boxes
col1, col2, col3 = st.columns(3)

with col1:
    source_cell = st.selectbox("Select Source Cell Type", source_labels)
with col2:
    target_cell = st.selectbox("Select Target Cell Type", target_labels)
with col3:
    dataset_option = st.selectbox("Select Dataset", list(datasets.keys()))

# Filter the data based on selections and weight
lr_df = datasets[dataset_option].uns['nichenet_lr_res']
filtered_df = lr_df[(lr_df['weight'] > 1.25)]
datasets[dataset_option].uns['nichenet_res_filtered'] = filtered_df

filtered_df_heatmap = lr_df[(lr_df['weight'] > 0.8)]

# Buttons for generating plots
#col1, col2 = st.columns(2)

if st.button("Generate Heatmap"):
    st.write("The heatmap visualizes the ligand-receptor interaction strengths between selected source and target cell types. Each cell represents the weight of interaction between a ligand and a receptor.")
    
    
    heatmap_df = filtered_df_heatmap[(filtered_df_heatmap['source'] == source_cell) & (filtered_df_heatmap['target'] == target_cell)]
    heatmap_df = heatmap_df.pivot_table(values='weight', index='ligand_complex', columns='receptor_complex', fill_value=0)

    fig, ax = plt.subplots(figsize=(40, 18))
    sns.heatmap(heatmap_df, annot=True, cmap='Blues', ax=ax, annot_kws={"size": 14})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=17)  
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=17)  
  
    ax.set_title('Ligand-Receptor Interaction Heatmap', fontsize=18)
    ax.set_xlabel('Receptor Complex')
    ax.set_ylabel('Ligand Complex')

    st.pyplot(fig)

if st.button("Generate Dotplot"):
    st.write("The dotplot illustrates the interactions between selected source and all target cell types.")
    
    dotplot = li.pl.dotplot(adata=datasets[dataset_option],
                            colour='weight',
                            size='weight',
                            source_labels=source_labels,
                            target_labels=target_labels,
                            figure_size=(20, 20),  
                            uns_key='nichenet_res_filtered')

    # Modify the ggplot object to rotate x-axis labels and increase label size
    dotplot += theme(
        axis_text_x=element_text(angle=45, hjust=1, size=12),
        axis_text_y=element_text(size=16)
    )

    with tempfile.NamedTemporaryFile(suffix=".png") as tmpfile:
        dotplot.save(tmpfile.name, width=20, height=20, dpi=300, limitsize=False)  
        col_center1, col_center2, col_center3 = st.columns([1, 2, 1])
        with col_center2:
            st.image(tmpfile.name, use_column_width=True)
