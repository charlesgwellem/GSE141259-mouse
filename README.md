# Exploratory Data Analysis of Mouse Dataset GSE141259

## **Overview**

This project involves the analysis of the **GSE141259** mouse dataset, focusing on exploratory data analysis (EDA) techniques. The analysis includes **data normalization**, **clustering**, **cell-cell communication** studies using **CellChat**, and **cell subtype characterization**. The goal is to gain insights into the cellular composition and interactions in the dataset through advanced computational techniques.

## **Project Structure**

### **Key Analyses Performed:**

1. **Data Normalization:**
   - The raw gene expression data was normalized to ensure that cells and genes are on a comparable scale. This step is essential for reliable downstream analysis, such as clustering and differential expression analysis.
   
2. **Clustering:**
   - **Unsupervised clustering** was performed to identify distinct cell populations in the dataset. Clustering allows us to categorize cells based on their gene expression profiles, highlighting the heterogeneity present in the sample.
   
3. **Cell-Cell Communication Analysis with CellChat:**
   - **CellChat**, a computational tool, was used to infer cell-cell communication networks. This analysis identifies the signaling pathways that cells use to interact with each other, revealing potential cellular crosstalk and communication patterns within the tissue.

4. **Cell Subtype Characterization:**
   - Identified clusters were further analyzed to **characterize cell subtypes** based on marker genes and functional annotations. This step helps in understanding the specific roles of different cell populations within the tissue.

### **Data Source:**
- **Accession Number:** [GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259)  
- The dataset contains gene expression profiles for mouse samples, which have been used in this study to explore cellular heterogeneity and communication.

---

## **Installation and Setup**

### **Required Packages**

This analysis was performed using **R** and the following R packages:

- **Seurat**: For data normalization, clustering, and visualization.
- **CellChat**: To infer cell-cell communication networks.
- **dplyr**: For data wrangling.
- **ggplot2**: For visualization.
- **patchwork**: For combining multiple plots.

You can install the required packages using the following code:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
BiocManager::install("CellChat")
```

### **Running the Analysis**

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/GSE141259-mouse-EDA.git
   cd GSE141259-mouse-EDA
   ```

2. **Run the RMarkdown File**:
   Open the `GSE141259_analysis.Rmd` file in RStudio and knit the document to generate the analysis report.

3. **Explore the Results**:
   The report will contain visualizations of the clustering, cell communication networks, and cell subtype characterizations.

---

## **Key Steps and Results**

### 1. **Data Normalization**
   The raw gene expression data from the GSE141259 dataset was normalized using **log normalization** to scale the gene expression values across cells.

   **Result:** 
   - The normalized data was used in subsequent clustering and communication analyses.

### 2. **Clustering Analysis**
   Using **Seurat**, clustering was performed based on the principal components derived from the normalized data. This step categorized the cells into distinct populations based on gene expression similarity.

   **Result:** 
   - The cells were grouped into clusters representing various cell types. These clusters were visualized using **UMAP** and **t-SNE** plots to capture the cellular diversity in the dataset.

### 3. **Cell-Cell Communication with CellChat**
   **CellChat** was employed to analyze intercellular communication networks. It inferred the signaling pathways active between the identified clusters, providing insights into how different cell types communicate with one another.

   **Result:**
   - A detailed communication network was established, showing the signaling interactions and pathways that connect the different cell types. Pathways such as **Wnt**, **Notch**, and **TNF** signaling were investigated in detail.

### 4. **Cell Subtype Characterization**
   Each cluster was characterized based on the expression of known **marker genes**. This step involved identifying specific cell subtypes such as immune cells, epithelial cells, and stromal cells based on the marker profiles.

   **Result:**
   - Subtypes of cells were assigned based on known markers, and their roles were inferred through functional enrichment analysis and literature review.

---

## **Visualizations**

- **UMAP/t-SNE Clustering**: Visualizes the clusters of cells in two dimensions, highlighting the separation between different cell types.
- **Cell Communication Network**: A graphical representation of the inferred cell-cell communication network, showing the active signaling pathways between clusters.
- **Heatmaps and Violin Plots**: Used to characterize and visualize the expression of marker genes across different cell types.

---

## **Conclusion**

This exploratory analysis of the **GSE141259-mouse** dataset reveals the complex cellular heterogeneity and cell-cell communication within the tissue. The use of **CellChat** provides insight into how cells interact through signaling pathways, while clustering and subtype characterization highlight the diversity of cell populations.

---

## **Session Info**

To view the R session information and the environment used for this analysis, run the following command:

```r
sessionInfo()
```
