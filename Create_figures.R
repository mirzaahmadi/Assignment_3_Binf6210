# Install required packages 
install.packages("ggVennDiagram")

# Load required libraries
library(ggVennDiagram)
library(ggplot2)



# Read data from plasmidfinder output
spades_plasmidfinder <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/spades_abricate_output/Spades_plasmidfinder_output.txt", header=TRUE)
flye_plasmidfinder <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/flye_abricate_output/Flye_plasmidfinder_output.txt", header=TRUE)
unicycler_plasmidfinder <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/unicycler_abricate_output/Unicycler_plasmidfinder_output.txt", header=TRUE)

# Read data from resfinder output
spades_resfinder <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/spades_abricate_output/Spades_resfinder_output.txt", header=TRUE)
flye_resfinder <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/flye_abricate_output/Flye_resfinder_output.txt", header=TRUE)
unicycler_resfinder <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/unicycler_abricate_output/Unicycler_resfinder_output.txt", header=TRUE)

# Read data from VFDB output
spades_VFDB <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/spades_abricate_output/Spades_vfdb_output.txt", header=TRUE)
flye_VFDB <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/flye_abricate_output/Flye_vfdb_output.txt", header=TRUE)
unicycler_VFDB <- read.delim("C:/Users/mamah/OneDrive - University of Guelph/Desktop/abricate_analysis_outputs/unicycler_abricate_output/Unicycler_vfdb_output.txt", header=TRUE)



# Extract unique gene lists for plasmids, resistance and virulence genes
spades_plasmids <- unique(spades_plasmidfinder$GENE)
flye_plasmids <- unique(flye_plasmidfinder$GENE)
unicycler_plasmids <- unique(unicycler_plasmidfinder$GENE)

spades_res_genes <- unique(spades_resfinder$GENE)
flye_res_genes <- unique(flye_resfinder$GENE)
unicycler_res_genes <- unique(unicycler_resfinder$GENE)

spades_virulence_genes <- unique(spades_VFDB$GENE)
flye_virulence_genes <- unique(flye_VFDB$GENE)
unicycler_virulence_genes <- unique(unicycler_VFDB$GENE)



# Create lists for Venn diagrams
plasmid_sets <- list(
  Spades = spades_plasmids,
  Flye = flye_plasmids,
  Unicycler = unicycler_plasmids
)

resistance_sets <- list(
  Spades = spades_res_genes,
  Flye = flye_res_genes,
  Unicycler = unicycler_res_genes
)

virulence_sets <- list(
  Spades = spades_virulence_genes,
  Flye = flye_virulence_genes,
  Unicycler = unicycler_virulence_genes
)



# Plot Venn diagram for plasmids
ggVennDiagram(plasmid_sets) + 
  scale_fill_gradient(low = "white", high = "blue") + 
  ggtitle("Plasmid Genes") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Larger, bold, centered title
    text = element_text(size = 20, face = "bold")  # Increases text size for categories
  )

# Plot Venn diagram for resistance genes
ggVennDiagram(resistance_sets) + 
  scale_fill_gradient(low = "white", high = "red") + 
  ggtitle("Resistance Genes") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  
    text = element_text(size = 20, face = "bold")  
  )

# Plot Venn diagram for virulence genes
ggVennDiagram(virulence_sets) + 
  scale_fill_gradient(low = "white", high = "green") + 
  ggtitle("Virulence Genes") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  
    text = element_text(size = 20, face = "bold")  
  )
