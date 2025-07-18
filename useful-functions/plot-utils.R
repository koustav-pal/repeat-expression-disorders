
output_prefix = function(){
    return(format(Sys.Date(),"%Y-%m-%d"))
}


# colour definitions 

region_colours <- c("exon" = "#EDF060", "intron" = "#F06449")

orientation_colours <- c("aligned" = "#0C6291", "opposed" = "#7E1946")

condition_colours <- c("wildtype" = "#ED254E", "mutant" = "#279AF1")
# colours from https://coolors.co/403f4c-ff7e6b-f9dc5c-3185fc-48e5c2
feature_type_colours <- c("AT-rich" = "#FF7E6B", "Gquad" = "#F9DC5C", "Retrotransposon" = "#3185FC")


molecule_colours <- c("#000000", "#677DB7", "#9CA3DB", "#9EC1A3",
	"#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
	"#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
	"#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
names(molecule_colours) <- c("c9orf72", "fus", "sod1", "tardbp", 
	"vcp_cyt_d0", "vcp_nuc_d0", "vcp_cyt_d3", "vcp_nuc_d3",
	"vcp_cyt_d7", "vcp_nuc_d7", "vcp_cyt_d14", "vcp_nuc_d14",
	"vcp_cyt_d22", "vcp_nuc_d22", "vcp_cyt_d35", "vcp_nuc_d35")

# https://coolors.co/8338ec-f71735-41ead4-3a86ff-ff9f1c
mol_name_colours <- c("#8338EC", "#F71735", "#41EAD4", "#3A86FF", "#FF9F1C")
names(mol_name_colours) <- c("c9orf72", "sod1", "fus", "tardbp", "vcp")