library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

my_graph <- grViz("
digraph workflow {
  graph [layout = dot, rankdir = TB, fontsize = 38]
  
  node [shape = rectangle, style = filled, fillcolor = LightGray, fontsize = 38, width = 3.5, height = 0.8]
  
  # === Nodes ===
  S1  [label = 'Collect run parameters']
  S2a [label = 'Merge paired-end reads', fillcolor = LightBlue]
  S2b [label = 'Merge raw fastq files', fillcolor = LightBlue]
  S3  [label = 'Filter reads by size']
  S4  [label = 'Demultiplex']
  S5  [label = 'Output demultiplexing report', fillcolor = White]
  S6a [label = 'Trim fwd primers from 5p end\n&\nTrim rev primers from 3p end', fillcolor = LightBlue]
  S6b [label = 'Trim fwd primers from 5p or 3p ends\n&\nTrim rev primers from opposite ends', fillcolor = LightBlue]
  S7  [label = 'Filter reads by size']
  S16 [label = 'Remove chimeras', fillcolor = Khaki]
  S8  [label = 'Cluster into ASVs']
  S9  [label = 'Filter ASVs by minimum read depth']
  S10 [label = 'Auto-trim ASVs', fillcolor = PaleGreen]
  S11 [label = 'Filter reads by size']
  S12 [label = 'Remove NUMTs\nand correct errors', fillcolor = PaleGreen]
  S13 [label = 'Cluster sample ASVs into run-wide ASVs']
  S14 [label = 'Identify run-wide ASVs']
  S15 [label = 'Denoise data']
  S17 [label = 'BIN match', fillcolor = PaleGreen]
  S18 [label = 'Reporting']
  
  edge [fontsize = 32]
  
  # === Initial pipeline ===
  S1 -> S2a [label = 'PE = Yes']
  S1 -> S2b [label = 'PE = No']
  S2a -> S3
  S2b -> S3
  S3 -> S4
  S4 -> S6a [label = 'PE = Yes']
  S4 -> S6b [label = 'PE = No']
  S6a -> S7
  S6b -> S7
  
  # === Chimera logic BEFORE clustering ===
  S7 -> S16 [label = 'amp_size >= 200 bp']
  S7 -> S8  [label = 'amp_size < 200 bp']
  S16 -> S8
  
  # === Continue pipeline ===
  S8 -> S9
  
  S9 -> S10 [label = 'marker = COI-5P\nAND\namp_size >= 500 bp']
  S10 -> S11
  S9 -> S11
  
  S11 -> S12 [label = 'marker = COI-5P\nAND\namp_size >= 200 bp']
  S11 -> S13
  S12 -> S13
  S13 -> S14 -> S15
  
  # === Bottom logic (simplified as requested) ===
  S15 -> S17 [label = 'marker = COI-5P\nAND\namp_size >= 300 bp']
  S15 -> S18
  S17 -> S18
  
  # Offshoot output
  S4 -> S5 [style = dashed, color = gray40]
  { rank = same; S4; S5 }
}
")

svg <- export_svg(my_graph)

writeLines(svg, "~/Downloads/workflow.svg")

rsvg_png("~/Downloads/workflow.svg", "~/Downloads/workflow.png", width = 2400, height = 3200)
