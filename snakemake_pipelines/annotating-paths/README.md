# Annotate the GFA files with P and W lines

## Basic Pipeline

- Align HPRC assemblies back to the graph
- Add BO and NO tags (initial tag) through the component decomposition of the GFA
- Sort the GAF file based on initial BO and NO tags
- Use the sorted GAF file to find the Allele Traversals in the graph
- Retag nodes with NO tag based on allele traversal
- Produce VCFs with the AT info (can be done in one of the previous steps)
