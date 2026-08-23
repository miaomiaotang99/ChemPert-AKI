############################################################
# Multi-database disease-associated gene integration
#
# Description:
#   Integrate disease-related genes from multiple databases
#   and identify high-confidence genes supported by
#   multiple independent sources.
#
# Databases:
#   - GeneCards
#   - GWAS Catalog
#   - MalaCards
#   - OMIM
#   - Open Targets
#
# Strategy:
#   Genes appearing in >=3 databases are considered
#   high-confidence disease-associated genes.
#
# Input:
#   Gene lists from different databases
#
# Output:
#   Integrated gene matrix
#   Overlap visualization
#   High-confidence gene list
#
# Author:
#   Weidong Huang
#
############################################################