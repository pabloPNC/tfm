import pandas as pd
import numpy as np

file_sheet_path = r".\tcga_data\gdc_sample_sheet.2024-06-11.tsv"
clinical_sheet_path = r".\tcga_data\clinical.cart.2024-06-11\clinical.tsv"

file_sheet = pd.read_csv(file_sheet_path, sep='\t')
clinical_sheet = pd.read_csv(clinical_sheet_path, sep='\t')

# Modify column names
file_sheet.rename(columns=lambda x: x.replace(
    " ", "_").lower(), inplace=True)


# Subset with gleason grades
clinical_sheet_filtered = clinical_sheet[[
    "case_id", "case_submitter_id", "primary_gleason_grade", "secondary_gleason_grade"]]

# Merge dfs
file_clinical_sheet = file_sheet.merge(
    clinical_sheet_filtered,
    how="left",
    left_on="case_id",
    right_on="case_submitter_id"
)

# Transform values to numeric
file_clinical_sheet = file_clinical_sheet.replace(
    to_replace=["'--", 'Pattern 2', 'Pattern 3', 'Pattern 4', 'Pattern 5'],
    value=[np.nan, 2, 3, 4, 5]
)

# Filter samples without gleason score
file_clinical_sheet = file_clinical_sheet[
    pd.notna(file_clinical_sheet["secondary_gleason_grade"]) |
    pd.notna(file_clinical_sheet["primary_gleason_grade"])
]

# Fix: set gleason_grades 'Solid tissue normal' to 0
file_clinical_sheet.loc[file_clinical_sheet['sample_type'] ==
                        "Solid Tissue Normal", 'primary_gleason_grade'] = 0
file_clinical_sheet.loc[file_clinical_sheet['sample_type'] ==
                        "Solid Tissue Normal", 'secondary_gleason_grade'] = 0

# Add gleason_score = primary_gleason_score + secondary_gleason_score
file_clinical_sheet['gleason_score'] = file_clinical_sheet['primary_gleason_grade'] + \
    file_clinical_sheet['secondary_gleason_grade']


# Delete duplicated rows
file_clinical_sheet = file_clinical_sheet.drop_duplicates()


# Export_filtered_data
file_clinical_sheet.to_csv("./file_clinical_sheet.csv")
