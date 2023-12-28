# HIV-5-leader-defect

all del.R:
- This script identifies all deletions from the input fasta file named "5PSD_master_040422_tilda_uniq.fas" (containing unique and aligned 5'L-defective sequences), and then outputs the results to "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx".

0405 all del together (for plot).R:
- This script plots all detected deletions using "5PSD_master_040522_hxb2.fas" (containing all unique 5'L sequences aligned to hxb2) and "0405_5psd_GagIntact.aa.fas" (containing all 5'L sequences that have intact Gag region), and outputs the plot in a pdf.

0425 all MH (consensus with spacer).R:
- This script detects all the microhomology stretches in the 5'L sequences using "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx" (containing the coordinates of all detected deletions in 5'L sequences when aligned to consensus sequence for HIV subtype-B) and "5PSD_master_040422_tilda_uniq.fas", and then outputs the results to "5PSD_master_042522_tilda_uniq_allDel_cleared.xlsx".
- have to manually remove the tilde

0419 all HP.R:
- This script detects all polymeric regions in the 5'L sequences using "5PSD_master_040422_tilda_uniq.fas" and "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx", and then outputs the results to "5PSD_master_041922_tilda_uniq_allHP.xlsx".
- have to manually change column names

0428 check real poly.R:
- This script applies several filters to get all desired conditions of polymeric regions based on the known info of the polymeric state of the 5'L region.
