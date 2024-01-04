# HIV-5L-defect

all del.R:
- This script identifies all deletions from the input fasta file named "5PSD_master_040422_tilda_uniq.fas" (containing unique and aligned 5'L-defective sequences), and then outputs the results to "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx".

0405 all del together (for plot).R:
- This script plots all detected deletions using "5PSD_master_040522_hxb2.fas" (containing all unique 5'L sequences aligned to hxb2) and "0405_5psd_GagIntact.aa.fas" (containing all 5'L sequences that have intact Gag region), and outputs the plot in a pdf.

0425 all MH (consensus with spacer).R:
- This script detects all the microhomology stretches at the deletion junctions in the 5'L sequences using "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx" (containing the coordinates of all detected deletions in 5'L sequences when aligned to consensus sequence for HIV subtype-B) and "5PSD_master_040422_tilda_uniq.fas", and then outputs the results to "5PSD_master_042522_tilda_uniq_allDel_cleared.xlsx".
- have to manually remove the tilde

0419 all HP.R:
- This script detects all polymeric regions in the 5'L sequences using "5PSD_master_040422_tilda_uniq.fas" and "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx", and then outputs the results to "5PSD_master_041922_tilda_uniq_allHP.xlsx".
- have to manually change column names

0428 check real poly.R:
- This script applies several filters to get all desired conditions of polymeric regions based on the known info of the polymeric state of the 5'L region.

0415 MH&del freq plots.R:
- This script maps the frequency of appearance of the 5'end, 3'end, and microhomology of each deletion junction using "5PSD_master_041422_tilda_uniq_consensus_WithSpacer_allMH_updated.xlsx" and "5PSD_master_041422_tilda_uniq_allDel_cleared.xlsx" and outputs the plots.

0429 find mH stretch everywhere.R
- This script detects all the microhomology stretches in the 5'L sequences (i.e. not just the ones at the deletion junctions) using "5PSD_master_041922_MH stretch counts.xlsx" and "ConsBandHXB2.fas" (HIV-1 subytype B consensus sequence), and outputs "5PSD_master_042522_tilda_uniq_STReverywhere_MH.xlsx".

0425 find MH stretch plot.R:
- This script plots the detected microhomology stretches using "5PSD_master_042922_tilda_uniq_STReverywhere.xlsx" (needs to have the first sheet with selected microhomology stretch of interest to map) and "PSD_master_042522_tilda_uniq_consensus_WithSpacer_allMH.xlsx" and "5PSD_master_041222_tilda_uniq_consensus.fas". It outputs the plot in pdf and outputs a table listing the selected microhomology stretch of interest "5PSD_master_042622_tilda_uniq_selected STR table.xlsx".

0502 find poly stretch everywhere.R:
- This script detects all the polymeric regions in the 5'L sequences (i.e. not just the ones at the deletion junctions) using "5PSD_master_042922_tilda_uniq_checkedHP.xlsx" (containing all the polymeric regions at the deletion junctions of 5'L sequences), "5PSD_master_050222_tilda_uniq_polySTRcount.xlsx" (containing the count of appearance of each unique polymeric region in 5'L sequences), "5PSD_master_050222_tilda_uniq_3endPolyBetw750to770.xlsx" (containing the counts of appearance of polymeric region of coordinates 750-770, which is identified as highly polymeric and variable), and "5PSD_master_040422_tilda_uniq.fas", and outputs "5PSD_master_050922_tilda_uniq_STReverywhere_poly.xlsx".
- 5PSD_master_050222_tilda_uniq_polySTRcount.xlsx is obtained by making a pivot table using the "true poly" in the "checkedHP.xlsx" file (CAUTION: remember to concatenate to get the "stretch" column and change the poly_len accordingly for the "corrected_poly_len" column, if applicable), so that we get all the poly stretches and their counts of occurrences. We copy and paste the pivot table into "polySTRcount.xlsx" file, and use the R script to scan through the HXB2 to see the occurrence of theses stretches everywhere. The chosen range of HXB2 coordinates should be corresponding to the range of consensus coordinats chosen previously for "finding MH everywhere". AT-table is made using the "3end to 5end.xlsx" file (which is made from filtering the true-poly with respect to their del junctions) is to avoid the null detection in HXB2 of certain stretches due to indels. 

0504 find poly stretch plot.R:
- This script plots the detected polymeric regions of the 5'L sequences using "5PSD_master_050922_tilda_uniq_STReverywhere_poly.xlsx" and "5PSD_master_042922_tilda_uniq_checkedHP.xlsx", and outputs the plot in a pdf.
- Using "5PSD_master_050922_tilda_uniq_STReverywhere_poly.xlsx" file, we check for the stretches that have AT region, combine them into one category called "AAAAATTTT region" and sum up their occurrences in the "refined" sheet (CAUTION: the 1st page is still the RAW data). We then get a plot in excel. Next, using this R script, we scan through the raw data of this excel file (CAUTION: check the order of the stretches in the raw data first, and make sure that the non-AT-region stretches [the ones that have >1 for occurrence everywhere] are on the top, and the AT-region stretches are at the bottom. Script might need to be changed accordingly regarding the number of rows that are non-AT-region stretches vs AT-region stretches.) For the "checkedHP.xlsx" file, make sure that the "all true poly" sheet is on the first page. Run the script and get the plot (NOTE: the AAAAATTTT regions are all on one row no matter for their occurrences in 5PSD or in HXB2).

0523 specific positions.R:
- This script scans for the presence of certain RNA features in 5'L sequences and identifies whether the RNA feature is complete, partial, or absent, using "5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx" and "5PSD_master_040422_tilda_uniq.fas", and outputs the results in "5PSD_master_0524_SpecPos.xlsx".
