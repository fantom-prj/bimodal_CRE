#!/bin/sh

cd /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/
for file in $(ls -d TEAD*/); do

	cd /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/BINDetect_compare_output/${file%/}/beds/

# Visualize the difference in footprints between two conditions for all accessible sites
	TOBIAS PlotAggregate \
		--TFBS "${file%/}_all.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/iPS_ATACorrect_test/iPS_corrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_corrected.bw" \
		--output "${file%/}_footprint_iPS_NSC.all.pdf" \
		--share_y both \
		--plot_boundaries \
		--signal-on-x
		
	TOBIAS PlotAggregate \
		--TFBS "${file%/}_all.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_corrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_corrected.bw" \
		--output "${file%/}_footprint_NSC_NRN.all.pdf" \
		--share_y both \
		--plot_boundaries \
		--signal-on-x
		
	TOBIAS PlotAggregate \
		--TFBS "${file%/}_all.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_corrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_corrected.bw" \
		--output "${file%/}_footprint_iPS_NRN.all.pdf" \
		--share_y both \
		--plot_boundaries \
		--signal-on-x


# Visualize the difference in footprints between two conditions exclusively for bound sites
	TOBIAS PlotAggregate \
		--TFBS "${file%/}_iPS_footprints_bound.bed" "${file%/}_NSC_footprints_bound.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/iPS_ATACorrect_test/iPS_corrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_corrected.bw" \
		--output "${file%/}_footprint_iPS_NSC.subset.pdf" \
		--share_y both \
		--plot_boundaries

	TOBIAS PlotAggregate \
	--TFBS "NSC_footprints_bound.bed" "${file%/}_NRN_footprints_bound.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_corrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_corrected.bw" \
		--output "${file%/}_footprint_NSC_NRN.subset.pdf" \
		--share_y both \
		--plot_boundaries

	TOBIAS PlotAggregate \
		--TFBS "${file%/}_iPS_footprints_bound.bed" "${file%/}_NRN_footprints_bound.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/iPS_ATACorrect_test/iPS_corrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_corrected.bw" \
		--output "${file%/}_footprint_iPS_NRN.subset.pdf" \
		--share_y both \
		--plot_boundaries


#Visualize the split of bound/unbound sites for one condition
	TOBIAS PlotAggregate \
		--TFBS "${file%/}_all.bed" "${file%/}_iPS_footprints_bound.bed" "${file%/}_iPS_footprints_unbound.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/iPS_ATACorrect_test/iPS_uncorrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/iPS_ATACorrect_test/iPS_expected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/iPS_ATACorrect_test/iPS_corrected.bw" \
		--output "${file%/}_footprint_iPS.pdf" \
		--share_y sites \
		--plot_boundaries

	TOBIAS PlotAggregate \
		--TFBS "${file%/}_all.bed" "${file%/}_NSC_footprints_bound.bed" "${file%/}_NSC_footprints_unbound.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_uncorrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_expected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NSC_ATACorrect_test/NSC_corrected.bw" \
		--output "${file%/}_footprint_NSC.pdf" \
		--share_y sites \
		--plot_boundaries

	TOBIAS PlotAggregate \
		--TFBS "${file%/}_all.bed" "${file%/}_NRN_footprints_bound.bed" "${file%/}_NRN_footprints_unbound.bed" \
		--signals "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_uncorrected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_expected.bw" "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/NRN_ATACorrect_test/NRN_corrected.bw" \
		--output "${file%/}_footprint_NRN.pdf" \
		--share_y sites \
		--plot_boundaries;
		
done




