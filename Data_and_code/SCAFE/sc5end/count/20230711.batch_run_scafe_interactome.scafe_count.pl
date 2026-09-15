#!/usr/bin/perl -w
use strict;

my $baseDir='/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/count';
my $solo_dir= '/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/solo';
my $annotate_dir= '/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/merge/annotate';
my $resource_dir= '/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/resource';
my $scafe_dir = '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/v1.0.1/scripts/';
my $tag = 'sc5end.iPSC_NSC_Neuron';
my $ctss_scope_bed_path = "$annotate_dir/$tag/bed/$tag.cluster.coord.bed.gz";
my $genome = "hg38.gencode_v39";
my $sbatch_dir = "$baseDir/00_sbatch/";
my $scope_info_hsh_ref = {};
system "mkdir -pm 755 $baseDir";
system "mkdir -pm 744 $sbatch_dir";
system "cp $0 $baseDir";

my $countRegion_hsh_ref = {
	'tCRE' => "$annotate_dir/$tag/bed/$tag.CRE.coord.bed.gz",
	'aCRE_original' => "$resource_dir/peaks.merged.all.bed",
	'aCRE_500nt' => "$resource_dir/peaks.merged.all.500nt.bed",
};

my $start_time_log_path = "$baseDir/00_start.time.log.txt";
my $finish_time_log_path = "$baseDir/00_finish.time.log.txt";
system"echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system"echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

foreach my $lib (qw/iPS1 iPS2 NSC1 NSC2 Neuron1 Neuron2/) {
	my $ctss_bed_path = "$solo_dir/$lib/bam_to_ctss/$lib/bed/$lib.CB.ctss.bed.gz";
	my $cellBarcode_list_path = "$solo_dir/$lib/count/$lib/matrix/barcodes.tsv";
	
	foreach my $countRegion (keys %{$countRegion_hsh_ref}) {
		my $countRegion_bed_path = $countRegion_hsh_ref->{$countRegion};
		my $outputPrefix = $lib;

		my $sbatch_sh_path = "$sbatch_dir/$outputPrefix.$countRegion.sbatch.cmd.sh";
		my $sbatch_stderr_path = "$sbatch_dir/$outputPrefix.$countRegion.sbatch.stderr.txt";
		my $sbatch_stdout_path = "$sbatch_dir/$outputPrefix.$countRegion.sbatch.stdout.txt";

		open (SBATCHCMD, ">", $sbatch_sh_path);
		print SBATCHCMD "#!/bin/bash\n";
		print SBATCHCMD "echo \"$lib is started at \$(date)\" >>$start_time_log_path\n";
		print SBATCHCMD "$scafe_dir/scafe.tool.sc.count \\\n";
		print SBATCHCMD "--overwrite=yes \\\n";
		print SBATCHCMD "--genome=$genome \\\n";
		print SBATCHCMD "--ctss_scope_bed_path=$ctss_scope_bed_path \\\n";
		print SBATCHCMD "--countRegion_bed_path=$countRegion_bed_path \\\n";
		print SBATCHCMD "--cellBarcode_list_path=$cellBarcode_list_path \\\n";
		print SBATCHCMD "--ctss_bed_path=$ctss_bed_path \\\n";
		print SBATCHCMD "--outputPrefix=$outputPrefix \\\n";
		print SBATCHCMD "--outDir=$baseDir/$countRegion\n";
		print SBATCHCMD "echo \"$lib is finished at \$(date)\" >>$finish_time_log_path\n";
		close SBATCHCMD;

		print "$lib $countRegion sbatch submitted.\n";
		system "sbatch -c 1 -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
	}
}

