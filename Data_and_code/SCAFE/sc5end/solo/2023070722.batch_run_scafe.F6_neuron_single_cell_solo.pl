#!/usr/bin/perl -w
use strict;

my $outDir = "/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/solo/";
my $in_lib_list_path = "/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/cellranger_dir.txt";
my $scafe_sc_solo_path = "/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/v1.0.1/scripts/scafe.workflow.sc.solo";
my $detect_TS_oligo = 'no';
my $genome = "hg38.gencode_v40";

system "mkdir -pm 755 $outDir";
system "cp $0 $outDir";

my $start_time_log_path = "$outDir/00_start.time.log.txt";
my $finish_time_log_path = "$outDir/00_finish.time.log.txt";
system"echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system"echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

my $lib_info_hsh_ref = {};

open CTSSLIST, "<", $in_lib_list_path;
while (<CTSSLIST>) {
	chomp;
	next if $_ =~ m/^#/;
	my ($libID, $cellranger_out_dir) = split /\t/;
	
	my $run_bam_path = $cellranger_out_dir.'/outs/possorted_genome_bam.bam';
	my $run_cellbarcode_path = $cellranger_out_dir.'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz';
	
	die "run_bam_path $run_bam_path does not exists for $libID\n" if not -s $run_bam_path;
	die "run_cellbarcode_path does not exists for $libID\n" if not -s $run_cellbarcode_path;

	$lib_info_hsh_ref->{$libID}{'run_bam_path'} = $run_bam_path;
	$lib_info_hsh_ref->{$libID}{'run_cellbarcode_path'} = $run_cellbarcode_path;
}
close CTSSLIST;

my $num_lib = keys %{$lib_info_hsh_ref};
print "$num_lib libraries read\n";

foreach my $libID (sort keys %{$lib_info_hsh_ref}) {
	my $run_bam_path = $lib_info_hsh_ref->{$libID}{'run_bam_path'};
	my $run_cellbarcode_path = $lib_info_hsh_ref->{$libID}{'run_cellbarcode_path'};
	my $run_tag = $libID;
	my $run_outDir = "$outDir/$run_tag/";
	my $force_rerun = 'no';
	my $max_thread = 10;
	my $cmd = join " ", (
		"perl $scafe_sc_solo_path",
		"--TSS_mode=softclip",
		"--run_bam_path=$run_bam_path",
		"--run_cellbarcode_path=$run_cellbarcode_path",
		"--max_thread=$max_thread",
		"--detect_TS_oligo=auto",
		"--genome=$genome",
		"--run_tag=$run_tag",
		"--run_outDir=$run_outDir",
	);

	my $sub_outDir = "$outDir/00_sbatch/";
	system "mkdir -pm 744 $sub_outDir";
	
	my $sbatch_sh_path = "$sub_outDir/$run_tag.sbatch.cmd.sh";
	my $sbatch_stderr_path = "$sub_outDir/$run_tag.sbatch.stderr.txt";
	my $sbatch_stdout_path = "$sub_outDir/$run_tag.sbatch.stdout.txt";

	open (SBATCHCMD, ">", $sbatch_sh_path);
	print SBATCHCMD "#!/bin/bash\n";
	print SBATCHCMD "echo \"$libID is started at \$(date)\" >>$start_time_log_path\n";
	print SBATCHCMD "$cmd;\n";
	print SBATCHCMD "echo \"$libID is finished at \$(date)\" >>$finish_time_log_path\n";
	close SBATCHCMD;

	print "$libID sbatch submitted.\n";
	system "sbatch -c $max_thread -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
}
