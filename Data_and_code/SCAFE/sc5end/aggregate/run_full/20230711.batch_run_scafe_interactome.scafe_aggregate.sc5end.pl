#!/usr/bin/perl -w
use strict;

my $baseDir='/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/aggregate/run_full';
my $run_script = '/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/script/00_run_full.sh';
my $ctss_list_path = '/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20230614/neuron/sc5end/aggregate/00_iPSC_NSC_Neuron.aggregate_list.tsv';
my $genome = "hg38.gencode_v39";
my $sbatch_dir = "$baseDir/00_sbatch/";
my $aggregate_list_dir = "$baseDir/00_aggregate_list/";
my $scope_info_hsh_ref = {};
system "mkdir -pm 755 $baseDir";
system "mkdir -pm 744 $sbatch_dir";
system "mkdir -pm 744 $aggregate_list_dir";
system "cp $0 $baseDir";

my $start_time_log_path = "$baseDir/00_start.time.log.txt";
my $finish_time_log_path = "$baseDir/00_finish.time.log.txt";
system"echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system"echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

open (INCTSSLIST, "<", $ctss_list_path);
while (<INCTSSLIST>) {
	chomp;
	my ($libID, $all_ctss_path, $ung_ctss_path) = split /\t/;
	my ($CAGE_type, $cell_type)  = split /\./, $libID;
	my $cell_type_scope = "sc5end.$cell_type";
	my $all_scope = "sc5end.iPSC_NSC_Neuron";
	foreach my $scopeID ($cell_type_scope, $all_scope) {
		$scope_info_hsh_ref->{$scopeID}{$libID} = [$all_ctss_path, $ung_ctss_path];
	}
}
close INCTSSLIST;

foreach my $scopeID (sort keys %{$scope_info_hsh_ref}) {
	my $lib_list_path = "$aggregate_list_dir/00_".$scopeID.".aggregate_list.tsv";
	open AGRLIST, ">", $lib_list_path;
	foreach my $libID (sort keys %{$scope_info_hsh_ref->{$scopeID}}) {
		my ($all_ctss_path, $ung_ctss_path) = @{$scope_info_hsh_ref->{$scopeID}{$libID}};
		print AGRLIST join "", (join "\t", ($libID, $all_ctss_path, $ung_ctss_path)), "\n";
	}
	close AGRLIST;

	my $outputPrefix = $scopeID;
	my $sbatch_sh_path = "$sbatch_dir/$outputPrefix.sbatch.cmd.sh";
	my $sbatch_stderr_path = "$sbatch_dir/$outputPrefix.sbatch.stderr.txt";
	my $sbatch_stdout_path = "$sbatch_dir/$outputPrefix.sbatch.stdout.txt";

	open (SBATCHCMD, ">", $sbatch_sh_path);
	print SBATCHCMD "#!/bin/bash\n";
	print SBATCHCMD "echo \"$scopeID is started at \$(date)\" >>$start_time_log_path\n";
	print SBATCHCMD "run_script=\"$run_script\"\n";
	print SBATCHCMD "lib_list_path=\"$lib_list_path\"\n";
	print SBATCHCMD "genome=\"$genome\"\n";
	print SBATCHCMD "baseDir=\"$baseDir\"\n";
	print SBATCHCMD "outputPrefix=\"$outputPrefix\"\n";
	print SBATCHCMD 'sh $run_script $lib_list_path $genome $baseDir $outputPrefix'."\n";
	print SBATCHCMD "echo \"$scopeID is finished at \$(date)\" >>$finish_time_log_path\n";
	close SBATCHCMD;

	print "$scopeID sbatch submitted.\n";
	system "sbatch -c 10 -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
}
