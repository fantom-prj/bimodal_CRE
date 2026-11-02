#!/usr/bin/perl -w
use strict;

my $outDir = "/analysisdata/fantom6/Interactome/single_cell_wallace/2026_10_12_Wall_RITM0039305_CAGE_iPSC_TEAD24_NFYAB_KD/SCAFE/bam_to_ctss";
my $run_scafe_script_path = '/home/yip/anaconda3/envs/scafe-env/bin/scafe.tool.bk.bam_to_ctss';
my $in_lib_list_path = "/analysisdata/fantom6/Interactome/single_cell_wallace/2026_10_12_Wall_RITM0039305_CAGE_iPSC_TEAD24_NFYAB_KD/SCAFE/bam_path.tsv";
my $genome = "hg38.gencode_v39";

my $sbatch_dir = "$outDir/00_sbatch/";
system "mkdir -pm 744 $sbatch_dir";

system "mkdir -pm 755 $outDir";
system "cp $0 $outDir";

my $start_time_log_path = "$outDir/00_start.time.log.txt";
my $finish_time_log_path = "$outDir/00_finish.time.log.txt";
system"echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system"echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

my $lib_info_hsh_ref = {};

open BAMLIST, "<", $in_lib_list_path;
while (<BAMLIST>) {
	chomp;
	next if $_ =~ m/^#/;
	my ($libID, $bamPath) = split /\t/;
	
	die "bamPath does not exists for $libID\n" if not -s $bamPath;

	$lib_info_hsh_ref->{$libID}{'bamPath'} = $bamPath;
}
close BAMLIST;

my $num_lib = keys %{$lib_info_hsh_ref};
print "$num_lib libraries read\n";

foreach my $libID (sort keys %{$lib_info_hsh_ref}) {
	my $bamPath = $lib_info_hsh_ref->{$libID}{'bamPath'};
	my $outputPrefix = $libID;
	my $max_thread = 5;
	my $cmd = join " ", (
		"$run_scafe_script_path",
		"--TSS_mode=softclip",
		"--bamPath=$bamPath",
		"--unencoded_G_upstrm_nt=3",
		"--max_thread=$max_thread",
		"--min_MAPQ=225",
		"--genome=$genome",
		"--max_softclip_length=3",
		"--outputPrefix=$outputPrefix",
		"--outDir=$outDir",
	);
		
	my $sbatch_sh_path = "$sbatch_dir/$outputPrefix.sbatch.cmd.sh";
	my $sbatch_stderr_path = "$sbatch_dir/$outputPrefix.sbatch.stderr.txt";
	my $sbatch_stdout_path = "$sbatch_dir/$outputPrefix.sbatch.stdout.txt";

	open (SBATCHCMD, ">", $sbatch_sh_path) or die $!;

	print SBATCHCMD "#!/bin/bash\n";
	print SBATCHCMD "set -euo pipefail\n\n";

	# Activate SCAFE conda environment
	print SBATCHCMD "export PATH=/anaconda3/envs/scafe-env/bin:\$PATH\n";
	print SBATCHCMD "export CONDA_PREFIX=/anaconda3/envs/scafe-env\n\n";

	# Run SCAFE
	print SBATCHCMD "echo \"$libID is started at \$(date)\" >> $start_time_log_path\n";
	print SBATCHCMD "$cmd\n";
	print SBATCHCMD "echo \"$libID is finished at \$(date)\" >> $finish_time_log_path\n";

	close SBATCHCMD;

	print "$libID sbatch submitted.\n";
	system "sbatch -c $max_thread -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
}
