#!/usr/bin/env perl

###############################################################################
# Nanopore-quality-filter.pl
#
# Author: Cyrill Hofer
# Created: 29.11.2023
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a standardized preprocessing
# workflow for raw Nanopore reads. It performs length/quality filtering, adapter
# trimming, and renaming, and produces NanoStat reports before and after
# processing.
#
# WORKFLOW STEPS
# -----------------------------------------------------------------------------
#  1) NanoStat report on raw reads
#  2) Filter reads by quality and length using chopper
#  3) Trim adapters using porechop_abi
#  4) Filter short reads again using chopper
#  5) Rename reads with BBMap rename.sh for consistent read IDs
#  6) NanoStat report on final cleaned reads
#  7) Cleanup intermediates; keep final output as: <prefix>.qf.fq.gz
#
# REPRODUCIBILITY / HPC DESIGN
# -----------------------------------------------------------------------------
# The script writes a self-contained bash script:
#    <prefix>.qualcheck.sh
#
# By default it ONLY generates the bash script. Use --run to execute immediately
# or submit the bash script to your scheduler (SLURM, etc.).
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#  - NanoStat     (from NanoStat / NanoPack)
#  - chopper
#  - porechop_abi
#  - rename.sh    (BBMap)
#  - gzip / gunzip (usually available)
#
# Optional directory overrides:
#  --nanostat_dir DIR    Directory containing NanoStat
#  --chopper_dir DIR     Directory containing chopper
#  --porechop_dir DIR    Directory containing porechop_abi
#  --bbmap_dir DIR       Directory containing rename.sh
#
###############################################################################

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use Cwd qw(getcwd);
use POSIX qw(strftime);

my $USAGE = <<"USAGE";

Usage:
  $0 -f reads.fastq.gz -p PREFIX -o OUTPUT_DIR [options]

Required:
  -f, --file FILE         Nanopore reads (fastq or fastq.gz)
  -p, --prefix PREFIX     Prefix for output files
  -o, --output DIR        Output directory (will be created)

Options:
  --ncpu INT              Threads (default: 40)
  --memory INT            Memory in GB for BBMap tools (default: 100; used as -Xmx100g)

  --min_len INT           Minimum read length (default: 500)
  --min_qual INT          Minimum read Q-score for first chopper pass (default: 10)

  --nanostat_dir DIR      Directory containing NanoStat (optional; else PATH)
  --chopper_dir DIR       Directory containing chopper (optional; else PATH)
  --porechop_dir DIR      Directory containing porechop_abi (optional; else PATH)
  --bbmap_dir DIR         Directory containing rename.sh (optional; else PATH)

  --run                   Execute generated workflow immediately
  --force                 Allow existing output folder (may overwrite files)
  -h, --help              Show help

Example:
  $0 -f raw.fq.gz -p sampleA -o sampleA.nanopore.qc --ncpu 64 --run

USAGE

# ----------------------------
# Args + defaults
# ----------------------------
my ($file,$prefix,$output,$help);
my ($NCPU,$MEMORY);
my ($min_len,$min_qual);
my ($nanostat_dir,$chopper_dir,$porechop_dir,$bbmap_dir);
my ($run,$force);

$NCPU     = 40;
$MEMORY   = 100;
$min_len  = 500;
$min_qual = 10;

GetOptions(
  "f|file=s"        => \$file,
  "p|prefix=s"      => \$prefix,
  "o|output=s"      => \$output,
  "ncpu=i"          => \$NCPU,
  "memory=i"        => \$MEMORY,

  "min_len=i"       => \$min_len,
  "min_qual=i"      => \$min_qual,

  "nanostat_dir=s"  => \$nanostat_dir,
  "chopper_dir=s"   => \$chopper_dir,
  "porechop_dir=s"  => \$porechop_dir,
  "bbmap_dir=s"     => \$bbmap_dir,

  "run!"            => \$run,
  "force!"          => \$force,
  "h|help"          => \$help,
) or die $USAGE;

die $USAGE if $help;
die "ERROR: Missing -f/--file\n$USAGE"   unless $file;
die "ERROR: Missing -p/--prefix\n$USAGE" unless $prefix;
die "ERROR: Missing -o/--output\n$USAGE" unless $output;

# ----------------------------
# Helpers
# ----------------------------
sub resolve_exec {
  my ($dir, $exe, $display_name) = @_;
  $display_name ||= $exe;

  if (defined $dir && $dir ne "") {
    $dir =~ s{/$}{};
    my $path = "$dir/$exe";
    die "ERROR: Cannot find executable for $display_name at: $path\n" unless -x $path;
    return $path;
  }

  my $which = `command -v $exe 2>/dev/null`;
  chomp $which;
  die "ERROR: '$exe' ($display_name) not found in PATH. Provide *_dir or load module/conda env.\n"
    unless $which;

  return $exe;
}

sub shell_quote {
  my ($s) = @_;
  $s =~ s/'/'"'"'/g;
  return "'$s'";
}

# ----------------------------
# Resolve tools
# ----------------------------
my $NANOSTAT = resolve_exec($nanostat_dir, "NanoStat",     "nanostat");
my $CHOPPER  = resolve_exec($chopper_dir,  "chopper",      "chopper");
my $PORECHOP = resolve_exec($porechop_dir, "porechop_abi", "porechop");
my $RENAME   = resolve_exec($bbmap_dir,    "rename.sh",    "bbmap");

# ----------------------------
# Output dir
# ----------------------------
if (-d $output && !$force) {
  die "ERROR: Output folder '$output' already exists. Use --force to proceed.\n";
}
make_path($output) unless -d $output;

# ----------------------------
# Derived file names
# ----------------------------
my $scriptfile = "$prefix.qualcheck.sh";

my $raw_report = "$output/$prefix.QualityReport.raw.txt";
my $qf_report  = "$output/$prefix.QualityReport.qf.txt";

my $step1_filt = "$output/$prefix.Filtered.fq.gz";
my $step2_trim = "$output/$prefix.Filtered.AdapterTrimmed.fq.gz";
my $step3_filt = "$output/$prefix.Filtered.AdapterTrimmed.Filtered.fq.gz";
my $step4_ren  = "$output/$prefix.Filtered.AdapterTrimmed.Filtered.Renamed.fq.gz";
my $final_qf   = "$output/$prefix.qf.fq.gz";

my $tmpdir     = "$output/$prefix.tmp";

my $BBMEM = "-Xmx${MEMORY}g";

# ----------------------------
# Generate workflow script
# ----------------------------
open(my $OUT, ">", $scriptfile) or die "ERROR: Cannot write '$scriptfile': $!\n";

my $date = strftime("%Y-%m-%d %H:%M:%S", localtime);
my $cwd  = getcwd();

print $OUT <<"SCRIPT";
#!/usr/bin/env bash
set -euo pipefail

# Generated by: $0
# Generated on: $date
# Working directory: $cwd
# Input: $file
# Prefix: $prefix
# Output: $output

echo "---------------------- [1/6] NanoStat on raw reads"
${\shell_quote($NANOSTAT)} --fastq ${\shell_quote($file)} -t $NCPU > ${\shell_quote($raw_report)}
echo "---------------------- DONE NanoStat raw"

echo "---------------------- [2/6] Filter by quality+length (chopper)"
# Using gunzip -c to stream. Works for .gz and (usually) also for plain files if gunzip supports it,
# but safest is gz input as intended.
gunzip -c ${\shell_quote($file)} | \\
  ${\shell_quote($CHOPPER)} -q $min_qual -l $min_len -t $NCPU | \\
  gzip > ${\shell_quote($step1_filt)}
echo "---------------------- DONE chopper pass 1"

echo "---------------------- [3/6] Adapter trimming (porechop_abi)"
rm -rf ${\shell_quote($tmpdir)}
${\shell_quote($PORECHOP)} -abi \\
  -i ${\shell_quote($step1_filt)} \\
  -o ${\shell_quote($step2_trim)} \\
  -tmp ${\shell_quote($tmpdir)} \\
  -t $NCPU \\
  --format fastq.gz
echo "---------------------- DONE porechop_abi"

echo "---------------------- [4/6] Filter by length again (chopper)"
gunzip -c ${\shell_quote($step2_trim)} | \\
  ${\shell_quote($CHOPPER)} -l $min_len -t $NCPU | \\
  gzip > ${\shell_quote($step3_filt)}
echo "---------------------- DONE chopper pass 2"

echo "---------------------- [5/6] Rename reads (BBMap rename.sh)"
${\shell_quote($RENAME)} $BBMEM \\
  in=${\shell_quote($step3_filt)} \\
  out=${\shell_quote($step4_ren)} \\
  prefix=${\shell_quote($prefix)} \\
  interleaved=false pigz=t unpigz=t threads=$NCPU
echo "---------------------- DONE rename"

echo "---------------------- [6/6] NanoStat on processed reads"
${\shell_quote($NANOSTAT)} --fastq ${\shell_quote($step4_ren)} -t $NCPU > ${\shell_quote($qf_report)}
echo "---------------------- DONE NanoStat processed"

echo "---------------------- Cleanup + finalize"
rm -f ${\shell_quote($step1_filt)} ${\shell_quote($step2_trim)} ${\shell_quote($step3_filt)}
mv -f ${\shell_quote($step4_ren)} ${\shell_quote($final_qf)}
rm -rf ${\shell_quote($tmpdir)}
echo "Final output: $final_qf"
echo "Reports: $raw_report ; $qf_report"
echo "---------------------- DONE"
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Generated workflow script: $scriptfile\n";
print "Output folder: $output\n";
print "Final output file will be: $final_qf\n";

# ----------------------------
# Optional execution
# ----------------------------
if ($run) {
  print "Executing workflow (bash $scriptfile)...\n";
  system("bash", $scriptfile) == 0 or die "ERROR: Workflow execution failed.\n";
} else {
  print "Not executed (default). Run with: bash $scriptfile\n";
  print "Or re-run this generator with --run\n";
}

exit 0;

__END__
