#!/usr/bin/env perl

###############################################################################
# Nanopore-16S-amplicon-quality-filter.pl
#
# Author: Cyrill Hofer
# Created: 29.11.2023
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a preprocessing workflow for
# raw Nanopore 16S amplicon reads. It performs length/quality filtering,
# adapter trimming, and a second strict filtering step tailored to the expected
# amplicon size range.
#
# WORKFLOW STEPS
# -----------------------------------------------------------------------------
#  1) NanoStat report on raw reads
#  2) Initial filter (broad): chopper quality + length range
#  3) Adapter trimming: porechop_abi
#  4) Final filter (strict): chopper quality + tighter length range
#  5) NanoStat report on final reads
#  6) Cleanup intermediates; keep final output as: <prefix>.qf.fq.gz
#
# REPRODUCIBILITY / HPC DESIGN
# -----------------------------------------------------------------------------
# This Perl script writes a self-contained bash script:
#    <prefix>.qualcheck.sh
#
# By default it ONLY generates the bash script. Use --run to execute immediately
# or submit the generated script to your HPC scheduler.
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#  - NanoStat (NanoPack)
#  - chopper
#  - porechop_abi
#  - gzip / gunzip
#
# Optional directory overrides:
#  --nanostat_dir DIR    Directory containing NanoStat
#  --chopper_dir DIR     Directory containing chopper
#  --porechop_dir DIR    Directory containing porechop_abi
#
# NOTES / FIXES vs ORIGINAL
# -----------------------------------------------------------------------------
# - Fixed a bug in the original script: options were incorrectly reused as "b|..."
#   and path assignment used "$$var" (double dereference). This version uses
#   correct option names and assignments.
# - Removed unused BBMap dependency (it was never used in this script).
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
  $0 -f reads.fq.gz -p PREFIX -o OUTPUT_DIR [options]

Required:
  -f, --file FILE         Nanopore 16S amplicon reads (FASTQ/FASTQ.GZ)
  -p, --prefix PREFIX     Prefix for output files
  -o, --output DIR        Output directory (will be created)

Options:
  --ncpu INT              Threads (default: 40)

  # First (broad) filter pass:
  --q1 INT                Min Q-score (default: 10)
  --minlen1 INT           Min length (default: 1000)
  --maxlen1 INT           Max length (default: 2000)

  # Second (strict) filter pass:
  --q2 INT                Min Q-score (default: 17)
  --minlen2 INT           Min length (default: 1200)
  --maxlen2 INT           Max length (default: 1600)

  --nanostat_dir DIR      Directory containing NanoStat (optional; else PATH)
  --chopper_dir DIR       Directory containing chopper (optional; else PATH)
  --porechop_dir DIR      Directory containing porechop_abi (optional; else PATH)

  --run                   Execute generated workflow immediately
  --force                 Allow existing output folder (may overwrite files)
  -h, --help              Show help

Example:
  $0 -f reads.fq.gz -p sample16S -o sample16S.qc --ncpu 48 --run

USAGE

# ----------------------------
# Args + defaults
# ----------------------------
my ($file,$prefix,$output,$help);
my ($NCPU);

my ($q1,$minlen1,$maxlen1);
my ($q2,$minlen2,$maxlen2);

my ($nanostat_dir,$chopper_dir,$porechop_dir);
my ($run,$force);

$NCPU = 40;

$q1 = 10;    $minlen1 = 1000; $maxlen1 = 2000;
$q2 = 17;    $minlen2 = 1200; $maxlen2 = 1600;

GetOptions(
  "f|file=s"        => \$file,
  "p|prefix=s"      => \$prefix,
  "o|output=s"      => \$output,
  "ncpu=i"          => \$NCPU,

  "q1=i"            => \$q1,
  "minlen1=i"       => \$minlen1,
  "maxlen1=i"       => \$maxlen1,

  "q2=i"            => \$q2,
  "minlen2=i"       => \$minlen2,
  "maxlen2=i"       => \$maxlen2,

  "nanostat_dir=s"  => \$nanostat_dir,
  "chopper_dir=s"   => \$chopper_dir,
  "porechop_dir=s"  => \$porechop_dir,

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

# ----------------------------
# Output dir
# ----------------------------
if (-d $output && !$force) {
  die "ERROR: Output folder '$output' already exists. Use --force to proceed.\n";
}
make_path($output) unless -d $output;

# ----------------------------
# Files
# ----------------------------
my $scriptfile   = "$prefix.qualcheck.sh";

my $raw_report   = "$output/$prefix.QualityReport.raw.txt";
my $qf_report    = "$output/$prefix.QualityReport.qf.txt";

my $step1_filt   = "$output/$prefix.Filtered.fq.gz";
my $step2_trim   = "$output/$prefix.Filtered.AdapterTrimmed.fq.gz";
my $step3_filt   = "$output/$prefix.Filtered.AdapterTrimmed.Filtered.fq.gz";

my $tmpdir       = "$output/$prefix.tmp";
my $final_qf     = "$output/$prefix.qf.fq.gz";

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

echo "---------------------- [1/5] NanoStat on raw reads"
${\shell_quote($NANOSTAT)} --fastq ${\shell_quote($file)} -t $NCPU > ${\shell_quote($raw_report)}
echo "---------------------- DONE NanoStat raw"

echo "---------------------- [2/5] Filter pass 1 (chopper): q>=$q1, len=$minlen1..$maxlen1"
gunzip -c ${\shell_quote($file)} | \\
  ${\shell_quote($CHOPPER)} -q $q1 --minlength $minlen1 --maxlength $maxlen1 -t $NCPU | \\
  gzip > ${\shell_quote($step1_filt)}
echo "---------------------- DONE chopper pass 1"

echo "---------------------- [3/5] Adapter trimming (porechop_abi)"
rm -rf ${\shell_quote($tmpdir)}
${\shell_quote($PORECHOP)} -abi \\
  -i ${\shell_quote($step1_filt)} \\
  -tmp ${\shell_quote($tmpdir)} \\
  -o ${\shell_quote($step2_trim)} \\
  -t $NCPU \\
  --format fastq.gz
echo "---------------------- DONE porechop_abi"

echo "---------------------- [4/5] Filter pass 2 (chopper): q>=$q2, len=$minlen2..$maxlen2"
gunzip -c ${\shell_quote($step2_trim)} | \\
  ${\shell_quote($CHOPPER)} -q $q2 --minlength $minlen2 --maxlength $maxlen2 -t $NCPU | \\
  gzip > ${\shell_quote($step3_filt)}
echo "---------------------- DONE chopper pass 2"

echo "---------------------- [5/5] NanoStat on final reads"
${\shell_quote($NANOSTAT)} --fastq ${\shell_quote($step3_filt)} -t $NCPU > ${\shell_quote($qf_report)}
echo "---------------------- DONE NanoStat final"

echo "---------------------- Cleanup + finalize"
rm -f ${\shell_quote($step1_filt)} ${\shell_quote($step2_trim)}
rm -rf ${\shell_quote($tmpdir)}
mv -f ${\shell_quote($step3_filt)} ${\shell_quote($final_qf)}
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
