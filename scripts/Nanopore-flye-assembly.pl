#!/usr/bin/env perl

###############################################################################
# Nanopore-metaflye-assembly.pl
#
# Author: Cyrill Hofer
# Created: 29.11.2023
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a standardized workflow for
# assembling nanopore metagenomic reads with metaFlye, filtering short contigs,
# renaming contigs to a stable prefix-based naming scheme, and generating basic
# assembly statistics.
#
# Workflow steps:
#   1) Assembly using Flye in metagenome mode:
#        flye --meta --nano-hq reads.fq.gz -o <prefix>.flye -t <NCPU>
#   2) Filter contigs shorter than a threshold (default: 2500 bp) using lenfilter
#   3) Rename contigs using BBMap's rename.sh for standardized contig IDs
#   4) Generate sequence statistics using seqstat
#   5) Cleanup intermediate files, keep raw contigs and final contigs
#
# DESIGN / REPRODUCIBILITY
# -----------------------------------------------------------------------------
# The script does not directly run the assembly commands. Instead, it writes a
# self-contained bash script:
#    <prefix>.flye.sh
#
# This improves reproducibility and makes the executed commands transparent.
# Use --run to execute immediately, or submit the generated script to an HPC
# scheduler (SLURM, etc.).
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#   - flye
#   - lenfilter
#   - rename.sh  (BBMap)
#   - seqstat
#
# Optional tool directory overrides:
#   --flye_dir <dir>       (directory containing flye)
#   --lenfilter_dir <dir>  (directory containing lenfilter)
#   --bbmap_dir <dir>      (directory containing rename.sh)
#   --seqstat_dir <dir>    (directory containing seqstat)
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
  $0 -f reads.fastq.gz -p PREFIX [options]

Required:
  -f, --file FILE         Nanopore reads (fastq or fastq.gz)
  -p, --prefix PREFIX     Prefix for output files/folder

Options:
  --ncpu INT              Threads (default: 40; Flye may cap internally)
  --min_len INT           Minimum contig length to keep (default: 2500)

  --flye_dir DIR          Directory containing flye (optional; else PATH)
  --lenfilter_dir DIR     Directory containing lenfilter (optional; else PATH)
  --bbmap_dir DIR         Directory containing rename.sh (optional; else PATH)
  --seqstat_dir DIR       Directory containing seqstat (optional; else PATH)

  --run                   Execute generated workflow immediately
  --force                 Allow existing output folder (will overwrite some files)
  -h, --help              Show help

Example:
  $0 -f reads.fq.gz -p sampleA --ncpu 64 --min_len 2500 --run

USAGE

# ----------------------------
# Args
# ----------------------------
my ($file,$prefix,$help,$NCPU,$min_len);
my ($flye_dir,$lenfilter_dir,$bbmap_dir,$seqstat_dir);
my ($run,$force);

$NCPU    = 40;
$min_len = 2500;

GetOptions(
  "f|file=s"         => \$file,
  "p|prefix=s"       => \$prefix,
  "ncpu=i"           => \$NCPU,
  "min_len=i"        => \$min_len,

  "flye_dir=s"       => \$flye_dir,
  "lenfilter_dir=s"  => \$lenfilter_dir,
  "bbmap_dir=s"      => \$bbmap_dir,
  "seqstat_dir=s"    => \$seqstat_dir,

  "run!"             => \$run,
  "force!"           => \$force,
  "h|help"           => \$help,
) or die $USAGE;

die $USAGE if $help;
die "ERROR: Missing -f/--file\n$USAGE"   unless $file;
die "ERROR: Missing -p/--prefix\n$USAGE" unless $prefix;

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
  die "ERROR: '$exe' ($display_name) not found in PATH. Provide --${display_name}_dir or load a module/conda env.\n"
    unless $which;

  return $exe; # call via PATH
}

sub shell_quote {
  my ($s) = @_;
  $s =~ s/'/'"'"'/g;
  return "'$s'";
}

# ----------------------------
# Resolve tools
# ----------------------------
my $FLYE      = resolve_exec($flye_dir,      "flye",      "flye");
my $LENFILTER = resolve_exec($lenfilter_dir, "lenfilter", "lenfilter");
my $RENAME    = resolve_exec($bbmap_dir,     "rename.sh", "bbmap");     # display_name used only for message
my $SEQSTAT   = resolve_exec($seqstat_dir,   "seqstat",   "seqstat");

# ----------------------------
# Output paths
# ----------------------------
my $outdir     = "$prefix.flye";
my $scriptfile = "$prefix.flye.sh";

if (-d $outdir && !$force) {
  die "ERROR: Output folder '$outdir' already exists. Use --force if you really want to proceed.\n";
}
make_path($outdir) unless -d $outdir;

# ----------------------------
# Generate workflow script
# ----------------------------
open(my $OUT, ">", $scriptfile) or die "ERROR: Cannot write output script '$scriptfile': $!\n";

my $date = strftime("%Y-%m-%d %H:%M:%S", localtime);
my $cwd  = getcwd();

my $q_file   = shell_quote($file);
my $q_prefix = shell_quote($prefix);
my $q_outdir = shell_quote($outdir);

my $raw_asm          = "$outdir/assembly.fasta";
my $filtered_fasta   = "$outdir/$prefix.flye.filtered.fasta";
my $final_contigs    = "$outdir/$prefix.contigs.fasta";
my $seqstat_report   = "$outdir/$prefix.contigs.seqstat.txt";
my $raw_contigs_keep = "$outdir/$prefix.raw_contigs.fasta";

my $q_raw_asm        = shell_quote($raw_asm);
my $q_filtered_fasta = shell_quote($filtered_fasta);
my $q_final_contigs  = shell_quote($final_contigs);
my $q_seqstat_report = shell_quote($seqstat_report);
my $q_raw_keep       = shell_quote($raw_contigs_keep);

print $OUT <<"SCRIPT";
#!/usr/bin/env bash
set -euo pipefail

# Generated by: $0
# Generated on: $date
# Working directory: $cwd
# Input reads: $file
# Prefix: $prefix
# Output dir: $outdir

echo "---------------------- [1/4] Running metaFlye assembly"
${\shell_quote($FLYE)} --meta --nano-hq $q_file -o $q_outdir -t ${NCPU}
echo "---------------------- DONE assembling"

if [[ ! -s $q_raw_asm ]]; then
  echo "ERROR: Expected assembly file not found or empty: $raw_asm" >&2
  exit 1
fi

echo "---------------------- [2/4] Filtering contigs < ${min_len} bp (lenfilter)"
cat $q_raw_asm | ${\shell_quote($LENFILTER)} -f fasta -me -len ${min_len} > $q_filtered_fasta
echo "---------------------- DONE length filtering"

echo "---------------------- [3/4] Renaming contigs (BBMap rename.sh)"
${\shell_quote($RENAME)} in=$q_filtered_fasta out=$q_final_contigs prefix=${\shell_quote("$prefix.contig")} interleaved=false threads=${NCPU}
echo "---------------------- DONE renaming"

echo "---------------------- [4/4] Generating contig statistics (seqstat)"
${\shell_quote($SEQSTAT)} $q_final_contigs > $q_seqstat_report
echo "---------------------- DONE seqstat"

echo "---------------------- Cleanup"
rm -f $q_filtered_fasta
mv -f $q_raw_asm $q_raw_keep
echo "---------------------- DONE"

echo "All done. Final contigs: $final_contigs"
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Generated workflow script: $scriptfile\n";
print "Output folder: $outdir\n";

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
