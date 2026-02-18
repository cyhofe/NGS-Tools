#!/usr/bin/env perl

###############################################################################
# Nanopore-racon-polishing.pl
#
# Author: Cyrill Hofer
# Created: 29.11.2023
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a standardized long-read-only
# polishing workflow for a Nanopore assembly using Racon.
#
# WORKFLOW STEPS
# -----------------------------------------------------------------------------
#  1) Map quality-filtered Nanopore reads to the assembly using minimap2
#     (SAM output with -a for Racon input)
#  2) Polish the assembly using racon based on read-to-contig alignments
#  3) Cleanup intermediate alignment file
#
# REPRODUCIBILITY / HPC DESIGN
# -----------------------------------------------------------------------------
# The script writes a self-contained bash script:
#   <prefix>.racon-polishing.sh
#
# By default it ONLY generates the bash script. Use --run to execute immediately,
# or submit the generated bash script to your HPC scheduler.
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#  - minimap2
#  - racon
#
# Optional directory overrides:
#  --minimap2_dir DIR   Directory containing minimap2
#  --racon_dir DIR      Directory containing racon
#
# NOTES
# -----------------------------------------------------------------------------
# - This runs a single Racon round. For multiple rounds, run iteratively or
#   extend the script (can be added if you want).
# - Mapping preset is map-ont (Nanopore reads).
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
  $0 -f1 assembly.fasta -f2 reads.fq.gz -p PREFIX -o OUTPUT_DIR [options]

Required:
  -f1 FILE              Flye (or other) assembly (FASTA)
  -f2 FILE              Quality-filtered Nanopore reads (FASTQ/FASTQ.GZ)
  -p PREFIX             Prefix for output files
  -o OUTPUT_DIR          Output directory (will be created)

Options:
  --ncpu INT             Threads (default: 40)

  --minimap2_dir DIR     Directory containing minimap2 (optional; else PATH)
  --racon_dir DIR        Directory containing racon (optional; else PATH)

  --racon_window INT     Racon -w window length (default: 500)
  --run                  Execute generated workflow immediately
  --force                Allow existing output folder (may overwrite files)
  -h, --help             Show help

Example:
  $0 -f1 assembly.fasta -f2 reads.fq.gz -p sampleA -o sampleA.racon --run

USAGE

# ----------------------------
# Args + defaults
# ----------------------------
my ($f1,$f2,$prefix,$output,$help);
my ($NCPU);
my ($minimap2_dir,$racon_dir);
my ($racon_window);
my ($run,$force);

$NCPU = 40;
$racon_window = 500;

GetOptions(
  "f1=s"             => \$f1,
  "f2=s"             => \$f2,
  "p=s"              => \$prefix,
  "o=s"              => \$output,
  "ncpu=i"           => \$NCPU,

  "minimap2_dir=s"   => \$minimap2_dir,
  "racon_dir=s"      => \$racon_dir,

  "racon_window=i"   => \$racon_window,

  "run!"             => \$run,
  "force!"           => \$force,
  "h|help"           => \$help,
) or die $USAGE;

die $USAGE if $help;

die "ERROR: Missing -f1\n$USAGE" unless $f1;
die "ERROR: Missing -f2\n$USAGE" unless $f2;
die "ERROR: Missing -p\n$USAGE"  unless $prefix;
die "ERROR: Missing -o\n$USAGE"  unless $output;

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
my $MINIMAP2 = resolve_exec($minimap2_dir, "minimap2", "minimap2");
my $RACON    = resolve_exec($racon_dir,    "racon",    "racon");

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
my $scriptfile = "$prefix.racon-polishing.sh";

my $sam      = "$output/$prefix.minimap2.sam";
my $racon_fa = "$output/$prefix.racon.fasta";

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
# Assembly: $f1
# Reads: $f2
# Prefix: $prefix
# Output: $output

echo "---------------------- [1/2] Map reads to assembly (minimap2)"
${\shell_quote($MINIMAP2)} -x map-ont -a -t $NCPU ${\shell_quote($f1)} ${\shell_quote($f2)} > ${\shell_quote($sam)}
echo "---------------------- DONE mapping"

echo "---------------------- [2/2] Polish with Racon"
${\shell_quote($RACON)} -m 8 -x -6 -g -8 -w $racon_window -t $NCPU \\
  ${\shell_quote($f2)} ${\shell_quote($sam)} ${\shell_quote($f1)} > ${\shell_quote($racon_fa)}
echo "---------------------- DONE polishing"

echo "---------------------- Cleanup"
rm -f ${\shell_quote($sam)}
echo "Final output: $racon_fa"
echo "---------------------- DONE"
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Generated workflow script: $scriptfile\n";
print "Output folder: $output\n";
print "Final output file will be: $racon_fa\n";

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
