#!/usr/bin/env perl

###############################################################################
# Nanopore-short-read-mapping.pl
#
# Author: Cyrill Hofer
# Created: 26.04.2024
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates and optionally executes a standardized workflow for
# mapping short-read (Illumina) sequencing data against a corrected long-read
# (Nanopore) assembly.
#
# The workflow performs:
#
#   1) Mapping of Illumina reads to the assembly using BBMap (bbwrap.sh)
#   2) Conversion of SAM → sorted BAM using samtools
#   3) Indexing of the sorted BAM file
#   4) Generation of mapping statistics (samtools flagstat)
#   5) Cleanup of intermediate files
#
# The script does NOT directly execute mapping commands itself. Instead,
# it generates a self-contained bash workflow script:
#
#     <prefix>.short-read-mapping.sh
#
# This ensures:
#   - reproducibility
#   - transparency of executed commands
#   - easier debugging
#   - compatibility with HPC job schedulers
#
# The generated bash script can be:
#   - executed automatically using --run
#   - submitted to a cluster scheduler (SLURM, etc.)
#   - executed manually
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH or provided via arguments):
#
#   - bbwrap.sh  (BBMap suite)
#   - samtools
#
# Default behavior:
#   Tools are resolved from the system PATH.
#
# Optional:
#   --bbmap_dir <dir>     Directory containing bbwrap.sh
#   --samtools_dir <dir>  Directory containing samtools
#
# This makes the script portable across:
#   - Conda environments
#   - HPC module systems
#   - Local installations
#
# DESIGN PRINCIPLES
# -----------------------------------------------------------------------------
# - No hard-coded personal paths
# - Explicit dependency checks
# - Safe quoting of file paths
# - HPC-ready (thread control via --ncpu)
# - Reproducible bash workflow generation
#
###############################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use Cwd qw(getcwd);
use POSIX qw(strftime);

############################
# Usage
############################

my $USAGE = <<"USAGE";

Usage:
  $0 -f1 assembly.fasta -f2 reads.fastq.gz -p PREFIX -o OUTPUT_DIR [options]

Required:
  -f1              Corrected long-read assembly (FASTA)
  -f2              Quality-filtered Illumina reads (FASTQ.GZ)
  -p               Prefix for output files
  -o               Output directory (will be created)

Optional:
  --ncpu INT       Number of threads (default: 40)
  --memory INT     Memory in GB for BBMap (default: 100)
  --bbmap_dir DIR  Directory containing bbwrap.sh
  --samtools_dir DIR Directory containing samtools
  --run            Execute generated workflow immediately
  -h, --help       Show this help

Example:
  $0 -f1 assembly.fasta -f2 reads.fq.gz -p sample1 -o sample1_out --ncpu 32 --run

USAGE

############################
# Parse arguments
############################

my ($f1,$f2,$prefix,$output,$help,$NCPU,$MEMORY,$bbmap_dir,$samtools_dir,$run);

$NCPU   = 40;
$MEMORY = 100;

GetOptions(
    "f1=s"           => \$f1,
    "f2=s"           => \$f2,
    "p=s"            => \$prefix,
    "o=s"            => \$output,
    "ncpu=i"         => \$NCPU,
    "memory=i"       => \$MEMORY,
    "bbmap_dir=s"    => \$bbmap_dir,
    "samtools_dir=s" => \$samtools_dir,
    "run!"           => \$run,
    "h|help"         => \$help,
);

die $USAGE if $help;
die "ERROR: Missing -f1\n$USAGE" unless $f1;
die "ERROR: Missing -f2\n$USAGE" unless $f2;
die "ERROR: Missing -p\n$USAGE" unless $prefix;
die "ERROR: Missing -o\n$USAGE" unless $output;

############################
# Dependency resolution
############################

sub resolve_exec {
    my ($dir, $exe) = @_;

    if (defined $dir && $dir ne "") {
        $dir =~ s{/$}{};
        my $path = "$dir/$exe";
        die "ERROR: Cannot find executable: $path\n" unless -x $path;
        return $path;
    }

    my $which = `command -v $exe 2>/dev/null`;
    chomp $which;

    die "ERROR: '$exe' not found in PATH.\n" unless $which;
    return $exe;
}

my $BBWRAP   = resolve_exec($bbmap_dir,   "bbwrap.sh");
my $SAMTOOLS = resolve_exec($samtools_dir,"samtools");

############################
# Prepare output directory
############################

die "ERROR: Output directory '$output' already exists.\n" if -d $output;
make_path($output);

############################
# Generate workflow script
############################

my $scriptfile = "$prefix.short-read-mapping.sh";
open(my $OUT, ">", $scriptfile) or die "Cannot write $scriptfile\n";

my $date = strftime("%Y-%m-%d %H:%M:%S", localtime);
my $cwd  = getcwd();

print $OUT <<"SCRIPT";
#!/usr/bin/env bash
set -euo pipefail

# Generated on: $date
# Working directory: $cwd

echo "Starting BBMap mapping..."
"$BBWRAP" ref="$f1" in="$f2" out="$output/$prefix.sam" \\
  threads=$NCPU -Xmx${MEMORY}G \\
  kfilter=31 subfilter=15 maxindel=80 overwrite=true \\
  covstats="$output/$prefix.CovStats.txt"

echo "Sorting BAM..."
"$SAMTOOLS" view -bShu "$output/$prefix.sam" | \\
"$SAMTOOLS" sort -@ $NCPU -o "$output/$prefix.sorted.bam" -

echo "Indexing BAM..."
"$SAMTOOLS" index -@ $NCPU "$output/$prefix.sorted.bam"

echo "Generating mapping statistics..."
"$SAMTOOLS" flagstat -@ $NCPU "$output/$prefix.sorted.bam" > "$output/$prefix.MappedReads.txt"

echo "Cleaning intermediate files..."
rm -f "$output/$prefix.sam"

echo "Workflow completed successfully."
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Workflow script generated: $scriptfile\n";

############################
# Optional execution
############################

if ($run) {
    print "Executing workflow...\n";
    system("bash $scriptfile") == 0
        or die "ERROR: Workflow execution failed.\n";
}

print "Done.\n";
