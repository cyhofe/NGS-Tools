#!/usr/bin/env perl

###############################################################################
# Hybrid-pilon-polishing.pl
#
# Author: Cyrill Hofer
# Created: 29.11.2023
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a standardized hybrid polishing
# workflow for correcting a long-read (Nanopore) assembly using short reads
# (Illumina) and Pilon.
#
# INPUTS
# -----------------------------------------------------------------------------
#  - A long-read assembly that has already been corrected (e.g., Flye + Medaka)
#  - Quality-controlled Illumina reads (paired-end) in FASTQ(.gz)
#
# WORKFLOW STEPS
# -----------------------------------------------------------------------------
#  1) Index the assembly with BWA
#  2) Map Illumina reads to the assembly (bwa mem)
#  3) Convert SAM -> BAM, sort BAM, and index BAM (samtools)
#  4) Produce mapping statistics (samtools flagstat + depth)
#  5) Run Pilon to correct bases supported by Illumina read evidence
#
# REPRODUCIBILITY / HPC DESIGN
# -----------------------------------------------------------------------------
# This script writes a self-contained bash script:
#    <prefix>.pilon-polishing.sh
#
# Benefits:
#   - reproducible and transparent commands
#   - easy debugging
#   - can be submitted to SLURM or run interactively
#
# By default this Perl script ONLY generates the workflow bash script.
# Use --run to execute immediately.
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#   - bwa
#   - samtools
#   - java
#   - pilon (either pilon.jar path or Pilon installed in PATH)
#
# This generator supports two ways to specify Pilon:
#   (A) --pilon_jar /path/to/pilon.jar  (recommended, most explicit)
#   (B) --pilon_dir DIR where pilon-*.jar is located (will use latest match)
#
# Optional tool directory overrides:
#   --bwa_dir DIR       (directory containing bwa)
#   --samtools_dir DIR  (directory containing samtools)
#   --java_dir DIR      (directory containing java)
#
###############################################################################

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use Cwd qw(getcwd abs_path);
use POSIX qw(strftime);

my $USAGE = <<"USAGE";

Usage:
  $0 -f1 assembly.fasta -f2 illumina.fastq.gz -p PREFIX -o OUTPUT_DIR [options]

Required:
  -f1 FILE              Medaka-corrected (or otherwise corrected) assembly (FASTA)
  -f2 FILE              Illumina reads (FASTQ/FASTQ.GZ). Use interleaved reads if using -p with bwa mem.
  -p PREFIX             Prefix for output files
  -o OUTPUT_DIR          Output directory (will be created)

Options:
  --ncpu INT             Threads (default: 40)
  --memory INT           Java memory in GB for Pilon (default: 100)

  --bwa_dir DIR          Directory containing bwa (optional; else PATH)
  --samtools_dir DIR     Directory containing samtools (optional; else PATH)
  --java_dir DIR         Directory containing java (optional; else PATH)

  --pilon_jar FILE       Path to pilon.jar (recommended)
  --pilon_dir DIR        Directory containing pilon*.jar (uses latest match)

  --run                  Execute generated workflow immediately
  --force                Allow existing output folder (may overwrite files)
  -h, --help             Show help

Example:
  $0 -f1 assembly.fasta -f2 reads.interleaved.fq.gz -p sampleA -o sampleA_pilon \\
     --ncpu 64 --memory 120 --pilon_jar /path/to/pilon.jar --run

USAGE

# ----------------------------
# Args
# ----------------------------
my ($f1,$f2,$prefix,$output,$help);
my ($NCPU,$MEMORY);
my ($bwa_dir,$samtools_dir,$java_dir);
my ($pilon_jar,$pilon_dir);
my ($run,$force);

$NCPU   = 40;
$MEMORY = 100;

GetOptions(
  "f1=s"          => \$f1,
  "f2=s"          => \$f2,
  "p=s"           => \$prefix,
  "o=s"           => \$output,
  "ncpu=i"        => \$NCPU,
  "memory=i"      => \$MEMORY,

  "bwa_dir=s"     => \$bwa_dir,
  "samtools_dir=s"=> \$samtools_dir,
  "java_dir=s"    => \$java_dir,

  "pilon_jar=s"   => \$pilon_jar,
  "pilon_dir=s"   => \$pilon_dir,

  "run!"          => \$run,
  "force!"        => \$force,
  "h|help"        => \$help,
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
  die "ERROR: '$exe' ($display_name) not found in PATH. Provide --${display_name}_dir or load a module/conda env.\n"
    unless $which;

  return $exe;
}

sub shell_quote {
  my ($s) = @_;
  $s =~ s/'/'"'"'/g;
  return "'$s'";
}

sub resolve_pilon_jar {
  my ($jar, $dir) = @_;

  if (defined $jar && $jar ne "") {
    die "ERROR: Pilon jar not found: $jar\n" unless -f $jar;
    return abs_path($jar);
  }

  if (defined $dir && $dir ne "") {
    $dir =~ s{/$}{};
    die "ERROR: Pilon dir not found: $dir\n" unless -d $dir;

    # Use newest pilon*.jar in that directory (lexicographic sort usually OK for versions;
    # to be safe, we sort by mtime here).
    opendir(my $DH, $dir) or die "ERROR: Cannot open pilon_dir $dir: $!\n";
    my @jars = grep { /^pilon.*\.jar$/ && -f "$dir/$_" } readdir($DH);
    closedir($DH);

    die "ERROR: No pilon*.jar found in --pilon_dir $dir\n" unless @jars;

    @jars = sort { (stat("$dir/$b"))[9] <=> (stat("$dir/$a"))[9] } @jars; # newest first
    return abs_path("$dir/$jars[0]");
  }

  die "ERROR: You must provide --pilon_jar /path/to/pilon.jar (recommended) OR --pilon_dir DIR\n";
}

# ----------------------------
# Resolve tools
# ----------------------------
my $BWA      = resolve_exec($bwa_dir,     "bwa",     "bwa");
my $SAMTOOLS = resolve_exec($samtools_dir,"samtools","samtools");
my $JAVA     = resolve_exec($java_dir,    "java",    "java");
my $PILON_JAR = resolve_pilon_jar($pilon_jar, $pilon_dir);

# ----------------------------
# Output dir
# ----------------------------
if (-d $output && !$force) {
  die "ERROR: Output directory '$output' already exists. Use --force to proceed.\n";
}
make_path($output) unless -d $output;

# ----------------------------
# Generate workflow script
# ----------------------------
my $scriptfile = "$prefix.pilon-polishing.sh";
open(my $OUT, ">", $scriptfile) or die "ERROR: Cannot write '$scriptfile': $!\n";

my $date = strftime("%Y-%m-%d %H:%M:%S", localtime);
my $cwd  = getcwd();

my $q_f1      = shell_quote($f1);
my $q_f2      = shell_quote($f2);
my $q_prefix  = shell_quote($prefix);
my $q_output  = shell_quote($output);

my $sam        = "$output/$prefix.sam";
my $bam        = "$output/$prefix.bam";
my $sorted_bam = "$output/$prefix.sorted.bam";
my $flagstat   = "$output/$prefix.MappedReads.txt";
my $depth      = "$output/$prefix.depth.txt";

my $q_sam        = shell_quote($sam);
my $q_bam        = shell_quote($bam);
my $q_sorted_bam = shell_quote($sorted_bam);
my $q_flagstat   = shell_quote($flagstat);
my $q_depth      = shell_quote($depth);

# Pilon output naming: --output writes files as <outprefix>.fasta etc in the working dir,
# so we set --outdir to keep everything inside $output.
my $pilon_out_prefix = "$prefix.pilon";
my $q_pilon_out_prefix = shell_quote($pilon_out_prefix);
my $q_pilon_outdir = shell_quote($output);

print $OUT <<"SCRIPT";
#!/usr/bin/env bash
set -euo pipefail

# Generated by: $0
# Generated on: $date
# Working directory: $cwd
# Assembly: $f1
# Reads: $f2
# Prefix: $prefix
# Output dir: $output

echo "---------------------- [1/6] Index assembly with bwa"
${\shell_quote($BWA)} index $q_f1
echo "---------------------- DONE indexing"

echo "---------------------- [2/6] Map reads with bwa mem"
# Note: This assumes reads are interleaved (bwa mem -p).
# If you have R1/R2 separately, modify this section accordingly.
${\shell_quote($BWA)} mem -t ${NCPU} $q_f1 -p $q_f2 > $q_sam
echo "---------------------- DONE mapping"

echo "---------------------- [3/6] Convert SAM -> BAM"
${\shell_quote($SAMTOOLS)} view -@ ${NCPU} -b $q_sam -o $q_bam
echo "---------------------- DONE converting"

echo "---------------------- [4/6] Sort + index BAM"
${\shell_quote($SAMTOOLS)} sort -@ ${NCPU} -o $q_sorted_bam $q_bam
${\shell_quote($SAMTOOLS)} index -@ ${NCPU} $q_sorted_bam
echo "---------------------- DONE sort/index"

echo "---------------------- [5/6] Mapping statistics"
${\shell_quote($SAMTOOLS)} flagstat -@ ${NCPU} $q_sorted_bam > $q_flagstat
${\shell_quote($SAMTOOLS)} depth -a $q_sorted_bam > $q_depth
echo "---------------------- DONE stats"

echo "---------------------- [6/6] Pilon polishing"
${\shell_quote($JAVA)} -Xmx${MEMORY}G -jar ${\shell_quote($PILON_JAR)} \\
  --genome $q_f1 \\
  --frags $q_sorted_bam \\
  --changes \\
  --output $q_pilon_out_prefix \\
  --outdir $q_pilon_outdir
echo "---------------------- DONE pilon"

echo "---------------------- Cleanup"
rm -f $q_sam $q_bam
# Keep sorted BAM + index + stats by default.
echo "---------------------- DONE"

echo "All done."
echo "Pilon output (FASTA) should be: $output/${pilon_out_prefix}.fasta"
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Generated workflow script: $scriptfile\n";
print "Output folder: $output\n";
print "Pilon jar: $PILON_JAR\n";

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
