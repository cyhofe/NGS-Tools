#!/usr/bin/env perl

###############################################################################
# Illumina-quality-filter.pl
#
# Author: Cyrill Hofer (modified from stamp-quality-check.pl)
# Created: 29.11.2023
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a standardized Illumina
# short-read preprocessing workflow (paired-end), producing a cleaned,
# interleaved FASTQ(.gz) suitable for downstream assembly / mapping workflows.
#
# WORKFLOW STEPS
# -----------------------------------------------------------------------------
#  1) Combine R1/R2 into a single interleaved FASTQ using BBMap reformat.sh
#  2) Run FastQC on the raw interleaved reads
#  3) Quality trimming using BBMap bbduk.sh (qtrim + trimq; includes tbo)
#  4) Rename reads using BBMap rename.sh to apply a consistent prefix
#  5) Remove vector contamination (e.g., PhiX, cloning vectors) using bbduk.sh
#  6) Trim Illumina adapters using bbduk.sh with a user-provided adapter FASTA
#  7) Post-trim de novo adapter detection using bbmerge.sh
#  8) Run FastQC on the final cleaned reads
#  9) Cleanup intermediates; keep final output as: <prefix>.qf.fq.gz
#
# REPRODUCIBILITY / HPC DESIGN
# -----------------------------------------------------------------------------
# This Perl script writes a self-contained bash script:
#   <prefix>.qualcheck.sh
#
# By default it ONLY generates the bash script. Use --run to execute immediately,
# or submit the generated bash script to your HPC scheduler.
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#  - reformat.sh   (BBMap)
#  - bbduk.sh      (BBMap)
#  - rename.sh     (BBMap)
#  - bbmerge.sh    (BBMap)
#  - fastqc
#
# Optional (only used if present AND --vecscreen is enabled):
#  - run_vecscreen (NCBI VecScreen wrapper; site-specific)
#  - lenseq, prettytable (site-specific helpers; optional)
#
# NOTE ON VEC-SCREEN / SITE-SPECIFIC HELPERS
# -----------------------------------------------------------------------------
# Your original workflow called "run_vecscreen", "lenseq", and "prettytable".
# Those are not standard tools on most systems. This publishable version:
#   - keeps the de novo adapter detection fasta
#   - optionally runs vecscreen if --vecscreen is given AND run_vecscreen exists
#   - does NOT require lenseq/prettytable; instead it prints the fasta headers
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
  $0 -f1 R1.fq.gz -f2 R2.fq.gz -p PREFIX -o OUTPUT_DIR [options]

Required:
  -f1, --file1 FILE        Forward reads (R1), fastq(.gz)
  -f2, --file2 FILE        Reverse reads (R2), fastq(.gz)
  -p,  --prefix PREFIX     Prefix for output files
  -o,  --output DIR        Output directory (will be created)

Required reference FASTAs (defaults can be overridden):
  -a, --adapter FILE       Adapter FASTA for trimming (REQUIRED)
  -v, --vectors FILE       Vector/contaminant FASTA (REQUIRED)

Options:
  --ncpu INT               Threads (default: 40)
  --memory INT             Java-like memory flag for BBTools (default: 100; used as -Xmx100g)

  --bbmap_dir DIR          Directory containing BBMap scripts (reformat.sh, bbduk.sh, rename.sh, bbmerge.sh)
                           (optional; otherwise resolved from PATH)
  --fastqc_dir DIR         Directory containing fastqc (optional; otherwise resolved from PATH)

  --trimq INT              Trimming quality threshold (default: 18)
  --qtrim STR              BBduk qtrim mode (default: rl)
  --k INT                  kmer length for bbduk ref trimming (default: 21)

  --vecscreen              If set, attempt to run run_vecscreen on adapters_after_trimming.fa (optional)
  --run                    Execute generated workflow immediately
  --force                  Allow existing output folder (may overwrite files)
  -h, --help               Show help

Defaults for reference FASTAs in your environment (can override):
  --adapter  /home/opt/progs/tools_cyrill/Dependencies/Illumina_SequencingAdapters.fasta
  --vectors  /home/opt/progs/tools_cyrill/Dependencies/Illumina_CloningVectors.fasta

Example:
  $0 -f1 sample_R1.fq.gz -f2 sample_R2.fq.gz -p sample -o sample.qc \\
     --ncpu 32 --run

USAGE

# ----------------------------
# Args + defaults
# ----------------------------
my ($file1,$file2,$prefix,$output,$help);
my ($NCPU,$MEMORY);
my ($adapterfile,$vectorfile);
my ($bbmap_dir,$fastqc_dir);
my ($trimq,$qtrim,$kmer);
my ($run,$force,$vecscreen);

$NCPU   = 40;
$MEMORY = 100;

# Keep your current defaults, but allow override
$adapterfile = "/home/opt/progs/tools_cyrill/Dependencies/Illumina_SequencingAdapters.fasta";
$vectorfile  = "/home/opt/progs/tools_cyrill/Dependencies/Illumina_CloningVectors.fasta";

$trimq = 18;
$qtrim = "rl";
$kmer  = 21;

GetOptions(
  "f1|file1=s"   => \$file1,
  "f2|file2=s"   => \$file2,
  "p|prefix=s"   => \$prefix,
  "o|output=s"   => \$output,

  "a|adapter=s"  => \$adapterfile,
  "v|vectors=s"  => \$vectorfile,

  "ncpu=i"       => \$NCPU,
  "memory=i"     => \$MEMORY,

  "bbmap_dir=s"  => \$bbmap_dir,
  "fastqc_dir=s" => \$fastqc_dir,

  "trimq=i"      => \$trimq,
  "qtrim=s"      => \$qtrim,
  "k=i"          => \$kmer,

  "vecscreen!"   => \$vecscreen,
  "run!"         => \$run,
  "force!"       => \$force,
  "h|help"       => \$help,
) or die $USAGE;

die $USAGE if $help;

die "ERROR: Missing -f1/--file1\n$USAGE" unless $file1;
die "ERROR: Missing -f2/--file2\n$USAGE" unless $file2;
die "ERROR: Missing -p/--prefix\n$USAGE" unless $prefix;
die "ERROR: Missing -o/--output\n$USAGE" unless $output;

die "ERROR: Adapter FASTA not found: $adapterfile\n" unless $adapterfile && -f $adapterfile;
die "ERROR: Vectors FASTA not found: $vectorfile\n" unless $vectorfile && -f $vectorfile;

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
  die "ERROR: '$exe' ($display_name) not found in PATH. Provide appropriate *_dir or load module/conda env.\n"
    unless $which;

  return $exe; # call via PATH
}

sub shell_quote {
  my ($s) = @_;
  $s =~ s/'/'"'"'/g;
  return "'$s'";
}

# Resolve BBTools scripts (they might be executable scripts, so -x should work when given a dir)
my $REFORMAT = resolve_exec($bbmap_dir, "reformat.sh", "bbmap");
my $BBDUK    = resolve_exec($bbmap_dir, "bbduk.sh",    "bbmap");
my $RENAME   = resolve_exec($bbmap_dir, "rename.sh",   "bbmap");
my $BBMERGE  = resolve_exec($bbmap_dir, "bbmerge.sh",  "bbmap");

my $FASTQC   = resolve_exec($fastqc_dir, "fastqc",     "fastqc");

# ----------------------------
# Output dir
# ----------------------------
if (-d $output && !$force) {
  die "ERROR: Output folder '$output' already exists. Use --force to proceed.\n";
}
make_path($output) unless -d $output;

# ----------------------------
# Paths
# ----------------------------
my $scriptfile = "$prefix.qualcheck.sh";

my $interleaved_raw  = "$output/$prefix.interleaved.fq.gz";
my $raw_fastqc_dir   = "$output/fastqc_raw";
my $trimmed          = "$output/$prefix.trimmed.fq.gz";
my $trimmed_renamed  = "$output/$prefix.trimmed.renamed.fq.gz";
my $no_vector        = "$output/$prefix.trimmed.renamed.novector.fq.gz";
my $no_adapt         = "$output/$prefix.trimmed.renamed.novector.noadapt.fq.gz";
my $final_qf         = "$output/$prefix.qf.fq.gz";

my $adapter_after_fa = "$output/adapters_after_trimming.fa";
my $processed_fastqc_dir = "$output/fastqc_qf";

my $qualtrim_stats   = "$output/qualtrim.stats.txt";
my $vector_stats     = "$output/vector.stats.txt";
my $adapter_stats    = "$output/adapter.stats.txt";

# BBTools memory flag
my $BBMEM = "-Xmx${MEMORY}g";

# ----------------------------
# Generate workflow bash script
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
# R1: $file1
# R2: $file2
# Prefix: $prefix
# Output: $output

mkdir -p ${\shell_quote($raw_fastqc_dir)} ${\shell_quote($processed_fastqc_dir)}

echo "---------------------- [1/8] Create interleaved reads (reformat.sh)"
${\shell_quote($REFORMAT)} $BBMEM in=${\shell_quote($file1)} in2=${\shell_quote($file2)} \\
  out=${\shell_quote($interleaved_raw)} interleaved=f pigz=t unpigz=t threads=$NCPU
echo "---------------------- DONE interleaving"

echo "---------------------- [2/8] FastQC on raw interleaved reads"
${\shell_quote($FASTQC)} -t $NCPU -o ${\shell_quote($raw_fastqc_dir)} ${\shell_quote($interleaved_raw)}
echo "---------------------- DONE FastQC raw"

echo "---------------------- [3/8] Quality trimming (bbduk.sh)"
${\shell_quote($BBDUK)} $BBMEM in=${\shell_quote($interleaved_raw)} out=${\shell_quote($trimmed)} \\
  qtrim=${\shell_quote($qtrim)} trimq=$trimq stats=${\shell_quote($qualtrim_stats)} \\
  pigz=t unpigz=t threads=$NCPU tbo
echo "---------------------- DONE quality trimming"

echo "---------------------- [4/8] Rename reads (rename.sh)"
${\shell_quote($RENAME)} $BBMEM in=${\shell_quote($trimmed)} out=${\shell_quote($trimmed_renamed)} \\
  prefix=${\shell_quote($prefix)} interleaved=true pigz=t unpigz=t threads=$NCPU
echo "---------------------- DONE renaming"

echo "---------------------- [5/8] Remove vector contamination (bbduk.sh)"
${\shell_quote($BBDUK)} $BBMEM in=${\shell_quote($trimmed_renamed)} out=${\shell_quote($no_vector)} \\
  k=$kmer ref=${\shell_quote($vectorfile)} ordered cardinality stats=${\shell_quote($vector_stats)} \\
  pigz=t unpigz=t threads=$NCPU tbo
echo "---------------------- DONE vector removal"

echo "---------------------- [6/8] Trim adapters (bbduk.sh)"
${\shell_quote($BBDUK)} $BBMEM in=${\shell_quote($no_vector)} out=${\shell_quote($no_adapt)} \\
  k=$kmer ref=${\shell_quote($adapterfile)} ordered cardinality stats=${\shell_quote($adapter_stats)} \\
  pigz=t unpigz=t threads=$NCPU tbo
echo "---------------------- DONE adapter trimming"

echo "---------------------- [7/8] De novo check for remaining adapters (bbmerge.sh)"
${\shell_quote($BBMERGE)} $BBMEM in=${\shell_quote($no_adapt)} out=${\shell_quote($adapter_after_fa)} reads=1m \\
  pigz=t unpigz=t threads=$NCPU
echo "---------------------- DONE de novo adapter identification"
echo "Adapter-like sequences saved to: $adapter_after_fa"
echo "FASTA headers:"
grep '^>' ${\shell_quote($adapter_after_fa)} || true

SCRIPT

# Optional vecscreen section (only if requested and tool exists)
if ($vecscreen) {
print $OUT <<"SCRIPT";
echo "---------------------- Optional: VecScreen on adapter-like sequences"
if command -v run_vecscreen >/dev/null 2>&1; then
  run_vecscreen -i ${\shell_quote($adapter_after_fa)} || true
  echo "VecScreen finished (see outputs in current directory, depending on your installation)."
else
  echo "run_vecscreen not found in PATH; skipping VecScreen."
fi
SCRIPT
}

print $OUT <<"SCRIPT";
echo "---------------------- [8/8] FastQC on cleaned reads"
${\shell_quote($FASTQC)} -t $NCPU -o ${\shell_quote($processed_fastqc_dir)} ${\shell_quote($no_adapt)}
echo "---------------------- DONE FastQC cleaned"

echo "---------------------- Cleanup + finalize"
rm -f ${\shell_quote($interleaved_raw)} ${\shell_quote($trimmed)} ${\shell_quote($trimmed_renamed)} ${\shell_quote($no_vector)}
mv -f ${\shell_quote($no_adapt)} ${\shell_quote($final_qf)}
echo "Final output: $final_qf"
echo "FastQC raw: $raw_fastqc_dir"
echo "FastQC cleaned: $processed_fastqc_dir"
echo "Stats: $qualtrim_stats, $vector_stats, $adapter_stats"

echo "---------------------- DONE"
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Generated workflow script: $scriptfile\n";
print "Output folder: $output\n";
print "Final output file will be: $final_qf\n";
print "FastQC outputs: $raw_fastqc_dir and $processed_fastqc_dir\n";

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
