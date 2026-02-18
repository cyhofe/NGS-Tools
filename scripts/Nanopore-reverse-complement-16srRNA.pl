#!/usr/bin/env perl

###############################################################################
# Nanopore-reverse-complement-16srRNA.pl
#
# Author: Cyrill Hofer
# Created: 22.04.2024
# Updated: 18.02.2026
#
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates (and optionally executes) a workflow for Nanopore 16S
# rRNA reads where read orientation can be mixed. It maps reads to a 16S
# reference database (e.g., SILVA) and splits reads by mapping strand:
#
#   - "+" strand reads are treated as "forward"
#   - "-" strand reads are treated as "reverse" and are reverse-complemented
#
# The final output is a single FASTQ.GZ containing:
#   forward reads + reverse-complemented reverse reads
# with standardized read headers via BBMap rename.sh.
#
# WORKFLOW STEPS
# -----------------------------------------------------------------------------
#  1) Map reads to reference with minimap2 (PAF output)
#  2) Identify reads that map uniquely to either + or - strand
#     (reads mapping to both strands are excluded as ambiguous)
#  3) Extract forward and reverse reads using seqtk subseq
#  4) Reverse-complement the reverse reads using seqtk seq -r
#  5) Combine forward + reverse-complemented reads
#  6) Rename reads using BBMap rename.sh
#  7) Run NanoStat on final output (optional but enabled by default)
#  8) Cleanup intermediate files
#
# REPRODUCIBILITY / HPC DESIGN
# -----------------------------------------------------------------------------
# This Perl script writes a self-contained bash script:
#   <prefix>.ReverseComplement.sh
#
# By default it ONLY generates the bash script. Use --run to execute immediately,
# or submit the bash script to your HPC scheduler.
#
# DEPENDENCIES
# -----------------------------------------------------------------------------
# Required executables (must be in PATH OR provided via *_dir options):
#  - minimap2
#  - seqtk
#  - rename.sh (BBMap)
#  - NanoStat (optional; can disable with --no_nanostat)
#
# Optional directory overrides:
#  --minimap2_dir DIR   directory containing minimap2
#  --seqtk_dir DIR      directory containing seqtk
#  --bbmap_dir DIR      directory containing rename.sh
#  --nanostat_dir DIR   directory containing NanoStat
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
  $0 -f reads.fq.gz -p PREFIX -o OUTPUT_DIR --ref reference.fasta [options]

Required:
  -f, --file FILE           Nanopore 16S reads (FASTQ/FASTQ.GZ)
  -p, --prefix PREFIX       Prefix for output files
  -o, --output DIR          Output directory (will be created)
  --ref, --reference FILE   Reference database (FASTA), e.g. SILVA

Options:
  --ncpu INT                Threads (default: 40)

  --minimap2_dir DIR        Directory containing minimap2 (optional; else PATH)
  --seqtk_dir DIR           Directory containing seqtk (optional; else PATH)
  --bbmap_dir DIR           Directory containing rename.sh (optional; else PATH)
  --nanostat_dir DIR        Directory containing NanoStat (optional; else PATH)

  --no_nanostat             Do not run NanoStat on final output

  --run                     Execute generated workflow immediately
  --force                   Allow existing output folder (may overwrite files)
  -h, --help                Show help

Example:
  $0 -f reads.fq.gz -p sample16S -o sample16S.rc --ref SILVA_138.1.fasta --run

USAGE

# ----------------------------
# Args + defaults
# ----------------------------
my ($file,$prefix,$output,$reference,$help);
my ($NCPU);
my ($minimap2_dir,$seqtk_dir,$bbmap_dir,$nanostat_dir);
my ($run,$force);
my $do_nanostat = 1;

$NCPU = 40;

GetOptions(
  "f|file=s"        => \$file,
  "p|prefix=s"      => \$prefix,
  "o|output=s"      => \$output,
  "ref|reference=s" => \$reference,
  "ncpu=i"          => \$NCPU,

  "minimap2_dir=s"  => \$minimap2_dir,
  "seqtk_dir=s"     => \$seqtk_dir,
  "bbmap_dir=s"     => \$bbmap_dir,
  "nanostat_dir=s"  => \$nanostat_dir,

  "no_nanostat!"    => sub { $do_nanostat = 0 },

  "run!"            => \$run,
  "force!"          => \$force,
  "h|help"          => \$help,
) or die $USAGE;

die $USAGE if $help;

die "ERROR: Missing -f/--file\n$USAGE"     unless $file;
die "ERROR: Missing -p/--prefix\n$USAGE"   unless $prefix;
die "ERROR: Missing -o/--output\n$USAGE"   unless $output;
die "ERROR: Missing --ref/--reference\n$USAGE" unless $reference;

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
my $MINIMAP2  = resolve_exec($minimap2_dir, "minimap2", "minimap2");
my $SEQTK     = resolve_exec($seqtk_dir,    "seqtk",    "seqtk");
my $RENAME    = resolve_exec($bbmap_dir,    "rename.sh","bbmap");
my $NANOSTAT  = $do_nanostat ? resolve_exec($nanostat_dir, "NanoStat", "nanostat") : "";

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
my $scriptfile = "$prefix.ReverseComplement.sh";

my $paf        = "$output/$prefix.paf";
my $ambig_ids  = "$output/$prefix.AmbiguousReads.txt";
my $fwd_ids    = "$output/$prefix.ForwardReads.txt";
my $rev_ids    = "$output/$prefix.ReverseReads.txt";

my $fwd_fq     = "$output/$prefix.ForwardReads.fq.gz";
my $rev_fq     = "$output/$prefix.ReverseReads.fq.gz";
my $rev_rc_fq  = "$output/$prefix.ReverseComplementedReads.fq.gz";
my $combined   = "$output/$prefix.RevComp.fq.gz";
my $final_rc   = "$output/$prefix.rc.fq.gz";
my $nanostat_out = "$output/$prefix.rc.NanoStat.txt";

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
# Reference: $reference
# Prefix: $prefix
# Output: $output

echo "---------------------- [1/7] Map reads to reference (minimap2 -> PAF)"
${\shell_quote($MINIMAP2)} -x map-ont ${\shell_quote($reference)} ${\shell_quote($file)} -t $NCPU > ${\shell_quote($paf)}
echo "---------------------- DONE mapping"

echo "---------------------- [2/7] Identify strand-specific reads"
# PAF fields: qname (col1), strand (col5). We exclude reads mapping to both + and - (ambiguous).
cut -f1,5 ${\shell_quote($paf)} | sort -u > ${\shell_quote("$output/$prefix.qname_strand.tsv")}

cut -f1 ${\shell_quote("$output/$prefix.qname_strand.tsv")} | sort | uniq -d > ${\shell_quote($ambig_ids)}

grep -vFf ${\shell_quote($ambig_ids)} ${\shell_quote("$output/$prefix.qname_strand.tsv")} | awk '\$2=="+"' | cut -f1 > ${\shell_quote($fwd_ids)}
grep -vFf ${\shell_quote($ambig_ids)} ${\shell_quote("$output/$prefix.qname_strand.tsv")} | awk '\$2=="-"' | cut -f1 > ${\shell_quote($rev_ids)}

echo "Forward IDs: \$(wc -l < ${\shell_quote($fwd_ids)} || echo 0)"
echo "Reverse IDs: \$(wc -l < ${\shell_quote($rev_ids)} || echo 0)"
echo "Ambiguous IDs excluded: \$(wc -l < ${\shell_quote($ambig_ids)} || echo 0)"
echo "---------------------- DONE fetching ID lists"

echo "---------------------- [3/7] Extract forward/reverse reads (seqtk subseq)"
${\shell_quote($SEQTK)} subseq ${\shell_quote($file)} ${\shell_quote($fwd_ids)} | gzip > ${\shell_quote($fwd_fq)}
${\shell_quote($SEQTK)} subseq ${\shell_quote($file)} ${\shell_quote($rev_ids)} | gzip > ${\shell_quote($rev_fq)}
echo "---------------------- DONE extracting reads"

echo "---------------------- [4/7] Reverse-complement reverse reads (seqtk seq -r)"
${\shell_quote($SEQTK)} seq -r ${\shell_quote($rev_fq)} | gzip > ${\shell_quote($rev_rc_fq)}
echo "---------------------- DONE reverse complementing"

echo "---------------------- [5/7] Combine forward + reverse-complemented reads"
cat ${\shell_quote($fwd_fq)} ${\shell_quote($rev_rc_fq)} > ${\shell_quote($combined)}
echo "---------------------- DONE combining"

echo "---------------------- [6/7] Rename reads (BBMap rename.sh)"
${\shell_quote($RENAME)} in=${\shell_quote($combined)} out=${\shell_quote($final_rc)} \\
  prefix=${\shell_quote($prefix)} interleaved=false pigz=t unpigz=t threads=$NCPU
echo "---------------------- DONE renaming"

SCRIPT

if ($do_nanostat) {
print $OUT <<"SCRIPT";
echo "---------------------- [7/7] NanoStat on final output"
${\shell_quote($NANOSTAT)} --fastq ${\shell_quote($final_rc)} -t $NCPU > ${\shell_quote($nanostat_out)}
echo "---------------------- DONE NanoStat"
SCRIPT
} else {
print $OUT <<"SCRIPT";
echo "---------------------- [7/7] NanoStat disabled (--no_nanostat)"
SCRIPT
}

print $OUT <<"SCRIPT";

echo "---------------------- Cleanup"
rm -f ${\shell_quote("$output/$prefix.qname_strand.tsv")}
rm -f ${\shell_quote($fwd_ids)} ${\shell_quote($rev_ids)} ${\shell_quote($ambig_ids)}
rm -f ${\shell_quote($fwd_fq)} ${\shell_quote($rev_fq)} ${\shell_quote($rev_rc_fq)} ${\shell_quote($paf)} ${\shell_quote($combined)}
echo "Final output: $final_rc"
echo "---------------------- DONE"
SCRIPT

close $OUT;
chmod 0755, $scriptfile;

print "Generated workflow script: $scriptfile\n";
print "Output folder: $output\n";
print "Final output file will be: $final_rc\n";

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
