#!/bin/bash
#SBATCH -J p3_gsize
#SBATCH -A pi-hohndor
#SBATCH -p batch
#SBATCH -t 04:00:00
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -o /ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1/logs/%x_%j.out
# Genome length, contig count and GC of each of the 5,438 representatives (refined, retained contigs),
# with CheckM2 completeness/contamination from 10_thermal/index/genome_index.tsv. 16 files in parallel.
W=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1
IDX=/ibex/scratch/projects/c2014/rob/metagenomics/10_thermal/index/genome_index.tsv
one() { awk -v g="$1" -v c="$3" -v k="$4" "/^>/{n++;next}{L+=length(\$0); gc+=gsub(/[GCgc]/,\"\")} END{printf \"%s\t%d\t%d\t%.4f\t%s\t%s\n\", g, L, n, gc/L, c, k}" "$2"; }
export -f one
tail -n +2 $IDX | awk -F"\t" "\$2==\"soil_representative\"{print \$1\"\t\"\$4\"\t\"\$7\"\t\"\$8}" > $W/results/gsize_targets.tsv
{ echo -e "genome_id\tlength_bp\tn_contigs\tgc\tcompleteness\tcontamination"
  xargs -P 16 -L 1 bash -c "one \$0 \$1 \$2 \$3" < $W/results/gsize_targets.tsv | sort; } > $W/results/rep_genome_size.tsv
wc -l $W/results/rep_genome_size.tsv; sha256sum $W/results/rep_genome_size.tsv
