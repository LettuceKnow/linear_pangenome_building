configfile: "config.yaml"

rule all:
    input:
        "output/multiqc/multiqc_report.html", #mapping stats genome
        "output/genome.pav.pdf",
        expand("output/filtering_assembly/{accession}/possible_contaminants.txt", accession=config["data"]),

        "output/busco/pangenome/short_summary.specific.{lineage}.pangenome.txt".format(lineage=config["busco"]), #pangenome BUSCO
        "output/pangenome_blobplot/pangenome.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png", #blobplot before final filtering of GC content and blastn contamination
        "output/pangenome.fa", #pangenome sequence
        "output/pangenome.gff", #pangenome annotation

        "output/pangenome_multiqc/multiqc_report.html", #mapping stats pangenome
        "output/pangenome.pav.pdf",

rule fastQC_before_trimming:
    input:
        lambda wildcards: config["data"][wildcards.accession][wildcards.sample]
    output:
        "output/fastQC_before_trimming/{accession}/{sample}_1_fastqc.html",
        "output/fastQC_before_trimming/{accession}/{sample}_2_fastqc.html",
    log:
        "output/logs/fastQC_before_trimming/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/fastQC_before_trimming/{accession}/{sample}.txt"
    conda:
        "envs/qualitycontrol.yaml"
    threads:
        5
    shell:
        "fastqc -t {threads} -o `dirname {output} | uniq` {input} &> {log}"

rule trimmomatic:
    input:
        lambda wildcards: config["data"][wildcards.accession][wildcards.sample]
    output:
        #forward=temp("output/trimmed_reads/{accession}/{sample}_1.fq.gz"),
        #backward=temp("output/trimmed_reads/{accession}/{sample}_2.fq.gz"),
        forward="output/trimmed_reads/{accession}/{sample}_1.fq.gz",
        backward="output/trimmed_reads/{accession}/{sample}_2.fq.gz",
        trimlog="output/trimmed_reads/{accession}/{sample}.trimlog"
    log:
        "output/logs/trimmomatic/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/trimmomatic/{accession}/{sample}.txt"
    conda:
        "envs/qualitycontrol.yaml"
    threads:
        5
    params:
        adapter=config["adapter"],
        illuminaclip="2:30:10",
        leading=3,
        trailing=3,
        slidingwindow="4:15",
        minlen=36
    shell:
        "trimmomatic PE -threads {threads} -trimlog {output.trimlog} {input} {output.forward} /dev/null {output.backward} /dev/null "
        "ILLUMINACLIP:{params.adapter}:{params.illuminaclip} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen} &> {log}"

rule fastQC_after_trimming:
    input:
        "output/trimmed_reads/{accession}/{sample}_1.fq.gz",
        "output/trimmed_reads/{accession}/{sample}_2.fq.gz",
    output:
        "output/fastQC_after_trimming/{accession}/{sample}_1_fastqc.html",
        "output/fastQC_after_trimming/{accession}/{sample}_2_fastqc.html",
    log:
        "output/logs/fastQC_after_trimming/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/fastQC_after_trimming/{accession}/{sample}.txt"
    conda:
        "envs/qualitycontrol.yaml"
    threads:
        5
    shell:
        "fastqc -t {threads} -o `dirname {output} | uniq` {input} &> {log}"

rule get_genome:
    input:
        config["genome"]
    output:
        fa="output/data/genome.fa", #sequence
        fai="output/data/genome.fa.fai", #sequence index
        tsv="output/data/genome.tsv", #sequence ID + length
        txt="output/data/genome.txt", #IDs
    log:
        "output/logs/get_genome.log"
    benchmark:
        "output/benchmarks/get_genome.txt"
    conda:
        "envs/alignment.yaml"
    shell:
        """
        (
        cp {input} {output.fa};
        samtools faidx {output.fa};
        awk 'BEGIN{{FS = OFS = \"\\t\";}} {{print $1,$2;}}' {output.fai} > {output.tsv};
        grep '^>' {input} | cut -d ' ' -f 1 | cut -c 2- > {output.txt};
        ) &> {log}
        """

rule correct_annotation:
    input:
        gff=config["annotation"],
    output:
        temporary("output/data/genome.gff.tmp"),
    log:
        "output/logs/correct_annotation.log"
    benchmark:
        "output/benchmarks/correct_annotation.txt"
    conda:
        "envs/agat.yaml"
    shell:
        "agat_convert_sp_gxf2gxf.pl -g {input} -o {output} &> {log}"

def genome_annotation_input(wildcards):
    if config["correct_annotation"] == True:
        return ["output/data/genome.gff.tmp"] #use AGAT to fix gff input
    else:
        return config["annotation"]

rule get_annotation:
    input:
        gff=genome_annotation_input,
        genome="output/data/genome.tsv",
    output:
        gff="output/data/genome.gff", #annotation
        sortedgff="output/data/genome.sorted.gff", #sorted annotation
        table="output/data/genome.id_to_type.tsv", #conversion table between ID and type in annotation
    log:
        "output/logs/get_annotation.log"
    benchmark:
        "output/benchmarks/get_annotation.txt"
    params:
        exclude=config["exclude"]
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        (
        exclude=$(echo {params.exclude} | sed 's/ /|/g');
        grep -vP \"^(${{exclude}})\" {input.gff} | awk 'BEGIN{{FS = OFS = \"\\t\"; print \"##gff-version 3\";}} $1!~/^#/' > {output.gff};
        bedtools sort -g {input.genome} -i {output.gff} > {output.sortedgff};
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"ID\",\"type\";}} /^[^#]/{{split($9,a,\";\"); for (i in a){{if (a[i]~/^ID/){{split(a[i],b,\"=\"); id=b[2];}}}}; print id,$3;}}' {output.sortedgff} > {output.table};
        ) &> {log}
        """

rule bwa_index:
    input:
        "output/data/genome.fa"
    output:
        "output/data/genome.fa.0123",
        "output/data/genome.fa.amb",
        "output/data/genome.fa.ann",
        "output/data/genome.fa.pac"
    log:
        "output/logs/bwa_index.log"
    benchmark:
        "output/benchmarks/bwa_index.txt"
    conda:
        "envs/alignment.yaml"
    shell:
        "bwa-mem2 index -p {input} {input} &> {log}"

rule bwa_mapping:
    input:
        index="output/data/genome.fa",
        index1="output/data/genome.fa.0123",
        index2="output/data/genome.fa.amb",
        index3="output/data/genome.fa.ann",
        index4="output/data/genome.fa.pac",
        forward="output/trimmed_reads/{accession}/{sample}_1.fq.gz",
        backward="output/trimmed_reads/{accession}/{sample}_2.fq.gz"
    output:
        #temp("output/alignment/{accession}/{sample}.bam")
        "output/alignment/{accession}/{sample}.bam"
    log:
        "output/logs/alignment/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/alignment/{accession}/{sample}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "(bwa-mem2 mem -t {threads} {input.index} {input.forward} {input.backward} | samtools sort -@ $(({threads}-1)) -o {output} -) &> {log}"

rule alignment_index:
    input:
        "output/alignment/{accession}/{sample}.bam"
    output:
        #temp("output/alignment/{accession}/{sample}.bam.bai")
        "output/alignment/{accession}/{sample}.bam.bai"
    log:
        "output/logs/alignment_index/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/alignment_index/{accession}/{sample}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools index -@ $(({threads}-1)) {input} > {output} 2> {log}"

rule alignment_stats:
    input:
        bam="output/alignment/{accession}/{sample}.bam",
        bai="output/alignment/{accession}/{sample}.bam.bai"
    output:
        "output/alignment/{accession}/{sample}.txt"
    log:
        "output/logs/alignment_stats/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/alignment_stats/{accession}/{sample}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools stats -@ $(({threads}-1)) {input.bam} > {output} 2> {log}"

def get_multiqc_input():
    """This function returns all fastqc, alignment stats and duplication stats"""
    combinations = []
    for accession in list(config["data"].keys()):
        combination=expand(f"{accession}/{{sample}}", sample=[x for x in config["data"][accession]])
        combinations.append(combination)
    combinations = [x for y in combinations for x in y]

    allfiles = []
    allfiles.append(expand("output/fastQC_before_trimming/{combination}_{nr}_fastqc.html", combination=combinations, nr=[1,2]))
    allfiles.append(expand("output/fastQC_after_trimming/{combination}_{nr}_fastqc.html", combination=combinations, nr=[1,2]))
    allfiles.append(expand("output/alignment/{combination}.txt", combination=combinations))

    return allfiles

rule multiqc:
    input:
        *get_multiqc_input()
    output:
        "output/multiqc/multiqc_report.html"
    log:
        "output/logs/multiqc.log"
    benchmark:
        "output/benchmarks/multiqc.txt"
    conda:
        "envs/qualitycontrol.yaml"
    shell:
        "multiqc -f -s -d `dirname {input} | uniq` -n {output} &> {log}"

rule alignment_merge:
    input:
        bam=lambda wildcards: expand("output/alignment/{{accession}}/{sample}.bam", sample=[x for x in config["data"][wildcards.accession]]),
        bai=lambda wildcards: expand("output/alignment/{{accession}}/{sample}.bam.bai", sample=[x for x in config["data"][wildcards.accession]]),
    output:
        #temp("output/alignment_merged/{accession}.bam")
        "output/alignment_merged/{accession}.bam"
    log:
        "output/logs/alignment_merge/{accession}.log"
    benchmark:
        "output/benchmarks/alignment_merge/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools merge -@ $(({threads}-1)) {output} {input.bam} 2> {log}"

rule alignment_merge_index:
    input:
        "output/alignment_merged/{accession}.bam"
    output:
        #temp("output/alignment_merged/{accession}.bam.bai")
        "output/alignment_merged/{accession}.bam.bai"
    log:
        "output/logs/alignment_merge_index/{accession}.log"
    benchmark:
        "output/benchmarks/alignment_merge_index/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools index -@ $(({threads}-1)) {input} 2> {log}"

rule calculate_breadth:
    input:
        genome="output/data/genome.tsv",
        annotation="output/data/genome.sorted.gff",
        bam="output/alignment_merged/{accession}.bam",
        bai="output/alignment_merged/{accession}.bam.bai",
    output:
        allcov="output/coverage/breadth/{accession}.all.bed",
        exoncov="output/coverage/breadth/{accession}.exons.tsv",
    log:
        "output/logs/coverage_breadth/{accession}.log"
    benchmark:
        "output/benchmarks/coverage_breadth/{accession}.txt"
    conda:
        "envs/bedtools.yaml"
    threads:
        5
    shell:
        """
        (
        bedtools coverage -g {input.genome} -sorted -bed -a {input.annotation} -b {input.bam} 1> {output.allcov};
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"parent\",\"number exons\",\"total length exons\",\"total length exons covered\",\"fraction length exons covered\";}} $3==\"exon\"{{split($9,a,\";\"); for (i in a){{if (a[i]~/^Parent/){{split(a[i],tmp,\"=\"); parent=tmp[2];}}}} num[parent]++; len[parent]+=$12; covered[parent]+=$11;}} END{{for (parent in len){{print parent,num[parent],len[parent],covered[parent],covered[parent]/len[parent];}}}}' {output.allcov} > {output.exoncov};
        ) &> {log}
        """

rule calculate_depth:
    input:
        genome="output/data/genome.tsv",
        annotation="output/data/genome.sorted.gff",
        bam="output/alignment_merged/{accession}.bam",
        bai="output/alignment_merged/{accession}.bam.bai",
    output:
        allcov="output/coverage/depth/{accession}.all.bed",
        exoncov="output/coverage/depth/{accession}.exons.tsv",
    log:
        "output/logs/coverage_depth/{accession}.log"
    benchmark:
        "output/benchmarks/coverage_depth/{accession}.txt"
    conda:
        "envs/bedtools.yaml"
    threads:
        5
    shell:
        """
        (
        bedtools coverage -g {input.genome} -mean -sorted -bed -a {input.annotation} -b {input.bam} 1> {output.allcov};
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"parent\",\"number exons\",\"total length exons\",\"mean depth exons\";}} $3==\"exon\"{{split($9,a,\";\"); for (i in a){{if (a[i]~/^Parent/){{split(a[i],tmp,\"=\"); parent=tmp[2];}}}} num[parent]++; len[parent]+=($5-$4+1); normdepth[parent]+=($10*($5-$4+1));}} END{{for (parent in len){{print parent,num[parent],len[parent],normdepth[parent]/len[parent];}}}}' {output.allcov} > {output.exoncov};
        ) &> {log}
        """

rule combine_breadth_depth:
    input:
        breadth="output/coverage/breadth/{accession}.exons.tsv",
        depth="output/coverage/depth/{accession}.exons.tsv",
    output:
        "output/coverage/combined/{accession}.tsv"
    log:
        "output/logs/combine_breadth_depth/{accession}.log"
    benchmark:
        "output/benchmarks/combine_breadth_depth/{accession}.txt"
    shell:
        """
        (
        awk 'BEGIN{{FS = OFS = \"\\t\";}} NR==FNR{{a[$1] = $4; next;}} NR!=1{{print $0,a[$1];}}' {input.depth} {input.breadth} > {output}
        ) &> {log}
        """

rule calculate_pav:
    input:
        pav=expand("output/coverage/combined/{accession}.tsv", accession=config["data"]),
        conversion="output/data/genome.id_to_type.tsv",
    output:
        depth="output/genome.depth.tsv",
        breadth="output/genome.breadth.tsv",
        plot="output/genome.pav.pdf",
    log:
        "output/logs/calculate_pav.log"
    benchmark:
        "output/benchmarks/calculate_pav.txt"
    conda:
        "envs/rbase.yaml"
    shell:
        "../scripts/combine_pav.R {output.depth} {output.breadth} {output.plot} {input.conversion} {input.pav} &> {log}"

rule get_unmapped_reads:
    input:
        bam="output/alignment_merged/{accession}.bam",
        bai="output/alignment_merged/{accession}.bam.bai",
    output:
        forward="output/reads_unmapped/{accession}_1.fq.gz",
        backward="output/reads_unmapped/{accession}_2.fq.gz",
    log:
        "output/logs/get_unmapped_reads/{accession}.log"
    benchmark:
        "output/benchmarks/get_unmapped_reads/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools fastq -@ $(({threads}-1)) -Nf4 -1 {output.forward} -2 {output.backward} -0 /dev/null -s /dev/null {input.bam} 2> {log}"

rule assembly_unmapped_reads:
    input:
        forward="output/reads_unmapped/{accession}_1.fq.gz",
        backward="output/reads_unmapped/{accession}_2.fq.gz",
    output:
        "output/assembly_unmapped/{accession}/final.contigs.fa"
    log:
        "output/logs/assembly_unmapped_reads/{accession}.log"
    benchmark:
        "output/benchmarks/assembly_unmapped_reads/{accession}.txt"
    conda:
        "envs/assembly.yaml"
    threads:
        10
    shell:
        "rm -rd `dirname {output}` && megahit -t {threads} -1 {input.forward} -2 {input.backward} -o `dirname {output}` &> {log}"

rule assembly_rename:
    input:
        "output/assembly_unmapped/{accession}/final.contigs.fa"
    output:
        "output/filtering_assembly/{accession}/contigs_filter_renamed.fa"
    log:
        "output/logs/assembly_rename/{accession}.log"
    benchmark:
        "output/benchmarks/assembly_rename/{accession}.txt"
    shell:
        "awk -vid={wildcards.accession} '/^>/{{num+=1; split($1,a,\">\"); $1=\">\"id\"_\"num\" orig=\"a[2];}} {{print;}}' {input} > {output} 2> {log}"

rule assembly_unmapped_stats:
    input:
        "output/filtering_assembly/{accession}/contigs_filter_renamed.fa"
    output:
        "output/filtering_assembly/{accession}/contigs_filter_renamed.stats.tsv"
    log:
        "output/logs/assembly_unmapped_stats/{accession}.log"
    benchmark:
        "output/benchmarks/assembly_unmapped_stats/{accession}.txt"
    conda:
        "envs/seq.yaml"
    shell:
        "seqkit fx2tab -igl {input} | cut -f 1,4- > {output} 2> {log}"

rule assembly_unmapped_filter_length_gc:
    input:
        assembly="output/filtering_assembly/{accession}/contigs_filter_renamed.fa",
        stats="output/filtering_assembly/{accession}/contigs_filter_renamed.stats.tsv",
    output:
        assembly="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.fa",
        to_keep="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.txt",
    log:
        "output/logs/assembly_filter_length_gc/{accession}.log"
    benchmark:
        "output/benchmarks/assembly_filter_length_gc/{accession}.txt"
    conda:
        "envs/seq.yaml"
    params:
        minlen=config["minlength"],
        maxgc=config["maxgc"],
    shell:
        """
        (
        awk 'BEGIN{{FS = OFS = \"\\t\"}} $2 > {params.minlen} && $3 < {params.maxgc} {{print $1;}}' {input.stats} > {output.to_keep};
        seqkit grep -f {output.to_keep} {input.assembly} > {output.assembly};
        ) 2> {log}
        """

rule run_kraken2:
    input:
        fasta="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.fa",
    output:
        out="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.out",
        report="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.report.txt",
    log:
        "output/logs/run_kraken2/{accession}.log"
    benchmark:
        "output/benchmarks/run_kraken2/{accession}.txt"
    params:
        db=config["kraken_db"]
    threads:
        10
    conda:
        "envs/kraken.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --output {output.out} --report {output.report} {input} > {output} 2> {log}"

rule visualise_kraken2:
    input:
        "output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.out",
    output:
        "output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.krona.html"
    log:
        "output/logs/visualise_kraken2/{accession}.log"
    benchmark:
        "output/benchmarks/visualise_kraken2/{accession}.txt"
    conda:
        "envs/kraken.yaml"
    shell:
        """
        (
        output={output};
        krona=${{output%.html}};
        cut -f 2,3 {input} > $krona;
        ktImportTaxonomy $krona -o {output};
        ) &> {log}
        """

rule filter_kraken2:
    input:
        fasta="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.fa",
        krakenout="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.out",
        krakenreport="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.report.txt",
        visualisation="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc.kraken2.krona.html", #not used computationally, but needed for visual inspection
    output:
        fasta="output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc_contamination.fa",
        contaminants="output/filtering_assembly/{accession}/possible_contaminants.txt",
    log:
        "output/logs/filter_kraken2/{accession}.log"
    benchmark:
        "output/benchmarks/filter_kraken2/{accession}.txt"
    params:
        taxid=config["taxid_filtering"]
    conda:
        "envs/kraken.yaml"
    shell:
        """
        (
        awk '$1>1 && $4=="S"' {input.krakenreport} > {output.contaminants};
        extract_kraken_reads.py -k {input.krakenout} -r {input.krakenreport} -s {input.fasta} --include-children -o {output.fasta} --exclude -t {params.taxid}
        ) &> {log}
        """

def get_current_pangenome(wildcards):
    n = int(wildcards.n)
    if n == 1:
        return config["genome"]
    else:
        return f"output/iterative_build/{n-1}/current_pangenome.fa"

rule iterative_build_pangenome:
    input:
        current_pangenome=get_current_pangenome,
        addition=lambda wildcards: "output/filtering_assembly/{accession}/contigs_filter_renamed_length_gc_contamination.fa".format(accession=[x for x in config["data"]][int(wildcards.n)-1])
    output:
        paf="output/iterative_build/{n}/addition_against_previous_pangenome.paf",
        mapped="output/iterative_build/{n}/mapped_sequences_against_previous_pangenome.txt",
        next_iteration="output/iterative_build/{n}/current_pangenome.fa",
    log:
        "output/logs/iterative_build_pangenome/{n}.log"
    benchmark:
        "output/benchmarks/iterative_build_pangenome/{n}.txt"
    conda:
        "envs/iterative_pangenome.yaml"
    threads:
        min(workflow.cores, 10)
    shell:
        """
        (
        minimap2 -cxasm5 --secondary=no -t {threads} {input.current_pangenome} {input.addition} > {output.paf};
        awk '$2<$7{{print $1;}} $2>=$7{{print $6;}}' {output.paf} | sort | uniq > {output.mapped};
        cat {input.current_pangenome} {input.addition} | seqkit grep -vf {output.mapped} > {output.next_iteration};
        ) 2> {log}
        """

rule pangenome_get_new_sequences:
    input:
        pan="output/iterative_build/{n}/current_pangenome.fa".format(n=len([x for x in config["data"]])),
        ref_ID="output/data/genome.txt", #IDs
    output:
        "output/pangenome_filtering/pangenome.new.fa",
    log:
        "output/logs/pangenome_get_new_sequences.log"
    benchmark:
        "output/benchmarks/pangenome_get_new_sequences.txt"
    conda:
        "envs/seq.yaml"
    shell:
        "seqkit grep -vf {input.ref_ID} {input.pan} > {output} 2> {log}"

rule pangenome_filter_duplicates:
    input:
        "output/pangenome_filtering/pangenome.new.fa",
    output:
        new="output/pangenome_filtering/pangenome.new.clustering.fa",
        clstr="output/pangenome_filtering/pangenome.new.clustering.fa.clstr",
    log:
        "output/logs/pangenome_filter_duplicates.log"
    benchmark:
        "output/benchmarks/pangenome_filter_duplicates.txt"
    conda:
        "envs/cd-hit.yaml"
    params:
        memory=10000, #Uses 10G memory
        similarity=0.9,
    threads:
        workflow.cores
    shell:
        "cd-hit-est -T {threads} -M {params.memory} -c {params.similarity} -i {input} -o {output.new} &> {log}"

rule pangenome_blastn:
    input:
        "output/pangenome_filtering/pangenome.new.clustering.fa",
    output:
        "output/pangenome_filtering/pangenome.new.clustering.blastn.out",
    log:
        "output/logs/pangenome_blastn.log"
    benchmark:
        "output/benchmarks/pangenome_blastn.txt"
    conda:
        "envs/busco.yaml"
    params:
        nt=config["nt"],
    threads:
        max(workflow.cores - 1, 10)
    shell:
        "blastn -query {input} -db {params.nt} -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -num_threads {threads} -evalue 1e-25 -out {output} &> {log}"

rule pangenome_new_bwa_index:
    input:
        "output/pangenome_filtering/pangenome.new.clustering.fa",
    output:
        "output/pangenome_filtering/pangenome.new.clustering.fa.0123",
        "output/pangenome_filtering/pangenome.new.clustering.fa.amb",
        "output/pangenome_filtering/pangenome.new.clustering.fa.ann",
        "output/pangenome_filtering/pangenome.new.clustering.fa.pac",
    log:
        "output/logs/pangenome_new_bwa_index.log"
    benchmark:
        "output/benchmarks/pangenome_new_bwa_index.txt"
    conda:
        "envs/alignment.yaml"
    shell:
        "bwa-mem2 index -p {input} {input} &> {log}"

rule pangenome_new_bwa_mapping:
    input:
        pangenome="output/pangenome_filtering/pangenome.new.clustering.fa",
        forward="output/reads_unmapped/{accession}_1.fq.gz",
        backward="output/reads_unmapped/{accession}_2.fq.gz",
        index1="output/pangenome_filtering/pangenome.new.clustering.fa.0123",
        index2="output/pangenome_filtering/pangenome.new.clustering.fa.amb",
        index3="output/pangenome_filtering/pangenome.new.clustering.fa.ann",
        index4="output/pangenome_filtering/pangenome.new.clustering.fa.pac",
    output:
        #temp("output/pangenome_filtering/alignment/{accession}.sorted.bam")
        "output/pangenome_filtering/alignment/{accession}.sorted.bam"
    log:
        "output/logs/pangenome_new_bwa_mapping/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_new_bwa_mapping/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "(bwa-mem2 mem -t {threads} {input.pangenome} {input.forward} {input.backward} | samtools sort -@ $(({threads}-1)) -o {output} -) &> {log}"

rule pangenome_new_bwa_mapping_index:
    input:
        "output/pangenome_filtering/alignment/{accession}.sorted.bam",
    output:
        "output/pangenome_filtering/alignment/{accession}.sorted.bam.bai",
    log:
        "output/logs/pangenome_new_bwa_mapping_index/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_new_bwa_mapping_index/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools index -@ $(({threads}-1)) {input} &> {log}"

rule pangenome_new_merge_map_back:
    input:
        bam=lambda wildcards: expand("output/pangenome_filtering/alignment/{accession}.sorted.bam", accession=config["data"]),
        bai=lambda wildcards: expand("output/pangenome_filtering/alignment/{accession}.sorted.bam.bai", accession=config["data"]),
    output:
        "output/pangenome_filtering/alignment.merged.sorted.bam",
    log:
        "output/logs/pangenome_new_merge_map_back.log"
    benchmark:
        "output/benchmarks/pangenome_new_merge_map_back.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools merge -@ $(({threads}-1)) {output} {input.bam} &> {log}"

rule pangenome_new_index_map_back:
    input:
        "output/pangenome_filtering/alignment.merged.sorted.bam",
    output:
        "output/pangenome_filtering/alignment.merged.sorted.bam.bai",
    log:
        "output/logs/pangenome_new_index_map_back.log"
    benchmark:
        "output/benchmarks/pangenome_new_index_map_back.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools index -@ $(({threads}-1)) {input} &> {log}"

rule pangenome_blobtools_create:
    input:
        fasta="output/pangenome_filtering/pangenome.new.clustering.fa",
        hits="output/pangenome_filtering/pangenome.new.clustering.blastn.out",
        bam="output/pangenome_filtering/alignment.merged.sorted.bam",
        bai="output/pangenome_filtering/alignment.merged.sorted.bam.bai",
    output:
        cov="output/pangenome_blobplot/pangenome.alignment.merged.sorted.bam.cov",
        json="output/pangenome_blobplot/pangenome.blobDB.json",
    log:
        "output/benchmarks/pangenome_blobtools_create.log"
    benchmark:
        "output/benchmarks/pangenome_blobtools_create.txt"
    conda:
        "envs/blobtools.yaml"
    shell:
        "blobtools create -i {input.fasta} -b {input.bam} -t {input.hits} -o $(echo {output.json} | sed 's/.blobDB.json//g') &> {log}"

rule pangenome_blobtools_view:
    input:
        json="output/pangenome_blobplot/pangenome.blobDB.json",
    output:
        "output/pangenome_blobplot/pangenome.blobDB.table.txt",
    log:
        "output/logs/pangenome_blobtools_view.log"
    benchmark:
        "output/benchmarks/pangenome_blobtools_view.txt"
    conda:
        "envs/blobtools.yaml"
    shell:
        "blobtools view -i {input.json} -o $(dirname {output})/ &> {log}"

rule pangenome_blobtools_plot:
    input:
        json="output/pangenome_blobplot/pangenome.blobDB.json",
        table="output/pangenome_blobplot/pangenome.blobDB.table.txt",
    output:
        file1="output/pangenome_blobplot/pangenome.blobDB.json.bestsum.phylum.p8.span.100.blobplot.stats.txt",
        file2="output/pangenome_blobplot/pangenome.blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.bam0.png",
        file3="output/pangenome_blobplot/pangenome.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png",
    log:
        "output/logs/pangenome_blobtools_plot.log"
    benchmark:
        "output/benchmarks/pangenome_blobtools_plot.txt"
    conda:
        "envs/blobtools.yaml"
    shell:
        "blobtools plot -i {input.json} -o $(dirname {output.file1})/ &> {log}"

rule pangenome_filter_blobplot:
    input:
        assembly="output/pangenome_filtering/pangenome.new.clustering.fa",
        table="output/pangenome_blobplot/pangenome.blobDB.table.txt",
    output:
        assembly="output/pangenome_filtering/pangenome.new.clustering.filtered.fa",
        to_keep="output/pangenome_filtering/pangenome.new.clustering.filtered.txt",
    log:
        "output/logs/pangenome_filter_blobplot.log"
    benchmark:
        "output/benchmarks/pangenome_filter_blobplot.txt"
    params:
        taxid=config["blastn_taxid"],
        mincov=config["mincov"],
    conda:
        "envs/seq.yaml"
    shell: #basically, filter out based on taxonomy and coverage
        """
        (
        awk 'BEGIN{{FS = OFS = \"\\t\";}} /^[^#]/ && ($6==\"{params.taxid}\" || $6==\"no-hit\") && $5 > {params.mincov}{{print $1;}}' {input.table} > {output.to_keep};
        seqkit grep -f {output.to_keep} {input.assembly} > {output.assembly};
        ) 2> {log}
        """

rule panproteome_filtering_prep:
    input:
        new="output/pangenome_filtering/pangenome.new.clustering.filtered.fa", #the sequences to add
        full_pan="output/iterative_build/{n}/current_pangenome.fa".format(n=len([x for x in config["data"]])), #for retrieving the reference genome
        ref_ID="output/data/genome.txt", #IDs that belong to reference genome
    output:
        full="output/panproteome_filtering/pangenome.fa",
        new="output/panproteome_filtering/pangenome.new.fa",
    log:
        "output/logs/panproteome_filtering_prep.log"
    benchmark:
        "output/benchmarks/panproteome_filtering_prep.txt"
    conda:
        "envs/seq.yaml"
    shell: #first get genome sequences (from pangenome, just to be sure); then add novel sequences
        """
        (
        seqkit grep -f {input.ref_ID} {input.full_pan} > {output.full};
        cat {input.new} >> {output.full};
        cp {input.new} {output.new};
        ) &> {log}
        """

rule panproteome_genome:
    input:
        "output/panproteome_filtering/pangenome.fa",
    output:
        fai="output/panproteome_filtering/pangenome.fa.fai",
        tsv="output/panproteome_filtering/pangenome.tsv",
    log:
        "output/logs/panproteome_genome.log"
    benchmark:
        "output/benchmarks/panproteome_genome.txt"
    conda:
        "envs/alignment.yaml"
    shell:
        """
        (
        samtools faidx {input};
        awk 'BEGIN{{FS = OFS = \"\\t\";}} {{print $1,$2;}}' {output.fai} > {output.tsv};
        ) &> {log}
        """

rule pangenome_repeat_database:
    input:
        "output/panproteome_filtering/pangenome.fa",
    output:
        "output/panproteome_filtering/pangenome_xdb.nhr",
        "output/panproteome_filtering/pangenome_xdb.nin",
        "output/panproteome_filtering/pangenome_xdb.nnd",
        "output/panproteome_filtering/pangenome_xdb.nni",
        "output/panproteome_filtering/pangenome_xdb.nog",
        "output/panproteome_filtering/pangenome_xdb.nsq",
        "output/panproteome_filtering/pangenome_xdb.translation",
    log:
        "output/logs/pangenome_repeat_database.log"
    benchmark:
        "output/benchmarks/pangenome_repeat_database.txt"
    container:
        "docker://harbor.containers.wurnet.nl/proxy-cache/dfam/tetools:1.5" #Using proxy to prevent too many requests to docker hub
    shell:
        "BuildDatabase -name $(echo {output} | cut -d '.' -f 1 | uniq) {input} &> {log}"

rule pangenome_repeat_modeler:
    input:
        xdb1="output/panproteome_filtering/pangenome_xdb.nhr",
        xdb2="output/panproteome_filtering/pangenome_xdb.nin",
        xdb3="output/panproteome_filtering/pangenome_xdb.nnd",
        xdb4="output/panproteome_filtering/pangenome_xdb.nni",
        xdb5="output/panproteome_filtering/pangenome_xdb.nog",
        xdb6="output/panproteome_filtering/pangenome_xdb.nsq",
        xdb7="output/panproteome_filtering/pangenome_xdb.translation",
    output:
        "output/panproteome_filtering/pangenome_xdb-families.fa"
    log:
        "output/logs/pangenome_repeat_modeler.log"
    benchmark:
        "output/benchmarks/pangenome_repeat_modeler.txt"
    container:
        "docker://harbor.containers.wurnet.nl/proxy-cache/dfam/tetools:1.5" #Using proxy to prevent too many requests to docker hub
    threads:
        workflow.cores - 1
    shell:
        "RepeatModeler -pa {threads} -database $(echo {input.xdb1} | sed 's/.translation//g') -LTRStruct &> {log}"

rule pangenome_repeat_masker:
    input:
        genome="output/panproteome_filtering/pangenome.new.fa",
        database="output/panproteome_filtering/pangenome_xdb-families.fa",
    output:
        "output/panproteome_filtering/pangenome.new.fa.masked"
    log:
        "output/logs/pangenome_repeat_masker.log"
    benchmark:
        "output/benchmarks/pangenome_repeat_masker.txt"
    container:
        "docker://harbor.containers.wurnet.nl/proxy-cache/dfam/tetools:1.5" #Using proxy to prevent too many requests to docker hub
    threads:
        workflow.cores - 1
    shell:
        "RepeatMasker -pa {threads} -lib {input.database} {input.genome} -xsmall -gff &> {log}"

rule pangenome_gene_prediction:
    input:
        genome="output/panproteome_filtering/pangenome.new.fa.masked"
    output:
        "output/panproteome_filtering/pangenome.new.gff"
    log:
        "output/logs/pangenome_gene_prediction.log"
    benchmark:
        "output/benchmarks/pangenome_gene_prediction.txt"
    conda:
        "envs/busco.yaml"
    params:
        augustus_config_path=config["augustus_config_path"],
        species=config["augustus_species"],
    shell:
        "AUGUSTUS_CONFIG_PATH={params.augustus_config_path} augustus --species={params.species} --gff3=on --softmasking=true --outfile={output} {input.genome} &> {log}"

rule panproteome_creation:
    input:
        gff="output/panproteome_filtering/pangenome.new.gff",
        fa="output/panproteome_filtering/pangenome.new.fa",
    output:
        "output/panproteome_filtering/pangenome.new.pep.fa",
    log:
        "output/logs/panproteome_creation.log"
    benchmark:
        "output/benchmarks/panproteome_creation.txt"
    conda:
        "envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl -p --cis --cfs -g {input.gff} -f {input.fa} -o {output} &> {log}"

rule panproteome_mmseqs2:
    input:
        "output/panproteome_filtering/pangenome.new.pep.fa",
    output:
        "output/panproteome_filtering/pangenome.new.pep.fa_tophit_aln",
        "output/panproteome_filtering/pangenome.new.pep.fa_tophit_report",
        "output/panproteome_filtering/pangenome.new.pep.fa_report",
        "output/panproteome_filtering/pangenome.new.pep.fa_lca.tsv",
    log:
        "output/logs/panproteome_mmseqs2.log"
    benchmark:
        "output/benchmarks/panproteome_mmseqs2.txt"
    conda:
        "envs/mmseqs2.yaml"
    params:
        db=config["mmseqs2"],
        tmp=config["mmseqs2_tmpdir"],
    threads:
        max(workflow.cores - 1, 1)
    shell:
        """
        (
        [ ! -d {params.tmp} ] && mkdir -p {params.tmp};
        rm -rd {params.tmp};
        mmseqs easy-taxonomy {input} {params.db} {input} {params.tmp} --threads {threads}
        ) &> {log}
        """

rule panproteome_contamination_visualisation:
    input:
        "output/panproteome_filtering/pangenome.new.pep.fa_lca.tsv",
    output:
        krona="output/panproteome_filtering/pangenome.new.pep.fa.krona",
        html="output/panproteome_filtering/pangenome.new.pep.fa.krona.html",
    log:
        "output/logs/panproteome_contamination_visualisation.log"
    benchmark:
        "output/benchmarks/panproteome_contamination_visualisation.txt"
    conda:
        "envs/kraken.yaml"
    shell:
        """
        (
        cut -f 1,2 {input} > {output.krona};
        ktImportTaxonomy {output.krona} -o {output.html}
        ) &> {log}
        """

rule panproteome_all_contamination:
    input:
        report="output/panproteome_filtering/pangenome.new.pep.fa_report",
        visualisation="output/panproteome_filtering/pangenome.new.pep.fa.krona.html",
    output:
        "output/panproteome_filtering/pangenome.new.pep.fa.taxid.tsv",
    log:
        "output/logs/panproteome_all_contamination.log"
    benchmark:
        "output/benchmarks/panproteome_all_contamination.txt"
    conda:
        "envs/blobtools.yaml"
    shell:
        "ete3 ncbiquery --info --search $(awk 'BEGIN{{FS = OFS = \"\\t\";}} FNR!=1{{printf \"%s \",$5;}}' {input.report} | sed 's/ $/\\n/g') > {output} 2> {log}"

rule panproteome_filter_Eukaryota:
    input:
        taxid="output/panproteome_filtering/pangenome.new.pep.fa.taxid.tsv",
        mmseqs2="output/panproteome_filtering/pangenome.new.pep.fa_lca.tsv",
        gff="output/panproteome_filtering/pangenome.new.gff",
        fa="output/panproteome_filtering/pangenome.new.fa",
    output:
        clean_mrna_id="output/panproteome_filtering/pangenome.new.pep.fa.Eukaryota_and_0.txt",
        clean_seq_id="output/panproteome_filtering/pangenome.new.Eukaryota_and_0.txt",
        clean_fa="output/pangenome.new.fa",
        clean_gff="output/pangenome.new.unchecked.gff",
    log:
        "output/logs/panproteome_filter_Eukaryota.log"
    benchmark:
        "output/benchmarks/panproteome_filter_Eukaryota.txt"
    conda:
        "envs/seq.yaml"
    shell:
        """
        (
        awk 'BEGIN{{FS = OFS = \"\\t\";}} FNR==NR{{taxid[$1]=$4; next;}} $2==0 || taxid[$2]~/,Eukaryota/{{print $1;}}' {input.taxid} {input.mmseqs2} | sort > {output.clean_mrna_id};
        awk 'BEGIN{{FS = OFS = \"\\t\";}} FNR==NR{{id[$1]=1; next;}} /^[^#]/ && $3=="transcript"{{split($9,a,\";\"); split(a[1],b,\"=\"); if (id[b[2]]==1){{print $1;}}}}' {output.clean_mrna_id} {input.gff} | sort | uniq > {output.clean_seq_id};
        awk 'FNR==NR && /^[^#]/{{annotation[$1]=1;}} FNR==NR{{next;}} /^>/{{id=substr($1,2); if (annotation[id]!=1){{print id;}}}}' {input.gff} {input.fa} >> {output.clean_seq_id};
        seqkit grep -f {output.clean_seq_id} {input.fa} > {output.clean_fa};
        awk 'BEGIN{{FS = OFS = \"\\t\";}} FNR==NR{{contamination[$1]=1; next;}} $3==\"transcript\"{{$3=\"mRNA\"}} /^[^#]/ && contamination[$1]==1{{print;}}' {output.clean_seq_id} {input.gff} > {output.clean_gff};
        ) &> {log}
        """

rule pangenome_final:
    input:
        genome_fa="output/data/genome.fa",
        new_fa="output/pangenome.new.fa",
    output:
        fa="output/pangenome.fa",
    log:
        "output/logs/pangenome_final.log"
    benchmark:
        "output/benchmarks/pangenome_final.txt"
    shell:
        """
        (
        cat {input.genome_fa} {input.new_fa} > {output.fa};
        ) &> {log}
        """

rule pangenome_genome:
    input:
        "output/pangenome.fa",
    output:
        fai="output/pangenome.fa.fai",
        tsv="output/pangenome.tsv",
    log:
        "output/logs/pangenome_genome.log"
    benchmark:
        "output/benchmarks/pangenome_genome.txt"
    conda:
        "envs/alignment.yaml"
    shell:
        """
        (
        samtools faidx {input};
        awk 'BEGIN{{FS = OFS = \"\\t\";}} {{print $1,$2;}}' {output.fai} > {output.tsv};
        ) &> {log}
        """

rule pangenome_fix_annotation:
    input:
        "output/pangenome.new.unchecked.gff"
    output:
        "output/pangenome.new.gff"
    log:
        "output/logs/pangenome_fix_annotation.log"
    benchmark:
        "output/benchmarks/pangenome_fix_annotation.txt"
    conda:
        "envs/agat.yaml"
    shell:
        "agat_convert_sp_gxf2gxf.pl -g {input} -o {output} &> {log}"

rule pangenome_annotation:
    input:
        genome_annotation="output/data/genome.sorted.gff",
        new_annotation="output/pangenome.new.gff",
        sorting="output/pangenome.tsv",
    output:
        gff="output/pangenome.gff",
        table="output/pangenome.id_to_type.tsv",
    log:
        "output/logs/pangenome_annotation.log"
    benchmark:
        "output/benchmarks/pangenome_annotation.txt"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        (
        cat {input.genome_annotation} {input.new_annotation} | grep -vE '^#' > {output.gff}.tmp;
        bedtools sort -g {input.sorting} -i {output.gff}.tmp > {output.gff};
        rm {output.gff}.tmp;
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"ID\",\"type\";}} /^[^#]/{{split($9,a,\";\"); for (i in a){{if (a[i]~/^ID/){{split(a[i],b,\"=\"); id=b[2];}}}}; print id,$3;}}' {output.gff} > {output.table};
        ) &> {log}
        """

rule busco_pangenome:
    input:
        "output/pangenome.fa"
    output:
        "output/busco/pangenome/short_summary.specific.{lineage}.pangenome.txt"
    log:
        "output/logs/busco_pangenome/{lineage}.log"
    benchmark:
        "output/benchmarks/busco_pangenome/{lineage}.txt"
    conda:
        "envs/busco.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "busco -f -i {input} -m geno --out_path $(dirname $(dirname {output})) -o $(basename $(dirname {output})) -l {wildcards.lineage} -c {threads} --tar &> {log}"

rule pangenome_bwa_index:
    input:
        "output/pangenome.fa"
    output:
        "output/pangenome.fa.0123",
        "output/pangenome.fa.amb",
        "output/pangenome.fa.ann",
        "output/pangenome.fa.pac"
    log:
        "output/logs/pangenome_bwa_index.log"
    benchmark:
        "output/benchmarks/pangenome_bwa_index.txt"
    conda:
        "envs/alignment.yaml"
    shell:
        "bwa-mem2 index -p {input} {input} &> {log}"

rule pangenome_bwa_mapping:
    input:
        index="output/pangenome.fa",
        index1="output/pangenome.fa.0123",
        index2="output/pangenome.fa.amb",
        index3="output/pangenome.fa.ann",
        index4="output/pangenome.fa.pac",
        forward="output/trimmed_reads/{accession}/{sample}_1.fq.gz",
        backward="output/trimmed_reads/{accession}/{sample}_2.fq.gz",
    output:
        #temp("output/pangenome_alignment/{accession}/{sample}.bam")
        "output/pangenome_alignment/{accession}/{sample}.bam"
    log:
        "output/logs/pangenome_alignment/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/pangenome_alignment/{accession}/{sample}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "(bwa-mem2 mem -t {threads} {input.index} {input.forward} {input.backward} | samtools sort -@ $(({threads}-1)) -o {output} -) &> {log}"

rule pangenome_alignment_index:
    input:
        "output/pangenome_alignment/{accession}/{sample}.bam"
    output:
        #temp("output/pangenome_alignment/{accession}/{sample}.bam.bai")
        "output/pangenome_alignment/{accession}/{sample}.bam.bai"
    log:
        "output/logs/pangenome_alignment_index/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/pangenome_alignment_index/{accession}/{sample}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools index -@ $(({threads}-1)) {input} > {output} 2> {log}"

rule pangenome_alignment_stats:
    input:
        bam="output/pangenome_alignment/{accession}/{sample}.bam",
        bai="output/pangenome_alignment/{accession}/{sample}.bam.bai"
    output:
        "output/pangenome_alignment/{accession}/{sample}.txt"
    log:
        "output/logs/pangenome_alignment_stats/{accession}/{sample}.log"
    benchmark:
        "output/benchmarks/pangenome_alignment_stats/{accession}/{sample}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools stats -@ $(({threads}-1)) {input.bam} > {output} 2> {log}"

def get_pangenome_multiqc_input():
    """This function returns all fastqc, alignment stats and duplication stats"""
    combinations = []
    for accession in list(config["data"].keys()):
        combination=expand(f"{accession}/{{sample}}", sample=[x for x in config["data"][accession]])
        combinations.append(combination)
    combinations = [x for y in combinations for x in y]

    allfiles = []
    allfiles.append(expand("output/fastQC_before_trimming/{combination}_{nr}_fastqc.html", combination=combinations, nr=[1,2]))
    allfiles.append(expand("output/fastQC_after_trimming/{combination}_{nr}_fastqc.html", combination=combinations, nr=[1,2]))
    allfiles.append(expand("output/pangenome_alignment/{combination}.txt", combination=combinations))

    return allfiles

rule pangenome_multiqc:
    input:
        *get_pangenome_multiqc_input()
    output:
        "output/pangenome_multiqc/multiqc_report.html"
    log:
        "output/logs/pangenome_multiqc.log"
    benchmark:
        "output/benchmarks/pangenome_multiqc.txt"
    conda:
        "envs/qualitycontrol.yaml"
    shell:
        "multiqc -f -s -d `dirname {input} | uniq` -n {output} &> {log}"

rule pangenome_alignment_merge:
    input:
        bam=lambda wildcards: expand("output/pangenome_alignment/{{accession}}/{sample}.bam", sample=[x for x in config["data"][wildcards.accession]]),
        bai=lambda wildcards: expand("output/pangenome_alignment/{{accession}}/{sample}.bam.bai", sample=[x for x in config["data"][wildcards.accession]]),
    output:
        #temp("output/pangenome_alignment_merged/{accession}.bam")
        "output/pangenome_alignment_merged/{accession}.bam"
    log:
        "output/logs/pangenome_alignment_merge/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_alignment_merge/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools merge -@ $(({threads}-1)) {output} {input.bam} 2> {log}"

rule pangenome_alignment_merge_index:
    input:
        "output/pangenome_alignment_merged/{accession}.bam"
    output:
        #temp("output/pangenome_alignment_merged/{accession}.bam.bai")
        "output/pangenome_alignment_merged/{accession}.bam.bai"
    log:
        "output/logs/pangenome_alignment_merge_index/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_alignment_merge_index/{accession}.txt"
    conda:
        "envs/alignment.yaml"
    threads:
        10
    shell:
        "samtools index -@ $(({threads}-1)) {input} 2> {log}"

rule pangenome_calculate_breadth:
    input:
        genome="output/pangenome.tsv",
        annotation="output/pangenome.gff",
        bam="output/pangenome_alignment_merged/{accession}.bam",
        bai="output/pangenome_alignment_merged/{accession}.bam.bai",
    output:
        allcov="output/pangenome_coverage/breadth/{accession}.all.bed",
        exoncov="output/pangenome_coverage/breadth/{accession}.exons.tsv",
    log:
        "output/logs/pangenome_calculate_breadth/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_calculate_breadth/{accession}.txt"
    conda:
        "envs/bedtools.yaml"
    threads:
        5
    shell:
        """
        (
        bedtools coverage -g {input.genome} -sorted -bed -a {input.annotation} -b {input.bam} 1> {output.allcov};
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"parent\",\"number exons\",\"total length exons\",\"total length exons covered\",\"fraction length exons covered\";}} $3==\"exon\"{{split($9,a,\";\"); for (i in a){{if (a[i]~/^Parent/){{split(a[i],tmp,\"=\"); parent=tmp[2];}}}} num[parent]++; len[parent]+=$12; covered[parent]+=$11;}} END{{for (parent in len){{print parent,num[parent],len[parent],covered[parent],covered[parent]/len[parent];}}}}' {output.allcov} > {output.exoncov};
        ) &> {log}
        """

rule pangenome_calculate_depth:
    input:
        genome="output/pangenome.tsv",
        annotation="output/pangenome.gff",
        bam="output/pangenome_alignment_merged/{accession}.bam",
        bai="output/pangenome_alignment_merged/{accession}.bam.bai",
    output:
        allcov="output/pangenome_coverage/depth/{accession}.all.bed",
        exoncov="output/pangenome_coverage/depth/{accession}.exons.tsv",
    log:
        "output/logs/pangenome_calculate_depth/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_calculate_depth/{accession}.txt"
    conda:
        "envs/bedtools.yaml"
    threads:
        5
    shell:
        """
        (
        bedtools coverage -g {input.genome} -mean -sorted -bed -a {input.annotation} -b {input.bam} 1> {output.allcov};
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"parent\",\"number exons\",\"total length exons\",\"mean depth exons\";}} $3==\"exon\"{{split($9,a,\";\"); for (i in a){{if (a[i]~/^Parent/){{split(a[i],tmp,\"=\"); parent=tmp[2];}}}} num[parent]++; len[parent]+=($5-$4+1); normdepth[parent]+=($10*($5-$4+1));}} END{{for (parent in len){{print parent,num[parent],len[parent],normdepth[parent]/len[parent];}}}}' {output.allcov} > {output.exoncov};
        ) &> {log}
        """

rule pangenome_combine_breadth_depth:
    input:
        breadth="output/pangenome_coverage/breadth/{accession}.exons.tsv",
        depth="output/pangenome_coverage/depth/{accession}.exons.tsv",
    output:
        "output/pangenome_coverage/combined/{accession}.tsv"
    log:
        "output/logs/pangenome_combine_breadth_depth/{accession}.log"
    benchmark:
        "output/benchmarks/pangenome_combine_breadth_depth/{accession}.txt"
    shell:
        """
        (
        awk 'BEGIN{{FS = OFS = \"\\t\";}} NR==FNR{{a[$1] = $4; next;}} NR!=1{{print $0,a[$1];}}' {input.depth} {input.breadth} > {output}
        ) &> {log}
        """

rule pangenome_calculate_pav:
    input:
        pav=expand("output/pangenome_coverage/combined/{accession}.tsv", accession=config["data"]),
        conversion="output/pangenome.id_to_type.tsv",
    output:
        depth="output/pangenome.depth.tsv",
        breadth="output/pangenome.breadth.tsv",
        plot="output/pangenome.pav.pdf",
    log:
        "output/logs/pangenome_calculate_pav.log"
    benchmark:
        "output/benchmarks/pangenome_calculate_pav.txt"
    conda:
        "envs/rbase.yaml"
    shell:
        "../scripts/combine_pav.R {output.depth} {output.breadth} {output.plot} {input.conversion} {input.pav} &> {log}"
