workdir: config["out_dir"]

def get_read1(wildcards):
    return config["reads"][wildcards.sample_id][0]

def get_read2(wildcards):
    return config["reads"][wildcards.sample_id][1]

if "fasta" in config:
    if config["fasta"].lower() != "no":
        file_type = "f"
        seq_type = "fa"
    else:
        file_type = "q"
        seq_type = "fq"
else:
    file_type = "q"
    seq_type = "fq"

def get_reads_Trinity(LR):
    if LR == "left":
        sample_index = 0
    elif LR == "right":
        sample_index = 1
    else:
        sys.exit() # replace with message
    file_list = []
    for sample_id in config["reads"]:
        file_list.append(config["reads"][sample_id][sample_index])
    return ",".join(file_list)


rule all:
    input:
        [
          #"trinity_out_dir/Trinity.fasta",
          #expand("rsem_trinity_counts/{sample_id}/{sample_id}.isoforms.results",sample_id=config["reads"].keys())
          "dammit_out/dammit.done"
          expand("trinity_abundance/{sample_id}/quant.sf.genes",sample_id = config["reads"].keys())
        ]

rule dammit:
    input:
        trinity_assembly = "trinity_out_dir/Trinity.fasta",
        db = config["reference_p"]
    output:
        done_flag = touch("dammit_out/dammit.done")
    threads : 20
    params:
        dammit_out_dir = "dammit_out"
    conda: "env/dammit.yaml"
    shell:
        """
        dammit databases --install --quick --busco-group metazoa &&
        dammit annotate  {input.trinity_assembly} --quick --user-databases {input.db} -e 1e-10 --output-dir {params.dammit_out_dir} --n_threads {threads}
        """

rule index_bowtie2_trinity:
    input:
        reference = "trinity_out_dir/Trinity.fasta"
    output:
        touch('bowtie2_trinity_index/index.done')
    params:
        index_basename="bowtie2_trinity_index/index"
    conda:
        "../env/assembly_map_quant.yml"
    shell:
        "bowtie2-build {input.reference} {params.index_basename}"

rule align_bowtie2_trinity:
    input:
        reference_flag = "bowtie2_trinity_index/index.done",
        reads1 = get_read1,
        reads2 = get_read2
    output:
        "bowtie2_trinity_alignments/{sample_id}.bam"
    threads: workflow.cores/len(config["reads"])
    benchmark:
        "benchmarks/{sample_id}/bowtie2.txt"
    params:
        index_basename=rules.index_bowtie2_trinity.params.index_basename,
        file_type = file_type
    conda:
        "../env/assembly_map_quant.yml"
    shell:
        "bowtie2 -p {threads} -x {params.index_basename} -{params.file_type} -1 {input.reads1} -2 {input.reads2} --no-mixed --no-discordant --gbar 100 --end-to-end -k 200 | samtools view -bS - > {output}"

rule rsem_prepare_reference_trinity:
    input:
        reference="trinity_out_dir/Trinity.fasta"
    output:
        touch('rsem_trinity_index/index.done')
    params:
        rsem_index_basename="rsem_trinity_index/index"

    conda:
        "../env/assembly_map_quant.yml"
    shell:
        "rsem-prepare-reference {input.reference} {params.rsem_index_basename}"

rule rsem_trinity:
    input:
        reference = "trinity_out_dir/Trinity.fasta",
        alignments="bowtie2_trinity_alignments/{sample_id}.bam",
        index_flag="rsem_trinity_index/index.done"
    params:
        rsem_index_basename=rules.rsem_prepare_reference_trinity.params.rsem_index_basename,
        rsem_out_dir="rsem_trinity_counts/{sample_id}/{sample_id}"
    output:
        "rsem_trinity_counts/{sample_id}/{sample_id}.isoforms.results",
        "rsem_trinity_counts/{sample_id}/{sample_id}.genes.results"
    conda:
        "../env/assembly_map_quant.yml"
    shell:
        "rsem-calculate-expression -q --no-bam-output --alignments --paired-end {input.alignments} {params.rsem_index_basename} {params.rsem_out_dir}"

rule abundance_est:
    output:
        "trinity_abundance/{sample_id}/quant.sf",
        "trinity_abundance/{sample_id}/quant.sf.genes"
    input:
        assembly = "trinity_out_dir/Trinity.fasta",
        reads1 = get_read1,
        reads2 = get_read2
    params:
        seq_type = seq_type,
        trinity_abundance_out_dir = "trinity_abundance/{sample_id}/"
    threads: 5
    conda:
        "../env/assembly_map_quant.yml"

    shell:
        "$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.assembly} --left {input.reads1} --right {input.reads2} --seqType {params.seq_type} --est_method salmon --output_dir {params.trinity_abundance_out_dir} --thread_count {threads} --trinity_mode --prep_reference"

rule trinity:
    output:
      "trinity_out_dir/Trinity.fasta"
    threads: 20
    params:
        seq_type = seq_type,
        left = get_reads_Trinity("left"),
        right = get_reads_Trinity("right")
    benchmark:
        "benchmarks/assembly.txt"
    conda:
        "../env/assembly_map_quant.yml"
    shell:
        "Trinity --seqType {params.seq_type} --max_memory 90G --left {params.left} --right {params.right} --CPU {threads}"
