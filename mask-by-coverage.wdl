version 1.0

task make_mask_file {
	input {
		File bam
		Int min_coverage

		# runtime attributes
		Int addldisk = 250
		Int cpu      = 16
		Int retries  = 1
		Int memory   = 32
		Int preempt  = 1
	}
	String basestem = basename(bam, ".bam")
	Int finalDiskSize = ceil(size(bam, "GB")) + addldisk
	
	command <<<
	set -eux pipefail
	cp ~{bam} .
	samtools sort -u ~{basestem}.bam > sorted_u_~{basestem}.bam
	bedtools genomecov -ibam sorted_u_~{basestem}.bam -bga | \
		awk '$4 < ~{min_coverage}' > \
		~{basestem}_below_~{min_coverage}x_coverage.bedgraph
	>>>

	runtime {
		cpu: cpu
		docker: "ashedpotatoes/sranwrp:1.0.2"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	output {
		File mask_file = glob("*coverage.bedgraph")[0]
	}

}


workflow MaskByCoverage {
	input {
		Array[File] bams
		Int min_coverage = 1
	}

	scatter(bam in bams) {
		call make_mask_file {
			input:
				bam = bam,
				min_coverage = min_coverage
		}
	}

	output {
		Array[File] mask_files = make_mask_file.mask_file
	}
}
