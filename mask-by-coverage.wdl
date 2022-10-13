version 1.0

task make_mask_file {
	input {
		File sam
		Int min_coverage

		# runtime attributes
		Int addldisk = 250
		Int cpu      = 16
		Int retries  = 1
		Int memory   = 32
		Int preempt  = 1
	}
	String basestem = basename(sam, ".sam")
	Int finalDiskSize = ceil(size(sam, "GB")) + addldisk
	
	command <<<
	set -eux pipefail
	cp ~{sam} .
	samtools sort -u ~{basestem}.sam > sorted_u_~{basestem}.sam
	bedtools genomecov -ibam sorted_u_~{basestem}.sam -bga | \
		awk '$4 < ~{min_coverage}' > \
		~{basestem}_below_~{min_coverage}x_coverage.bga
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
		File mask_file = glob("*coverage.bga")[0]
	}

}


workflow MaskByCoverage {
	input {
		Array[File] sams
		Int min_coverage = 1
	}

	scatter(sam in sams) {
		call make_mask_file {
			input:
				sam = sam,
				min_coverage = min_coverage
		}
	}

	output {
		Array[File] mask_files = make_mask_file.mask_file
	}
}
