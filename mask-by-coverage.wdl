version 1.0

task make_mask_file {
	input {
		File sam
		Int min_coverage = 1

		# runtime attributes
		Int addldisk = 250
		Int cpu      = 16
		Int retries  = 1
		Int memory   = 32
		Int preempt  = 1
	}
	String basestem = basename(sam, ".sam")

	command <<<
	samtools sort -u ~{sam} > sorted_~{basestem}.sam
	bedtools genomecov -ibam sorted_u_SAMEA2534421.sam -bga | awk '$4 < ~{min_coverage}' > low_coverage.bga
	>>>

	runtime {
		cpu: cpu
		docker: "ashedpotatoes/iqbal-unofficial-clockwork-mirror:latest"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	output {
		File mask_file = "low_coverage.bga"
	}

}

#task bedtools_maskfasta {
#	input {
#		File mask_file
#		Array[File] fastas
#
#		# runtime attributes
#		Int addldisk = 250
#		Int cpu      = 16
#		Int retries  = 1
#		Int memory   = 32
#		Int preempt  = 1
#	}
#
#	command <<<
#	bedtools maskfasta
#	>>>
#
#	runtime {
#		cpu: cpu
#		docker: "ashedpotatoes/iqbal-unofficial-clockwork-mirror:latest"
#		disks: "local-disk " + finalDiskSize + " HDD"
#		maxRetries: "${retries}"
#		memory: "${memory} GB"
#		preemptible: "${preempt}"
#	}
#
#	output {
#		File mask_file = "low_coverage.bga"
#	}
#}

task seqtk_mask {
	input {
		
	}

	command <<<
	seqtk seq -M reg.bed in.fa > out.fa
	>>>
}







biocontainers/seqtk:v1.3-1-deb_cv1


workflow MaskByCoverage {
	input:
		Array[Array[File]] fastqs
		Array[File] sam

	scatter(zip(samples, fastqs))
		call depth {
			input:
				sam = sam
		}
	
}
