#!/usr/bin/env nextflow

params.inputFile = ''
params.cutoff = ''

if (!params.inputFile || !params.cutoff) {
    println "Please provide inputFile and cutoff as parameters"
    System.exit(1)
}

inputFasta = file(params.inputFile)

process filterSequences {
    input:
    path inputFasta
    val cutoff

    output:
    path "output.txt"

    script:
        """
        #!/usr/bin/env python3
        import sys
        from Bio import SeqIO

        # Get input parameters from Nextflow
        fasta_file = '$inputFasta'
        cutoff = $cutoff

        # Iterate over input sequences
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                # Write sequences longer than cutoff to output file
                if len(record.seq) > cutoff:
                    with open('output.txt', 'a') as out:
                        out.write(f'>{record.id}\\n{record.seq}\\n')
        """
}

workflow {
    filterSequences(inputFasta, params.cutoff)
}
