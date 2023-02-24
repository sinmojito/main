#!/usr/bin/env nextflow

// Define input parameters
params.inputFile = null // Path to input FASTA file
params.cutoff = null    // Sequence length cutoff

// Define command line parameters
inputFileVal = file(args.inputFile)
cutoffVal = args.cutoff.toInteger()

// Check if input parameters have been provided
if (!inputFileVal || !cutoffVal) {
    println "Please provide inputFile and cutoff as parameters"
    System.exit(1)
}

// Create path channel to input FASTA file
inputFasta = file(inputFileVal)

// Define process to filter sequences longer than cutoff
process filterSequences {
    input:
    path inputFasta    // Input FASTA file
    val cutoff         // Sequence length cutoff

    output:
    path "output.txt"  // Filtered output file

    script:
        """
        #!/usr/bin/env python3
        from Bio import SeqIO

        fasta_file = '${inputFasta}'
        cutoff = "${cutoff}"
        # Iterate over input sequences
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                # Write sequences longer than cutoff to output file
                if len(record.seq) > cutoff:
                    with open('output.txt', 'a') as out:
                        out.write(f'>{record.id}\\n{record.seq}\\n')
        """

}

// Define a workflow that includes the process
workflow {
    filterSequences(inputFasta, cutoffVal)
}
