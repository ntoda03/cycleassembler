
///////////////////////////////////////////////////////////////////////////////
/*                                                                           
            read_fastq Get fastqs into a channel                   

 inputs:
 reads - Path to input data (must be surrounded with quotes), can be gzipped
          Can include a semicolon separate list of files.
          For multiple files use the * wildcard.
          Indicate read pairs forward and reverse with {1,2}

 output:
 reads_ch - a reads channel containing all the files found
                                                                             */
///////////////////////////////////////////////////////////////////////////////

def read_fastq(reads){

    file_list = reads.tokenize(';')
    file_list.each{ if( ! file(it) ){ exit 1, "Error: $it does not exist!"} }

    reads_ch = Channel.fromFilePairs( file_list, checkIfExists: true, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\n" }
        //.map { get_meta_format(it) }

    return reads_ch
}
