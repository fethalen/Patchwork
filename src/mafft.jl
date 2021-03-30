# Wrapper for running MAFFT from the command-line.

""" 
    mafft(input, flags)

Runs mafft on `input` of the type `MultipleSequenceAlignment`. May include
optional `flags` such as `["--thread", 4]`.
"""
function mafft(input::MultipleSequenceAlignment, flags=[])::MultipleSequenceAlignment
    input_fasta = mktemp_fasta(input, removehyphens=true)
    output, io = mktemp()
    mafft_cmd = pipeline(`mafft $flags $input_fasta`, stdout=output)
    run(mafft_cmd)
    close(io)
    return readmsa(output, '@')
end

""" 
    mafft_linsi(input, flags)

Runs mafft with an accurate option (L-INS-i) on `input` of the type
`MultipleSequenceAlignment`. May include optional `flags` such as
`["--thread", 4]`. `flags` may not include `"--maxiterate"` and or 
`"--localpair"`, since these options are already set.
"""
function mafft_linsi(input::MultipleSequenceAlignment, flags=[])
    input_fasta = mktemp_fasta(input, removehyphens=true)
    output, io = mktemp()
    mafft_cmd = pipeline(`mafft --maxiterate 1000 --localpair $flags $input_fasta`, 
                         stdout=output)
    run(mafft_cmd)
    close(io)
    return readmsa(output, '@')
end
