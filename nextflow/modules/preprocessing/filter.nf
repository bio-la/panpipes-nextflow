    #!/usr/bin/env nextflow
nextflow.enable.dsl=2
process filter {

    input:
    path input_from_A

    output:
    path "outputB.txt", emit: resultB

    script:
    """
    # Simulate processing
    cp $input_from_A outputB.txt
    """
}