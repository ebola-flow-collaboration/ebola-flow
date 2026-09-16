process TREEANNOTATOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beast2:2.7.7--haf2fa61_0' :
        'quay.io/biocontainers/beast2:2.7.7--haf2fa61_0' }"

    input:
    tuple val(meta), path(trees)

    output:
    tuple val(meta), path("*.tree"), emit: tree
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args     ?: ''
    prefix       = task.ext.prefix   ?: "${meta.id}"
    def burnin   = task.ext.burnin   != null ? "-burnin ${task.ext.burnin}" : '-burnin 0'
    def topology = task.ext.topology ?: 'CCD0Approx'
    def low_mem  = task.ext.lowMem   ? '-lowMem true' : ''
    // Size the JVM heap to the task's allocated memory (leaving ~4 GB headroom).
    // The BEAST treeannotator launcher shell wrapper HARDCODES -Xmx8g and ignores
    // both the BEAST_XMX environment variable and any -Xmx passed on the command
    // line (the latter causes a parse error). The only reliable way to set the
    // heap is to bypass the wrapper and invoke `java` directly with our own -Xmx,
    // using the same TreeAnnotatorLauncher class the wrapper would have used.
    def avail_mb = task.memory ? (task.memory.toMega() as long) : 8192
    def xmx_mb   = Math.max(1024, avail_mb - 4096)
    """
    # Resolve the real BEAST install dir from the treeannotator wrapper.
    BEAST_BIN=\$(readlink -f "\$(command -v treeannotator)")
    BEAST_HOME=\$(dirname "\$(dirname "\${BEAST_BIN}")")

    # Prefer the JRE bundled with BEAST; fall back to the PATH java.
    JAVA="\${BEAST_HOME}/jre/bin/java"
    if [ ! -x "\${JAVA}" ]; then
        JAVA="java"
    fi

    "\${JAVA}" -Dlauncher.wait.for.exit=true -Xss256m -Xmx${xmx_mb}m \\
        -Djava.library.path="\${BEAST_HOME}/lib" -Duser.language=en \\
        -cp "\${BEAST_HOME}/lib/launcher.jar" \\
        beast.pkgmgmt.launcher.TreeAnnotatorLauncher \\
            ${burnin} \\
            -topology ${topology} \\
            ${low_mem} \\
            ${args} \\
            ${trees} \\
            ${prefix}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        treeannotator: \$(treeannotator -help 2>&1 | head -n1 | sed 's/.*v//' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        treeannotator: \$(treeannotator -help 2>&1 | head -n1 | sed 's/.*v//' || echo 'unknown')
    END_VERSIONS
    """
}
