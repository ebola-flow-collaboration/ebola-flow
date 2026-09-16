process CHECK_NEW_SUBMISSIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/curl:8.5.0' :
        'quay.io/biocontainers/curl:8.5.0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("has_submissions.txt"), emit: proceed, optional: true
    tuple val(meta), path("README.md")          , emit: no_submissions_report, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: [:]
    def organism = args.organism ?: meta.organism ?: 'ebola-bdbv'
    def today    = new Date().format('yyyy-MM-dd')
    """
    #!/usr/bin/env python3

    import json
    import sys
    import urllib.request
    from datetime import date
    from platform import python_version

    organism = '${organism}'
    today_str = '${today}'
    today_date = date.fromisoformat(today_str)

    # Fetch today's submissions
    url_today = (
        f"https://lapis.pathoplexus.org/{organism}/sample/details"
        f"?submittedDate={today_str}&versionStatus=LATEST_VERSION"
    )
    with urllib.request.urlopen(url_today) as resp:
        today_data = json.loads(resp.read().decode())

    count = len(today_data.get("data", []))

    if count > 0:
        with open("has_submissions.txt", "w") as out:
            out.write(f"{count} new submissions found\\n")
    else:
        # Fetch all submissions to find the most recent date
        url_all = (
            f"https://lapis.pathoplexus.org/{organism}/sample/details"
            f"?versionStatus=LATEST_VERSION"
        )
        with urllib.request.urlopen(url_all) as resp:
            all_data = json.loads(resp.read().decode())

        records = all_data.get("data", [])
        submitted_dates = [
            rec.get("submittedDate", "")
            for rec in records
            if rec.get("submittedDate")
        ]

        if submitted_dates:
            last_date_str = max(submitted_dates)
            last_date = date.fromisoformat(last_date_str)
            days_ago = (today_date - last_date).days
            last_date_display = last_date.strftime("%B %d, %Y")
        else:
            last_date_str = "unknown"
            days_ago = "unknown"
            last_date_display = "unknown"

        lines = [
            "# No New Ebola Bundibugyo Virus (EBDBV) Submissions",
            "",
            f"**Date checked:** {today_str}",
            "",
            "No new Ebola Bundibugyo Virus (EBDBV) sequences were submitted to "
            f"[Pathoplexus](https://pathoplexus.org/ebola-bdbv/search) on {today_str}.",
            "",
            f"**Last submission date:** {last_date_display}",
            "",
            f"**Days since last submission:** {days_ago}",
            "",
            "The pipeline did not proceed as there was no new data to process.",
            "",
        ]
        with open("README.md", "w") as f:
            f.write("\\n".join(lines))

    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"CHECK_NEW_SUBMISSIONS":\\n')
        f.write(f'    python: {python_version()}\\n')

    sys.exit(0)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    from platform import python_version

    with open("has_submissions.txt", "w") as f:
        f.write("stub\\n")

    with open('versions.yml', 'w') as f:
        f.write('"CHECK_NEW_SUBMISSIONS":\\n')
        f.write(f'    python: {python_version()}\\n')
    """
}
