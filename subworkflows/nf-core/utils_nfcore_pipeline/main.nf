/*
    subworkflows/nf-core/utils_nfcore_pipeline/main.nf
    Minimal local stand-in for the nf-core utils subworkflow. Provides
    softwareVersionsToYAML used to collate versions.yml channels.

    For a full nf-core setup, install the real subworkflow with:
        nf-core subworkflows install utils_nfcore_pipeline
    and remove this file.
*/

// Collate a channel of versions.yml file paths into unique YAML lines.
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
        .unique()
        .map { version_file -> version_file.text }
        .unique()
}
