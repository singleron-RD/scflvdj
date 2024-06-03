def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "scrna/stats": {"fn": "*scrna.*stats.json"},
        "scrna/umi_count": {
            "fn": "*scrna.umi_count.json",
        },
        "scrna/saturation": {
            "fn": "*scrna.saturation.json",
        },
        "scrna/median_gene": {
            "fn": "*scrna.median_gene.json",
        },
        "scflvdj/stats": {"fn": "*scflvdj.*stats.json"},
        "scflvdj/gene": {"fn": "*scflvdj.gene.json"},
        "scflvdj/count": {"fn": "*scflvdj.count.json"},
        "scflvdj/meta": {"fn": "*scflvdj.meta.json"},
    }
    config.update_dict(config.sp, sgr_search_patterns)
