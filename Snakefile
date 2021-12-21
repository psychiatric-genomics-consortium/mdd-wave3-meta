configfile: "config.yaml"
import os
import json
import gzip
import itertools
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider

HTTP = HTTPRemoteProvider()
DBox_dist = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["analyst"])
DBox_dist_public = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["public"])

wildcard_constraints:
    version="[\d.]+"
    
# current European ancestries analysis
# analysis version format: v3.[PGC Cohorts Count].[Other Cohorts Count].[Revision]
analysis_version_eur = ["3.49.24.09"]
analysis_version = analysis_version_eur

# current East Asian ancestries analysis
# analysis version format: v3.[PGC Cohorts Count].[Other Cohorts Count].[Revision]
analysis_version_eas = ["3.00.04"]

include: "rules/resources.smk"
include: "rules/ldsc.smk"
include: "rules/vcf.smk"
include: "rules/meta.smk"
include: "rules/meta_ricopili.smk"
include: "rules/meta_carpa.smk"
include: "rules/meta_mtag.smk"
include: "rules/meta_metal.smk"
include: "rules/meta_gsem.smk"
include: "rules/metaqc.smk"
include: "rules/meta_post.smk"
include: "rules/cojo.smk"
include: "rules/distribution.smk"
include: "rules/finemapping.smk"
include: "rules/open_gwas.smk"
include: "rules/reports.smk"

