configfile: "config_23andMe.yaml"
import os
import json
import gzip
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider

HTTP = HTTPRemoteProvider()
DBox_dist = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["analyst"])
DBox_dist_public = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["public"])

wildcard_constraints:
    version="[\d.]+"

include: "rules/ldsc.smk"
include: "rules/meta.smk"
include: "rules/distribution.smk"
include: "rules/vcf.smk"
include: "rules/reports.smk"
include: "rules/finemapping.smk"