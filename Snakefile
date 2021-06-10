configfile: "config_23andMe.yaml"
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

include: "rules/resources.smk"
include: "rules/ldsc.smk"
include: "rules/meta.smk"
include: "rules/meta_carpa.smk"
include: "rules/metaqc.smk"
include: "rules/distribution.smk"
include: "rules/vcf.smk"
include: "rules/reports.smk"
include: "rules/finemapping.smk"