configfile: "config.yaml"
import os
import json
import gzip
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider

HTTP = HTTPRemoteProvider()
DBox_dist = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["analyst"])
DBox_dist_no23am = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["public"])

include: "rules/meta.smk"
include: "rules/vcf.smk"
include: "rules/reports.smk"
include: "rules/twas.smk"
