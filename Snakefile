configfile: "config.yaml"
import os
import json
import gzip
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider

HTTP = HTTPRemoteProvider()
DBox_dist = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["full"])
DBox_dist_no23am = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["no23am"])



include: "rules/meta.smk"
include: "rules/vcf.smk"