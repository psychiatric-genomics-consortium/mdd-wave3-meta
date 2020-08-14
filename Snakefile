configfile: "config.yaml"

from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
DBox_dist = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["full"])
DBox_dist_no23am = DropboxRemoteProvider(oauth2_access_token=config["remote"]["distribution"]["no23am"])

include: "rules/meta.smk"