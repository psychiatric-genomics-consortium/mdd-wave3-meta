configfile: "config.yaml"

from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
DBox = DropboxRemoteProvider(oauth2_access_token=config["remote"]["dropbox"]["full"])
DBox_no23am = DropboxRemoteProvider(oauth2_access_token=config["remote"]["dropbox"]["no23am"])

include: "rules/meta.smk"