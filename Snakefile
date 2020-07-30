configfile: "config.yaml"

from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
DBox = DropboxRemoteProvider(oauth2_access_token=config["remote"]["dropbox"])

include: "rules/meta.smk"